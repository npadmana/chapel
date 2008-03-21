/*** normalize
 ***
 *** This pass and function normalizes parsed and scope-resolved AST.
 ***/

#include "astutil.h"
#include "build.h"
#include "expr.h"
#include "passes.h"
#include "runtime.h"
#include "stmt.h"
#include "stringutil.h"
#include "symbol.h"
#include "symscope.h"

bool normalized = false;
Vec<const char*> usedConfigParams;

static void change_method_into_constructor(FnSymbol* fn);
static void normalize_returns(FnSymbol* fn);
static void call_constructor_for_class(SymExpr* se);
static void hack_resolve_types(Expr* expr);
static void apply_getters_setters(FnSymbol* fn);
static void insert_call_temps(CallExpr* call);
static void fix_user_assign(CallExpr* call);
static void fix_def_expr(VarSymbol* var);
static void tag_global(FnSymbol* fn);
static void fixup_array_formals(FnSymbol* fn);
static void clone_parameterized_primitive_methods(FnSymbol* fn);
static void fixup_query_formals(FnSymbol* fn);
static void checkConfigParams();

static void
checkUseBeforeDefs() {
  forv_Vec(FnSymbol, fn, gFns) {
    if (fn->defPoint->parentSymbol) {
      ModuleSymbol* mod = fn->getModule();
      Vec<const char*> undeclared;
      Vec<Symbol*> undefined;
      Vec<BaseAST*> asts;
      Vec<Symbol*> defined;
      collect_asts_postorder(&asts, fn);
      forv_Vec(BaseAST, ast, asts) {
        if (CallExpr* call = toCallExpr(ast)) {
          if (call->isPrimitive(PRIMITIVE_MOVE))
            defined.set_add(toSymExpr(call->get(1))->var);
        } else if (DefExpr* def = toDefExpr(ast)) {
          if (isArgSymbol(def->sym))
            defined.set_add(def->sym);
        } else if (SymExpr* sym = toSymExpr(ast)) {
          CallExpr* call = toCallExpr(sym->parentExpr);
          if (call && call->isPrimitive(PRIMITIVE_MOVE) && call->get(1) == sym)
            continue;
          if (toModuleSymbol(sym->var))
            USR_FATAL_CONT(sym, "illegal use of module '%s'", sym->var->name);
          if ((!call || call->baseExpr != sym) && sym->unresolved) {
            if (!undeclared.set_in(sym->unresolved)) {
              if (!toFnSymbol(fn->defPoint->parentSymbol)) {
                USR_FATAL_CONT(sym, "'%s' undeclared (first use this function)",
                               sym->unresolved);
                undeclared.set_add(sym->unresolved);
              }
            }
          }
          if (isVarSymbol(sym->var) || isArgSymbol(sym->var)) {
            if (sym->var->defPoint &&
                (sym->var->defPoint->parentSymbol == fn ||
                 (sym->var->defPoint->parentSymbol == mod && mod->initFn == fn))) {
              if (!defined.set_in(sym->var) && !undefined.set_in(sym->var)) {
                if (strcmp(sym->var->name, "this")) {
                  USR_FATAL_CONT(sym, "'%s' used before defined (first used here)", sym->var->name);
                  undefined.set_add(sym->var);
                }
              }
            }
          }
        }
      }
    }
  }
}


static void
flattenGlobalFunctions() {
  forv_Vec(ModuleSymbol, mod, allModules) {
    for_alist(expr, mod->initFn->body->body) {
      if (DefExpr* def = toDefExpr(expr))
        if ((toVarSymbol(def->sym) && !def->sym->isCompilerTemp) ||
            toTypeSymbol(def->sym) ||
            toFnSymbol(def->sym))
          if (!((!strncmp("_anon_record", def->sym->name, 12)) ||
                (!strncmp("_forallexpr", def->sym->name, 11)) ||
                (!strncmp("_let_fn", def->sym->name, 7)) ||
                (!strncmp("_if_fn", def->sym->name, 6)) ||
                (!strncmp("_reduce_scan", def->sym->name, 12)) ||
                (!strncmp("_forif_fn", def->sym->name, 9))))
            mod->block->insertAtTail(def->remove());
    }
  }
}


static void insertHeapAllocate(CallExpr* move, bool global = false) {
  currentLineno = move->lineno;
  currentFilename = move->filename;
  VarSymbol* tmp = new VarSymbol("_tmp");
  tmp->isCompilerTemp = true;
  move->insertBefore(new DefExpr(tmp));
  move->insertBefore(new CallExpr(PRIMITIVE_MOVE, tmp, move->get(2)->remove()));
  if (global) {
    if (move->get(1)->isConstant())
      move->insertAtTail(new CallExpr("_heapAllocConstGlobal", tmp));
    else 
      move->insertAtTail(new CallExpr("_heapAllocGlobal", tmp));
  } else
    move->insertAtTail(new CallExpr("_heapAlloc", tmp));
}


static void insertHeapAccess(SymExpr* se, bool global = false) {
  currentLineno = se->lineno;
  currentFilename = se->filename;
  CallExpr* call =
    new CallExpr((global && se->var->isConstant())
                 ? "_heapAccessConstGlobal"
                 : "_heapAccess");
  Expr* stmt = se->getStmtExpr();
  if (!stmt) {
    se->replace(call);
    call->insertAtTail(se);
  } else {
    VarSymbol* tmp = new VarSymbol("_tmp");
    tmp->isCompilerTemp = true;
    tmp->isExprTemp = true;
    stmt->insertBefore(new DefExpr(tmp));
    stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, tmp, call));
    se->replace(new SymExpr(tmp)); call->insertAtTail(se);
  }
}


//
// return true if expr is in function's body
//
static bool inBody(FnSymbol* fn, Expr* expr) {
  if (expr == fn->body)
    return true;
  if (expr->parentExpr)
    return inBody(fn, expr->parentExpr);
  if (expr->parentSymbol && expr->parentSymbol != fn)
    return inBody(fn, expr->parentSymbol->defPoint);
  return false;
}


//
// change local variables and arguments into _heap classes if they
// should be allocated on the heap; do this for local variables and
// arguments accessed in begin- or on-statements
//
static void heapAllocateLocals() {
  Vec<Symbol*> heapSet; // set of symbols that should be heap allocated

  compute_sym_uses();

  // for each nested begin and on function
  forv_Vec(FnSymbol, fn, gFns) {
    if (fn->hasPragma("begin") || fn->hasPragma("on")) {

      Vec<BaseAST*> asts;
      collect_asts(&asts, fn);

      // collect defs of all local variables
      Vec<DefExpr*> defSet;
      forv_Vec(BaseAST, ast, asts) {
        if (DefExpr* def = toDefExpr(ast)) {
          defSet.set_add(def);
        }
      }

      // collect symbols that should be heap allocated
      forv_Vec(BaseAST, ast, asts) {
        if (SymExpr* se = toSymExpr(ast)) {

          // collect arguments
          if (isArgSymbol(se->var))
            heapSet.set_add(se->var);

          // collect symbols defined in outer functions
          if (VarSymbol* var = toVarSymbol(se->var))
            if (var->defPoint &&
                !defSet.set_in(var->defPoint) &&
                isFnSymbol(var->defPoint->parentSymbol) &&
                !var->isParam &&
                !var->hasPragma("private"))
              heapSet.set_add(var);
        }
      }
    }
  }

  // for each symbol that should be heap allocated
  forv_Vec(Symbol, sym, heapSet) {

    // allocate local variables 'v' on heap
    //   replace first definition 'v = ...' with 'v = _heap(...)'
    //   replace all other accesses 'v' with 'v._val'
    if (isVarSymbol(sym)) {

      // disable for index of parameter for loops
      bool disable = false;
      forv_Vec(SymExpr, se, sym->uses) {
        if (CallExpr* parent = toCallExpr(se->parentExpr))
          if (parent->isPrimitive(PRIMITIVE_LOOP_PARAM))
            disable = true;
      }
      if (disable)
        continue;

      bool first = true;
      forv_Vec(SymExpr, se, sym->defs) {

        // ack!! this is troublesome: we're assuming that the first
        // definition in the defs vector is the initial definition

        if (first) {
          CallExpr* move = toCallExpr(se->parentExpr);
          INT_ASSERT(move);
          INT_ASSERT(move->isPrimitive(PRIMITIVE_MOVE));
          insertHeapAllocate(move);
          first = false;
        } else {
          insertHeapAccess(se);
        }
      }
      forv_Vec(SymExpr, se, sym->uses) {
        insertHeapAccess(se);
      }
    }

    // replace arguments 'a' with heap-allocated temporaries
    //   insert definition at beginning of function 't = _heap(a)'
    //   replace all accesses 'a' with 't._val'
    if (ArgSymbol* arg = toArgSymbol(sym)) {
      FnSymbol* fn = sym->getFunction();
      VarSymbol* tmp = new VarSymbol("_tmp");
      tmp->isCompilerTemp = true;
      fn->insertAtHead(new CallExpr(PRIMITIVE_MOVE, tmp, new CallExpr("_heapAlloc", sym)));
      fn->insertAtHead(new DefExpr(tmp));
      forv_Vec(SymExpr, se, sym->defs) {
        se->var = tmp;
        insertHeapAccess(se);
      }
      forv_Vec(SymExpr, se, sym->uses) {

        // disable insertion of heap access on arguments accessed
        // outside of the function's body
        if (!inBody(fn, se))
          continue;

        se->var = tmp;
        insertHeapAccess(se);
      }

      // write back to arguments of out or inout intent
      if (arg->intent == INTENT_OUT || arg->intent == INTENT_INOUT) {
        SymExpr* se = new SymExpr(tmp);
        fn->insertBeforeReturnAfterLabel(new CallExpr(PRIMITIVE_MOVE, arg, new CallExpr(PRIMITIVE_GET_REF, se)));
        insertHeapAccess(se);
      }
    }
  }
}


//
// change global variables into _heap classes
//
// for constants of scalar type, use replicated scalars
//
// insert PRIMITIVE_COMM_ALL call to make it so that _heap classes on
// each locale point to the same data or, for replicated scalars, to
// make it so that each replicated scalar contains the same value
//
static void heapAllocateGlobals() {
  if (fLocal)
    return;

  compute_sym_uses();

  //
  // collect all global variables less parameters and compiler temps
  //
  forv_Vec(BaseAST, ast, gAsts) {
    if (DefExpr* def = toDefExpr(ast)) {
      if (VarSymbol* var = toVarSymbol(def->sym)) {
        if (var->defPoint && toModuleSymbol(var->defPoint->parentSymbol)) {

          //
          // do not change parameters and private global variables
          //
          if (var->isParam || var->hasPragma("private"))
            continue;

          bool first = true;
          forv_Vec(SymExpr, se, var->defs) {

            // ack!! this is troublesome: we're assuming that the first
            // definition in the defs vector is the initial definition

            if (first) {
              CallExpr* move = toCallExpr(se->parentExpr);
              INT_ASSERT(move);
              INT_ASSERT(move->isPrimitive(PRIMITIVE_MOVE));
              insertHeapAllocate(move, true);
              first = false;
            } else {
              insertHeapAccess(se, true);
            }
          }
          forv_Vec(SymExpr, se, var->uses) {
            insertHeapAccess(se, true);
          }

        }
      }
    }
  }
}


void normalize(void) {
  normalize(theProgram);
  normalized = true;
  checkUseBeforeDefs();
  flattenGlobalFunctions();
  heapAllocateGlobals();
  heapAllocateLocals();
}

void normalize(BaseAST* base) {
  Vec<BaseAST*> asts;

  asts.clear();
  collect_asts( &asts, base);
  forv_Vec(BaseAST, ast, asts) {
    if (FnSymbol* fn = toFnSymbol(ast)) {
      currentLineno = fn->lineno;
      currentFilename = fn->filename;
      if (!fn->hasPragma("type constructor") &&
          !fn->hasPragma("default constructor") &&
          !fn->isWrapper)
        fixup_array_formals(fn);
      clone_parameterized_primitive_methods(fn);
      fixup_query_formals(fn);
      change_method_into_constructor(fn);
    }
  }

  asts.clear();
  collect_asts(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    if (FnSymbol* fn = toFnSymbol(ast)) {
      if (fn->noParens && !fn->isMethod)
        fn->visible = false;
      normalize_returns(fn);
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    currentLineno = ast->lineno;
    currentFilename = ast->filename;
    if (FnSymbol* a = toFnSymbol(ast))
      if (!a->defSetGet)
        apply_getters_setters(a);
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    currentLineno = ast->lineno;
    currentFilename = ast->filename;
    if (SymExpr* a = toSymExpr(ast)) {
      call_constructor_for_class(a);
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    currentLineno = ast->lineno;
    currentFilename = ast->filename;
    if (DefExpr* a = toDefExpr(ast)) {
      if (VarSymbol* var = toVarSymbol(a->sym))
        if (toFnSymbol(a->parentSymbol))
          fix_def_expr(var);
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    if (SymExpr* sym = toSymExpr(ast)) {
      if (sym == sym->getStmtExpr() && isFnSymbol(sym->parentSymbol)) {
        CallExpr* call = new CallExpr("_statementLevelSymbol");
        sym->insertBefore(call);
        call->insertAtTail(sym->remove());
      }
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    currentLineno = ast->lineno;
    currentFilename = ast->filename;
    if (CallExpr* a = toCallExpr(ast)) {
      insert_call_temps(a);
      fix_user_assign(a);
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    if (FnSymbol *fn = toFnSymbol(ast)) {
      tag_global(fn);
    }
  }

  asts.clear();
  collect_asts_postorder(&asts, base);
  forv_Vec(BaseAST, ast, asts) {
    currentLineno = ast->lineno;
    currentFilename = ast->filename;
    if (Expr* a = toExpr(ast)) {
      hack_resolve_types(a);
    }
  }
  checkConfigParams();
}


static void normalize_returns(FnSymbol* fn) {
  Vec<BaseAST*> asts;
  Vec<CallExpr*> rets;
  collect_asts(&asts, fn);
  forv_Vec(BaseAST, ast, asts) {
    if (CallExpr* returnStmt = toCallExpr(ast)) {
      if (returnStmt->isPrimitive(PRIMITIVE_RETURN) ||
          returnStmt->isPrimitive(PRIMITIVE_YIELD))
        if (returnStmt->parentSymbol == fn) // not in a nested function
          rets.add(returnStmt);
    }
  }
  if (rets.n == 0) {
    if (fn->fnTag == FN_ITERATOR)
      USR_FATAL(fn, "iterator does not yield or return a value");
    fn->insertAtTail(new CallExpr(PRIMITIVE_RETURN, gVoid));
    return;
  }
  if (rets.n == 1) {
    CallExpr* ret = rets.v[0];
    if (ret == fn->body->body.last() && toSymExpr(ret->get(1)))
      return;
  }
  SymExpr* retSym = toSymExpr(rets.v[0]->get(1));
  bool returns_void = retSym && retSym->var == gVoid;
  LabelSymbol* label = new LabelSymbol(astr("_end_", fn->name));
  fn->insertAtTail(new DefExpr(label));
  VarSymbol* retval = NULL;
  if (returns_void) {
    fn->insertAtTail(new CallExpr(PRIMITIVE_RETURN, gVoid));
  } else {
    retval = new VarSymbol(astr("_ret_", fn->name), fn->retType);
    retval->isCompilerTemp = true;
    if (fn->retTag == RET_PARAM)
      retval->isParam = true;
    if (fn->retTag == RET_TYPE)
      retval->isTypeVariable = true;
    if (fn->retExprType && fn->retTag != RET_VAR) {
      fn->insertAtHead(new CallExpr(PRIMITIVE_MOVE, retval, new CallExpr(PRIMITIVE_INIT, fn->retExprType->copy())));
      fn->addPragma("specified return type");
    }
    fn->insertAtHead(new DefExpr(retval));
    fn->insertAtTail(new CallExpr(PRIMITIVE_RETURN, retval));
  }
  bool label_is_used = false;
  forv_Vec(CallExpr, ret, rets) {
    if (retval) {
      Expr* ret_expr = ret->get(1);
      ret_expr->remove();
      if (fn->retTag == RET_VAR)
        ret->insertBefore(new CallExpr(PRIMITIVE_MOVE, retval, new CallExpr(PRIMITIVE_SET_REF, ret_expr)));
      else if (fn->retExprType)
        ret->insertBefore(new CallExpr(PRIMITIVE_MOVE, retval, new CallExpr("=", retval, ret_expr)));
      else
        ret->insertBefore(new CallExpr(PRIMITIVE_MOVE, retval, new CallExpr(PRIMITIVE_GET_REF, ret_expr)));
    }
    if (fn->fnTag == FN_ITERATOR) {
      if (!retval)
        INT_FATAL(ret, "unexpected case");
      if (ret->isPrimitive(PRIMITIVE_RETURN)) {
        ret->insertAfter(new GotoStmt(GOTO_NORMAL, label));
        label_is_used = true;
      }
      ret->replace(new CallExpr(PRIMITIVE_YIELD, retval));
    } else if (ret->next != label->defPoint) {
      ret->replace(new GotoStmt(GOTO_NORMAL, label));
      label_is_used = true;
    } else {
      ret->remove();
    }
  }
  if (!label_is_used)
    label->defPoint->remove();
}


static void call_constructor_for_class(SymExpr* se) {
  if (TypeSymbol* ts = toTypeSymbol(se->var)) {
    CallExpr* call = toCallExpr(se->parentExpr);
    if (call && call->baseExpr == se) {
      if (ClassType* ct = toClassType(ts->type)) {
        CallExpr* parent = toCallExpr(call->parentExpr);
        CallExpr* parentParent = NULL;
        if (parent)
          parentParent = toCallExpr(parent->parentExpr);
        if (parent && parent->isPrimitive(PRIMITIVE_NEW)) {
          if (!ct->defaultConstructor)
            INT_FATAL(call, "class type has no default constructor");
          se->replace(new SymExpr(ct->defaultConstructor->name));
          parent->replace(call->remove());
        } else if (parentParent && parentParent->isPrimitive(PRIMITIVE_NEW) &&
                   call->partialTag == true) {
          if (!ct->defaultConstructor)
            INT_FATAL(call, "class type has no default constructor");
          se->replace(new SymExpr(ct->defaultConstructor->name));
          parentParent->replace(parent->remove());
        } else {
          if (!ct->defaultTypeConstructor)
            INT_FATAL(call, "class type has no default type constructor");
          se->replace(new SymExpr(ct->defaultTypeConstructor->name));
        }
      }
    } //else if (isClassType(ts->type)) {
      //se->replace(new CallExpr(se));
    //}
  }
}


static void apply_getters_setters(FnSymbol* fn) {
  // Most generally:
  //   x.f(a) --> f(_mt, x)(a)
  // which is the same as
  //   call(call(. x "f") a) --> call(call(f _mt x) a)
  // Also:
  //   x.f --> f(_mt, x)
  // Note:
  //   call(call or )( indicates partial
  Vec<BaseAST*> asts;
  collect_asts_postorder(&asts, fn);
  forv_Vec(BaseAST, ast, asts) {
    if (CallExpr* call = toCallExpr(ast)) {
      currentLineno = call->lineno;
      currentFilename = call->filename;
      if (call->getFunction() != fn) // in a nested function, handle
                                     // later, because it may be a
                                     // getter or a setter
        continue;
      if (call->isNamed(".")) { // handle getter
        CallExpr* getter = NULL;
        if (SymExpr* symExpr = toSymExpr(call->get(2))) {
          if (VarSymbol* var = toVarSymbol(symExpr->var)) {
            if (var->immediate->const_kind == CONST_KIND_STRING) {
              getter = new CallExpr(var->immediate->v_string,
                                    gMethodToken, call->get(1)->remove());
            }
          } else if (TypeSymbol* type = toTypeSymbol(symExpr->var)) {
            getter = new CallExpr(type,
                                  gMethodToken, call->get(1)->remove());
          }
        }
        if (!getter)
          INT_FATAL(call, "invalid dot call expression");
        getter->methodTag = true;
        call->replace(getter);
        if (CallExpr* parent = toCallExpr(getter->parentExpr))
          if (parent->baseExpr == getter)
            getter->partialTag = true;
      }
    }
  }
}


static void insert_call_temps(CallExpr* call) {
  if (!call->parentExpr || !call->getStmtExpr())
    return;

  if (call == call->getStmtExpr())
    return;
  
  if (toDefExpr(call->parentExpr))
    return;

  if (call->partialTag)
    return;

  if (call->isPrimitive(PRIMITIVE_TUPLE_EXPAND) ||
      call->isPrimitive(PRIMITIVE_GET_MEMBER_VALUE))
    return;

  if (CallExpr* parentCall = toCallExpr(call->parentExpr))
    if (parentCall->isPrimitive(PRIMITIVE_MOVE))
      return;

  Expr* stmt = call->getStmtExpr();
  VarSymbol* tmp = new VarSymbol("_tmp", dtUnknown);
  tmp->isCompilerTemp = true;
  tmp->isExprTemp = true;
  tmp->canParam = true;
  tmp->canType = true;
  call->replace(new SymExpr(tmp));
  stmt->insertBefore(new DefExpr(tmp));
  stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, tmp, call));
}


static void fix_user_assign(CallExpr* call) {
  if (!call->parentExpr ||
      call->getStmtExpr() == call->parentExpr ||
      !call->isNamed("="))
    return;
  CallExpr* move = new CallExpr(PRIMITIVE_MOVE, call->get(1)->copy());
  call->replace(move);
  move->insertAtTail(call);
}

//
// fix_def_expr removes DefExpr::exprType and DefExpr::init from a
//   variable's def expression, normalizing the AST with primitive
//   moves, calls to _copy, _init, and _cast, and assignments.
//
static void
fix_def_expr(VarSymbol* var) {
  Expr* type = var->defPoint->exprType;
  Expr* init = var->defPoint->init;
  Expr* stmt = var->defPoint; // insertion point
  VarSymbol* constTemp = var; // temp for constants

  if (!type && !init)
    return; // already fixed

  //
  // handle var ... : ... => ...;
  //
  if (var->isUserAlias) {
    CallExpr* partial;
    VarSymbol* arrTemp = new VarSymbol("_arrTmp");
    arrTemp->isCompilerTemp = true;
    stmt->insertBefore(new DefExpr(arrTemp));
    stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, arrTemp, init->remove()));
    if (!type) {
      stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, var, arrTemp));
    } else {
      partial = new CallExpr("reindex", gMethodToken, arrTemp);
      partial->partialTag = true;
      partial->methodTag = true;
      stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, var, new CallExpr(partial, type->remove())));
    }
    return;
  }

  //
  // insert temporary for constants to assist constant checking
  //
  if (var->isConst) {
    constTemp = new VarSymbol("_constTmp");
    constTemp->isCompilerTemp = true;
    stmt->insertBefore(new DefExpr(constTemp));
    stmt->insertAfter(new CallExpr(PRIMITIVE_MOVE, var, constTemp));
  }

  //
  // insert code to initialize config variable from the command line
  //
  if (var->isConfig) {
    if (!var->isParam) {
      Expr* noop = new CallExpr(PRIMITIVE_NOOP);
      ModuleSymbol* module = var->getModule();
      CallExpr* strToValExpr =
        new CallExpr("_cast",
                     new CallExpr(PRIMITIVE_TYPEOF, constTemp),
                     new CallExpr(primitives_map.get("_config_get_value"),
                                  new_StringSymbol(var->name),
                                  new_StringSymbol(module->name)));
      stmt->insertAfter(
        new CondStmt(
          new CallExpr("!",
            new CallExpr(primitives_map.get("_config_has_value"),
                         new_StringSymbol(var->name),
                         new_StringSymbol(module->name))),
          noop,
          new CallExpr(PRIMITIVE_MOVE, constTemp, strToValExpr)));
      strToValExpr->filename = astr("<command line setting of '", var->name,
                                    "'>");
      strToValExpr->lineno = 0;
                     

      stmt = noop; // insert regular definition code in then block
    } else {
      if (const char* value = configParamMap.get(astr(var->name))) {
        usedConfigParams.add(astr(var->name));
        if (SymExpr* symExpr = toSymExpr(init)) {
          if (VarSymbol* varSymbol = toVarSymbol(symExpr->var)) {
            if (varSymbol->immediate) {
              Immediate* imm;
              if (varSymbol->immediate->const_kind == CONST_KIND_STRING) {
                imm = new Immediate(value);
              } else {
                imm = new Immediate(*varSymbol->immediate);
                convert_string_to_immediate(value, imm);
              }
              init->replace(new SymExpr(new_ImmediateSymbol(imm)));
              init = var->defPoint->init;
            }
          } else if (EnumSymbol* sym = toEnumSymbol(symExpr->var)) {
            if (EnumType* et = toEnumType(sym->type)) {
              for_enums(constant, et) {
                if (!strcmp(constant->sym->name, value)) {
                  init->replace(new SymExpr(constant->sym));
                  init = var->defPoint->init;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  if (type) {

    //
    // use cast for parameters to avoid multiple parameter assignments
    //
    if (init && var->isParam) {
      stmt->insertAfter(
        new CallExpr(PRIMITIVE_MOVE, var,
          new CallExpr("_cast", type->remove(), init->remove())));
      return;
    }

    //
    // initialize variable based on specified type and then assign it
    // the initialization expression if it exists
    //
    VarSymbol* typeTemp = new VarSymbol("_typeTmp");
    typeTemp->isCompilerTemp = true;
    stmt->insertBefore(new DefExpr(typeTemp));
    CallExpr* callType = toCallExpr(type);
    if (callType && callType->isNamed("_build_sparse_subdomain_type"))
      stmt->insertBefore(
        new CallExpr(PRIMITIVE_MOVE, typeTemp, type->remove()));
    else
      stmt->insertBefore(
        new CallExpr(PRIMITIVE_MOVE, typeTemp,
          new CallExpr(PRIMITIVE_INIT, type->remove())));
    if (init) {
      VarSymbol* initTemp = new VarSymbol("_tmp");
      initTemp->isCompilerTemp = true;
      initTemp->canParam = true;
      stmt->insertBefore(new DefExpr(initTemp));
      stmt->insertBefore(new CallExpr(PRIMITIVE_MOVE, initTemp, init->remove()));
      stmt->insertAfter(new CallExpr(PRIMITIVE_MOVE, constTemp, typeTemp));
      stmt->insertAfter(
        new CallExpr(PRIMITIVE_MOVE, typeTemp,
          new CallExpr("=", typeTemp, initTemp)));
    } else {
      if (constTemp->isTypeVariable)
        stmt->insertAfter(new CallExpr(PRIMITIVE_MOVE, constTemp, new CallExpr(PRIMITIVE_TYPEOF, typeTemp)));
      else
        stmt->insertAfter(new CallExpr(PRIMITIVE_MOVE, constTemp, typeTemp));
    }

  } else {

    //
    // initialize untyped variable with initialization expression
    //
    stmt->insertAfter(
      new CallExpr(PRIMITIVE_MOVE, constTemp,
        new CallExpr("_copy", init->remove())));

  }
}


static void checkConfigParams() {
  bool anyBadConfigParams = false;
  Vec<const char*> configParamSetNames;
  configParamMap.get_keys(configParamSetNames);
  forv_Vec(const char, name, configParamSetNames) {
    if (!usedConfigParams.in(name)) {
      USR_FATAL_CONT("Trying to set unrecognized config param '%s' via -s flag", name);
      anyBadConfigParams = true;
    }
  }
  if (anyBadConfigParams) {
    USR_STOP();
  }
}


static void hack_resolve_types(Expr* expr) {
  if (DefExpr* def = toDefExpr(expr)) {
    if (ArgSymbol* arg = toArgSymbol(def->sym)) {
      if (arg->type == dtUnknown || arg->type == dtAny) {
        if (!arg->isTypeVariable && !arg->typeExpr && arg->defaultExpr) {
          SymExpr* se = NULL;
          if (arg->defaultExpr->body.length() == 1)
            se = toSymExpr(arg->defaultExpr->body.tail);
          if (!se || se->var != gNil) {
            arg->typeExpr = arg->defaultExpr->copy();
            FnSymbol* fn = def->getFunction();
            insert_help(arg->typeExpr, NULL, arg, fn->argScope);
          }
        }
        if (arg->typeExpr && arg->typeExpr->body.length() == 1) {
          Type* type = arg->typeExpr->body.only()->typeInfo();
          if (type != dtUnknown && type != dtAny) {
            arg->type = type;
            arg->typeExpr->remove();
          }
        }
      }
    }
  }
}


static void tag_global(FnSymbol* fn) {
  if (!fn->global && !fn->isWrapper) {
    if (ClassType* ct = toClassType(fn->_this)) {
      if (ct->classTag == CLASS_CLASS &&
          !ct->symbol->hasPragma("ref") &&
          !ct->symbol->hasPragma("domain") &&
          !ct->symbol->hasPragma("array")) {
        fn->global = true;
      }
    }
    if (fn->global) {
      fn->parentScope->removeVisibleFunction(fn);
      theProgram->block->blkScope->addVisibleFunction(fn);
      if (toFnSymbol(fn->defPoint->parentSymbol)) {
        ModuleSymbol* mod = fn->getModule();
        Expr* def = fn->defPoint;
        CallExpr* noop = new CallExpr(PRIMITIVE_NOOP);
        def->insertBefore(noop);
        fn->visiblePoint = noop;
        def->remove();
        mod->block->insertAtTail(def);
      }
    }
  }
}


static void fixup_array_formals(FnSymbol* fn) {
  Vec<BaseAST*> asts;
  collect_top_asts(&asts, fn);
  Vec<BaseAST*> all_asts;
  collect_asts(&all_asts, fn);
  forv_Vec(BaseAST, ast, asts) {
    if (CallExpr* call = toCallExpr(ast)) {
      if (call->isNamed("_build_array_type")) {
        SymExpr* sym = toSymExpr(call->get(1));
        DefExpr* def = toDefExpr(call->get(1));
        ArgSymbol* arg = toArgSymbol(call->parentSymbol);
        if (call->numActuals() == 1)
          if (!arg || !arg->typeExpr || arg->typeExpr->body.tail != call)
            USR_FATAL(call, "array declaration has no element type");
        if (def || (sym && sym->var == gNil) || call->numActuals() == 1) {
          if (!arg || !arg->typeExpr || arg->typeExpr->body.tail != call)
            USR_FATAL(call, "array with empty or queried domain can "
                      "only be used as a formal argument type");
          arg->typeExpr->replace(new BlockStmt(new SymExpr(chpl_array), BLOCK_SCOPELESS));
          if (!fn->where) {
            fn->where = new BlockStmt(new SymExpr(gTrue));
            insert_help(fn->where, NULL, fn, fn->argScope);
          }
          Expr* expr = fn->where->body.tail;
          if (call->numActuals() == 2)
            expr->replace(new CallExpr("&", expr->copy(),
                            new CallExpr("==", call->get(2)->remove(),
                              new CallExpr(".", arg, new_StringSymbol("eltType")))));
          if (def) {
            forv_Vec(BaseAST, ast, all_asts) {
              if (SymExpr* sym = toSymExpr(ast)) {
                if (sym->var == def->sym)
                  sym->replace(new CallExpr(".", arg, new_StringSymbol("_dom")));
              }
            }
          } else if (!sym || sym->var != gNil) {
            VarSymbol* tmp = new VarSymbol(astr("_reindex_", arg->name));
            forv_Vec(BaseAST, ast, all_asts) {
              if (SymExpr* sym = toSymExpr(ast)) {
                if (sym->var == arg)
                  sym->var = tmp;
              }
            }
            fn->insertAtHead(
              new CondStmt(
                new CallExpr("!=", dtNil->symbol, arg),
                new CallExpr(PRIMITIVE_MOVE, tmp,
                  new CallExpr(new CallExpr(".", arg,
                                 new_StringSymbol("reindex")),
                               call->get(1)->copy())),
                new CallExpr(PRIMITIVE_MOVE, tmp, gNil)));
            fn->insertAtHead(new DefExpr(tmp));
          }
        } else {  //// DUPLICATED CODE ABOVE AND BELOW
          if (ArgSymbol* arg = toArgSymbol(call->parentSymbol)) {
            if (arg->typeExpr && arg->typeExpr->body.tail == call) {
              arg->typeExpr->replace(new BlockStmt(new SymExpr(chpl_array), BLOCK_SCOPELESS));
              VarSymbol* tmp = new VarSymbol(astr("_reindex_", arg->name));
              forv_Vec(BaseAST, ast, all_asts) {
                if (SymExpr* sym = toSymExpr(ast)) {
                  if (sym->var == arg)
                    sym->var = tmp;
                }
              }
              fn->insertAtHead(
                new CondStmt(
                  new CallExpr("!=", dtNil->symbol, arg),
                  new CallExpr(PRIMITIVE_MOVE, tmp,
                    new CallExpr(new CallExpr(".", arg,
                                   new_StringSymbol("reindex")),
                                 call->get(1)->copy())),
                  new CallExpr(PRIMITIVE_MOVE, tmp, gNil)));
              fn->insertAtHead(new DefExpr(tmp));
            }
          }
        }
      }
    }
  }
}


static void clone_parameterized_primitive_methods(FnSymbol* fn) {
  if (toArgSymbol(fn->_this)) {
    if (fn->_this->type == dtInt[INT_SIZE_32]) {
      for (int i=INT_SIZE_1; i<INT_SIZE_NUM; i++) {
        if (dtInt[i] && i != INT_SIZE_32) {
          FnSymbol* nfn = fn->copy();
          nfn->_this->type = dtInt[i];
          fn->defPoint->insertBefore(new DefExpr(nfn));
        }
      }
    }
    if (fn->_this->type == dtUInt[INT_SIZE_32]) {
      for (int i=INT_SIZE_1; i<INT_SIZE_NUM; i++) {
        if (dtUInt[i] && i != INT_SIZE_32) {
          FnSymbol* nfn = fn->copy();
          nfn->_this->type = dtUInt[i];
          fn->defPoint->insertBefore(new DefExpr(nfn));
        }
      }
    }
    if (fn->_this->type == dtReal[FLOAT_SIZE_64]) {
      for (int i=FLOAT_SIZE_16; i<FLOAT_SIZE_NUM; i++) {
        if (dtReal[i] && i != FLOAT_SIZE_64) {
          FnSymbol* nfn = fn->copy();
          nfn->_this->type = dtReal[i];
          fn->defPoint->insertBefore(new DefExpr(nfn));
        }
      }
    }
    if (fn->_this->type == dtImag[FLOAT_SIZE_64]) {
      for (int i=FLOAT_SIZE_16; i<FLOAT_SIZE_NUM; i++) {
        if (dtImag[i] && i != FLOAT_SIZE_64) {
          FnSymbol* nfn = fn->copy();
          nfn->_this->type = dtImag[i];
          fn->defPoint->insertBefore(new DefExpr(nfn));
        }
      }
    }
    if (fn->_this->type == dtComplex[COMPLEX_SIZE_128]) {
      for (int i=COMPLEX_SIZE_32; i<COMPLEX_SIZE_NUM; i++) {
        if (dtComplex[i] && i != COMPLEX_SIZE_128) {
          FnSymbol* nfn = fn->copy();
          nfn->_this->type = dtComplex[i];
          fn->defPoint->insertBefore(new DefExpr(nfn));
        }
      }
    }
  }
}


static void
clone_for_parameterized_primitive_formals(FnSymbol* fn,
                                          DefExpr* def,
                                          int width) {
  ASTMap map;
  FnSymbol* newfn = fn->copy(&map);
  DefExpr* newdef = toDefExpr(map.get(def));
  Symbol* newsym = newdef->sym;
  newdef->replace(new SymExpr(new_IntSymbol(width)));
  Vec<BaseAST*> asts;
  map.get_values(asts);
  forv_Vec(BaseAST, ast, asts) {
    if (SymExpr* se = toSymExpr(ast))
      if (se->var == newsym)
        se->var = new_IntSymbol(width);
  }
  fn->defPoint->insertAfter(new DefExpr(newfn));
  fixup_query_formals(newfn);
}

static void
replace_query_uses(ArgSymbol* formal, DefExpr* def, ArgSymbol* arg,
                   Vec<BaseAST*>& asts) {
  if (!arg->isTypeVariable && arg->intent != INTENT_PARAM)
    USR_FATAL(def, "query variable is not type or parameter: %s", arg->name);
  forv_Vec(BaseAST, ast, asts) {
    if (SymExpr* se = toSymExpr(ast)) {
      if (se->var == def->sym) {
        if (formal->variableExpr) {
          CallExpr* parent = toCallExpr(se->parentExpr);
          if (!parent || parent->numActuals() != 1)
            USR_FATAL(se, "illegal access to query type or parameter");
          se->replace(new SymExpr(formal));
          parent->replace(se);
          se->replace(new CallExpr(".", parent, new_StringSymbol(arg->name)));
        } else {
          se->replace(new CallExpr(".", formal, new_StringSymbol(arg->name)));
        }
      }
    }
  }
}

static void
add_to_where_clause(ArgSymbol* formal, Expr* expr, ArgSymbol* arg) {
  if (!arg->isTypeVariable && arg->intent != INTENT_PARAM)
    USR_FATAL(expr, "type actual is not type or parameter");
  FnSymbol* fn = formal->defPoint->getFunction();
  if (!fn->where) {
    fn->where = new BlockStmt(new SymExpr(gTrue));
    insert_help(fn->where, NULL, fn, fn->argScope);
  }
  Expr* where = fn->where->body.tail;
  Expr* clause;
  if (formal->variableExpr)
    clause = new CallExpr(PRIMITIVE_TUPLE_AND_EXPAND, formal,
                          new_StringSymbol(arg->name), expr->copy());
  else
    clause = new CallExpr("==", expr->copy(),
               new CallExpr(".", formal, new_StringSymbol(arg->name)));
  where->replace(new CallExpr("&", where->copy(), clause));
}

static void
fixup_query_formals(FnSymbol* fn) {
  for_formals(formal, fn) {
    if (!formal->typeExpr)
      continue;
    if (DefExpr* def = toDefExpr(formal->typeExpr->body.tail)) {
      Vec<BaseAST*> asts;
      collect_asts(&asts, fn);
      forv_Vec(BaseAST, ast, asts) {
        if (SymExpr* se = toSymExpr(ast)) {
          if (se->var == def->sym) {
            se->replace(new CallExpr(PRIMITIVE_TYPEOF, formal));
          }
        }
      }
      formal->typeExpr->remove();
      formal->type = dtAny;
    } else if (CallExpr* call = toCallExpr(formal->typeExpr->body.tail)) {
      // clone query primitive types
      if (call->numActuals() == 1) {
        if (DefExpr* def = toDefExpr(call->get(1))) {
          if (call->isNamed("int") || call->isNamed("uint")) {
            for( int i=INT_SIZE_1; i<INT_SIZE_NUM; i++)
              if (dtInt[i])
                clone_for_parameterized_primitive_formals(fn, def,
                                                          get_width(dtInt[i]));
            fn->defPoint->remove();
            return;
          } else if (call->isNamed("real") || call->isNamed("imag")) {
            for( int i=FLOAT_SIZE_16; i<FLOAT_SIZE_NUM; i++)
              if (dtReal[i])
                clone_for_parameterized_primitive_formals(fn, def,
                                                          get_width(dtReal[i]));
            fn->defPoint->remove();
            return;
          } else if (call->isNamed("complex")) {
            for( int i=COMPLEX_SIZE_32; i<COMPLEX_SIZE_NUM; i++)
              if (dtComplex[i])
                clone_for_parameterized_primitive_formals(fn, def,
                                                          get_width(dtComplex[i]));
            fn->defPoint->remove();
            return;
          }
        }
      }
      bool queried = false;
      for_actuals(actual, call) {
        if (toDefExpr(actual))
          queried = true;
        if (NamedExpr* named = toNamedExpr(actual))
          if (toDefExpr(named->actual))
            queried = true;
      }
      if (queried) {
        Vec<BaseAST*> asts;
        collect_asts(&asts, fn);
        SymExpr* base = toSymExpr(call->baseExpr);
        if (!base)
          USR_FATAL(base, "illegal queried type expression");
        TypeSymbol* ts = toTypeSymbol(base->var);
        if (!ts)
          USR_FATAL(base, "illegal queried type expression");
        Vec<ArgSymbol*> args;
        for_formals(arg, ts->type->defaultConstructor) {
          args.add(arg);
        }
        for_actuals(actual, call) {
          if (NamedExpr* named = toNamedExpr(actual)) {
            for (int i = 0; i < args.n; i++) {
              if (args.v[i]) {
                if (!strcmp(named->name, args.v[i]->name)) {
                  if (DefExpr* def = toDefExpr(named->actual)) {
                    replace_query_uses(formal, def, args.v[i], asts);
                  } else {
                    add_to_where_clause(formal, named->actual, args.v[i]);
                  }
                  args.v[i] = NULL;
                  break;
                }
              }
            }
          }
        }
        for_actuals(actual, call) {
          if (!toNamedExpr(actual)) {
            for (int i = 0; i < args.n; i++) {
              if (args.v[i]) {
                if (DefExpr* def = toDefExpr(actual)) {
                  replace_query_uses(formal, def, args.v[i], asts);
                } else {
                  add_to_where_clause(formal, actual, args.v[i]);
                }
                args.v[i] = NULL;
                break;
              }
            }
          }
        }
        formal->typeExpr->remove();
        formal->type = ts->type;
      }
    }
  }
}


static void change_method_into_constructor(FnSymbol* fn) {
  if (fn->numFormals() <= 1)
    return;
  if (fn->getFormal(1)->type != dtMethodToken)
    return;
  if (fn->getFormal(2)->type == dtUnknown)
    INT_FATAL(fn, "this argument has unknown type");
  if (strcmp(fn->getFormal(2)->type->symbol->name, fn->name))
    return;
  ClassType* ct = toClassType(fn->getFormal(2)->type);
  if (!ct)
    INT_FATAL(fn, "constructor on non-class type");
  fn->_this = new VarSymbol("this", ct);
  fn->insertAtHead(new CallExpr(PRIMITIVE_MOVE, fn->_this, new CallExpr(ct->defaultConstructor)));
  fn->insertAtHead(new DefExpr(fn->_this));
  fn->insertAtTail(new CallExpr(PRIMITIVE_RETURN, new SymExpr(fn->_this)));
  ASTMap map;
  map.put(fn->getFormal(2), fn->_this);
  fn->formals.get(2)->remove();
  fn->formals.get(1)->remove();
  update_symbols(fn, &map);
  fn->defPoint->remove();
  fn->name = astr(astr("_construct_", fn->name));
  ct->symbol->defPoint->insertBefore(fn->defPoint);
  fn->retType = ct;
  if (ct->defaultConstructor->visible) {
    Expr* before = ct->defaultConstructor->defPoint->prev;
    ct->defaultConstructor->defPoint->remove();
    ct->defaultConstructor->visible = false;
    before->insertAfter(ct->defaultConstructor->defPoint);
  }
}
