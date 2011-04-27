config param noRefCount = false;

var privatizeLock$: sync int;

pragma "privatized class"
proc _isPrivatized(value) param
  return !_local & ((_privatization & value.dsiSupportsPrivatization()) | value.dsiRequiresPrivatization());

proc _newPrivatizedClass(value) {
  privatizeLock$.writeEF(true);
  var n = __primitive("chpl_numPrivatizedClasses");
  var hereID = here.uid;
  const privatizeData = value.dsiGetPrivatizeData();
  on Realms(0) do
    _newPrivatizedClassHelp(value, value, n, hereID, privatizeData);

  proc _newPrivatizedClassHelp(parentValue, originalValue, n, hereID, privatizeData) {
    var newValue = originalValue;
    if hereID != here.uid {
      newValue = parentValue.dsiPrivatize(privatizeData);
      __primitive("chpl_newPrivatizedClass", newValue);
      newValue.pid = n;
    } else {
      __primitive("chpl_newPrivatizedClass", newValue);
      newValue.pid = n;
    }
    cobegin {
      if chpl_localeTree.left then
        on chpl_localeTree.left do
          _newPrivatizedClassHelp(newValue, originalValue, n, hereID, privatizeData);
      if chpl_localeTree.right then
        on chpl_localeTree.right do
          _newPrivatizedClassHelp(newValue, originalValue, n, hereID, privatizeData);
    }
  }

  privatizeLock$.readFE();
  return n;
}

proc _reprivatize(value) {
  var pid = value.pid;
  var hereID = here.uid;
  const reprivatizeData = value.dsiGetReprivatizeData();
  on Realms(0) do
    _reprivatizeHelp(value, value, pid, hereID, reprivatizeData);
  proc _reprivatizeHelp(parentValue, originalValue, pid, hereID, reprivatizeData) {
    var newValue = originalValue;
    if hereID != here.uid {
      newValue = chpl_getPrivatizedCopy(newValue.type, pid);
      newValue.dsiReprivatize(parentValue, reprivatizeData);
    }
    cobegin {
      if chpl_localeTree.left then
        on chpl_localeTree.left do
          _reprivatizeHelp(newValue, originalValue, pid, hereID, reprivatizeData);
      if chpl_localeTree.right then
        on chpl_localeTree.right do
          _reprivatizeHelp(newValue, originalValue, pid, hereID, reprivatizeData);
    }
  }
}

//
// Take a rank and value and check that the value is a rank-tuple or not a
// tuple. If the value is not a tuple and expand is true, copy the value into
// a rank-tuple. If the value is a scalar and rank is 1, copy it into a 1-tuple.
//
proc _makeIndexTuple(param rank, t: _tuple, param expand: bool=false) where rank == t.size {
  return t;
}

proc _makeIndexTuple(param rank, t: _tuple, param expand: bool=false) where rank != t.size {
  compilerError("index rank must match domain rank");
}

proc _makeIndexTuple(param rank, val:integral, param expand: bool=false) {
  if expand || rank == 1 {
    var t: rank*val.type;
    for param i in 1..rank do
      t(i) = val;
    return t;
  } else {
    compilerWarning(typeToString(val.type));
    compilerError("index rank must match domain rank");
    return val;
  }
}

proc _newArray(value) {
  if _isPrivatized(value) then
    return new _array(_newPrivatizedClass(value), value);
  else
    return new _array(value, value);
}

proc _newDomain(value) {
  if _isPrivatized(value) then
    return new _domain(_newPrivatizedClass(value), value);
  else
    return new _domain(value, value);
}

proc _getDomain(value) {
  if _isPrivatized(value) then
    return new _domain(value.pid, value);
  else
    return new _domain(value, value);
}

proc _newDistribution(value) {
  if _isPrivatized(value) then
    return new _distribution(_newPrivatizedClass(value), value);
  else
    return new _distribution(value, value);
}

proc _getDistribution(value) {
  if _isPrivatized(value) then
    return new _distribution(value.pid, value);
  else
    return new _distribution(value, value);
}


//
// Support for domain types
//
pragma "has runtime type"
proc chpl__buildDomainRuntimeType(d: _distribution, param rank: int,
                                 type idxType = int(32),
                                 param stridable: bool = false) type
  return _newDomain(d.newRectangularDom(rank, idxType, stridable));

pragma "has runtime type"
proc chpl__buildDomainRuntimeType(d: _distribution, type idxType,
                                  param parSafe: bool = true) type
  return _newDomain(d.newAssociativeDom(idxType, parSafe));

pragma "has runtime type"
proc chpl__buildDomainRuntimeType(d: _distribution, type idxType,
                                  param parSafe: bool = true) type
 where idxType == _OpaqueIndex
  return _newDomain(d.newOpaqueDom(idxType, parSafe));

// This function has no 'has runtime type' pragma since the idxType of
// opaque domains is _OpaqueIndex, not opaque.  This function is
// essentially a wrapper around the function that actually builds up
// the runtime type.
proc chpl__buildDomainRuntimeType(d: _distribution, type idxType) type
 where idxType == opaque
  return chpl__buildDomainRuntimeType(d, _OpaqueIndex);

pragma "has runtime type"
proc chpl__buildSparseDomainRuntimeType(d: _distribution, dom: domain) type
  return _newDomain(d.newSparseDom(dom.rank, dom._value.idxType, dom));

proc chpl__convertValueToRuntimeType(dom: domain) type
 where dom._value:BaseRectangularDom
  return chpl__buildDomainRuntimeType(dom.dist, dom._value.rank,
                            dom._value.idxType, dom._value.stridable);

proc chpl__convertValueToRuntimeType(dom: domain) type
 where dom._value:BaseAssociativeDom
  return chpl__buildDomainRuntimeType(dom.dist, dom._value.idxType, dom._value.parSafe);

proc chpl__convertValueToRuntimeType(dom: domain) type
 where dom._value:BaseOpaqueDom
  return chpl__buildDomainRuntimeType(dom.dist, dom._value.idxType);

proc chpl__convertValueToRuntimeType(dom: domain) type
 where dom._value:BaseSparseDom
  return chpl__buildSparseDomainRuntimeType(dom.dist, dom._value.parentDom);

proc chpl__convertValueToRuntimeType(dom: domain) type {
  compilerError("the global domain class of each domain map implementation must be a subclass of BaseRectangularDom, BaseAssociativeDom, BaseOpaqueDom, or BaseSparseDom", 0);
  return 0; // dummy
}

//
// Support for array types
//
pragma "has runtime type"
proc chpl__buildArrayRuntimeType(dom: domain, type eltType) type
  return dom.buildArray(eltType);

proc chpl__convertValueToRuntimeType(arr: []) type
  return chpl__buildArrayRuntimeType(arr.domain, arr.eltType);

proc chpl__getDomainFromArrayType(type arrayType) {
  var A: arrayType;
  pragma "no copy" var D = A.domain;
  pragma "dont disable remote value forwarding"
  proc help() {
    D._value._domCnt$ += 1;
  }
  if !noRefCount then
    help();
  return D;
}

//
// Support for subdomain types
//
// Note the domain of a subdomain is not yet part of the runtime type
//
proc chpl__buildSubDomainType(dom: domain) type
  return chpl__convertValueToRuntimeType(dom);

//
// Support for domain expressions, e.g., [1..3, 1..3]
//
proc chpl__buildDomainExpr(x: domain)
  return x;

proc chpl__buildDomainExpr(ranges: range(?) ...?rank) {
  for param i in 2..rank {
    if ranges(1).idxType != ranges(i).idxType then
      compilerError("domain has mixed dimensional type");
    if ranges(i).boundedType != BoundedRangeType.bounded then
      compilerError("domain has dimension of unbounded range");
  }
  var d: domain(rank, ranges(1).idxType, chpl__anyStridable(ranges));
  d.setIndices(ranges);
  return d;
}

//
// Support for distributed domain expression, e.g., [1..3, 1..3] distributed Dist
//
proc chpl__distributed(d: _distribution, dom: domain) {
  if isRectangularDom(dom) {
    var distDom: domain(dom.rank, dom._value.idxType, dom._value.stridable) dmapped d = dom;
    return distDom;
  } else {
    var distDom: domain(dom._value.idxType) dmapped d = dom;
    return distDom;
  }
}

proc chpl__distributed(d: _distribution, ranges: range(?) ...?rank) {
  return chpl__distributed(d, chpl__buildDomainExpr((...ranges)));
}

proc chpl__isRectangularDomType(type domainType) param {
  var dom: domainType;
  return isRectangularDom(dom);
}

proc chpl__isSparseDomType(type domainType) param {
  var dom: domainType;
  return isSparseDom(dom);
}

proc chpl__distributed(d: _distribution, type domainType) type {
  if chpl__isRectangularDomType(domainType) {
    var dom: domainType;
    return chpl__buildDomainRuntimeType(d, dom._value.rank, dom._value.idxType,
                                        dom._value.stridable);
  } else if chpl__isSparseDomType(domainType) {
    //
    // this "no auto destroy" pragma is necessary as of 1/20 because
    // otherwise the parentDom gets destroyed in the sparse case; see
    // sparse/bradc/CSR/sparse.chpl as an example
    //
    pragma "no auto destroy" var dom: domainType;
    return chpl__buildSparseDomainRuntimeType(d, dom._value.parentDom);
  } else {
    var dom: domainType;
    return chpl__buildDomainRuntimeType(d, dom._value.idxType, dom._value.parSafe);
  }
}

//
// Support for index types
//
proc chpl__buildIndexType(param rank: int, type idxType) type where rank == 1 {
  var x: idxType;
  return x.type;
}

proc chpl__buildIndexType(param rank: int, type idxType) type where rank > 1 {
  var x: rank*idxType;
  return x.type;
}

proc chpl__buildIndexType(param rank: int) type
  return chpl__buildIndexType(rank, int);

proc chpl__buildIndexType(d: domain) type
  return chpl__buildIndexType(d.rank, d._value.idxType);

proc chpl__buildIndexType(type idxType) type where idxType == opaque
  return _OpaqueIndex;

proc isRectangularDom(d: domain) param {
  proc isRectangularDomClass(dc: BaseRectangularDom) param return true;
  proc isRectangularDomClass(dc) param return false;
  return isRectangularDomClass(d._value);
}

proc isRectangularArr(a: []) param return isRectangularDom(a.domain);

proc isIrregularDom(d: domain) param {
  return isSparseDom(d) || isAssociativeDom(d) || isOpaqueDom(d);
}

proc isIrregularArr(a: []) param return isIrregularDom(a.domain);

proc isAssociativeDom(d: domain) param {
  proc isAssociativeDomClass(dc: BaseAssociativeDom) param return true;
  proc isAssociativeDomClass(dc) param return false;
  return isAssociativeDomClass(d._value);
}

proc isAssociativeArr(a: []) param return isAssociativeDom(a.domain);

proc isEnumDom(d: domain) param {
  return isAssociativeDom(d) && _isEnumeratedType(d._value.idxType);
}

proc isEnumArr(a: []) param return isEnumDom(a.domain);

proc isOpaqueDom(d: domain) param {
  proc isOpaqueDomClass(dc: BaseOpaqueDom) param return true;
  proc isOpaqueDomClass(dc) param return false;
  return isOpaqueDomClass(d._value);
}

proc isSparseDom(d: domain) param {
  proc isSparseDomClass(dc: BaseSparseDom) param return true;
  proc isSparseDomClass(dc) param return false;
  return isSparseDomClass(d._value);
}

proc isSparseArr(a: []) param return isSparseDom(a.domain);

//
// Support for distributions
//
pragma "syntactic distribution"
record dmap { }

proc chpl__buildDistType(type t) type where t: BaseDist {
  var x: t;
  var y = _newDistribution(x);
  return y.type;
}

proc chpl__buildDistType(type t) {
  compilerError("illegal domain map type specifier - must be a subclass of BaseDist");
}

proc chpl__buildDistValue(x) where x: BaseDist {
  return _newDistribution(x);
}

proc chpl__buildDistValue(x) {
  compilerError("illegal domain map value specifier - must be a subclass of BaseDist");
}

//
// Distribution wrapper record
//
pragma "distribution"
record _distribution {
  var _value;
  var _valueType;

  proc _distribution(_value, _valueType) { }

  pragma "inline"
  proc _value {
    if _isPrivatized(_valueType) {
      return chpl_getPrivatizedCopy(_valueType.type, _value);
    } else {
      return _value;
    }
  }

  proc ~_distribution() {
    if !_isPrivatized(_valueType) {
      on _value {
        var cnt = _value.destroyDist();
        if cnt == 0 {
          if !noRefCount {
            _value.dsiDestroyDistClass();
            delete _value;
          }
        }
      }
    }
  }

  proc clone() {
    return _newDistribution(_value.dsiClone());
  }

  proc newRectangularDom(param rank: int, type idxType, param stridable: bool) {
    var x = _value.dsiNewRectangularDom(rank, idxType, stridable);
    if x.linksDistribution() {
      var cnt = _value._distCnt$;
      _value._doms.append(x);
      if !noRefCount then
        _value._distCnt$ = cnt + 1;
      else
        _value._distCnt$ = cnt;
    }
    return x;
  }

  proc newAssociativeDom(type idxType, param parSafe: bool=true) {
    var x = _value.dsiNewAssociativeDom(idxType, parSafe);
    if x.linksDistribution() {
      var cnt = _value._distCnt$;
      _value._doms.append(x);
      if !noRefCount then
        _value._distCnt$ = cnt + 1;
      else
        _value._distCnt$ = cnt;
    }
    return x;
  }

  proc newAssociativeDom(type idxType, param parSafe: bool=true)
  where _isEnumeratedType(idxType) {
    var x = _value.dsiNewAssociativeDom(idxType, parSafe);
    if x.linksDistribution() {
      var cnt = _value._distCnt$;
      _value._doms.append(x);
      if !noRefCount then
        _value._distCnt$ = cnt + 1;
      else
        _value._distCnt$ = cnt;
    }
    const enumTuple = _enum_enumerate(idxType);
    for param i in 1..enumTuple.size do
      x.dsiAdd(enumTuple(i));
    return x;
  }

  proc newOpaqueDom(type idxType, param parSafe: bool=true) {
    var x = _value.dsiNewOpaqueDom(idxType, parSafe);
    if x.linksDistribution() {
      var cnt = _value._distCnt$;
      _value._doms.append(x);
      if !noRefCount then
        _value._distCnt$ = cnt + 1;
      else
        _value._distCnt$ = cnt;
    }
    return x;
  }

  proc newSparseDom(param rank: int, type idxType, dom: domain) {
    var x = _value.dsiNewSparseDom(rank, idxType, dom);
    if x.linksDistribution() {
      var cnt = _value._distCnt$;
      _value._doms.append(x);
      if !noRefCount then
        _value._distCnt$ = cnt + 1;
      else
        _value._distCnt$ = cnt;
    }
    return x;
  }

  proc idxToLocale(ind) return _value.dsiIndexLocale(ind);

  proc writeThis(x: Writer) {
    _value.writeThis(x);
  }

  proc displayRepresentation() { _value.dsiDisplayRepresentation(); }
}  // record _distribution


//
// Domain wrapper record
//
pragma "domain"
pragma "has runtime type"
record _domain {
  var _value;     // stores domain class, may be privatized
  var _valueType; // stores type of privatized domains
  var _promotionType: index(rank, _value.idxType);

  pragma "inline"
  proc _value {
    if _isPrivatized(_valueType) {
      return chpl_getPrivatizedCopy(_valueType.type, _value);
    } else {
      return _value;
    }
  }

  proc ~_domain () {
    if !_isPrivatized(_valueType) {
      on _value {
        var cnt = _value.destroyDom();
        if !noRefCount then
          if cnt == 0 then
            delete _value;
      }
    }
  }

  proc dist return _getDistribution(_value.dist);

  proc rank param {
    if isRectangularDom(this) || isSparseDom(this) then
      return _value.rank;
    else
      return 1;
  }

  proc stridable param where isRectangularDom(this) {
    return _value.stridable;
  }

  proc stridable param where isSparseDom(this) {
    compilerError("sparse domains do not currently support .stridable");
  }

  proc stridable param where isOpaqueDom(this) {
    compilerError("opaque domains do not support .stridable");  
  }

  proc stridable param where isEnumDom(this) {
    compilerError("enumerated domains do not support .stridable");  
  }

  proc stridable param where isAssociativeDom(this) {
    compilerError("associative domains do not support .stridable");  
  }

  pragma "inline"
  proc these() {
    return _value.these();
  }

  // see comments for the same method in _array
  //
  proc this(d: domain) where d.rank == rank
    return this((...d.getIndices()));

  proc this(ranges: range(?) ...rank) {
    param stridable = _value.stridable || chpl__anyStridable(ranges);
    var r: rank*range(_value.idxType,
                      BoundedRangeType.bounded,
                      stridable);

    for param i in 1..rank {
      r(i) = _value.dsiDim(i)(ranges(i));
    }
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType, stridable, r);
    if !noRefCount then
      if d.linksDistribution() then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
  }

  proc this(args ...rank) where _validRankChangeArgs(args, _value.idxType) {
    var ranges = _getRankChangeRanges(args);
    param newRank = ranges.size, stridable = chpl__anyStridable(ranges);
    var newRanges: newRank*range(idxType=_value.idxType, stridable=stridable);
    var newDistVal = _value.dist.dsiCreateRankChangeDist(newRank, args);
    var newDist = _getNewDist(newDistVal);
    var j = 1;
    var makeEmpty = false;

    for param i in 1..rank {
      if !isCollapsedDimension(args(i)) {
        newRanges(j) = dim(i)(args(i));
        j += 1;
      } else {
        if !dim(i).member(args(i)) then
          makeEmpty = true;
      }
    }
    if makeEmpty {
      for param i in 1..newRank {
        newRanges(i) = 1..0;
      }
    }
    var d = [(...newRanges)] dmapped newDist;
    return d;
  }

  proc dims() return _value.dsiDims();

  proc dim(d : int) return _value.dsiDim(d);

  proc dim(param d : int) return _value.dsiDim(d);

  iter dimIter(param d, ind) {
    for i in _value.dimIter(d, ind) do yield i;
  }

  proc buildArray(type eltType) {
    var x = _value.dsiBuildArray(eltType);
    pragma "dont disable remote value forwarding"
    proc help() {
      var cnt = _value._domCnt$; // lock
      _value._arrs.append(x);
      if !noRefCount then
        _value._domCnt$ = cnt + 1; // unlock
      else
        _value._domCnt$ = cnt;    // unlock
    }
    help();
    return _newArray(x);
  }

  proc clear() {
    _value.dsiClear();
  }

  proc create() {
    if _value.idxType != _OpaqueIndex then
      compilerError("domain.create() only applies to opaque domains");
    return _value.dsiCreate();
  }

  proc add(i) {
    _value.dsiAdd(i);
  }

  proc remove(i) {
    _value.dsiRemove(i);
  }

  proc numIndices return _value.dsiNumIndices;
  proc low return _value.dsiLow;
  proc high return _value.dsiHigh;
  proc stride return _value.dsiStride;
  proc alignment return _value.dsiAlignment;
  proc first return _value.dsiFirst;
  proc last return _value.dsiLast;
  proc alignedLow return _value.dsiAlignedLow;
  proc alignedHigh return _value.dsiAlignedHigh;

  proc member(i) {
    if isRectangularDom(this) then
      return _value.dsiMember(_makeIndexTuple(rank, i));
    else
      return _value.dsiMember(i);
  }

  // 1/5/10: do we want to support order() and position()?
  proc indexOrder(i) return _value.dsiIndexOrder(_makeIndexTuple(rank, i));

  proc position(i) {
    var ind = _makeIndexTuple(rank, i), pos: rank*_value.idxType;
    for d in 1..rank do
      pos(d) = _value.dsiDim(d).indexOrder(ind(d));
    return pos;
  }

  proc expand(off: rank*_value.idxType) where !isRectangularDom(this) {
    if isAssociativeDom(this) then
      compilerError("expand not supported on associative domains");
    else if isOpaqueDom(this) then
      compilerError("expand not supported on opaque domains");
    else if isSparseDom(this) then
      compilerError("expand not supported on sparse domains");
    else
      compilerError("expand not supported on this domain type");
  }
  proc expand(off: _value.idxType ...rank) return expand(off);
  proc expand(off: rank*_value.idxType) {
    var ranges = dims();
    for i in 1..rank do {
      ranges(i) = ranges(i).expand(off(i));
      if (ranges(i).low > ranges(i).high) {
        halt("***Error: Degenerate dimension created in dimension ", i, "***");
      }
    }

    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
  }
  proc expand(off: _value.idxType) where rank > 1 {
    var ranges = dims();
    for i in 1..rank do
      ranges(i) = dim(i).expand(off);
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
      d.dist._distCnt$ += 1;
    return _newDomain(d);
  }

  proc exterior(off: rank*_value.idxType) where !isRectangularDom(this) {
    if isAssociativeDom(this) then
      compilerError("exterior not supported on associative domains");
    else if isOpaqueDom(this) then
      compilerError("exterior not supported on opaque domains");
    else if isSparseDom(this) then
      compilerError("exterior not supported on sparse domains");
    else
      compilerError("exterior not supported on this domain type");
  }
  proc exterior(off: _value.idxType ...rank) return exterior(off);
  proc exterior(off: rank*_value.idxType) {
    var ranges = dims();
    for i in 1..rank do
      ranges(i) = dim(i).exterior(off(i));
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
   }
                  
  proc interior(off: rank*_value.idxType) where !isRectangularDom(this) {
    if isAssociativeDom(this) then
      compilerError("interior not supported on associative domains");
    else if isOpaqueDom(this) then
      compilerError("interior not supported on opaque domains");
    else if isSparseDom(this) then
      compilerError("interior not supported on sparse domains");
    else
      compilerError("interior not supported on this domain type");
  }
  proc interior(off: _value.idxType ...rank) return interior(off);
  proc interior(off: rank*_value.idxType) {
    var ranges = dims();
    for i in 1..rank do {
      if ((off(i) > 0) && (dim(i).high+1-off(i) < dim(i).low) ||
          (off(i) < 0) && (dim(i).low-1-off(i) > dim(i).high)) {
        halt("***Error: Argument to 'interior' function out of range in dimension ", i, "***");
      } 
      ranges(i) = _value.dsiDim(i).interior(off(i));
    }
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
  }

  //
  // NOTE: We eventually want to support translate on other domain types
  //
  proc translate(off) where !isRectangularDom(this) {
    if isAssociativeDom(this) then
      compilerError("translate not supported on associative domains");
    else if isOpaqueDom(this) then
      compilerError("translate not supported on opaque domains");
    else if isSparseDom(this) then
      compilerError("translate not supported on sparse domains");
    else
      compilerError("translate not supported on this domain type");
  }
  //
  // Notice that the type of the offset does not have to match the
  // index type.  This is handled in the range.translate().
  //
  proc translate(off: ?t ...rank) return translate(off);
  proc translate(off) where isTuple(off) {
    if off.size != rank then
      compilerError("the domain and offset arguments of translate() must be of the same rank");
    var ranges = dims();
    for i in 1..rank do
      ranges(i) = _value.dsiDim(i).translate(off(i));
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
   }

  //
  // intended for internal use only:
  //
  proc chpl__unTranslate(off: _value.idxType ...rank) return chpl__unTranslate(off);
  proc chpl__unTranslate(off: rank*_value.idxType) {
    var ranges = dims();
    for i in 1..rank do
      ranges(i) = dim(i).chpl__unTranslate(off(i));
    var d = _value.dsiBuildRectangularDom(rank, _value.idxType,
                                         _value.stridable, ranges);
    if !noRefCount then
      if (d.linksDistribution()) then
        d.dist._distCnt$ += 1;
    return _newDomain(d);
  }

  proc setIndices(x) {
    _value.dsiSetIndices(x);
    if _isPrivatized(_valueType) {
      _reprivatize(_value);
    }
  }

  proc getIndices()
    return _value.dsiGetIndices();

  proc writeThis(f: Writer) {
    _value.dsiSerialWrite(f);
  }

  proc localSlice(r: range(?)... rank) {
    return _value.dsiLocalSlice(chpl__anyStridable(r), r);
  }

  // associative array interface

  iter sorted() {
    for i in _value.dsiSorted() {
      yield i;
    }
  }

  proc displayRepresentation() { _value.dsiDisplayRepresentation(); }
}  // record _domain

proc _getNewDist(value) {
  return new dmap(value);
}

proc +(d: domain, i: index(d)) {
  if isRectangularDom(d) then
    compilerError("Cannot add indices to a rectangular domain");
  else
    compilerError("Cannot add indices to this domain type");
}

proc +(i, d: domain) where i: index(d) {
  if isRectangularDom(d) then
    compilerError("Cannot add indices to a rectangular domain");
  else
    compilerError("Cannot add indices to this domain type");
}

proc +(d: domain, i: index(d)) where isIrregularDom(d) {
  d.add(i);
  return d;
}

proc +(i, d: domain) where i:index(d) && isIrregularDom(d) {
  d.add(i);
  return d;
}

proc +(d1: domain, d2: domain) where
                                 (d1.type == d2.type) &&
                                 (isIrregularDom(d1) && isIrregularDom(d2)) {
  var d3: d1.type;
  // These should eventually become forall loops
  for e in d1 do d3.add(e);
  for e in d2 do d3.add(e);
  return d3;
}

proc +(d1: domain, d2: domain) {
  if (isRectangularDom(d1) || isRectangularDom(d2)) then
    compilerError("Cannot add indices to a rectangular domain");
  else
    compilerError("Cannot add indices to this domain type");
}

proc -(d: domain, i: index(d)) {
  if isRectangularDom(d) then
    compilerError("Cannot remove indices from a rectangular domain");
  else
    compilerError("Cannot remove indices from this domain type");
}

proc -(d: domain, i: index(d)) where isIrregularDom(d) {
  d.remove(i);
  return d;
}

proc -(d1: domain, d2: domain) where
                                 (d1.type == d2.type) &&
                                 (isIrregularDom(d1) && isIrregularDom(d2)) {
  var d3: d1.type;
  // These should eventually become forall loops
  for e in d1 do d3.add(e);
  for e in d2 do d3.remove(e);
  return d3;
}

proc -(d1: domain, d2: domain) {
  if (isRectangularDom(d1) || isRectangularDom(d2)) then
    compilerError("Cannot remove indices from a rectangular domain");
  else
    compilerError("Cannot remove indices from this domain type");
}

pragma "inline" proc ==(d1: domain, d2: domain) where isRectangularDom(d1) &&
                                                      isRectangularDom(d2) {
  if d1._value.rank != d2._value.rank then return false;
  for param i in 1..d1._value.rank do
    if (d1.dim(i) != d2.dim(i)) then return false;
  return true;
}

pragma "inline" proc !=(d1: domain, d2: domain) where isRectangularDom(d1) &&
                                                      isRectangularDom(d2) {
  if d1._value.rank != d2._value.rank then return true;
  for param i in 1..d1._value.rank do
    if (d1.dim(i) != d2.dim(i)) then return true;
  return false;
}

pragma "inline" proc ==(d1: domain, d2: domain) where (isAssociativeDom(d1) &&
                                                       isAssociativeDom(d2)) {
  if d1.numIndices != d2.numIndices then return false;
  for idx in d1 do
    if !d2.member(idx) then return false;
  return true;
}

pragma "inline" proc !=(d1: domain, d2: domain) where (isAssociativeDom(d1) &&
                                                       isAssociativeDom(d2)) {
  if d1.numIndices != d2.numIndices then return true;
  for idx in d1 do
    if !d2.member(idx) then return true;
  return false;
}

pragma "inline" proc ==(d1: domain, d2: domain) where (isSparseDom(d1) &&
                                                       isSparseDom(d2)) {
  if d1.numIndices != d2.numIndices then return false;
  if d1._value.parentDom != d2._value.parentDom then return false;
  for idx in d1 do
    if !d2.member(idx) then return false;
  return true;
}

pragma "inline" proc !=(d1: domain, d2: domain) where (isSparseDom(d1) &&
                                                       isSparseDom(d2)) {
  if d1.numIndices != d2.numIndices then return true;
  if d1._value.parentDom != d2._value.parentDom then return true;
  for idx in d1 do
    if !d2.member(idx) then return true;
  return false;
}

// any combinations not handled by the above

pragma "inline" proc ==(d1: domain, d2: domain) param {
  return false;
}

pragma "inline" proc !=(d1: domain, d2: domain) param {
  return true;
}


//
// Array wrapper record
//
pragma "array"
pragma "has runtime type"
record _array {
  var _value;     // stores array class, may be privatized
  var _valueType; // stores type of privatized arrays
  var _promotionType: _value.eltType;

  pragma "inline"
  proc _value {
    if _isPrivatized(_valueType) {
      return chpl_getPrivatizedCopy(_valueType.type, _value);
    } else {
      return _value;
    }
  }

  proc ~_array() {
    if !_isPrivatized(_valueType) {
      on _value {
        var cnt = _value.destroyArr();
        if !noRefCount then
          if cnt == 0 then
            delete _value;
      }
    }
  }

  proc eltType type return _value.eltType;
  proc _dom return _getDomain(_value.dom);
  proc rank param return this.domain.rank;

  pragma "inline"
  proc this(i: rank*_value.dom.idxType) var {
    if isRectangularArr(this) || isSparseArr(this) then
      return _value.dsiAccess(i);
    else
      return _value.dsiAccess(i(1));
  }

  pragma "inline"
  proc this(i: _value.dom.idxType ...rank) var
    return this(i);

  //
  // requires dense domain implementation that returns a tuple of
  // ranges via the getIndices() method; domain indexing is difficult
  // in the domain case because it has to be implemented on a
  // domain-by-domain basis; this is not terribly difficult in the
  // dense case because we can represent a domain by a tuple of
  // ranges, but in the sparse case, is there a general
  // representation?
  //
  proc this(d: domain) where d.rank == rank
    return this((...d.getIndices()));

  proc checkSlice(ranges: range(?) ...rank) {
    for param i in 1.._value.dom.rank do
      if !_value.dom.dsiDim(i).boundsCheck(ranges(i)) then
        halt("array slice out of bounds in dimension ", i, ": ", ranges(i));
  }

  proc this(ranges: range(?) ...rank) {
    if boundsChecking then
      checkSlice((... ranges));
    var d = _dom((...ranges));
    var a = _value.dsiSlice(d._value);
    a._arrAlias = _value;
    pragma "dont disable remote value forwarding"
    proc help() {
      d._value._domCnt$ += 1;
      a._arrAlias._arrCnt$ += 1;
    }
    if !noRefCount then
      help();
    return _newArray(a);
  }

  proc this(args ...rank) where _validRankChangeArgs(args, _value.dom.idxType) {
    if boundsChecking then
      checkRankChange(args);
    var ranges = _getRankChangeRanges(args);
    param rank = ranges.size, stridable = chpl__anyStridable(ranges);
    var d = _dom((...args));
    if !noRefCount then
      d._value._domCnt$ += 1;
    var a = _value.dsiRankChange(d._value, rank, stridable, args);
    a._arrAlias = _value;
    if !noRefCount then
      a._arrAlias._arrCnt$ += 1;
    return _newArray(a);
  }

  proc checkRankChange(args) {
    for param i in 1..args.size do
      if !_value.dom.dsiDim(i).boundsCheck(args(i)) then
        halt("array slice out of bounds in dimension ", i, ": ", args(i));
  }

  // Special cases of local slices for DefaultRectangularArrs because
  // we can't take an alias of the ddata class within that class
  proc localSlice(r: range(?)... rank) where _value.type: DefaultRectangularArr {
    if boundsChecking then
      checkSlice((...r));
    var dom = _dom((...r));
    return chpl__localSliceDefaultArithArrHelp(dom);
  }

  proc localSlice(d: domain) where _value.type: DefaultRectangularArr {
    if boundsChecking then
      checkSlice((...d.getIndices()));

    return chpl__localSliceDefaultArithArrHelp(d);
  }

  proc chpl__localSliceDefaultArithArrHelp(d: domain) {
    if (_value.locale != here) then
      halt("Attempting to take a local slice of an array on locale ",
           _value.locale.id, " from locale ", here.id);
    var A => this(d);
    return A;
  }

  proc localSlice(r: range(?)... rank) {
    if boundsChecking then
      checkSlice((...r));
    return _value.dsiLocalSlice(r);
  }

  proc localSlice(d: domain) {
    return localSlice(d.getIndices());
  }

  pragma "inline"
  proc these() var {
    return _value.these();
  }

  // 1/5/10: do we need this since it always returns domain.numIndices?
  proc numElements return _value.dom.dsiNumIndices;

  proc newAlias() {
    var x = _value.dsiReindex(_value.dom);
    x._arrAlias = _value;
    pragma "dont disable remote value forwarding"
    proc help() {
      _value.dom._domCnt$ += 1;
      x._arrAlias._arrCnt$ += 1;
    }
    if !noRefCount then
      help();
    return _newArray(x);
  }

  proc reindex(d: domain) {
    if rank != d.rank then
      compilerError("illegal implicit rank change");
    for param i in 1..rank do
      if d.dim(i).length != _value.dom.dsiDim(i).length then
        halt("extent in dimension ", i, " does not match actual");

    var newDist = new dmap(_value.dom.dist.dsiCreateReindexDist(d.dims(),
                                                                _value.dom.dsiDims()));
    var newDom = [(...d.dims())] dmapped newDist;
    var x = _value.dsiReindex(newDom._value);
    x._arrAlias = _value;
    pragma "dont disable remote value forwarding"
    proc help() {
      newDom._value._domCnt$ += 1;
      x._arrAlias._arrCnt$ += 1;
    }
    if !noRefCount then
      help();
    return _newArray(x);
  }

  proc writeThis(f: Writer) {
    _value.dsiSerialWrite(f);
  }

  // sparse array interface

  proc IRV var {
    return _value.IRV;
  }

  // associative array interface

  iter sorted() {
    for i in _value.dsiSorted() {
      yield i;
    }
  }

  proc displayRepresentation() { _value.dsiDisplayRepresentation(); }
}  // record _array

//
// Helper functions
//

proc isCollapsedDimension(r: range(?e,?b,?s,?a)) param return false;
proc isCollapsedDimension(r) param return true;


// computes || reduction over stridable of ranges
proc chpl__anyStridable(ranges, param d: int = 1) param {
  for param i in 1..ranges.size do
    if ranges(i).stridable then
      return true;
  return false;
}

// given a tuple args, returns true if the tuple contains only
// integers and ranges; that is, it is a valid argument list for rank
// change
proc _validRankChangeArgs(args, type idxType) param {
  proc _validRankChangeArg(type idxType, r: range(?)) param return true;
  proc _validRankChangeArg(type idxType, i: idxType) param return true;
  proc _validRankChangeArg(type idxType, x) param return false;

  proc help(param dim: int) param {
    if !_validRankChangeArg(idxType, args(dim)) then
      return false;
    else if dim < args.size then
      return help(dim+1);
    else
      return true;
  }

  return help(1);
}

proc chpl__isRange(r: range(?)) param return true;
proc chpl__isRange(r) param return false;

proc _getRankChangeRanges(args) {
  proc _tupleize(x) {
    var y: 1*x.type;
    y(1) = x;
    return y;
  }
  proc collectRanges(param dim: int) {
    if dim > args.size then
      compilerError("domain slice requires a range in at least one dimension");
    if chpl__isRange(args(dim)) then
      return collectRanges(dim+1, _tupleize(args(dim)));
    else
      return collectRanges(dim+1);
  }
  proc collectRanges(param dim: int, x: _tuple) {
    if dim > args.size {
      return x;
    } else if dim < args.size {
      if chpl__isRange(args(dim)) then
        return collectRanges(dim+1, ((...x), args(dim)));
      else
        return collectRanges(dim+1, x);
    } else {
      if chpl__isRange(args(dim)) then
        return ((...x), args(dim));
      else
        return x;
    }
  }
  return collectRanges(1);
}

//
// Support for += and -= over domains
//
proc chpl__isDomain(x: domain) param return true;
proc chpl__isDomain(x) param return false;

proc chpl__isArray(x: []) param return true;
proc chpl__isArray(x) param return false;

//
// Assignment of domains and arrays
//
proc =(a: _distribution, b: _distribution) {
  if a._value == nil {
    return chpl__autoCopy(b.clone());
  } else if a._value._doms.length == 0 {
    if a._value.type != b._value.type then
      compilerError("type mismatch in distribution assignment");
    a._value.dsiAssign(b._value);
    if _isPrivatized(a._value) then
      _reprivatize(a._value);
  } else {
    halt("assignment to distributions with declared domains is not yet supported");
  }
  return a;
}

proc =(a: domain, b: domain) {
  if !isIrregularDom(a) && !isIrregularDom(b) {
    for e in a._value._arrs do {
      on e do e.dsiReallocate(b);
    }
    a.setIndices(b.getIndices());
  } else {
    //
    // BLC: It's tempting to do a clear + add here, but because
    // we need to preserve array values that are in the intersection
    // between the old and new index sets, we use the following
    // instead.
    //
    // TODO: These should eventually become forall loops, hence the
    // warning
    //
    // NOTE: For the current implementation of associative domains,
    // the domain iteration is parallelized, but modification
    // of the underlying data structures (in particular, the _resize()
    // operation on the table) is not thread-safe.  Something more
    // intelligent will likely be needed before it is worth it to
    // parallelize whole-domain assignment for associative arrays.
    //
    compilerWarning("whole-domain assignment has been serialized (see note in $CHPL_HOME/STATUS)");
    for i in a._value.dsiIndsIterSafeForRemoving() {
      if !b.member(i) {
        a.remove(i);
      }
    }
    for i in b {
      if !a.member(i) {
        a.add(i);
      }
    }
  }
  return a;
}

proc =(a: domain, b: _tuple) {
  for ind in 1..b.size {
    a.add(b(ind));
  }
  return a;
}

proc =(d: domain, r: range(?)) {
  d = [r];
  return d;
}

//
// Return true if t is a tuple of ranges that is legal to assign to
// rectangular domain d
//
proc chpl__isLegalRectTupDomAssign(d, t) param {
  proc isRangeTuple(a) param {
    proc isRange(r: range(?e,?b,?s)) param return true;
    proc isRange(r) param return false;
    proc peelArgs(first, rest...) param {
      return if rest.size > 1 then
               isRange(first) && peelArgs((...rest))
             else
               isRange(first) && isRange(rest(1));
    }
    proc peelArgs(first) param return isRange(first);

    return if !isTuple(a) then false else peelArgs((...a));
  }

  proc strideSafe(d, rt, param dim: int=1) param {
    return if dim == d.rank then
             d.dim(dim).stridable || !rt(dim).stridable
           else
             (d.dim(dim).stridable || !rt(dim).stridable) && strideSafe(d, rt, dim+1);
  }
  return isRangeTuple(t) && d.rank == t.size && strideSafe(d, t);
}

proc =(d: domain, rt: _tuple) where chpl__isLegalRectTupDomAssign(d, rt) {
  d = [(...rt)];
  return d;
}

proc =(a: domain, b) {  // b is iteratable
  a._value.clearForIteratableAssign();
  for ind in b {
    a.add(ind);
  }
  return a;
}

pragma "inline" proc =(a: [], b : []) where (a._value.canCopyFromHost && b._value.canCopyFromHost) {
  if a.rank != b.rank then
    compilerError("rank mismatch in array assignment");
  compilerError("GPU to GPU transfers not yet implemented");
}

pragma "inline" proc =(a: [], b : []) where (a._value.canCopyFromDevice && b._value.canCopyFromHost) {
  if a.rank != b.rank then
    compilerError("rank mismatch in array assignment");
  __primitive("copy_gpu_to_host", 
              a._value.data, b._value.data, b._value.eltType, b._value.size);
  return a;
}

pragma "inline" proc =(a: [], b : []) where (a._value.canCopyFromHost && b._value.canCopyFromDevice) {
  if a.rank != b.rank then
    compilerError("rank mismatch in array assignment");
  __primitive("copy_host_to_gpu", 
              a._value.data, b._value.data, b._value.eltType, b._value.size);
  return a;
}

proc chpl__serializeAssignment(a: [], b) param {
  if a.rank != 1 && chpl__isRange(b) then
    return true;

  // Sparse and Opaque arrays do not yet support parallel iteration.  We
  // could let them fall through, but then we get multiple warnings for a
  // single assignment statement which feels like overkill
  //
  if ((!isRectangularArr(a) && !isAssociativeArr(a) && !isSparseArr(a)) ||
      (chpl__isArray(b) &&
       !isRectangularArr(b) && !isAssociativeArr(b) && !isSparseArr(b))) then
    return true;
  return false;
}

pragma "inline" proc =(a: [], b) {
  if (chpl__isArray(b) || chpl__isDomain(b)) && a.rank != b.rank then
    compilerError("rank mismatch in array assignment");
  if chpl__isArray(b) && b._value == nil then
    return a;
  if chpl__serializeAssignment(a, b) {
    compilerWarning("whole array assignment has been serialized (see note in $CHPL_HOME/STATUS)");
    for (aa,bb) in (a,b) do
      aa = bb;
  } else if chpl__tryToken { // try to parallelize using leader and follower iterators
    forall (aa,bb) in (a,b) do
      aa = bb;
  } else {
    compilerWarning("whole array assignment has been serialized (see note in $CHPL_HOME/STATUS)");
    for (aa,bb) in (a,b) do
      aa = bb;
  }
  return a;
}

proc =(a: [], b: _tuple) where isEnumArr(a) || isRectangularArr(a) {
  if isEnumArr(a) {
    if b.size != a.numElements then
      halt("tuple array initializer size mismatch");
    for (i,j) in (chpl_enumerate(index(a.domain)), 1..) {
      a(i) = b(j);
    }
  } else {
    proc chpl__tupleInit(j, param rank: int, b: _tuple) {
      const stride = a.domain.dim(a.rank-rank+1).stride,
            start = a.domain.dim(a.rank-rank+1).first;

      if rank == 1 {
        for param i in 1..b.size {
          j(a.rank-rank+1) = start + (i-1)*stride;
          a(j) = b(i);
        }
      } else {
        for param i in 1..b.size {
          j(a.rank-rank+1) = start + (i-1)*stride;
          chpl__tupleInit(j, rank-1, b(i));
        }
      }
    }
    var j: a.rank*int;
    chpl__tupleInit(j, a.rank, b);
  }
  return a;
}

proc _desync(type t) where t: _syncvar || t: _singlevar {
  var x: t;
  return x.value;
}

proc _desync(type t) {
  var x: t;
  return x;
}

proc =(a: [], b: _desync(a.eltType)) {
  if isRectangularArr(a) {
    forall e in a do
      e = b;
  } else {
    compilerWarning("Whole array assignment has been serialized (see note in $CHPL_HOME/STATUS)");
    for e in a do
      e = b;
  }
  return a;
}

proc by(a: domain, b) {
  var r: a.rank*range(a._value.idxType,
                    BoundedRangeType.bounded,
                    true);
  var t = _makeIndexTuple(a.rank, b, expand=true);
  for param i in 1..a.rank do
    r(i) = a.dim(i) by t(i);
  var d = a._value.dsiBuildRectangularDom(a.rank, a._value.idxType, true, r);
  if !noRefCount then
    if (d.linksDistribution()) then
      d.dist._distCnt$ += 1;
  return _newDomain(d);
}

//
// index for all opaque domains
//
class _OpaqueIndex { }

//
// Swap operators for arrays and domains
//
pragma "inline" proc _chpl_swap(x: [], y: []) {
  for (i,j) in (x.domain, y.domain) do
    x(i) <=> y(j);
}

pragma "inline" proc _chpl_swap(x: domain, y: domain) {
  const t = y;
  y = x;
  x = t;
}

//
// reshape function
//
proc reshape(A: [], D: domain) {
  var B: [D] A.eltType;
  for (i,a) in (D,A) do
    B(i) = a;
  return B;
}

iter linearize(Xs) {
  for x in Xs do yield x;
}

//
// module support for iterators
//
proc iteratorIndex(ic: _iteratorClass) {
  ic.advance();
  return ic.getValue();
}

pragma "expand tuples with values"
proc iteratorIndex(t: _tuple) {
  pragma "expand tuples with values"
  proc iteratorIndexHelp(t: _tuple, param dim: int) {
    if dim == t.size then
      return _build_tuple_always_allow_ref(iteratorIndex(t(dim)));
    else
      return _build_tuple_always_allow_ref(iteratorIndex(t(dim)),
                                           (...iteratorIndexHelp(t, dim+1)));
  }

  return iteratorIndexHelp(t, 1);
}

proc iteratorIndexType(x) type {
  pragma "no copy" var ic = _getIterator(x);
  pragma "no copy" var i = iteratorIndex(ic);
  _freeIterator(ic);
  return i.type;
}

proc _iteratorRecord.writeThis(f: Writer) {
  var first: bool = true;
  for e in this {
    if !first then
      f.write(" ");
    else
      first = false;
    f.write(e);
  }
}

proc =(ic: _iteratorRecord, xs) {
  for (e, x) in (ic, xs) do
    e = x;
  return ic;
}

proc =(ic: _iteratorRecord, x: iteratorIndexType(ic)) {
  for e in ic do
    e = x;
  return ic;
}

pragma "inline" proc _getIterator(x) {
  return _getIterator(x.these());
}

pragma "inline" proc _getIterator(ic: _iteratorClass)
  return ic;

pragma "inline" proc _getIterator(x: _tuple) {
  pragma "inline" proc _getIteratorHelp(x: _tuple, param dim: int) {
    if dim == x.size then
      return tuple(_getIterator(x(dim)));
    else
      return (_getIterator(x(dim)), (..._getIteratorHelp(x, dim+1)));
  }
  if x.size == 1 then
    return _getIterator(x(1));
  else
    return _getIteratorHelp(x, 1);
}

proc _getIterator(type t) {
  compilerError("cannot iterate over a type");
}

proc _checkIterator(type t) {
  compilerError("cannot iterate over a type");
}

pragma "inline" proc _checkIterator(x) {
  return x;
}

pragma "inline"
proc _freeIterator(ic: _iteratorClass) {
  __primitive("chpl_free", ic);
}

pragma "inline"
proc _freeIterator(x: _tuple) {
  for param i in 1..x.size do
    _freeIterator(x(i));
}

pragma "inline"
pragma "no implicit copy"
proc _toLeader(iterator: _iteratorClass)
  return chpl__autoCopy(__primitive("to leader", iterator));

pragma "inline"
proc _toLeader(ir: _iteratorRecord) {
  pragma "no copy" var ic = _getIterator(ir);
  pragma "no copy" var leader = _toLeader(ic);
  _freeIterator(ic);
  return leader;
}

pragma "inline"
proc _toLeader(x: _tuple)
  return _toLeader(x(1));

pragma "inline"
proc _toLeader(x)
  return _toLeader(x.these());

//
// returns lead entity
//
proc chpl__lead(x: _tuple) return chpl__lead(x(1));
proc chpl__lead(x) return x;

//
// return true if any iterator supports fast followers
//
proc chpl__staticFastFollowCheck(x) param {
  pragma "no copy" const lead = chpl__lead(x);
  if chpl__isDomain(lead) || chpl__isArray(lead) then
    return chpl__staticFastFollowCheck(x, lead);
  else
    return false;
}  

proc chpl__staticFastFollowCheck(x: _tuple, lead, param dim = 1) param {
  if x.size == dim then
    return chpl__staticFastFollowCheck(x(dim), lead);
  else
    return chpl__staticFastFollowCheck(x(dim), lead) || chpl__staticFastFollowCheck(x, lead, dim+1);
}

proc chpl__staticFastFollowCheck(x, lead) param {
  return false;
}

proc chpl__staticFastFollowCheck(x: [], lead) param {
  return x._value.dsiStaticFastFollowCheck(lead._value.type);
}

//
// return true if all iterators that support fast followers can use
// their fast followers
//
proc chpl__dynamicFastFollowCheck(x) {
  return chpl__dynamicFastFollowCheck(x, chpl__lead(x));
}

proc chpl__dynamicFastFollowCheck(x: _tuple, lead, param dim = 1) {
  if x.size == dim then
    return chpl__dynamicFastFollowCheck(x(dim), lead);
  else
    return chpl__dynamicFastFollowCheck(x(dim), lead) && chpl__dynamicFastFollowCheck(x, lead, dim+1);
}

proc chpl__dynamicFastFollowCheck(x, lead) {
  return true;
}

proc chpl__dynamicFastFollowCheck(x: [], lead) {
  if chpl__staticFastFollowCheck(x, lead) then
    return x._value.dsiDynamicFastFollowCheck(lead);
  else
    return false;
}

pragma "inline"
pragma "no implicit copy"
proc _toFollower(iterator: _iteratorClass, leaderIndex)
  return chpl__autoCopy(__primitive("to follower", iterator, leaderIndex));

pragma "inline"
proc _toFollower(ir: _iteratorRecord, leaderIndex) {
  pragma "no copy" var ic = _getIterator(ir);
  pragma "no copy" var follower = _toFollower(ic, leaderIndex);
  _freeIterator(ic);
  return follower;
}

pragma "inline"
proc _toFollower(x, leaderIndex) {
  return _toFollower(x.these(), leaderIndex);
}

pragma "inline"
proc _toFollowerHelp(x: _tuple, leaderIndex, param dim: int) {
  if dim == x.size-1 then
    return (_toFollower(x(dim), leaderIndex),
            _toFollower(x(dim+1), leaderIndex));
  else
    return (_toFollower(x(dim), leaderIndex),
            (..._toFollowerHelp(x, leaderIndex, dim+1)));
}

pragma "inline"
proc _toFollower(x: _tuple, leaderIndex) {
  return _toFollowerHelp(x, leaderIndex, 1);
}

pragma "inline"
pragma "no implicit copy"
proc _toFastFollower(iterator: _iteratorClass, leaderIndex, fast: bool) {
  return chpl__autoCopy(__primitive("to follower", iterator, leaderIndex, true));
}

pragma "inline"
proc _toFastFollower(ir: _iteratorRecord, leaderIndex, fast: bool) {
  pragma "no copy" var ic = _getIterator(ir);
  pragma "no copy" var follower = _toFastFollower(ic, leaderIndex, fast=true);
  _freeIterator(ic);
  return follower;
}

pragma "inline"
pragma "no implicit copy"
proc _toFastFollower(iterator: _iteratorClass, leaderIndex) {
  return _toFollower(iterator, leaderIndex);
}

pragma "inline"
proc _toFastFollower(ir: _iteratorRecord, leaderIndex) {
  return _toFollower(ir, leaderIndex);
}

pragma "inline"
proc _toFastFollower(x, leaderIndex) {
  if chpl__staticFastFollowCheck(x) then
    return _toFastFollower(x.these(), leaderIndex, fast=true);
  else
    return _toFollower(x.these(), leaderIndex);
}

pragma "inline"
proc _toFastFollowerHelp(x: _tuple, leaderIndex, param dim: int) {
  if dim == x.size-1 then
    return (_toFastFollower(x(dim), leaderIndex),
            _toFastFollower(x(dim+1), leaderIndex));
  else
    return (_toFastFollower(x(dim), leaderIndex),
            (..._toFastFollowerHelp(x, leaderIndex, dim+1)));
}

pragma "inline"
proc _toFastFollower(x: _tuple, leaderIndex) {
  return _toFastFollowerHelp(x, leaderIndex, 1);
}

proc chpl__initCopy(a: _distribution) {
  pragma "no copy" var b = chpl__autoCopy(a.clone());
  return b;
}

proc chpl__initCopy(a: domain) {
  var b: a.type;
  if isRectangularDom(a) && isRectangularDom(b) {
    b.setIndices(a.getIndices());
  } else {
    // TODO: These should eventually become forall loops, hence the
    // warning
    //
    // NOTE: See above note regarding associative domains
    //
    compilerWarning("whole-domain assignment has been serialized (see note in $CHPL_HOME/STATUS)");
    for i in a do
      b.add(i);
  }
  return b;
}

proc chpl__initCopy(a: []) {
  var b : [a._dom] a.eltType;
  b = a;
  return b;
}

proc chpl__initCopy(ir: _iteratorRecord) {
  iter _ir_copy_recursive(ir) {
    for e in ir do
      yield chpl__initCopy(e);
  }

  pragma "no copy" var irc = _ir_copy_recursive(ir);

  var i = 1, size = 4;
  pragma "insert auto destroy" var D = [1..size];

  // note that _getIterator is called in order to copy the iterator
  // class since for arrays we need to iterate once to get the
  // element type (at least for now); this also means that if this
  // iterator has side effects, we will see them; a better way to
  // handle this may be to get the static type (not initialize the
  // array) and use a primitive to set the array's element; that may
  // also handle skyline arrays
  var A: [D] iteratorIndexType(irc);

  for e in irc {
    pragma "no copy" pragma "insert auto destroy" var ee = e;
    if i > size {
      size = size * 2;
      D = [1..size];
    }
    A(i) = ee;
    i = i + 1;
  }
  D = [1..i-1];
  return A;
}
