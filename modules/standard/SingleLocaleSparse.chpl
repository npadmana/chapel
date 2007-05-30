use Search;

class SingleLocaleSparseDomain: BaseSparseArithmeticDomain {
  param rank : int;
  type dim_type;
  var parentDom: BaseArithmeticDomain;
  var nnz = 0;  // intention is that user might specify this to avoid reallocs
  //  type ind_type = rank*dim_type;

  var nnzDomSize = nnz;
  var nnzDom = [1..nnzDomSize];

  var indices: [nnzDom] index(rank);

  def initialize() {
    nnz = 0;
  }

  def numIndices return nnz;

  def getIndices() return 0;
  def setIndices(x);

  def buildArray(type eltType)
    return SingleLocaleSparseArray(eltType, rank, dim_type, dom=this);

  def buildEmptyDomain()
    return SingleLocaleSparseDomain(rank=rank, dim_type=dim_type, 
                                    parentDom = BaseArithmeticDomain());

  iterator ault() {
    for i in 1..nnz {
      yield indices(i);
    }
  }

  def dim(d : int) {
    return parentDom.bbox(d);
  }

  def find(ind: rank*dim_type) {
    return BinarySearch(indices, ind, 1, nnz);
  }

  def member(ind: rank*dim_type) {
    const (found, loc) = find(ind);
    return found;
  }

  def add(ind: rank*dim_type) {
    // find position in nnzDom to insert new index
    const (found, insertPt) = find(ind);

    // if the index already existed, then return
    if (found) then return;

    // increment number of nonzeroes
    nnz += 1;

    // double nnzDom if we've outgrown it; grab current size otherwise
    var oldNNZDomSize = nnzDomSize;
    if (nnz > nnzDomSize) {
      nnzDomSize = if (nnzDomSize) then 2*nnzDomSize else 1;

      nnzDom = [1..nnzDomSize];
    }

    // shift indices up
    for i in [insertPt..nnz) by -1 {
      indices(i+1) = indices(i);
    }

    indices(insertPt) = ind;

    // shift all of the arrays up and initialize nonzeroes if
    // necessary 
    //
    // BLC: Note: if arithmetic arrays had a user-settable
    // initialization value, we could set it to be the IRV and skip
    // this second initialization of any new values in the array.
    // we could also eliminate the oldNNZDomSize variable
    for a in _arrs {
      a.sparseShiftArray(insertPt..nnz-1, oldNNZDomSize+1..nnzDomSize);
    }
  }

  iterator dimIter(param d, ind) {
    if (d != rank-1) {
      compilerError("dimIter() not supported on sparse domains for dimensions other than the last");
    }
    halt("dimIter() not yet implemented for sparse domains");
    yield indices(1);
  }
}


class SingleLocaleSparseArray: BaseArray {
  type eltType;
  param rank : int;
  type dim_type;

  var dom : SingleLocaleSparseDomain(rank=rank, dim_type=dim_type);
  var data: [dom.nnzDom] eltType;
  var irv: eltType;

  //  def this(ind: dim_type ... 1) var where rank == 1
  //    return this(ind);

  def this(ind: rank*dim_type) {
    // make sure we're in the dense bounding box
    if boundsChecking then
      if !((dom.parentDom).member(ind)) then
        halt("array index out of bounds: ", ind);

    // lookup the index and return the data or IRV
    const (found, loc) = dom.find(ind);
    return if (found) then data(loc) else irv;
  }


  def =this(ind: rank*dim_type, val:eltType) {
    // make sure we're in the dense bounding box
    if boundsChecking then
      if !((dom.parentDom).member(ind)) then
        halt("array index out of bounds: ", ind);

    // lookup the index and return the data or IRV
    const (found, loc) = dom.find(ind);
    if found then
      data(loc) = val;
    else
      halt("attempting to assign a 'zero' value in a sparse array: ", ind);
  }

  def IRV var {
    return irv;
  }

  def sparseShiftArray(shiftrange, initrange) {
    for i in initrange {
      data(i) = irv;
    }
    for i in shiftrange by -1 {
      data(i+1) = data(i);
    }
    data(shiftrange.low) = irv;
  }
}


def SingleLocaleSparseDomain.writeThis(f: Writer) {
  f.writeln("[");
  if (nnz >= 1) {
    var prevInd = indices(1);
    write(" ", prevInd);
    for i in 2..nnz {
      if (prevInd(1) != indices(i)(1)) {
        writeln();
      }
      prevInd = indices(i);
      write(" ", prevInd);
    }
  }
  f.writeln("\n]");
}


def SingleLocaleSparseArray.writeThis(f: Writer) {
  if (dom.nnz >= 1) {
    var prevInd = dom.indices(1);
    write(data(1));
    for i in 2..dom.nnz {
      if (prevInd(1) != dom.indices(i)(1)) {
        writeln();
      } else {
        write(" ");
      }
      prevInd = dom.indices(i);
      write(data(i));
    }
    writeln();
  }
}
