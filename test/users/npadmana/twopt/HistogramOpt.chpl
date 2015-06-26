module HistogramOpt {

  class UniformBins {
    param dim : int;
    param nparallel : int;

    // Locks
    var locks$ : [0.. #nparallel] sync bool;
    var curid : atomic int;

    var nbins : dim*int;
    var Dhist : domain(dim);
    var lo, hi,dx, invdx : [1..dim] real;
    var arr : [0.. #nparallel][Dhist] real;


    proc UniformBins(param nparallel : int, param dim : int, nbins : dim*int, limits : dim*(real,real)) {
      var dd : dim*range;
      this.nbins = nbins;
      for param ii in 1..dim {
        lo[ii] = limits(ii)(1);
        hi[ii] = limits(ii)(2);
        dx[ii] = (hi[ii]-lo[ii])/nbins(ii);
        dd(ii) = 0.. #nbins(ii);
      }
      invdx = 1.0/dx;
      Dhist = {(...dd)};
      [arr1 in arr] arr1 = 0.0;
      curid.write(0);
    }

    proc lock() : int {
      var tid = curid.fetchAdd(1) % nparallel;
      locks$[tid] = true;
      return tid;
    }

    proc unlock(tid : int) {
      locks$[tid];
    }

    proc reset() {
      [arr1 in arr] arr1 = 0.0;
    }

    proc bins(idim : int) {
      var bounds : [0.. #nbins(idim)+1] real;
      for i in {0.. #(nbins(idim)+1)} do bounds[i] = lo[idim] + dx[idim]*i;
      return bounds;
    }

    // Note that this operation is not thread-safe by itself.
    proc add(tid : int, x : dim*real, w : real) {
      for param ii in 1..dim do
        if ((x(ii) < lo[ii]) | (x(ii) >= hi[ii])) then return;

      var pos : dim*int;
      for param ii in 1..dim do pos(ii) = ((x(ii)-lo(ii))*invdx[ii]) : int;
      arr[tid][pos] += w;
    }

    // Note that this operation is not thread-safe
    proc set(tid : int, ndx : dim*int, val : real) {
      arr[tid][ndx] = val;
    }

    // combine
    proc combine() {
      [arr1 in arr[1..]] arr[0] += arr1;
    }

    proc this(ndx) : real {
      return arr[0][ndx];
    }

  } // UniformBins

  proc writeHist(ff : channel, hh : UniformBins, fmt : string = "%20.14er ") {
    // Dump out values
    for xx in hh.bins(1) do ff.writef("%12.4dr",xx); 
    ff.writeln();
    for xx in hh.bins(2) do ff.writef("%12.4dr",xx); 
    ff.writeln("\n##");
    for ii in hh.Dhist.dim(1) {
      for jj in hh.Dhist.dim(2) {
        ff.writef(fmt, hh[(ii,jj)]);
      }
      ff.writeln();
    }
  }


} // End module Histogram
