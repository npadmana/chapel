use MPI;
use HistogramOpt;
use Random;
use Time;
use Memory;

// Test flags
config const isTest=false;
config const isPerf=false;

// Performance test parameters
config const nParticles=10000;

// Input/Output filenames
config const fn1 = "test.dat";
config const fn2 = "test.dat";
config const pairfn = "test-DD.dat";

// Dimension constanst
param NDIM=3;
param NEXT=2; // w, r2
param NTOT=NDIM+NEXT;
const W=NDIM;
const R2=NDIM+1;
const ParticleAttrib = 0.. #NTOT;
const DimSpace = 0.. #NDIM;

// Tree parameters
config const minpart=500;

// Histogram parameters
config const smax=200.0;
config const smax2=smax**2;
config const nmubins=5;
config const nsbins=5;
config param nParHist : int = 100;

// Testing variables
var nspawn : atomic uint;
nspawn.write(0);

proc main() {
  if memTrack then startVerboseMem();
  doPairs();
  MPI.finalize();
  if memTrack {
    stopVerboseMem();
    printMemAllocStats();
  }
}


class Particle3D {
  var npart : int = 0;
  var Darr : domain(2);
  var Dpart : domain(1);
  var arr : [Darr] real;

  var _tmp : [Dpart] real;
  var _n1, _ndx : [Dpart] int;

  proc Particle3D(npart1 : int, random : bool = false) {
    npart = npart1;
    Darr = {ParticleAttrib, 0.. #npart};
    Dpart = {0.. #npart};
    if random {
      var rng = new RandomStream();
      var x, y, z : real;
      for ii in Dpart {
        x = rng.getNext()*1000.0; y = rng.getNext()*1000.0; z = rng.getNext()*1000.0;
        arr[0,ii] = x; arr[1, ii] = y; arr[2, ii] = z;
        arr[3,ii] = 1.0;
        arr[4,ii] = x**2 + y**2 + z**2;
      }
      delete rng;
    }
  }

  proc this(ii : int) : NTOT*real {
    // ??? Hardcoded hack!
    return (arr[0,ii],arr[1,ii],arr[2,ii],arr[3,ii],arr[4,ii]);
  }


  // Reorder all the arrays according to _ndx (we assume you have set this)
  proc reorder(idom : domain(1)) {

    for kk in ParticleAttrib {
      forall ii in idom {
        _tmp[_ndx[ii]] = arr[kk,ii];
      }
      forall ii in idom {
        arr[kk,ii] = _tmp[ii];
      }
    }
  }

  // Shuffles the elements of the array... 
  // Note that the random number generator cannot cover the space
  // of all permutations, but that's fine for what we need.
  proc shuffle() {
    // Fill in _ndx
    forall ii in Dpart {
      _ndx[ii] = ii;
    }

    // Set the random number generator
    var rng = new RandomStream(41);
    var jj : int;
    for ii in 0..(npart-2) {
      jj = (rng.getNext()*(npart-ii)):int + ii;
      _ndx[jj] <=> _ndx[ii];
    }
    delete rng;

    reorder(Dpart);
  }

  
}

// Splits the domain into n pieces. Assumes the domain is contiguous
proc splitDomain(indom : domain(1), n : int) {
  var nel = indom.size;
  var lo = indom.low;
  var hi = indom.high;
  var nsplit = nel/n;
  var ret : [0.. #n] domain(1);
  for ii in 0.. #(n-1) {
    ret[ii] = {lo..(lo+nsplit-1)};
    lo+=nsplit;
  }
  ret[n-1] = {lo..hi};
  return ret;
}

  


proc countLines(fn : string) : int {
  var ff = open(fn, iomode.r);
  var ipart = 0;
  for iff in ff.lines() do ipart +=1;
  ff.close();
  return ipart;
}

proc readFile(fn : string) : Particle3D  {
  var npart = countLines(fn);
  var pp = new Particle3D(npart);

  var ff = openreader(fn);
  var ipart = 0;
  var x,y,z,w,r2 : real;
  while (ff.read(x,y,z,w)) {
    r2 = x**2 + y**2 + z**2;
    pp.arr[..,ipart] = [x,y,z,w,r2];
    ipart += 1;
  }
  ff.close();

  return pp;
}

proc syncParticles(root : int, ref pp : Particle3D) {
  var np : int;
  if (root == MPI.Rank) {
    np = pp.npart;
  } 
  MPI_Bcast(c_ptrTo(np) : c_void_ptr, 1, MPI_LONG, root : c_int, MPI_COMM_WORLD);
  if (root != MPI.Rank) then pp = new Particle3D(np);
  var nn : c_int = (NTOT*np) : c_int;
  MPI_Bcast(c_ptrTo(pp.arr[0,0]) : c_void_ptr, nn, MPI_DOUBLE, root : c_int, MPI_COMM_WORLD);
}

proc splitOn(pp : Particle3D, idom : domain(1), splitDim : int, xsplit : real) : int {
  // Setup for a prefix scan
  forall ii in idom {
    if (pp.arr[splitDim,ii] < xsplit) {
      pp._n1[ii] = 1;
    } else {
      pp._n1[ii] = 0;
    }
  }

  var lnpart : int = 0;
  for ii in idom do lnpart += pp._n1[ii];
  var li = idom.low;
  var ri = li+lnpart;
  for ii in idom {
    if (pp._n1[ii] == 1) {
      pp._ndx[ii] = li;
      li += 1;
    } else {
      pp._ndx[ii] = ri;
      ri += 1;
    }
  }

  pp.reorder(idom);

  return lnpart;
}
      
class KDNode {
  var lo, hi,npart,id : int;
  var dom : domain(1);
  var xcen : [DimSpace]real;
  var rcell : real;
  var left, right : KDNode;

  proc isLeaf() : bool {
    return (left==nil) && (right==nil);
  }

  proc ~KDNode() {
    if left then delete left;
    if right then delete right;
  }

}


proc BuildTree(pp : Particle3D, dom : domain(1), id : int) : KDNode  {
  var me : KDNode = new KDNode();
  me.lo = dom.low;
  me.hi = dom.high;
  me.dom = dom;
  me.id = id;
  me.npart = (me.hi-me.lo)+1;

  //  work out xcen and vantage point radius
  var pmin, pmax : [DimSpace] real;
  for idim in DimSpace {
    pmin[idim] = min reduce pp.arr[idim,me.lo..me.hi];
    pmax[idim] = max reduce pp.arr[idim,me.lo..me.hi];
  }
  me.xcen = (pmax+pmin)/2.0;
  me.rcell = 0.0;
  var r1 : real;
  for ii in me.dom {
    r1 = 0.0;
    for idim in DimSpace {
      r1 += (pp.arr[idim,ii]-me.xcen[idim])**2;
    }
    r1 = sqrt(r1);
    if (r1 > me.rcell) then me.rcell = r1;
  }

  // Continue to split
  if (me.npart <= minpart) then return me; // Don't split further.

  // Find dimension to split on
  var dx : [DimSpace] real;
  dx = pmax - pmin; 
  var splitDim = 0;
  for idim in DimSpace {
    if (dx[idim] > dx[splitDim]) then splitDim=idim;
  }

  // Split
  var lnpart = splitOn(pp,me.dom,splitDim, me.xcen[splitDim]);
  me.left  = BuildTree(pp, {me.lo..me.lo+lnpart-1}, 2*id+1);
  me.right = BuildTree(pp, {me.lo+lnpart..me.hi}, 2*id+2);
  return me;
}

record Job {
  var d1, d2 : domain(1);
}

proc TreeWalk(node1 : KDNode, node2 : KDNode, jobs : []Job, addjob : bool, ref count : int) {
  // Compute the distance between node1 and node2
  var rr = sqrt (+ reduce(node1.xcen - node2.xcen)**2);
  var rmin = rr - (node1.rcell+node2.rcell);

  // If distance is greater than all cases
  if (rmin > smax) then return;

  // If both nodes are leaves
  if (node1.isLeaf() & node2.isLeaf()) {
    if addjob {
      jobs[count] = new Job(node1.dom, node2.dom);
    }
    count += 1;
    return;
  }

  // If one node is a leaf 
  if (node1.isLeaf()) {
    TreeWalk(node1, node2.left, jobs, addjob, count);
    TreeWalk(node1, node2.right,jobs, addjob, count);
    return;
  }
  if (node2.isLeaf()) {
    TreeWalk(node1.left, node2, jobs, addjob, count);
    TreeWalk(node1.right, node2, jobs, addjob, count);
    return;
  }

  // Split the larger case;
  if (node1.npart > node2.npart) {
    TreeWalk(node1.left, node2, jobs, addjob, count);
    TreeWalk(node1.right, node2, jobs, addjob, count);
    return;
  } else {
    TreeWalk(node1, node2.left, jobs, addjob, count);
    TreeWalk(node1, node2.right, jobs, addjob, count);
    return;
  }

}


proc TreeAccumulate(hh : UniformBins, p1 : Particle3D, p2 : Particle3D, node1 : KDNode, node2 : KDNode, scale : real=1) {
  var nspawn : int = 0;
  // Define a record of jobs
  var Djobs : domain(1);
  var joblist : [Djobs] Job;

  TreeWalk(node1,node2,joblist, false, nspawn);
  //writef("%i jobs found...\n",nspawn);
  Djobs = {0.. #nspawn};
  nspawn=0;
  TreeWalk(node1,node2,joblist, true, nspawn);
  //writef("%i jobs queued...\n",nspawn);


  coforall itask in 0.. #nParHist {
    var i1 = itask;
    while (i1 < nspawn) {
      smuAccumulate(itask, hh, p1, p2, joblist[i1].d1, joblist[i1].d2,scale);
      i1 += nParHist;
    }
  }

}
  

// The basic pair counter
proc smuAccumulate(tid : int, hh : UniformBins, p1,p2 : Particle3D, d1,d2 : domain(1), scale : real) {
  var x1,y1,z1,w1,r2 : real;
  var sl, s2, l1, s1, l2, mu, wprod : real;
  
  for ii in d1 { // Loop over first set of particles
   
    (x1,y1,z1,w1,r2) = p1[ii];

    for jj in d2 { // Second set of particles
      mu=2*(p2.arr[0,jj]*x1 + p2.arr[1,jj]*y1 + p2.arr[2,jj]*z1);
      sl = r2 - p2.arr[R2,jj];
      l1 = r2 + p2.arr[R2,jj];
      s2 = l1 - mu;
      l2 = l1 + mu;
      if ((s2 < smax2) && (s2 > 1.0e-20)) {
        wprod = scale * w1 * p2.arr[W,jj];
        s1 = sqrt(s2);
        mu = sl/(s1*sqrt(l2));
        if (mu < 0) then mu = -mu;

        hh.add(tid,(s1,mu),wprod);
      }
    }
  }
}

proc doPairs() {
  // Informational
  if (!isTest & (MPI.Rank==0)) then writeln("Running on ", MPI.Size," ranks");

  var tt : Timer;

  // Read in the file
  var pp1, pp2 : Particle3D;
  if (MPI.Rank == 0) {
    tt.clear(); tt.start();
    if isPerf {
      pp1 = new Particle3D(nParticles, true);
      pp2 = new Particle3D(nParticles, true);
    } else {
      pp1 = readFile(fn1);
      pp2 = readFile(fn2);
      if !isTest {
        writef("Read in %i lines from file %s \n", pp1.npart, fn1);
        writef("Read in %i lines from file %s \n", pp2.npart, fn2);
      }
    }
    tt.stop();
    if !isTest {
      writef("Time to read : %r \n", tt.elapsed());
    }

    // Shuffle the particles -- this is not used here, but will be usef for the multilocale version
    tt.clear(); tt.start();
    pp1.shuffle();
    pp2.shuffle();
    tt.stop();
    if !isTest then writef("Time to shuffle : %r \n",tt.elapsed());
  }

  syncParticles(0, pp1);
  syncParticles(0, pp2);

  // Split
  var split1 = splitDomain(pp1.Dpart,MPI.Size);
  var split2 = splitDomain(pp2.Dpart,MPI.Size);


  // Build the tree
  tt.clear(); tt.start();
  var root1 : [0.. #MPI.Size] KDNode;
  var root2 : [0.. #MPI.Size] KDNode;
  forall ii in 0.. #MPI.Size {
    root1[ii] = BuildTree(pp1, split1[ii], 0);
    root2[ii] = BuildTree(pp2, split2[ii], 0);
  }
  tt.stop();
  if !isTest then writef("Time to build trees : %r \n", tt.elapsed());
  

  // Set up the histogram
  var hh = new UniformBins(2,nParHist,(nsbins,nmubins), ((0.0,smax),(0.0,1.0+1.e-10)));

  // Do the paircounts with a tree
  hh.reset();
  tt.clear(); tt.start();
  var ctr = 0;
  for iroot1 in root1 {
    for iroot2 in root2 {
      if (ctr%MPI.Size == MPI.Rank) then TreeAccumulate(hh,pp1,pp2,iroot1,iroot2,1.0);
      ctr += 1;
    }
  }
  hh.combine();
  if (MPI.Rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, c_ptrTo(hh.arr[0][0,0]) : c_void_ptr,
        (nsbins*nmubins):c_int, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(c_ptrTo(hh.arr[0][0,0]) : c_void_ptr, MPI_IN_PLACE, 
        (nsbins*nmubins):c_int, MPI_DOUBLE,MPI_SUM, 0, MPI_COMM_WORLD);
  }
  tt.stop();

  if (MPI.Rank == 0) {
    if (!isTest) {
      writef("Time to tree paircount : %r \n", tt.elapsed());
      if !isPerf {
        var ff = openwriter("%s.tree".format(pairfn));
        writeHist(ff,hh);
        ff.close();
      }
    } else {
      hh.set(0,(0,0),0.0);
      writeHist(stdout,hh,"%20.5er ");
    }
  }

  //
  // clean up
  //
  delete pp1;
  delete pp2;
  forall ii in 0.. #MPI.Size {
    delete root1[ii];
    delete root2[ii];
  }
  delete hh;
}

