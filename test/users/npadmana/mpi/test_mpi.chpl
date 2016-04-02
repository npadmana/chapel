/* A collection of tests/demonstrations of the C-MPI bindings.

These tests are an implementation of the examples here :
  https://computing.llnl.gov/tutorials/mpi/

Each test is in a separate procedure, each of which ends with an
 MPI_Barrier (this is especially important since, in some cases, 
 we only use a subset of the ranks).

*/

/* This initializes MPI by default, you can turn
that off with the autoInitMPI parameter
 */
use MPI; 
use C_MPI; // Include the C-API, to reduce verbosity of the code.
use Random;

const requiredSize = 4;

proc main() {
  if requiredSize != worldSize {
    writef("Please run with at least %i ranks..\n",requiredSize);
    MPI_Abort(MPI_COMM_WORLD, 10);
  }

  hello();
  point2point();
  ring();
  pi();
  test_scatter();

  MPI_Finalize();
}


proc hello() {
  /* Simple test of MPI initialization */
  writef("This is rank %i of %i processes saying Hello, Chapel!\n",worldRank, worldSize);
  MPI_Barrier(MPI_COMM_WORLD);
}

/* Simple blocking point-to-point communication */
proc point2point() {

  if worldRank < 2 {
    var recvpi : real = 0.0,
        sendpi : real = 3.1415926;
    var stat : MPI_Status;

    if worldRank==0 {
      MPI_Send(sendpi, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
      MPI_Recv(recvpi, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, stat);
    } else {
      MPI_Recv(recvpi, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, stat);
      MPI_Send(sendpi, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    var count = stat.getCount(MPI_DOUBLE);

    writef("Task %i: Received %i pi=%r from task %i with tag %i\n",
        worldRank, count, recvpi, stat.MPI_SOURCE, stat.MPI_TAG);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/* Non-blocking communication in a ring */
proc ring() {

  var left = mod(worldRank-1, worldSize);
  var right = mod(worldRank+1, worldSize);
  var toleft : c_int = 1;
  var toright : c_int = 2;
  var fromleft = toright,
      fromright = toleft;


  var buf : [1..2]int(32);
  var requests : [1..4]MPI_Request;
  var status : [1..4]MPI_Status;

  MPI_Irecv(buf[1], 1, MPI_INT, left, fromleft, MPI_COMM_WORLD, requests[1]);
  MPI_Irecv(buf[2], 1, MPI_INT, right, fromright, MPI_COMM_WORLD, requests[2]);

  MPI_Isend(worldRank, 1, MPI_INT, left, toleft, MPI_COMM_WORLD, requests[3]);
  MPI_Isend(worldRank, 1, MPI_INT, right, toright, MPI_COMM_WORLD, requests[4]);

  MPI_Waitall(4, requests, status);

  writef("Rank %i recieved %i from the left, and %i from the right\n",worldRank, buf[1], buf[2]);
  MPI_Barrier(MPI_COMM_WORLD);
}

/* Compute pi -- test MPI_Reduce */
proc pi() {
  const N=10000;
  var x : [0.. #N] real;
  var y : [0.. #N] real;

  fillRandom(x); fillRandom(y);
  var sum = 0.0;
  forall (x1,y1) in zip(x,y) with (+ reduce sum) {
    if (x1*x1+y1*y1) <= 1 then sum += 1;
  }
  sum /= N;

  var summ : real;
  MPI_Reduce(sum, summ, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if worldRank==0 {
    summ *= 4/worldSize;
    writef("I estimate pi to be %r\n",summ);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/* Test scatters. */
proc test_scatter() {
  var arr : [{0..3,0..3}]real(32);
  forall (i,j) in arr.domain {
    arr[i,j] = (i*4 + j):real(32);
  }
  var recbuf : [0..3]real(32);

  // Use MPI Scatter
  {
    MPI_Scatter(arr[0,0], 4, MPI_FLOAT, recbuf[0], 4, MPI_FLOAT, 0, MPI_COMM_WORLD);
    writef("Rank %i :",worldRank); writeln(recbuf);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Use a contiguous data type
  {
    arr += 1:real(32);
    var rowtype : MPI_Datatype;
    var stat : MPI_Status;

    MPI_Type_contiguous(4, MPI_FLOAT, rowtype);
    MPI_Type_commit(rowtype);

    if worldRank==0 {
      for irank in 0.. #worldSize do MPI_Send(arr[irank,0], 1, rowtype, irank, 1, MPI_COMM_WORLD);
    }
    MPI_Recv(recbuf[0], 4, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, stat);
    writef("Rank %i :",worldRank); writeln(recbuf);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(rowtype);
  }

  // Use a vector data dtype - send columns
  {
    var coltype : MPI_Datatype;
    var stat : MPI_Status;

    MPI_Type_vector(4, 1, 4, MPI_FLOAT, coltype);
    MPI_Type_commit(coltype);

    if worldRank==0 {
      for irank in 0.. #worldSize do MPI_Send(arr[0,irank], 1, coltype, irank, 1, MPI_COMM_WORLD);
    }
    MPI_Recv(recbuf[0], 4, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, stat);
    writef("Rank %i :",worldRank); writeln(recbuf);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(coltype);
  }

}
