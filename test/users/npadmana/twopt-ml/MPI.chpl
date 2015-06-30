/* A very basic MPI wrapper, to enable simple 
   interaction with MPI codes.

   This is a workaround to run chapel in SPMD mode; it will 
   think that's running only on a single locale. This turned out
   to be necessary to deal with very significant slowdowns in Chapel
   when running in multilocale mode.

   Nikhil Padmanabhan, Feb 14, 2015
*/

module MPI {

  use "mpi.h";
  use SysCTypes;

  // External types
  extern type MPI_Comm; // Opaque type for MPI communicator 
  extern type MPI_Datatype; 
  extern type MPI_Op;

  // MPI types
  extern const MPI_BYTE : MPI_Datatype;
  extern const MPI_INT : MPI_Datatype;
  extern const MPI_LONG : MPI_Datatype;
  extern const MPI_FLOAT : MPI_Datatype;
  extern const MPI_DOUBLE : MPI_Datatype;

  // MPI Ops
  extern const MPI_SUM : MPI_Op;
  extern const MPI_IN_PLACE : c_void_ptr;

  // Define MPI_COMM_WORLD
  extern const MPI_COMM_WORLD : MPI_Comm;

  // The MPI functions we need
  extern proc MPI_Init(argc, argv) : c_int;
  extern proc MPI_Initialized(ref flag : c_int) : c_int;
  extern proc MPI_Finalize() : c_int;
  extern proc MPI_Barrier(comm : MPI_Comm) : c_int;
  extern proc MPI_Comm_rank(comm : MPI_Comm, ref rank : c_int) : c_int;
  extern proc MPI_Comm_size(comm : MPI_Comm, ref size : c_int) : c_int;
  extern proc MPI_Bcast(buffer : c_void_ptr, count : c_int, datatype : MPI_Datatype,
      root : c_int, comm : MPI_Comm) : c_int;
  extern proc MPI_Reduce(const sendbuf : c_void_ptr, recvbuf : c_void_ptr, 
      count : c_int, datatype : MPI_Datatype, op : MPI_Op, 
      root : c_int, comm : MPI_Comm) : c_int;

  // Did we initialize MPI?
  var isInit : bool = false; 


  // Convenience wrappers
  proc finalize() {
    if isInit {
      MPI_Finalize();
    }
  }

  proc rank() : int {
    var _rank : c_int;
    MPI_Comm_rank(MPI_COMM_WORLD,_rank);
    return _rank;
  }
  proc size() : int {
    var _size : c_int;
    MPI_Comm_size(MPI_COMM_WORLD,_size);
    return _size;
  }


  // Module initialization will set this up
  var _flag : c_int;
  MPI_Initialized(_flag);
  if (_flag == 0) {
    isInit = true;
    MPI_Init(0, 0); // Send NULLS in.
  }

  var Rank = rank();
  var Size = size();


}
