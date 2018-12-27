#include <stdio.h>
#include <iostream>
using namespace std;
#include <parallel.h>


parallel::parallel() {

   initialized = false;

}


parallel::parallel(int argc, char **argv) {

   int err;
   err = this->initialize(argc, argv);

}


parallel::parallel(int comm) {

   int err;
   err = this->initialize(comm);

}



parallel::~parallel() {

   int err;
   err = this->finalize();

}



int parallel::initialize(int argc, char **argv) {

   int i, err;

   err  = MPI_Init(&argc, &argv);
   err  = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
   err  = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   err  = MPI_Comm_rank(comm, &rank);
   root = 0;
   ismaster    = false;
   initialized = true;
   if (rank==root) {
      ismaster = true;
      printf("parallel: initialized...\n");
   }

   return err;

}



int parallel::initialize(int comm) {

   int err, code;
   int isinit = 0;

   err  = MPI_Initialized(&isinit);
   if (isinit==0) {
      err = MPI_Abort(MPI_COMM_WORLD, code);
   }

   err  = MPI_Comm_dup(comm, &this->comm);
   err  = MPI_Comm_size(comm, &nprocs);
   err  = MPI_Comm_rank(comm, &rank);
   root = 0;
   ismaster    = false;
   initialized = true;
   if (rank==root) {
      ismaster = true;
      printf("parallel: initialized...\n");
   }

   return err;

}



int parallel::finalize() {

   int err;
   if (!initialized) printf("parallel: not initialized...\n");
   if (rank==root)   printf("parallel: finalized...\n");
   err = MPI_SUCCESS;
   if (initialized)  err = MPI_Finalize();
   return err;

}



int parallel::error(char *string) {

   int err, code;
   if (rank==root) printf("parallel: %s\n",string);
   err = MPI_Abort(comm, code);
   return err;

}



int parallel::sync() {

   int err;
   err = MPI_Barrier(comm);
   return err;

}



int parallel::sync(int comm) {

   int err;
   err = MPI_Barrier(comm);
   return err;

}



int parallel::bcast(void *buf, int cnt, int vtype, int root, int comm) {

   int err;
   if (root==-1) root = this->root;
   if (comm==-1) comm = this->comm;
   
   err = MPI_Bcast(buf,cnt,vtype,root,comm);

}


