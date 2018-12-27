#include <stdio.h>
#include <iostream>
//#include <parallel.h>
using namespace std;

#include <mpi.h>

class parallel {

   private:
   MPI_Comm comm;
   int nprocs, rank, root;
   bool ismaster;


   public:
   parallel();
   //parallel(int, char**);
   ~parallel();
   int ini(int, char**);
   int fin();

};


parallel::parallel() {


}


/*
parallel::parallel(int argc, char **argv) {

   int err;
   err = this->ini(&argc, &argv);

}
*/



parallel::~parallel() {

   int err;
   err = this->fin();

}


/*
int parallel::ini(int argc, char** argv) {

   int err;

int i;
   printf(" number of args = %d\n", argc);
   for (i=0;i<argc;i++){
     printf(" %d argument  = %s\n", i, argv[i]);
   }
printf("step 1\n");
   err  = MPI_Init(&argc, &argv);
   comm = MPI_COMM_WORLD;
printf("step 2\n");
   err  = MPI_Comm_size(comm, &nprocs);
printf("step 3\n");
   err  = MPI_Comm_rank(comm, &rank);
printf("step 4\n");
   root = 0;
   ismaster = false;
   if (rank==root) ismaster = true;

}
*/

int parallel::ini(int argc, char **argv) {

   int i, err; //, nprocs, rank;
//   MPI_Group grp;
//   MPI_Comm comm;

   printf(" number of args = %d\n", argc);
   for (i=0;i<argc;i++){
     printf(" %d argument  = %s\n", i, argv[i]);
   }

printf("step 0: %d \n",nprocs);
printf("step 1\n");
   err = MPI_Init(&argc, &argv);
printf("step 2: %d \n",err);
   err = MPI_Comm_dup(MPI_COMM_WORLD, &comm);
   //err = MPI_Comm_group(MPI_COMM_WORLD, &grp);
   //err = MPI_Comm_create(MPI_COMM_WORLD, grp, &comm);
printf("step 3: %d \n",err);
   err  = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
printf("step 4: %d \n",err);
   err  = MPI_Comm_rank(comm, &rank);
printf("step 5: %d \n",err);
   root = 0;
   ismaster = false;
   if (rank==root) ismaster = true;

   return err;

}

int parallel::fin() {

   int err;
   err = MPI_Finalize();
   return err;

}



int main(int argc, char **argv) {

   int ierr;

printf("main 0\n");
   parallel *par = new parallel();

printf("main 1\n");
   par->ini(argc, argv);
printf("main 2\n");

   par->fin();
printf("main 3\n");

}
