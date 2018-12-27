#include <mpi.h>

class parallel {

   private:


   public:

   MPI_Comm comm;
   int nprocs, rank, root;
   bool ismaster, initialized;

   parallel();
   ~parallel();
   parallel(int, char**);
   parallel(int);
   int check_init();
   int initialize(int, char**);
   int initialize(int);
   int finalize();
   int error(char *);
   int sync(int comm=-1);

};
