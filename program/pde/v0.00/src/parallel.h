#include <mpi.h>

class parallel {

   private:


   public:
   MPI_Comm comm;
   int nprocs, rank, root;
   bool ismaster, initialized;
   parallel();
   parallel(int, char**);
   parallel(int);
   ~parallel();
   int initialize(int, char**);
   int initialize(int);
   int finalize();
   int error(char *);
   int sync();
   int sync(int);
   int bcast(void *, int, int, int root=-1, int comm=-1);

};
