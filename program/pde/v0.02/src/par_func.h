#include <mpi.h>
#include <parallel.h>

class par_func {

   private:
   parallel *par;


   public:
   bool initialized;
   par_func();
   ~par_func();
   par_func(parallel *);
   int bcast(void *, int, int, int root=-1, int comm=-1);

};
