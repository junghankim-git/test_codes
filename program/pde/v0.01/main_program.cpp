#include <stdio.h>
#include <iostream>
using namespace std;
//#include <netcdf.h>
#include <parallel.h>
/*

   abstract: main program

   history log:
      date         name             contents
      20170518     junghan kim      first written

*/
int main(int argc, char **argv) {

   int err, val;

   // parallel *par;
   parallel *gpar = new parallel(argc, argv);
   parallel  *par = new parallel(gpar->comm);

   if (par->ismaster) {
      val = 1;
   } else {
      val = 0;
   }

   printf("before: rank %d: %d\n",par->rank,val);
   fflush(stdout);
   err = par->sync();
   err = par->bcast(&val,1,MPI_INT);
   err = par->sync();
   printf("after : rank %d: %d\n",par->rank,val);
   fflush(stdout);

   delete par, gpar;

   return 0

}
