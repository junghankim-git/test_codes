#include <stdio.h>
#include <iostream>
using namespace std;
//#include <netcdf.h>
//#include <parallel.h>
#include <carray.h>
#include <par_func.h>
/*

   abstract: main program

   history log:
      date         name             contents
      20170518     junghan kim      first written

*/
int main(int argc, char **argv) {

   int i, j, k, err, n;

   parallel *par = new parallel(argc, argv);
   par_func *pf  = new par_func(par);

   n = 2;
   carray *arr;
   int **val = arr->create(2,2);

   for (i=0;i<n;i++) {
      for (j=0;j<n;j++) {
         if (par->ismaster) {
            val[i][j] = 1;
         } else {
            val[i][j] = 0;
         }
      }
   }

   printf("before: rank %d: %d\n",par->rank,val[n-1][n-1]);
   cout << flush;
   err = par->sync();
   err = pf ->bcast(val[0],n*n,MPI_INT);
   err = par->sync();
   printf("after : rank %d: %d\n",par->rank,val[n-1][n-1]);
   cout << flush;

   delete arr, par, pf;

   return 0;

}




