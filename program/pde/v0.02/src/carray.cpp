#include <carray.h>

carray::carray() {

}



carray::~carray() {

}



int *carray::create(int n) {

   int *ptr  = new int[n];
   return ptr;

}


//template <typename T>
int **carray::create(int n, int m) {

   int i;
   int **ptr  = new int*[n];
   int  *pool = new int[n*m];
   for (i=0;i<n;i++) {
     ptr[i] = &pool[i*m];
   }
   return ptr;

}


/*
template <typename T>
void delete_1d(T *ptr) {

   delete ptr;

}
*/


/*
template <typename T>
void delete_2d(T **ptr) {

   int i, n;
   n = sizeof(T)/sizeof(T[0]);
   for (i=0;i<n;i++) delete []ptr[i];
   delete []ptr;

}
*/
