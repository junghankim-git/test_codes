#include <stdio.h>
#include <iostream>
using namespace std;
#include <par_func.h>
/*

   abstract: class for mpi functions.

   history log:
     date         name             contents
     20170518     junghan kim      first written

*/
par_func::par_func() {

   initialized = false;

}



par_func::~par_func() {

   initialized = false;
}



par_func::par_func(parallel *par) {

   int err;
   err = 0;
   //this->par = new parallel(par->comm);
   this->par = par;
   err = this->par->check_init();
   initialized = true;

}



int par_func::bcast(void *buf, int cnt, int vtype, int root, int comm) {

   int err;
   if (root==-1) root = this->par->root;
   if (comm==-1) comm = this->par->comm;
   
   err = MPI_Bcast(buf,cnt,vtype,root,comm);

}
