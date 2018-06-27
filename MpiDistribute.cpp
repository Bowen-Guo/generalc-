#include <iostream>
#include "mpi.h"
#include "math.h"
#include "MpiDistribute.h" 
#include "omp.h"
#include "limits.h"
using namespace std;


// ns: number of jobs in total
// ns_me: number of jobs at this thread
// is_me: integer array of (ns_me), which stores the index in the total job
MpiDistribute::MpiDistribute(int ns, int rank, int nsize){
   _ns = ns;
   _rank = rank;
   _nsize = nsize;
   
   // ns less than node number
   if (ns<nsize){
      if (rank == 0) 
         cout << "WARNING: # of data points "<< ns << " is less than # of nodes"<< nsize << endl;
      if (rank < ns){
         _ns_me=1;  
         _is_me = allocate_array_i(_ns_me);
         _is_me[0] = rank;
      }
      else{
         _ns_me=0;
         _is_me = nullptr;
      }
      return;
   }
 
   // one node case
   if (nsize == 1){
     _ns_me=ns;
     _is_me =allocate_array_i(_ns_me);
     for (int is=0; is<_ns_me; is++)
        _is_me[is]=is;
      return;
   }
   
   // Otherwise
   int ns_node = floor(ns/nsize); 
   int ns_left = ns - ns_node * nsize;
   
   //if ( rank < (nsize-ns_left) ) {
   //   _ns_me = ns_node;
   //   _is_me = allocate_array_i(_ns_me);
   //   for (int is=0;is<_ns_me;is++)
   //      _is_me[is] = is * nsize + rank;
   //}
   //else{
   //   _ns_me = ns_node + 1;
   //   _is_me = allocate_array_i(_ns_me);
   //   for (int is=0;is<_ns_me-1;is++)
   //      _is_me[is] = is * nsize + rank;
   //   _is_me[_ns_me-1] = (_ns_me-1)*nsize+rank-(nsize-ns_left);
   //}

   if (rank < ns_left){
      _ns_me = ns_node + 1;
      _is_me = allocate_array_i(_ns_me);
      for (int i=0; i<_ns_me; i++)
         _is_me[i]=  rank*_ns_me + i;
   }
   else{
      _ns_me = ns_node;
      _is_me = allocate_array_i(_ns_me);
      for (int i=0; i<_ns_me; i++)
         _is_me[i]=  ns_left + rank * ns_node + i;
   }
     
}

MpiDistribute::MpiDistribute(const MpiDistribute & md){
  _rank = md._rank;
  _nsize = md._nsize;
  _ns = md._ns;
  _ns_me = md._ns_me;
  _is_me = copy(_ns_me,md._is_me);
}

MpiDistribute& MpiDistribute::operator=(const MpiDistribute & md){
  _rank = md._rank;
  _nsize = md._nsize;
  _ns = md._ns;
  _ns_me = md._ns_me;
  _is_me = copy(_ns_me,md._is_me);
}


int * MpiDistribute::getIndex(){
   return _is_me;
}

int MpiDistribute::getNjob(){
   return _ns_me;
}

float * MpiDistribute::getArray(float * vin){
   float * vout = allocate_array(_ns_me);
   for (int i=0;i<_ns_me;i++)
      vout[i] = vin[_is_me[i]];
   return vout;
}

int * MpiDistribute::getArray(int * vin){
   int * vout = allocate_array_i(_ns_me);
   for (int i=0;i<_ns_me;i++)
      vout[i] = vin[_is_me[i]];
   return vout;
}


void MpiDistribute::minTPM(float *** t, float *** p, float *** pc, int *** m, int *** mc,
                          int n3, int n2, int n1){

   //int nn = n3*n2*n1;
   LLINT nn = (LLINT) n3*n2*n1;

   //struct mfi{
   //  float v; // value
   //  int r;   //rank
   //};
   struct mfi *in  = new struct mfi [nn];
   struct mfi *out  = new struct mfi [nn];

   #pragma omp parallel for
   for (int i=0;i<n3;i++){ 
      for (int j=0;j<n2;j++){ 
         for (int k=0;k<n1;k++){
            LLINT ind = (LLINT) i*n1*n2+j*n1+k;
            in[ind].v = t[i][j][k];
            in[ind].r = _rank;
         }
      }
   }
 
   //MPI_Barrier(MPI_COMM_WORLD);
   
   // Collect the minimum value to the master node 
   //MPI_Allreduce(in, out, nn, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);
   MPI_Allreduce_MINLOC_Array(in, out, nn);
   
   delete [] in; 
   
   if (_rank == 0){
      #pragma omp parallel for
      for (int i=0;i<n3;i++){ 
         for (int j=0;j<n2;j++){ 
            for (int k=0;k<n1;k++){
               LLINT ind = (LLINT) i*n1*n2+j*n1+k;
               t[i][j][k] = out[ind].v;
            }
         }
      }
   }
  
   //float *** pc = allocate_array(n3,n2,n1); // copy of p
   //int *** mc = allocate_array_i(n3,n2,n1);

   #pragma omp parallel for
   for (int i=0;i<n3;i++){ 
      for (int j=0;j<n2;j++){ 
         for (int k=0;k<n1;k++){
            LLINT ind = (LLINT) i*n1*n2+j*n1+k;
            int indp = out[ind].r;
            if (_rank == indp){ // this rank has the min t at [i,j,k]
               pc[i][j][k] = p[i][j][k];  // Record this p value 
               mc[i][j][k] = m[i][j][k];  // Record this p value 
            } 
         } 
      }
   }
   
   delete [] out; 

   if (_rank == 0){
      #pragma omp parallel for
      for (int i=0;i<n3;i++){ 
         for (int j=0;j<n2;j++){ 
            for (int k=0;k<n1;k++){
               p[i][j][k] = 0.0;
               m[i][j][k] = 0;
            }
         }
      }
   } 
  
   //MPI_Barrier(MPI_COMM_WORLD);
   
   //MPI_Reduce(pc[0][0], p[0][0], nn, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce_Array(pc,p,nn);
   MPI_Reduce_Array(mc,m,nn); 
}

void MPI_Bcast_Array(float *** x, int n3, int n2, int n1){
   LLINT n321 = (LLINT) n3* n2 *n1;
   LLINT sendsum = n321;
   int sendlimit = floor(INT_MAX/sizeof(float));
   float * ptr = & x[0][0][0];
   
   while (sendsum>sendlimit){
      MPI_Bcast(ptr, sendlimit, MPI_FLOAT, 0, MPI_COMM_WORLD);
      sendsum -= sendlimit;
      ptr += sendlimit;
      MPI_Barrier(MPI_COMM_WORLD);
   }
   MPI_Bcast(ptr, sendsum, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void MPI_Allreduce_MINLOC_Array(mfi * in, mfi * out, LLINT n){

   LLINT sendsum = n;
   int sendlimit = floor(INT_MAX/sizeof(mfi));
   mfi * ptr_in = in;
   mfi * ptr_out = out;
   while (sendsum>sendlimit){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(ptr_in, ptr_out, sendlimit, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);
      sendsum -= sendlimit;
      ptr_in += sendlimit;
      ptr_out += sendlimit;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(ptr_in, ptr_out, sendsum, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);
}

void MPI_Reduce_Array(float *** v_s, float *** v_m, LLINT n){

   // v_s is at the slave node
   // v_m is at the master node
   LLINT sendsum = n;
   int sendlimit = floor(INT_MAX/sizeof(float));
   float * ptr_s = v_s[0][0];
   float * ptr_m = v_m[0][0];
   while (sendsum>sendlimit){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(ptr_s, ptr_m, sendlimit, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      sendsum -= sendlimit;
      ptr_s += sendlimit;
      ptr_m += sendlimit;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Reduce(ptr_s, ptr_m, sendsum, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
}

void MPI_Reduce_Array(int *** v_s, int *** v_m, LLINT n){

   // v_s is at the slave node
   // v_m is at the master node
   LLINT sendsum = n;
   int sendlimit = floor(INT_MAX/sizeof(int));
   int * ptr_s = v_s[0][0];
   int * ptr_m = v_m[0][0];
   while (sendsum>sendlimit){
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(ptr_s, ptr_m, sendlimit, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      sendsum -= sendlimit;
      ptr_s += sendlimit;
      ptr_m += sendlimit;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Reduce(ptr_s, ptr_m, sendsum, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

MpiDistribute::~MpiDistribute(){
   free1(_is_me);
}



