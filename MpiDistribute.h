#ifndef MPIDISTRIBUTE_H
#define MPIDISTRIBUTE_H

#include "allocate_array.h"

struct mfi{
   float v;
   int r;
};

class MpiDistribute{
public:
   MpiDistribute(int ns, int rank, int nsize);
   MpiDistribute(const MpiDistribute & md);
   MpiDistribute& operator=(const MpiDistribute & md);   
   int * getIndex();
   float * getArray(float * vin); 
   int * getArray(int * vin); 
   ~MpiDistribute();
   int getNjob(); 
   void minTPM(float *** t, float *** p,float *** pc, int *** m, int *** mc, int n3, int n2, int n1);   

private:
  int _rank, _nsize,_ns;
  int _ns_me; 
  int * _is_me;
};

void MPI_Bcast_Array(float *** x, int n3, int n2, int n1);
void MPI_Allreduce_MINLOC_Array(mfi * in, mfi * out, LLINT n);
void MPI_Reduce_Array(float *** v_s, float *** v_m, LLINT n);
void MPI_Reduce_Array(int *** v_s, int *** v_m, LLINT n);


#endif


