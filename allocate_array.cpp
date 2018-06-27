/* Allocate 3D arrays */
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "allocate_array.h"
#include <stdio.h>
//#include "omp.h"
using namespace std;


//******************** allocate array ******************


//******************** 4D array ******************
// 4D array float
float **** allocate_array(int n4,int n3,int n2,int n1){

   // Allocate array with continous memory
   float **** xxxx = new float *** [n4];
   float *** xxx = new float ** [(LLINT) n4*n3];
   float ** xx = new float * [(LLINT) n4*n3*n2];
   float * x = new float [(LLINT) n4*n3*n2*n1]();
   //#pragma omp parallel for
   for (int i=0;i<n4*n3*n2;i++)
      xx[i] = &x[(LLINT) i*n1];
   //#pragma omp parallel for
   for (int i=0;i<n4*n3;i++)
      xxx[i] = &xx[(LLINT) i*n2];
   //#pragma omp parallel for
   for (int i=0;i<n4;i++)
      xxxx[i] = &xxx[(LLINT) i*n3];

   return xxxx;

}



//******************** 3D array ******************
// 3D array integer
int *** allocate_array_i(int n3,int n2,int n1){
   // Allocate array with continous memory
   int *** xxx = new int ** [n3];
   int ** xx = new int * [(LLINT) n3*n2];
   int * x = new int [(LLINT) n3*n2*n1]();

   //#pragma omp parallel for
   for (int i=0;i<n3*n2;i++)
      xx[i] = &x[(LLINT) i*n1];
   //#pragma omp parallel for
   for (int i=0;i<n3;i++)
      xxx[i] = &xx[(LLINT) i*n2];
   return xxx;
}

// 3D array short
short *** allocate_array_s(int n3,int n2,int n1){

   // Allocate array with continous memory
   short *** xxx = new short ** [n3];
   short ** xx = new short * [(LLINT) n3*n2];
   short * x = new short [(LLINT) n3*n2*n1]();

   //#pragma omp parallel for
   for (int i=0;i<n3*n2;i++)
      xx[i] = &x[(LLINT) i*n1];
   //#pragma omp parallel for
   for (int i=0;i<n3;i++)
      xxx[i] = &xx[(LLINT) i*n2];
   return xxx;
}

// 3D array float
float *** allocate_array(int n3,int n2,int n1){

   // Allocate array with continous memory
   float *** xxx = new float ** [n3];
   float ** xx = new float * [(LLINT) n3*n2];
   float * x = new float [(LLINT) n3*n2*n1]();

   //#pragma omp parallel for
   for (int i=0;i<n3*n2;i++)
      xx[i] = &x[(LLINT) i*n1];
   //#pragma omp parallel for
   for (int i=0;i<n3;i++)
      xxx[i] = &xx[(LLINT) i*n2];
   return xxx;
}

//******************** 2D array ******************
// 2D array short
short ** allocate_array_s(int n2,int n1){

   // Allocate array with continous memory
   short ** xx = new short *[n2];
   xx[0] = new short[(LLINT) n2*n1]();
   //#pragma omp parallel for
   for (int i=1;i<n2;i++)
      xx[i] = xx[i-1]+n1;
   return xx;
}

// 2D array integer
int ** allocate_array_i(int n2,int n1){

   // Allocate array with continous memory
   int ** xx = new int *[n2];
   xx[0] = new int[(LLINT) n2*n1]();
   //#pragma omp parallel for
   for (int i=1;i<n2;i++)
      xx[i] = xx[i-1]+n1;
   return xx;
}

// 2D array float
float ** allocate_array(int n2,int n1){

   // Allocate array with continous memory
   float ** xx = new float *[n2];
   xx[0] = new float[(LLINT) n2*n1]();
   //#pragma omp parallel for
   for (int i=1;i<n2;i++)
      xx[i] = xx[i-1]+n1;
   return xx;
}

// 2D array double
double ** allocate_array_d(int n2,int n1){

   // Allocate array with continous memory
   double ** xx = new double *[n2];
   xx[0] = new double[(LLINT) n2*n1]();
   //#pragma omp parallel for
   for (int i=1;i<n2;i++)
      xx[i] = xx[i-1]+n1;
   return xx;
}


//******************** 1D array ******************
// 1D array float
float * allocate_array(int n1){

   float * array = new float[n1]();
   return array;
}

// 1D array double
double * allocate_array_d(int n1){
   double * array = new double[n1]();
   return array;
}

// 1D array short
short * allocate_array_s(int n1){
   short * array = new short[n1]();
   return array;
}

// 1D array int
int * allocate_array_i(int n1){
   int * array = new int[n1]();
   return array;
}


//******************** copy array ******************


//******************** copy 3D array ******************
short *** copy(int n3,int n2,int n1, short * vin){

   short *** vout = allocate_array_s(n3,n2,n1);
   memcpy(vout[0][0], vin, (LLINT) n3*n2*n1*sizeof(short));
   return vout;
}

float *** copy(int n3,int n2,int n1, float * vin){

   float *** vout = allocate_array(n3,n2,n1);
   memcpy(vout[0][0], vin, (LLINT) n3*n2*n1*sizeof(float));
   return vout;
}

void copy(int n3,int n2,int n1, float * vin, float * vout){
   memcpy(vout, vin, (LLINT) n3*n2*n1*sizeof(float));
}

//******************** copy 2D array ******************
void copy(int n2,int n1, float * vin, float * vout){
   memcpy(vout, vin, (LLINT) n2*n1*sizeof(float));
}

float ** copy(int n2,int n1, float * vin){

   float ** vout = allocate_array(n2,n1);
   memcpy(vout[0], vin, (LLINT) n2*n1*sizeof(float));
   return vout;
}

int ** copy(int n2,int n1, int * vin){

   int ** vout = allocate_array_i(n2,n1);
   memcpy(vout[0], vin, (LLINT) n2*n1*sizeof(int));
   return vout;
}

//******************** copy 2D array ******************

int * copy(int n1, int * vin){

   int *vout = allocate_array_i(n1);
   memcpy(vout, vin, n1*sizeof(int));
   return vout;
}

float * copy(int n1, float * vin){

   float *vout = allocate_array(n1);
   memcpy(vout, vin, n1*sizeof(float));
   return vout;
}



//******************** copy 3D vector into array ******************
float *** copy_va(Dim3 x)
{
   int n1 = x[0][0].size();
   int n2 = x[0].size();
   int n3 = x.size();

   float *** y = allocate_array(n3,n2,n1);

   //#pragma omp parallel for
   for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
         for (int i1=0; i1<n1; i1++){
            y[i3][i2][i1] = x[i3][i2][i1];
         }
      }
   }
   return y;
}

//******************** copy 2D vector into array ******************
float ** copy_va(Dim2 x)
{
   int n1 = x[0].size();
   int n2 = x.size();

   float ** y = allocate_array(n2,n1);

   //#pragma omp parallel for
   for (int i2=0; i2<n2; i2++){
      for (int i1=0; i1<n1; i1++){
         y[i2][i1] = x[i2][i1];
      }
   }
   return y;
}

//******************** copy 1D vector into short array ******************
short * copy_va(Dims1 x)
{
  int n = x.size();
  short *y = allocate_array_s(n);
  //#pragma omp parallel for
  for (int i=0;i<n;i++)
     y[i] = x[i];
  return y;
}

int * copy_va(Dimi1 x)
{
  int n = x.size();
  int *y = allocate_array_i(n);
  //#pragma omp parallel for
  for (int i=0;i<n;i++)
     y[i] = x[i];
  return y;
}
//******************** copy 1D vector into float array ******************
float * copy_va(Dim1 x)
{
   int n1 = x.size();

   float * y = allocate_array(n1);

   //#pragma omp parallel for
   for (int i1=0; i1<n1; i1++){
      y[i1] = x[i1];
   }
   return y;
}


//******************** copy 4D array into vector ******************
void copy_av(float **** x, Dim4 &y)
{
   int n1 = y[0][0][0].size();
   int n2 = y[0][0].size();
   int n3 = y[0].size();
   int n4 = y.size();

   //#pragma omp parallel for
   for (int i4=0; i4<n4; i4++){
      for (int i3=0; i3<n3; i3++){
         for (int i2=0; i2<n2; i2++){
            for (int i1=0; i1<n1; i1++){
               y[i4][i3][i2][i1] = x[i4][i3][i2][i1];
            }
         }
      }
   }
}

//******************** copy 2D array into vector ******************
void copy_av(float ** x,Dim2 &y)
{
   int n1 = y[0].size();
   int n2 = y.size();

   //#pragma omp parallel for
   for (int i2=0; i2<n2; i2++){
      for (int i1=0; i1<n1; i1++){
         y[i2][i1] = x[i2][i1];
      }
   }
}

//******************** copy 1D array into vector ******************
void copy_av(float * x,Dim1 &y)
{
   int n = y.size();
   //#pragma omp parallel for
   for (int i=0; i<n; i++)
     y[i] = x[i];
}

void copy_av(int * x,Dimi1 &y)
{
   int n = y.size();
   //#pragma omp parallel for
   for (int i=0; i<n; i++)
     y[i] = x[i];
}


//3D copy vector to array
float *** copy_va_fill(Dim3 x, int n1)
{
   //int n1 = x[0][0].size();
   int n2 = x[0].size();
   int n3 = x.size();

   float *** y = allocate_array(n3,n2,n1);

   #pragma omp parallel for
   for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
         if (x[i3][i2].empty()){
            for (int i1=0; i1<n1; i1++)
               y[i3][i2][i1] = 0;
         }
         else{
            for (int i1=0; i1<n1; i1++)
               y[i3][i2][i1] = x[i3][i2][i1];
         }
      }
   }
   return y;
}

// Copy shiftedarray
void arraycopy(short *src,int src_pos,short *dest,int dest_pos,int n){
   for (int i=0;i<src_pos;i++)
      src ++;
   for (int i=0;i<dest_pos;i++)
      dest ++;
   memcpy(dest,src,n*sizeof(short));
}


//******************** free array******************

// 4D array ***************************************
// float
void free4(float ****p)
{
   if (p!=nullptr) {
      delete [] p[0][0][0];
      delete [] p[0][0];
      delete [] p[0];
      delete [] p;
   }
}

// 3D array ***************************************
// float
void free3(float ***p)
{
   if (p!=nullptr) {
      delete [] p[0][0];
      delete [] p[0];
      delete [] p;
   }
}

// double
void free3(double ***p)
{
   if (p!=nullptr) {
      delete [] p[0][0];
      delete [] p[0];
      delete [] p;
   }
}

// integer
void free3(int ***p)
{
   if (p!=nullptr) {
      delete [] p[0][0];
      delete [] p[0];
      delete [] p;
   }
}

// short
void free3(short ***p)
{
   //if (p!=nullptr) {
      delete [] p[0][0];
      delete [] p[0];
      delete [] p;
   //}
}

// 2D array ***************************************
// float
void free2(float **p)
{
   if (p!=nullptr) {
      delete [] p[0];
      delete [] p;
   }
}

// double
void free2(double **p)
{
   if (p!=nullptr) {
      delete [] p[0];
      delete [] p;
   }
}

// integer
void free2(int **p)
{
   if (p!=nullptr) {
      delete [] p[0];
      delete [] p;
   }
}

// short
void free2(short **p)
{
   if (p!=nullptr) {
      delete [] p[0];
      delete [] p;
   }
}

// 1D array ***************************************
// float
void free1(float *p)
{
   if (p!=nullptr) delete [] p;
}

// double
void free1(double *p)
{
   if (p!=nullptr) delete [] p;
}

// integer
void free1(int *p)
{
   if (p!=nullptr) delete [] p;
}

// short
void free1(short *p)
{
   if (p!=nullptr) delete [] p;
}
