#include <iostream>
#include "stdlib.h"
#include "allocate_array.h"
#include <time.h>
#include "math.h"
#include "Const.h"
//#include "omp.h"
#include <cstring>
#include "assert.h"
#include "float.h"

using namespace std;


// Sum the 3D array
float sum(float *** f, int n3, int n2, int n1){
   double sss = 0.0;
   #pragma omp parallel for
   for (int i=0;i<n3;i++){
      for (int j=0;j<n2;j++){
         for (int k=0;k<n1;k++){
            sss += f[i][j][k];
         }
      }
   }
   return (float) sss;
}

double sum_d(float *** f, int n3, int n2, int n1){
   double sss = 0.0;
   #pragma omp parallel for
   for (int i=0;i<n3;i++){
      for (int j=0;j<n2;j++){
         for (int k=0;k<n1;k++){
            sss += f[i][j][k];
         }
      }
   }
   return sss;
}

float * randfloat(int n1){
   srand(time(0));
   float *r = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0; i<n1; i++)
      r[i] = float(rand()%10000)/9999.0;
   return r;
}

float *** randfloat(int n3, int n2, int n1){
   srand(time(0));
   float ***r = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for (int i=0; i<n3; i++){
      for (int j=0; j<n2; j++){
         for (int k=0; k<n1; k++)
            r[i][j][k] = float(rand()%10000)/9999.0;
      }
   }
   return r;
}

// Real number multiple an array
//1D
float * mul(int c, float * f1,int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0;i<n1;i++)
      out[i] = c * f1[i];
   return out;
}

float * mul(float c, float * f1,int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0;i<n1;i++)
      out[i] = c * f1[i];
   return out;
}

void mul(float c, float * f1, float *f2, int n1){
   #pragma omp parallel for
   for (int i=0;i<n1;i++)
      f2[i] = c * f1[i];
}


float * mul(float *f1, float * f2,int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0;i<n1;i++)
      out[i] = f2[i] * f1[i];
   return out;
}

void mul(float c, Dim1 &f1, Dim1 &f2){
   int n1 = f1.size();
   #pragma omp parallel for
   for (int j=0;j<n1;j++)
      f2[j] = c * f1[j];
}


// 2D
void mul(float c, float ** f1, float ** f2, int n2, int n1){
   #pragma omp parallel for
   for (int i=0;i<n2;i++){
      for (int j=0;j<n1;j++)
         f2[i][j] = c * f1[i][j];
   }
}

void mul(float c, Dim2 &f1, Dim2 &f2){
   int n2 = f1.size();
   int n1 = f1[0].size();

   #pragma omp parallel for
   for (int i=0;i<n2;i++){
      for (int j=0;j<n1;j++)
         f2[i][j] = c * f1[i][j];
   }
}

//3D
float *** mul(float f, float *** f1, int n3, int n2, int n1){
   float *** out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            out[i][j][k] = f1[i][j][k] * f  ;
      }
   }
return out;
}

float *** mul(float *** f1,float f, int n3, int n2, int n1){
   float *** out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            out[i][j][k] = f1[i][j][k] * f  ;
      }
   }
return out;
}


void mul(float *** vin1, float *** vin2, float *** vout,int n3, int n2, int n1){
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            vout[i][j][k] = vin1[i][j][k] * vin2[i][j][k];
      }
   }
}

float *** mul(float *** vin1, float *** vin2, int n3, int n2, int n1){
   float *** vout = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            vout[i][j][k] = vin1[i][j][k] * vin2[i][j][k];
      }
   }
   return vout;
}

// div

float * div(float c, float * f1,int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0;i<n1;i++)
      out[i] = c / f1[i];
   return out;
}

float *** div(float *** f1, float *** f2, int n3, int n2, int n1){

   float ***out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for (int i=0;i<n3;i++){
      for (int j=0;j<n2;j++){
         for (int k=0;k<n1;k++)
            out[i][j][k] = f1[i][j][k]/ f2[i][j][k];
      }
   }
   return out;
}

void div(float a, float *** vin, float *** vout, int n3, int n2, int n1){

   #pragma omp parallel for
   for (int i=0;i<n3;i++){
      for (int j=0;j<n2;j++){
         for (int k=0;k<n1;k++)
            vout[i][j][k] = a/vin[i][j][k];
      }
   }
}

// neg
float * neg(float * f1, int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for(int i =0; i<n1; i++)
      out[i] = -f1[i];
return out;
}

float *** neg(float *** f1, int n3, int n2, int n1){
   float ***out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for (int i=0;i<n3;i++){
      for (int j=0;j<n2;j++){
         for (int k=0;k<n1;k++)
            out[i][j][k] = -f1[i][j][k];
      }
   }
return out;
}


// add
void add(float *** f1, float f2, float *** fout, int n3, int n2, int n1){
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            fout[i][j][k] = f1[i][j][k] + f2  ;
      }
   }
}

float *** add(float *** f1, float *** f2, int n3, int n2, int n1){
   float ***out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            out[i][j][k] = f1[i][j][k] + f2[i][j][k]  ;
      }
   }
return out;
}


// f1-f2

float *** sub(float *** f1, float  f2, int n3, int n2, int n1){
   float ***out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            out[i][j][k] = (double) f1[i][j][k] - (double) f2   ;
      }
   }
return out;
}

float *** sub(float *** f1, double  f2, int n3, int n2, int n1){
   float ***out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for(int i =0; i<n3; i++){
      for(int j =0; j<n2; j++){
         for(int k =0; k<n1; k++)
            out[i][j][k] = f1[i][j][k] - f2  ;
      }
   }
return out;
}


float * sub(float * f1, float * f2, int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for(int i =0; i<n1; i++)
      out[i] = f1[i] - f2[i];
return out;
}

// f1-f2
float * sub(float * f1, float  f2, int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for(int i =0; i<n1; i++)
      out[i] = f1[i] - f2;
return out;
}


// f1^f2
float * pow(float * f1, float * f2, int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for(int i =0; i<n1; i++)
      out[i] = pow(f1[i],  f2[i]);
return out;
}

// f1^f2
float * pow(float * f1, float  f2, int n1){
   float *out = allocate_array(n1);
   #pragma omp parallel for
   for(int i =0; i<n1; i++)
      out[i] = pow(f1[i],  f2);
return out;
}

// clip ( rmin< r <rmax)
float * clip(float rxmin, float rxmax, float * rx,int n1) {
   float * out = allocate_array(n1);
   #pragma omp parallel for
   for (int i=0; i<n1; i++){
      out[i] = (rx[i]<rxmin)?rxmin:(rx[i]>rxmax)?rxmax:rx[i];
  }
  return out;
}

//  fill the 3D array with a number
float *** fill(float f1, int n3, int n2, int n1){
   float *** out = allocate_array(n3,n2,n1);
   #pragma omp parallel for
   for (int i = 0; i<n3; i++){
      for (int j = 0; j<n2; j++){
         for (int k = 0; k<n1; k++)
            out[i][j][k] = f1;
      }
   }
   return out;
}

int *** fill(int f1, int n3, int n2, int n1){
   int *** out = allocate_array_i(n3,n2,n1);
   #pragma omp parallel for
   for (int i = 0; i<n3; i++){
      for (int j = 0; j<n2; j++){
         for (int k = 0; k<n1; k++)
            out[i][j][k] = f1;
      }
   }
   return out;
}

//  fill the 2D array with a number
float ** fill(float f1, int n2, int n1){
   float ** out = allocate_array(n2,n1);
   #pragma omp parallel for
     for (int j = 0; j<n2; j++){
        for (int k = 0; k<n1; k++)
           out[j][k] = f1;
     }
   return out;
}

// fillfloat fill the array with a float number
float * fill(float f1, int n){
   float *out = allocate_array(n);
   #pragma omp parallel for
   for (int i = 0; i<n; i++){
            out[i] = f1;
   }
   return out;
}

// fillfloat fill the array with a float number
int * fill(int f1, int n){
   int * out = allocate_array_i(n);
   #pragma omp parallel for
   for (int i = 0; i<n; i++){
            out[i] = f1;
   }
   return out;
}

void fill(short s1, short *** sss, int n3, int n2, int n1){
   #pragma omp parallel for
   for (int i = 0; i<n3; i++){
      for (int j = 0; j<n2; j++){
         for (int k = 0; k<n1; k++)
            sss[i][j][k] = s1;
      }
   }
}

void fill(float s1, float *** sss, int n3, int n2, int n1){
   #pragma omp parallel for
   for (int i = 0; i<n3; i++){
      for (int j = 0; j<n2; j++){
         for (int k = 0; k<n1; k++)
            sss[i][j][k] = s1;
      }
   }
}

void fill(float s1, float ** ss, int n2, int n1){
#pragma omp parallel for
  for (int i = 0; i<n2; i++){
     for (int j = 0; j<n1; j++)
        ss[i][j] = s1;
  }
}




// Performs a binary search in a monotonic array of values
int binarySearch(float * a, float x, int i, int n) {
   int nm1 = n-1;
   int low = 0;
   int high = nm1;
   bool increasing = n<2 || a[0]<a[1];
   if (i<n) {
     high = (0<=i)?i:-(i+1);
     low = high-1;
     int step = 1;
     if (increasing) {
       for (; 0<low && x<a[low]; low-=step,step+=step)
         high = low;
       for (; high<nm1 && a[high]<x; high+=step,step+=step)
         low = high;
     } else {
       for (; 0<low && x>a[low]; low-=step,step+=step)
         high = low;
       for (; high<nm1 && a[high]>x; high+=step,step+=step)
         low = high;
     }
     if (low<0) low = 0;
     if (high>nm1) high = nm1;
   }
   if (increasing) {
     while (low<=high) {
       int mid = (low+high)>>1;
       float amid = a[mid];
       if (amid<x)
         low = mid+1;
       else if (amid>x)
         high = mid-1;
       else
         return mid;
     }
   } else {
     while (low<=high) {
       int mid = (low+high)>>1;
       float amid = a[mid];
       if (amid>x)
         low = mid+1;
       else if (amid<x)
         high = mid-1;
       else
         return mid;
     }
   }
   return -(low+1);
}

int binarySearch(double * a, double x, int i, int n) {
   int nm1 = n-1;
   int low = 0;
   int high = nm1;
   bool increasing = n<2 || a[0]<a[1];
   if (i<n) {
     high = (0<=i)?i:-(i+1);
     low = high-1;
     int step = 1;
     if (increasing) {
       for (; 0<low && x<a[low]; low-=step,step+=step)
         high = low;
       for (; high<nm1 && a[high]<x; high+=step,step+=step)
         low = high;
     } else {
       for (; 0<low && x>a[low]; low-=step,step+=step)
         high = low;
       for (; high<nm1 && a[high]>x; high+=step,step+=step)
         low = high;
     }
     if (low<0) low = 0;
     if (high>nm1) high = nm1;
   }
   if (increasing) {
     while (low<=high) {
       int mid = (low+high)>>1;
       float amid = a[mid];
       if (amid<x)
         low = mid+1;
       else if (amid>x)
         high = mid-1;
       else
         return mid;
     }
   } else {
     while (low<=high) {
       int mid = (low+high)>>1;
       float amid = a[mid];
       if (amid>x)
         low = mid+1;
       else if (amid<x)
         high = mid-1;
       else
         return mid;
     }
   }
   return -(low+1);
}



float min(float v1, float v2){

if (v1>v2)
   v1 = v2;
   return v1;
}

float max(float v1, float v2){

if (v1<v2)
   v1 = v2;
   return v1;
}

float max(float *p, int n){
   float v = -FLT_MAX;
   for (int i=0;i<n;i++){
      if (p[i]>v) v = p[i];
   }
   return v;
}

float min(float *p, int n){
   float v = FLT_MAX;
   for (int i=0;i<n;i++){
      if (p[i] < v) v = p[i];
   }
   return v;
}

double toRadians(double angdeg){
   angdeg = angdeg/180.0*PI;
   return angdeg;
}

float toDegrees(float angdeg){
   angdeg = angdeg/PI*180.0;
   return angdeg;
}

void shuffle(float * f, int n) {
  srand(314159); // constant seed for consistency
  float ii;
  for (int i=n-1; i>0; --i) {
    int j = rand() % (i+1);
    //ii = i1[i]; i1[i] = i1[j]; i1[j] = ii;
    //ii = i2[i]; i2[i] = i2[j]; i2[j] = ii;
    //ii = i3[i]; i3[i] = i3[j]; i3[j] = ii;
    ii = f[i];  f[i] = f[j];  f[j] = ii;
  }
}

void shuffle(int * f, int n) {
  srand(314159); // constant seed for consistency
  int ii;
  for (int i=n-1; i>0; --i) {
    int j = rand() % (i+1);
    //ii = i1[i]; i1[i] = i1[j]; i1[j] = ii;
    //ii = i2[i]; i2[i] = i2[j]; i2[j] = ii;
    //ii = i3[i]; i3[i] = i3[j]; i3[j] = ii;
    ii = f[i];  f[i] = f[j];  f[j] = ii;
  }
}

void clip(float *** p, float cmin, float cmax, int n3, int n2, int n1){
  assert(cmin<=cmax);
  #pragma omp parallel for
  for (int i3=0;i3<n3;i3++){
     for (int i2=0;i2<n2;i2++){
        for (int i1=0;i1<n1;i1++){
           if (p[i3][i2][i1]>cmax) p[i3][i2][i1] = cmax;
           if (p[i3][i2][i1]<cmin) p[i3][i2][i1] = cmin;
        }
     }
  }
}
