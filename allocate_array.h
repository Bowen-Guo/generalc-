#ifndef ALLOCATE_ARRAY_H
#define ALLOCATE_ARRAY_H

#include <vector>

using namespace std;

float **** allocate_array(int n4,int n3,int n2,int n1);
float ***  allocate_array(int n3,int n2,int n1);
short ***  allocate_array_s(int n3,int n2,int n1);
int ***    allocate_array_i(int n3,int n2,int n1);
float **   allocate_array(int n2,int n1);
int **     allocate_array_i(int n2,int n1);
short **   allocate_array_s(int n2,int n1);
double **  allocate_array_d(int n2,int n1);
float *    allocate_array(int n1);
double *   allocate_array_d(int n1);
short *    allocate_array_s(int n1);
int *    allocate_array_i(int n1);

typedef vector<float> Dim1;
typedef vector<Dim1> Dim2;
typedef vector<Dim2> Dim3;
typedef vector<Dim3> Dim4;

typedef vector<int> Dimi1;
typedef vector<Dimi1> Dimi2;
typedef vector<Dimi2> Dimi3;
typedef vector<Dimi3> Dimi4;

typedef vector<short> Dims1;
typedef vector<Dims1> Dims2;
typedef vector<Dims2> Dims3;
typedef vector<Dims3> Dims4;

typedef long long int LLINT;

// Copy array
float *** copy(int n3,int n2,int n1, float *f);
short *** copy(int n3,int n2,int n1, short *f);
void copy(int n3,int n2,int n1, float * vin, float * vout);
void copy(int n2,int n1, float * vin, float * vout);

float ** copy(int n2,int n1, float *f);
int ** copy(int n2,int n1, int *f);
int * copy(int n1, int * vin);
float * copy(int n1, float * vin);

void copy_av(float **** x, Dim4  &y);
void copy_av(float ** x, Dim2 &y);
void copy_av(float * x,Dim1 &y);
void copy_av(int * x,Dimi1 &y);
float * copy_va(Dim1 x);
int * copy_va(Dimi1 x);

float ** copy_va(Dim2 x);
float *** copy_va(Dim3 x);
float *** copy_va_fill(Dim3 x, int n1);
short * copy_va(Dims1 x);
void arraycopy(short *src,int src_pos,short *dest,int dest_pos,int n);

// free array
void free4(float ****p);
void free3(double ***p);
void free3(float ***p);
void free3(int ***p);
void free3(short ***p);
void free2(double **p);
void free2(float **p);
void free2(int **p);
void free2(short **p);
void free1(double *p);
void free1(float *p);
void free1(int *p);
void free1(short *p);




#endif
