#ifndef ARRAYMATH_H
#define ARRAYMATH_H

#include "allocate_array.h"

float * randfloat(int n1);
float *** randfloat(int n3,int n2, int n1);

// add
void add(float *** f1, float f2, float *** fout, int n3, int n2, int n1);
float *** add(float *** f1,float  *** f2, int n3, int n2, int n1);
float sum(float *** f, int n3, int n2, int n1);
double sum_d(float *** f, int n3, int n2, int n1);

//multiply
float *** mul(float *** vin1, float *** vin2, int n3, int n2, int n1);
float *** mul(float f, float *** f1,int n3, int n2, int n1);
float *** mul(float *** f1,float f, int n3, int n2, int n1);
void *** mul(float *** vin1, float *** vin2, float *** vout,int n3, int n2, int n1);
void  mul(float c, float ** f1 ,float **f2, int n2, int n1);
void  mul(float c, Dim2 &f1 , Dim2 &f2);
void  mul(float c, Dim1 &f1 , Dim1 &f2);
float * mul(int c, float * f1,int n1);
float * mul(float c, float * f1,int n1);
void * mul(float c, float * f1 ,float *f2, int n1);
float * mul(float *f1, float * f2,int n1);

// divide
float * div(float c, float * f1,int n1);
float *** div(float ***f1, float *** f2,int n3, int n2, int n1);
void div(float a, float *** vin, float *** vout, int n3, int n2, int n1);

// others
float * clip(float rxmin, float rxmax, float * rx,int n1);
float * neg(float * f1, int n1);
float *** neg(float *** f1, int n3, int n2, int n1);
float * sub(float * f1,float * f2, int n1);
float * sub(float * f1,float  f2, int n1);
float *** sub(float *** f1,float  f2, int n3, int n2, int n1);
float *** sub(float *** f1,double  f2, int n3, int n2, int n1);
float * pow(float * f1, float * f2,int n1);
float * pow(float * f1, float  f2,int n1);

float *** fill(float f1, int n3, int n2, int n1);
int *** fill(int f1, int n3, int n2, int n1);
float ** fill(float f1, int n2, int n1);


float * fill(float f1, int n);
int * fill(int f1, int n);


void fill(short s, short *** sss, int n3, int n2, int n1);
void fill(float s, float *** sss, int n3, int n2, int n1);
void fill(float s1, float ** ss, int n2, int n1);

int binarySearch(float * a, float x, int i, int n);
int binarySearch(double * a, double x, int i, int n);
float min(float v1, float v2);
float max(float v1, float v2);

float max(float *p, int n);
float min(float *p, int n);

double toRadians (double x);
float toDegrees (float x);
void shuffle(float *f, int n);
void shuffle(int *f, int n);
void clip(float *** p, float cmin, float cmax, int n3, int n2, int n1);

#endif
