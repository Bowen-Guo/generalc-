#ifndef INOUT_H
#define INOUT_H

void readdata_1d (char *, int *x,    int n);
void readdata_1d (char *, float *x,  int n);
void readdata_1d (char *, double *x, int n);
void readdata_2d (char *f, float ** x, int n2, int n1);
void readdata_3d (char *, float ***x,int n3,int n2,int n1);
void readdata_3d (char *, int ***x,int n3,int n2,int n1);

void writedata_1d (char *, int *x,    int n);
void writedata_1d (char *, float *x,  int n);
void writedata_1d (char *, double *x, int n);
void writedata_2d (char *, float **x,int n2,int n1);
void writedata_2d (char *f, int ** x, int n2,int n1);
void writedata_3d (char *, float ***x,int nx,int ny,int nz);
void writedata_3d (char *f, int ***x, int nx,int ny, int nz);
#endif
