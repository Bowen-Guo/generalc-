#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "inout.h"
using namespace std;


//******************** write array into disk ******************
void writedata_1d (char *f, float *x, int n1) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *) x, n1 * sizeof(float));
   out.close();
}

void writedata_1d (char *f, int *x, int n) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *)x, n * sizeof(int));
   out.close();
}

void writedata_1d (char *f, double *x, int n) {
   ofstream out (f, ios::out|ios::binary);
   out.write ((char *)x, n * sizeof(double));
   out.close();
}

void writedata_2d (char *f, float ** x, int n2,int n1) {
   ofstream out (f, ios::out|ios::binary);
   for (int i2 = 0; i2 < n2; i2++){
         out.write ((char *) x[i2], n1 * sizeof(float));
      }
   out.close();
}

void writedata_2d (char *f, int ** x, int n2,int n1) {
   ofstream out (f, ios::out|ios::binary);
   for (int i2 = 0; i2 < n2; i2++){
         out.write ((char *) x[i2], n1 * sizeof(int));
      }
   out.close();
}


void writedata_3d (char *f, float ***x, int nx,int ny, int nz) {
   ofstream out (f, ios::out|ios::binary);
   for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy <ny; iy++){
         out.write ((char *) x[ix][iy], nz * sizeof(float));
      }
   }
   out.close();
}

void writedata_3d (char *f, int ***x, int nx,int ny, int nz) {
   ofstream out (f, ios::out|ios::binary);
   for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy <ny; iy++){
         out.write ((char *) x[ix][iy], nz * sizeof(int));
      }
   }
   out.close();
}

//******************** read array from disk ******************
void readdata_3d (char *f, float *** x, int nx,int ny,int nz) {
   ifstream in (f, ios::in|ios::binary);
    if(!in){
      cerr<<"Open error !" << endl;
      exit(1);
   }

   for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy <ny; iy++){
         in.read ((char *) x[ix][iy], nz * sizeof(float));
      }
   }
   in.close();
}

void readdata_3d (char *f, int *** x, int nx,int ny,int nz) {
   ifstream in (f, ios::in|ios::binary);
    if(!in){
      cerr<<"Open error !" << endl;
      exit(1);
   }

   for (int ix = 0; ix < nx; ix++){
      for (int iy = 0; iy <ny; iy++){
         in.read ((char *) x[ix][iy], nz * sizeof(int));
      }
   }
   in.close();
}

void readdata_2d (char *f, float ** x, int n2, int n1) {
   ifstream in (f, ios::in|ios::binary);
    if(!in){
      cerr<<"Open error !" << endl;
      exit(1);
   }

   for (int i2 = 0; i2 < n2; i2++)
     in.read ((char *) x[i2], n1 * sizeof(float));

   in.close();
}


void readdata_1d (char *f, int *x, int n) {
   ifstream in (f, ios::in|ios::binary);
   in.read ((char *)x, n * sizeof(int));
   in.close();
}

void readdata_1d (char *f, double *x, int n) {
   ifstream in (f, ios::in|ios::binary);
   in.read ((char *)x, n * sizeof(double));
   in.close();
}
