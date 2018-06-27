#ifndef BLOCK_H
#define BLOCK_H

#include "mpi.h"

class block{

public:
  block(int nprocs, int rank, int n3, int n2, int n1);
  ~block();
  int * coords;
  int n1_b, n2_b, n3_b;
  int x2_b, x3_b;
  void mpiReadModel(char * fin, float *** vin);
  void mpiWriteModel(char * fout, float *** vout);
  void setCoords();
  void setBlockInfo();
  void divideModel();

private:
  int _ndims;
  int * _dims;
  int * _periods;
  int _nprocs, _rank;
  int _n1, _n2, _n3;// whole model size
  MPI_Comm COMM_CART;
};

#endif
