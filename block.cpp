

#include <cmath>
#include <fstream>
#include "mpi.h"
#include "block.h"

// head file for debug
#include <iostream>

using namespace std;

block::block(int nprocs, int rank, int n3, int n2, int n1){
  _nprocs = nprocs;
  _rank = rank;
  _n1 = n1;
  _n2 = n2;
  _n3 = n3;
  _ndims = 2;
}

/*
 Compute the number of segments in each dimension
 */
void block::divideModel(){
  _dims = new int[_ndims]();
  int n2_d = int (sqrt((float) _nprocs));
  float n3_d = _nprocs/n2_d;

  while ((floor(n3_d) != ceil(n3_d)) && (n2_d>0)){
   // n3_w not integer
   n3_d = _nprocs/n2_d;
   n2_d--;
  }
  _dims[0] = n2_d;
  _dims[1] = n3_d;
};

void block::setCoords(){
  _periods = new int[_ndims]();
  for (int i=0;i<_ndims;i++)
     _periods[i] = false;
  MPI_Cart_create(MPI_COMM_WORLD, _ndims, _dims, _periods, false,
                  &COMM_CART);
  coords = new int[_ndims]();
  MPI_Cart_coords(COMM_CART, _rank, _ndims, coords);
}

void block::setBlockInfo()
{
  // dims[0] -> n2; dims[1] ->n3
  int n2_t = floor(_n2/_dims[0]); //temporary n2
  int n2_left = _n2 - n2_t * _dims[0];
  if (coords[0]<n2_left){
    n2_b = n2_t + 1;
    x2_b = n2_b * coords[0];
  }
  else{
    n2_b = n2_t;
    x2_b = (coords[0]-n2_left)*n2_b
            + (n2_b+1)*n2_left;
  }

  int n3_t = floor(_n3/_dims[1]); //temporary n3
  int n3_left = _n3 - n3_t * _dims[1];
  if (coords[1]<n3_left){
    n3_b = n3_t + 1;
    x3_b = n3_b * coords[1];
  }
  else{
    n3_b = n3_t;
    x3_b = (coords[1]-n3_left)*n3_b
            + (n3_b+1)*n3_left;
  }
  n1_b = _n1;
}

void block::mpiReadModel(char * fin, float *** vin){
  ifstream in (fin, ios::in|ios::binary);
  if(!in){
    cerr << "Open error !" << endl;
    exit(1);
  }
  in.close();

  // cout << "n3_b = " << n3_b << endl;
  // cout << "n2_b = " << n2_b << endl;
  // cout << "n1_b = " << n1_b << endl;
  // cout << "x3_b = " << x3_b << endl;
  // cout << "x2_b = " << x2_b << endl;
  // cout << "_n1 = " << _n1 << endl;
  // cout << "_n2 = " << _n2 << endl;

  MPI_File fh;
  MPI_Offset disp;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, fin,
               MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  //LLINT fileSize = n2*n3;
  for (int i3=0;i3<n3_b;i3++){
   int x3 = i3 + x3_b;
   for (int i2=0;i2<n2_b;i2++){
     //cout << "i3 = " << i3 << endl;
     //cout << "i2 = " << i2 << endl;

     int x2 = i2 + x2_b;
     disp = x3 * _n1 * _n2 + _n1 * x2;
     disp = disp * sizeof(float);
     MPI_File_read_at(fh, disp, vin[i3][i2], n1_b, MPI_FLOAT,
                     &status);
   }
  }
  cout <<"Finish reading " << endl;
  MPI_File_close(&fh);
}

void block::mpiWriteModel(char * fout, float *** vin){
  MPI_File fh;
  MPI_Offset disp;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD, fout,
               MPI_MODE_RDWR|MPI_MODE_CREATE, MPI_INFO_NULL,&fh);

  // cout << "n3_b = " << n3_b << endl;
  // cout << "n2_b = " << n2_b << endl;
  // cout << "n1_b = " << n1_b << endl;

  //LLINT fileSize = n2*n3;
  for (int i3=0;i3<n3_b;i3++){
   int x3 = i3 + x3_b;
   for (int i2=0;i2<n2_b;i2++){
     int x2 = i2 + x2_b;
     disp = x3 * _n1 * _n2 + _n1 * x2;
     disp = disp * sizeof(float);
     MPI_File_write_at(fh, disp, vin[i3][i2], _n1, MPI_FLOAT,
                     &status);
   }
  }
  MPI_File_close(&fh);
}

block::~block(){
  if (coords!=nullptr)
    delete [] coords;
  if (_dims!=nullptr)
    delete [] _dims;
  if (_periods!=nullptr)
    delete [] _periods;
}
