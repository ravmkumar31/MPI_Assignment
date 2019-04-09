#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>
#include <math.h>

using namespace std;

float genA (int row, int col) {
  if (row > col)
    return 1.0;
  else
    return 0.0;
}

float genx0 (int i) {
  return 1.0;
}


void checkx (int iter, long i, float xval) {
  if (iter == 1) {
    float shouldbe = i;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99 )
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }

  if (iter == 2) {
    float shouldbe =(i-1)*i/2;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99)
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }
}

//perform dense y=Ax on an n \times n matrix
void matmul(float*A, float*x, float*y, long n) {
  
  for (long row = 0; row<n; ++row) {
    float sum = 0;
    
    for (long col = 0; col<n ; ++col) {
      //sum += x[col] *A[row][col]
       sum += x[col] * A[row*n+col];
    }

       y[row] = sum;
  }
}

int main (int argc, char* argv[]) {

  if (argc < 3) {
    std::cout<<"usage: "<<argv[0]<<" <n> <iteration>"<<std::endl;
  }


  //initialize data
  
   MPI_Init(&argc,&argv);

   long n = atol(argv[1]);

   int iteration = atoi(argv[2]);
   bool check = true;

   int Mpi_rank,valnp;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&Mpi_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&valnp);

   int p = sqrt(valnp);
   long div_val = n/p; 
   
   int Rdiv = Mpi_rank/p,Cdiv = Mpi_rank%p;
   
   MPI_Comm Rcomm;
   int Rrank;
   MPI_Comm_split(MPI_COMM_WORLD, Rdiv, Mpi_rank, &Rcomm);
   MPI_Comm_rank(Rcomm,&Rrank);

   MPI_Comm Ccomm;
   int Crank;
   MPI_Comm_split(MPI_COMM_WORLD, Cdiv, Mpi_rank, &Ccomm);
   MPI_Comm_rank(Ccomm,&Crank);

   float* arr = new float[div_val*div_val];
   long Rstart = (Rdiv*div_val),colstart = (Cdiv*div_val);
   long Rend = Rstart+div_val,colend = colstart+div_val;

  for (long row = Rstart,Rset=0; row<Rend; row++,Rset++) {
    for (long col= colstart,Cset=0; col<colend; col++,Cset++) {
      arr[(Rset*div_val)+Cset] = genA(row, col);
    }
  }
  float* x = new float[div_val];

  for (long i=0; i<div_val; i++)
    {
     x[i] = genx0(i);
    }
  float* y = new float[div_val];
  for (long i=0; i<div_val; i++)
    y[i] = 0.0;
 
   double start = MPI_Wtime(); 
   for (int ti = 0; ti<iteration; ti++) 
   {
      matmul(arr,x,y,div_val);
      MPI_Reduce(y,x,div_val,MPI_FLOAT,MPI_SUM,Rdiv,Rcomm);
      MPI_Bcast(x,div_val,MPI_FLOAT,Cdiv, Ccomm);
  if (check){
    for (long i = colstart,p=0; i<colend; ++i,++p){
        checkx (ti+1, i, x[p]);
             }
          }
  }
  
  if(Mpi_rank == 0){
        double end = MPI_Wtime(); 
        cerr<<end-start<<endl;
   }
 
  MPI_Comm_free(&Rcomm);
  MPI_Comm_free(&Ccomm);

  delete[] arr;
  delete[] x;
  delete[] y;

  MPI_Finalize();
 
  return 0;
}
