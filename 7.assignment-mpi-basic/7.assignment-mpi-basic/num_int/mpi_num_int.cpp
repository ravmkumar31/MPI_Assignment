#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>
#include <chrono>

using namespace std;
#ifdef __cplusplus
extern "C" {
#endif

float f1(float x, int intensity);
float f2(float x, int intensity);
float f3(float x, int intensity);
float f4(float x, int intensity);

#ifdef __cplusplus
}
#endif

void mpi_numint(int f,float a,float b,int intensity,int n ,int mpiRank,int mpisize)
{
  float function_eval;
  float y = 0;
  float mul = (b-a)/n;
  float compute_numeric;
  for(int x=(mpiRank)*n/mpisize;x<(mpiRank+1)*n/mpisize;x++)
    {
      switch(f)
      {
        case 1:
          y += f1(a + (x+0.5)*mul,intensity);
          break;
        case 2:
          y += f2(a + (x+0.5)*mul,intensity);
          break;
        case 3:
          y += f3(a + (x+0.5)*mul,intensity);
          break;
        case 4:
          y += f4(a + (x+0.5)*mul,intensity);
          break;
      }
    }
      y =  mul* (y);
      MPI_Reduce(&y, &compute_numeric, 1,MPI_FLOAT,MPI_SUM, 0, MPI_COMM_WORLD);
  if(mpiRank==0)
    cout<<compute_numeric;
}
  
int main (int argc, char* argv[]) {
  if (argc < 6) {
    std::cerr<<"usage: "<<argv[0]<<" <function_id> <a> <b> <n> <intensity>"<<std::endl;
    return 0;
  }

  int function_id = atoi(argv[1]);
  float a = atof(argv[2]);
  float b = atof(argv[3]);
  int n = atoi(argv[4]);
  int intensity = atoi(argv[5]);
  int mpiRank,mpisize;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  std::chrono::time_point<std::chrono::system_clock> start_time; 
  if(mpiRank ==0)
    start_time = std::chrono::system_clock::now();
  mpi_numint(function_id,a,b,intensity,n,mpiRank,mpisize);
  MPI_Finalize();
  if(mpiRank==0)
  {
    std::chrono::time_point<std::chrono::system_clock> end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> total_time = end_time-start_time;
    std::cerr<<total_time.count()<<std::endl;
  }
  return 0;
}