#include <mpi.h>
#include <unistd.h>
#include <iostream>

int main(int argc, char*argv[]) {
int rank , size;
char hostname[1024];
   hostname[1023] = '\0';
   gethostname(hostname, 1023);
MPI_Init(&argc,& argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
printf(" i am process %d out of %d. i am running on %s. ",rank, size,hostname);
MPI_Finalize();

}