#ifndef GENERIC_MPI_H_
#define GENERIC_MPI_H_

#include <mpi.h>

extern MPI_Datatype RMPI_ADOUBLE;

void RMPI_Init(int* argc, char** argv[]);
void RMPI_Finalize();


int RMPI_Send(void* buf,
              int count,
              MPI_Datatype datatype,
              int dest,
              int tag,
              MPI_Comm comm);

int RMPI_Recv(void* buf,
              int count,
              MPI_Datatype datatype,
              int src,
              int tag,
              MPI_Comm comm,
              MPI_Status* status);

int RMPI_Send_ind(void* buf,
              int count,
              MPI_Datatype datatype,
              int dest,
              int tag,
              MPI_Comm comm);

int RMPI_Recv_ind(void* buf,
              int count,
              MPI_Datatype datatype,
              int src,
              int tag,
              MPI_Comm comm,
              MPI_Status* status);

#endif // GENERIC_MPI_H_
