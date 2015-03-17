#ifndef GENERIC_MPI_TENSOR_H_
#define GENERIC_MPI_TENSOR_H_

int generic_mpi_tensor(short tag,
                      int dep,
                      int indep,
                      const double* basepoint,
                      int order);

#endif // GENERIC_MPI_TENSOR_H_
