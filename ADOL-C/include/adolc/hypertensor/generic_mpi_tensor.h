#ifndef GENERIC_MPI_TENSOR_H_
#define GENERIC_MPI_TENSOR_H_

int generic_mpi_tensor(short tag,
                      int dep,
                      int indep,
                      const double* basepoint,
                      int order,
                      int** nnz_p,
                      unsigned int**** indices_p,
                      double*** values_p);

#endif // GENERIC_MPI_TENSOR_H_
