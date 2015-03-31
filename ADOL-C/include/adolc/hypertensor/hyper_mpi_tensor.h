#ifndef HYPER_MPI_TENSOR_H_
#define HYPER_MPI_TENSOR_H_

int hyper_mpi_tensor(short tag,
                    int dep,
                    int indep,
                    const double* basepoint,
                    int order,
                    int** nnz_p,
                    unsigned int**** indices_p,
                    double*** values_p);

#endif // HYPER_MPI_TENSOR_H_
