#ifndef HYPER_MPI_TENSOR_H_
#define HYPER_MPI_TENSOR_H_

int hyper_mpi_tensor(short tag,
                    int dep,
                    int indep,
                    const double* basepoint,
                    int order);

#endif // HYPER_MPI_TENSOR_H_
