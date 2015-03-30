#ifndef HYPER_MPI_REVERSE_H_
#define HYPER_MPI_REVERSE_H_

#include <vector>
#include <map>
#include <adolc/adolc.h>

int hyper_mpi_reverse(short tag,
                  std::vector<locint>& hyper_index,
                  std::vector<double>& hyper_value,
                  std::map<locint, HyperDerivative<locint> >& global_gd);


void hyper_mpi_forward(std::map<locint, HyperDerivative<locint> >& global_gd);

#endif // HYPER_MPI_REVERSE_H_
