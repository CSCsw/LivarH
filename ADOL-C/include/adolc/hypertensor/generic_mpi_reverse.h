#ifndef GENERIC_MPI_REVERSE_H_
#define GENERIC_MPI_REVERSE_H_

#include <vector>
#include <adolc/adolc.h>
#include "generic_derivative.h"

int generic_mpi_reverse(short tag,
                        int order,
                        std::vector<locint>& hyper_index,
                        std::vector<double>& hyper_value,
                        std::map<locint, std::set<locint> >& live_set,
                        std::map<locint, GenericDerivative<locint> >& generic_derivative);

int generic_mpi_forward(int order,
                        std::map<locint, std::set<locint> >& live_set,
                        std::map<locint, GenericDerivative<locint> >& generic_derivative);
#endif // GENERIC_MPI_REVERSE_H_
