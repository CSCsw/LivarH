#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/hyper_derivative.h>
#include <adolc/hypertensor/hyper_tape.h>
#include <adolc/hypertensor/hyper_mpi_reverse.h>

// options[0] = 1, first order;
//            = 2, second order;
//            = 3, third order;
int hyper_mpi_tensor(short tag,
                    int dep,
                    int indep,
                    const double* basepoint,
                    int* nnz,
                    unsigned int** rind,
                    unsigned int** cind,
                    unsigned int** xind,
                    double** values,
                    int* options) {
  std::map<locint, HyperDerivative<locint> > global_gd;
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::cout << "In hyper_mpi_tensor" << std::endl;
  hyper_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);
//  for(const locint& x: hyper_index) {
//    std::cout << x << std::endl;
//  }
  hyper_mpi_reverse(tag, hyper_index, hyper_value, global_gd);

}
