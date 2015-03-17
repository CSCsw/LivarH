#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/generic_tape.h>
#include <adolc/hypertensor/opencomb.h>
#include <adolc/hypertensor/generic_derivative_table.h>
#include <adolc/hypertensor/generic_mpi_reverse.h>

int generic_mpi_tensor(short tag,
                       int dep,
                       int indep,
                       const double* basepoint,
                       int d) {
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::map<locint, GenericDerivative<locint> > generic_derivative;
  std::cout << "In generic_mpi_tensor" << std::endl;
  generic_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);
  special_derivative_table();
  generic_mpi_reverse(tag, d, hyper_index, hyper_value, generic_derivative);
  std::map<locint, GenericDerivative<locint> >::iterator iter;
  iter = generic_derivative.begin();
  while(iter != generic_derivative.end() ){
    std::cout << "Derivative For : " << iter->first << std::endl;
    iter->second.debug();
    iter++;
  }
}
