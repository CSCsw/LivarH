#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/hyper_tape.h>
#include <adolc/hypertensor/opencomb.h>
#include <adolc/hypertensor/generic_derivative_table.h>
#include <adolc/hypertensor/generic_reverse.h>

int generic_tensor(short tag,
                  int dep,
                  int indep,
                  const double* basepoint,
                  int d) {

  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  GenericDerivative<locint> generic_derivative(d);
  std::cout << "In generic_tensor" << std::endl;
  hyper_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);
/*
  std::cout << "index: " << std::endl;
  for(const locint& loc : hyper_index) {
    std::cout << loc << std::endl;
  }
  std::cout << "VAlues" << std::endl;
  for(const double& v : hyper_value) {
    std::cout << v << std::endl;
  }
*/
  special_derivative_table();
  generic_reverse(tag, d, hyper_index, hyper_value, generic_derivative);
  generic_derivative.debug();
}
