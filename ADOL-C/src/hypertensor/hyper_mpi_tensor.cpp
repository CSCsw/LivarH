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
                    int ndep,
                    int nindep,
                    const double* basepoint,
                    int order) {
  std::map<locint, HyperDerivative<locint> > global_gd;
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::cout << "In hyper_mpi_tensor" << std::endl;
  hyper_tape(tag, ndep, nindep, basepoint, ind_map, hyper_index, hyper_value);
//  for(const locint& x: hyper_index) {
//    std::cout << x << std::endl;
//  }
  hyper_mpi_reverse(tag, hyper_index, hyper_value, global_gd);

  typename std::map<locint, HyperDerivative<locint> >::iterator iter;
  locint dep;
/*
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    dep = iter->first;
    std::cout << "Derivatives for Dep " << dep << " : " << std::endl;
    global_gd[dep].debug();
    ++iter;
  }
*/
  hyper_mpi_forward(global_gd);
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    dep = iter->first;
    std::cout << "Derivatives for Dep " << dep << " : " << std::endl;
    global_gd[dep].debug();
    ++iter;
  }
}
