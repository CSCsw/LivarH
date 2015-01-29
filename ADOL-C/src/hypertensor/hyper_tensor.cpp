#include <vector>
#include <map>
#include <iostream>

#include <adolc/adolc.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/HyperGraphMap.h>
#include <adolc/hypertensor/hyper_tape.h>
#include <adolc/hypertensor/hyper_third_reverse.h>

int hyper_tensor(short tag,
                 int dep,
                 int indep,
                 const double* basepoint,
                 int* nnz,
                 unsigned int** rind,
                 unsigned int** cind,
                 unsigned int** xind,
                 double** values,
                 int* optinos) {
  VectorGraph<locint>* adjoints = new VectorGraphMap<locint>();
  MatrixGraph<locint>* hessian = new MatrixGraphMap<locint>();
  HyperGraph<locint>* tensor = new HyperGraphMap<locint>();
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::cout << "In hyper_tensor" << std::endl;
  hyper_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);
//  hyper_third_reverse(tag, hyper_graph, hyper_index, hyper_value);
  for(const locint& x: hyper_index) {
    std::cout << x << std::endl;
  }
  hyper_third_reverse(tag, hyper_index, hyper_value, adjoints, hessian, tensor);

// result check
  bool has_next;
  locint x;
  locint y;
  double w;
  has_next = adjoints->reset();
  while(has_next) {
    has_next = adjoints->get_next(x, w);
    std::cout << "A[" << x << "]=" << w << std::endl;
  }
  has_next = hessian->reset();
  while(has_next) {
    has_next = hessian->get_next(x, y, w);
    std::cout << "H[" << x << "," << y << "]=" << w << std::endl;
  }

  delete adjoints;
  delete hessian;
  delete tensor;
}
