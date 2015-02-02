#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/HyperGraphMap.h>
#include <adolc/hypertensor/hyper_tape.h>
#include <adolc/hypertensor/hyper_third_reverse.h>


// options[0] = 1, first order;
//            = 2, second order;
//            = 3, third order;
int hyper_tensor(short tag,
                 int dep,
                 int indep,
                 const double* basepoint,
                 int* nnz,
                 unsigned int** rind,
                 unsigned int** cind,
                 unsigned int** xind,
                 double** values,
                 int* options) {
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
  std::cout << "|A|=" << adjoints->get_size() << std::endl;
  while(has_next) {
    has_next = adjoints->get_next(x, w);
    std::cout << "A[" << x << "]=" << w << std::endl;
  }
  has_next = hessian->reset();
  std::cout << "|H|" << hessian->get_size() << std::endl;
  while(has_next) {
    has_next = hessian->get_next(x, y, w);
    std::cout << "H[" << x << "," << y << "]=" << w << std::endl;
  }
  if (options[0] == 1) {
    int n = adjoints->get_size();
    (*nnz) = n;
    free(*rind);
    (*rind) = (unsigned int*) malloc(sizeof(unsigned int) * n);
    free(*values);
    (*values) = (double*) malloc(sizeof(double) * n);
    has_next = adjoints->reset();
    int i = 0;
    while(has_next) {
      has_next = adjoints->get_next(x, w);
      (*rind)[i] = x;
      (*values)[i] = w;
      ++i;
    }
  } else if (options[0] == 2) {
    int n = hessian->get_size();
    (*nnz) = n;
    free(*rind);
    (*rind) = (unsigned int*) malloc(sizeof(unsigned int) * n);
    free(*cind);
    (*cind) = (unsigned int*) malloc(sizeof(unsigned int) * n);
    free(*values);
    (*values) = (double*) malloc(sizeof(double) * n);
    int i = 0;
    has_next = hessian->reset();
    while(has_next) {
      has_next = hessian->get_next(x, y, w);
      (*rind)[i] = x;
      (*cind)[i] = y;
      (*values)[i] = w;
      ++i;
    }
  }
  delete adjoints;
  delete hessian;
  delete tensor;
}
