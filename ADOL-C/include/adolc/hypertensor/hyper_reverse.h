#ifndef HYPER_REVERSE_H_
#define HYPER_REVERSE_H_

#include <vector>
#include <adolc/adolc.h>

int hyper_reverse(short tag,
                  std::vector<locint>& hyper_index,
                  std::vector<double>& hyper_value,
                  VectorGraph<locint>* adjoints,
                  MatrixGraph<locint>* hessian,
                  HyperGraph<locint>* tensor);


#endif // HYPER_REVERSE_H_
