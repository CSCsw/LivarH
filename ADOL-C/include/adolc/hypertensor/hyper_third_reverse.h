#ifndef HYPER_THIRD_REVERSE_H_
#define HYPER_THIRD_REVERSE_H_

#include <vector>
#include <adolc/adolc.h>

int hyper_third_reverse(short tag,
                        std::vector<locint>& hyper_index,
                        std::vector<double>& hyper_value,
                        VectorGraph<locint>* adjoints,
                        MatrixGraph<locint>* hessian,
                        HyperGraph<locint>* tensor);


#endif // HYPER_THRID_REVERSE_H_
