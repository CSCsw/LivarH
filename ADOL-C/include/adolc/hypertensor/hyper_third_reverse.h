#include <vector>
#include <map>

#include <adolc/adolc.h>

int hyper_third_reverse(short tag,
                        std::vector<locint>& hyper_index,
                        std::vector<double>& hyper_value,
                        VectorGraph<locint>* adjoints,
                        MatrixGraph<locint>* hessian,
                        HyperGraph<locint>* tensor);

