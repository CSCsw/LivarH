#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>

int hyper_third(DerivativeInfo<locint>& info,
                VectorGraph<locint>* adjoints,
                MatrixGraph<locint>* hessian,
                HyperGraph<locint>* tensor);

int hyper_hessian(DerivativeInfo<locint>& info,
                  VectorGraph<locint>* adjoints,
                  MatrixGraph<locint>* hessian);

int hyper_adjoints(DerivativeInfo<locint>& info,
                   VectorGraph<locint>* adjoints);
