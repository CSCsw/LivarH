#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>

int hyper_third(DerivativeInfo<locint>& info,
                VectorGraph<locint>* adjoints,
                MatrixGraph<locint>* hessian,
                HyperGraph<locint>* tensor,
                double w,
                VectorGraph<locint>* r,
                MatrixGraph<locint>* e);

int hyper_hessian(DerivativeInfo<locint>& info,
                  VectorGraph<locint>* adjoints,
                  MatrixGraph<locint>* hessian,
                  double w,
                  VectorGraph<locint>* r);

int hyper_adjoints(DerivativeInfo<locint>& info,
                   VectorGraph<locint>* adjoints,
                   double w);
