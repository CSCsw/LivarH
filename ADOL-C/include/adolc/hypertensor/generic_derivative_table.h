#ifndef GENERIC_DERIVATIVE_TABLE_H_
#define GENERIC_DERIVATIVE_TABLE_H_

#include "adolc/adolc.h"
#include "adolc/hypertensor/hyper_common.h"
#include "adolc/hypertensor/generic_derivative.h"
void populate_derivative_table(int order,
                               DerivativeInfo<locint>& info,
                               GenericDerivative<locint>& local_gd);

#endif // GENERIC_DERIVATIVE_TALBE_H_
