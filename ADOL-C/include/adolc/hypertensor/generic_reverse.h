#ifndef GENERIC_REVERSE_H_
#define GENERIC_REVERSE_H_

#include <vector>
#include <adolc/adolc.h>
#include "generic_derivative.h"

int generic_reverse(short tag,
                    int order,
                    std::vector<locint>& hyper_index,
                    std::vector<double>& hyper_value,
                    GenericDerivative<locint>& generic_derivative);

#endif // GENERIC_REVERSE_H_
