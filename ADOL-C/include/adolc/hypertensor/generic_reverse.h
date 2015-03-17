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

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      std::set<locint>& live_set,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd);

double binary_comb_coefficient(int max_level, int* dx, int* dy);

double unary_comb_coefficient(int max_level, int* dx);

void generate_binary_tuples(int curr_level,
                            int max_level,
                            int curr_order,
                            int max_order,
                            double cw,
                            typename GenericDerivative<locint>::iterator& prev,
                            int* dx,
                            int* dy,
                            double* ssw, int order, locint x, locint y);

void generate_unary_tuples(int curr_level,
                           int max_level,
                           int curr_order,
                           int max_order,
                           double cw,
                           typename GenericDerivative<locint>::iterator& prev,
                           int* dx,
                           double* sw);

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      std::set<locint>& live_set,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd);


#endif // GENERIC_REVERSE_H_
