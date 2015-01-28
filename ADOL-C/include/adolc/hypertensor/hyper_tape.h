#ifndef __HYPER_TAPE_H__
#define __HYPER_TAPE_H__

#include <vector>
#include <map>
#include <limits.h>
#include <adolc/adolc.h>

int hyper_tape(short tag,
               int dep,
               int indep,
               const double* basepoint,
               std::map<locint, locint>& ind_map,
               std::vector<locint>& hyper_index,
               std::vector<double>& hyper_value);


#endif // __HYPER_TAPE_H__
