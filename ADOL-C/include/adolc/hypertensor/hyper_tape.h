#ifndef __HYPER_TAPE_H__
#define __HYPER_TAPE_H__

#include <vector>
#include <map>

int hyper_tape(short tag,
               int dep,
               int indep,
               const double* basepoint,
               std::map<int, int> ind_map,
               std::vector<int>& hyper_index,
               std::vector<double>& hyper_value);


#endif // __HYPER_TAPE_H__
