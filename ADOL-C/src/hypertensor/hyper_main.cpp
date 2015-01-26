#include <vector>
#include <map>
#include <stdio.h>
#include <adolc/adolc.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/HyperGraphMap.h>

int hyper_tensor(short tag,
                 int dep,
                 int indep,
                 const double* basepoint,
                 int* nnz,
                 unsigned int** rind,
                 unsigned int** cind,
                 unsigned int** xind,
                 double** values,
                 int* optinos) {
  HyperGraph* hyper_graph = new HyperGraphMap();
  std::vector<int> hyper_index;
  std::vector<double> hyper_value;
  std::map<int, int> ind_map;
  std::cout << "In hyper_tensor" << std::endl;
//  hyper_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);
//  hyper_third_reverse(tag, hyper_graph, hyper_index, hyper_value);
  delete hyper_graph;
}
