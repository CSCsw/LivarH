#include <map>
#include <iostream>

#include "VectorGraph.h"
#include "VectorGraphMap.h"
#include "MatrixGraph.h"
#include "MatrixGraphMap.h"
#include "HyperGraph.h"
#include "HyperGraphMap.h"

int main() {
  HyperGraphMap<int>* hyper_graph = new HyperGraphMap<int>();
  hyper_graph->increase(1,1,1,1.0);
  hyper_graph->increase(1,1,2,2.0);
  hyper_graph->increase(1,2,1,3.0);
  hyper_graph->increase(1,1,3,4.0);
  hyper_graph->increase(2,1,1,5.0);
  MatrixGraphMap<int>* matrix_graph = (MatrixGraphMap<int>*)hyper_graph->get_and_erase(1);
  VectorGraphMap<int>* vector_graph = (VectorGraphMap<int>*)matrix_graph->get_and_erase(1);
  std::cout << vector_graph->get(1) << std::endl;
  std::cout << vector_graph->get(3) << std::endl;
  std::cout << vector_graph->data.size() << std::endl;
  std::cout << matrix_graph->data.size() << std::endl;
  std::cout << hyper_graph->data.size() << std::endl;
}
