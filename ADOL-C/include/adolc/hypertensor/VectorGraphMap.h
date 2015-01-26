#ifndef __VECTOR_GRAPH_MAP_H__
#define __VECTOR_GRAPH_MAP_H__

#include <map>

#include "VectorGraph.h"

template <typename T>
class VectorGraphMap : public VectorGraph {
 public:
  VectorGraphMap();
  VectorGraphMap(std::map<T, double>& source);
  VectorGraphMap(std::map<T, double>&& source);

  ~VectorGraphMap();

  void increase(T x, double v);
  double get_and_erase(T x);
  double get(T x);
 
// private:
  std::map<T, double> data;
};


#endif // __VECTOR_GRAPH_MAP_H__
