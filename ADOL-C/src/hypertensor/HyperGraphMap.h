#ifndef __HYPER_GRAPH_MAP_H__
#define __HYPER_GRAPH_MAP_H__

#include <map>

#include "HyperGraph.h"

template <typename T> class MatrixGraph;

template <typename T>
class HyperGraphMap : public HyperGraph<T> {
 public:
  HyperGraphMap();
  ~HyperGraphMap();
  void increase(T x, T y, T z, double v);
  MatrixGraph<T>* get_and_erase(T x);
  MatrixGraph<T>* get(T x);

// private:
  std::map<T, std::map<T, std::map<T, double> > > data;
};


#endif // __HYPER_GRAPH_MAP_H__
