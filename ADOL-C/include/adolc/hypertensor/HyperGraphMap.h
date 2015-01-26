#ifndef __HYPER_GRAPH_MAP_H__
#define __HYPER_GRAPH_MAP_H__

#include <map>

#include "MatrixGraph.h"
#include "MatrixGraphMap.h"
#include "HyperGraph.h"


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
template <typename T>
HyperGraphMap<T>::HyperGraphMap() {
  data.clear();
}

template <typename T>
HyperGraphMap<T>::~HyperGraphMap() {
  data.clear();
}

template <typename T>
void HyperGraphMap<T>::increase(T x, T y, T z, double v) {
  data[x][y][z]+=v;
}

template <typename T>
MatrixGraph<T>* HyperGraphMap<T>::get_and_erase(T x) {
  MatrixGraph<T>* ret = new MatrixGraphMap<T>(std::move(data[x]));
  data.erase(x);
  return ret;
}

template <typename T>
MatrixGraph<T>* HyperGraphMap<T>::get(T x) {
  MatrixGraph<T>* ret = new MatrixGraphMap<T>(data[x]);
  return ret;
}



#endif // __HYPER_GRAPH_MAP_H__
