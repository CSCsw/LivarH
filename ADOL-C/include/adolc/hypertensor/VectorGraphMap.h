#ifndef __VECTOR_GRAPH_MAP_H__
#define __VECTOR_GRAPH_MAP_H__

#include <map>

#include "VectorGraph.h"

template <typename T>
class VectorGraphMap : public VectorGraph<T> {
 public:
  VectorGraphMap();
  VectorGraphMap(std::map<T, double>& source);
  VectorGraphMap(std::map<T, double>&& source);

  ~VectorGraphMap();

  void increase(T x, double v);
  double get_and_erase(T x);
  double get(T x);
  bool reset();
  bool get_next(T& x, double& w); 
  int get_size();

  void debug();

// private:
  std::map<T, double> data;
  typename std::map<T, double>::iterator iter;
};

template <typename T>
VectorGraphMap<T>::VectorGraphMap() {
  data.clear();
}

// L-value c-tor
template <typename T>
VectorGraphMap<T>::VectorGraphMap(std::map<T, double>& source) {
  std::cout << "VectorGraph (L-ctor)" << std::endl;
  data = source;
}

// R-value c-tor
template <typename T>
VectorGraphMap<T>::VectorGraphMap(std::map<T, double>&& source) {
//  std::cout << "VectorGraph (R-ctor)" << std::endl;
  data = std::move(source);
}

// D-tor
template <typename T>
VectorGraphMap<T>::~VectorGraphMap() {
  data.clear();
}


template <typename T>
void VectorGraphMap<T>::increase(T x, double v) {
  if (v == 0.0) {return;}
  data[x]+=v;
}

template <typename T>
double VectorGraphMap<T>::get_and_erase(T x) {
  double ret = data[x];
  data.erase(x);
  return ret;
}

template <typename T>
double VectorGraphMap<T>::get(T x) {
  return data[x];
}

template <typename T>
bool VectorGraphMap<T>::reset() {
  iter = data.begin();
  if (iter == data.end()) {
    return false;
  }
  return true;
}

template <typename T>
bool VectorGraphMap<T>::get_next(T& x, double& w) {
  x = iter->first;
  w = iter->second;
  ++iter;
  if (iter == data.end()) {
    return false;
  }
  return true;
}

template <typename T>
int VectorGraphMap<T>::get_size() {
  return data.size();
}

template <typename T>
void VectorGraphMap<T>::debug() {
  typename std::map<T, double>::iterator t_iter;
  t_iter = data.begin();
  while(t_iter != data.end()) {
    std::cout << "A[" << t_iter->first << "]=" << t_iter->second << std::endl;
    ++t_iter;
  }
}

#endif // __VECTOR_GRAPH_MAP_H__
