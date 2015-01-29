#ifndef __MATRIX_GRAPH_MAP_H__
#define __MATRIX_GRAPH_MAP_H__

#include <map>

#include "VectorGraph.h"
#include "VectorGraphMap.h"
#include "MatrixGraph.h"

template <typename T>
class MatrixGraphMap : public MatrixGraph<T> {
 public:
  MatrixGraphMap();
  MatrixGraphMap(std::map<T, std::map<T, double> >& source);
  MatrixGraphMap(std::map<T, std::map<T, double> >&& source);

  ~MatrixGraphMap();
  void increase(T x, T y, double v);
  VectorGraph<T>* get_and_erase(T x);
  VectorGraph<T>* get(T x);
  bool reset();
  bool get_next(T& x, T& y, double& w);
// private:
  std::map<T, std::map<T, double> > data;
  typename std::map<T, std::map<T, double> >::iterator iter;
  typename std::map<T, double>::iterator iter2;
};

template <typename T>
MatrixGraphMap<T>::MatrixGraphMap() {
  data.clear();
}

// L-value ctor, copy assignment
template <typename T>
MatrixGraphMap<T>::MatrixGraphMap(std::map<T, std::map<T, double>>& source) {
//  std::cout << "MatrixGraphMap (L-ctor)" << std::endl;
  data = source;
}

// R-value c-tor, move assignment
template <typename T>
MatrixGraphMap<T>::MatrixGraphMap(std::map<T, std::map<T, double>>&& source) {
//  std::cout << "MatrixGraphMap (R-ctor)" << std::endl;
  data = std::move(source);
}

template <typename T>
MatrixGraphMap<T>::~MatrixGraphMap() {
  data.clear();
}

template <typename T>
void MatrixGraphMap<T>::increase(T x, T y, double v) {
  data[x][y]+=v;
}

template <typename T>
VectorGraph<T>* MatrixGraphMap<T>::get_and_erase(T x) {
  VectorGraph<T>* ret = new VectorGraphMap<T>(std::move(data[x]));
  data.erase(x);
  return ret;
}

template <typename T>
VectorGraph<T>* MatrixGraphMap<T>::get(T x) {
  VectorGraph<T>* ret = new VectorGraphMap<T>(data[x]);
  return ret;
}

template <typename T>
bool MatrixGraphMap<T>::reset() {
  iter = data.begin();
  while(iter != data.end()) {
    iter2 = iter->second.begin();
    if (iter2 != iter->second.end()) {
      return true;
    }
    ++iter;
  }
  return false;
}

template <typename T>
bool MatrixGraphMap<T>::get_next(T& x, T& y, double& w) {
  x = iter->first;
  y = iter2->first;
  w = iter2->second;
  ++iter2;
  while(iter != data.end()) {
    if (iter2 != iter->second.end()) {
      return true;
    }
    ++iter;
    iter2 = iter->second.begin();
  }
  return false;
}

#endif // __MATRIX_GRAPH_MAP_H__
