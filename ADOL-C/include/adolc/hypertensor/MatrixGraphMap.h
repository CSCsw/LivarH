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

// private:
  std::map<T, std::map<T, double> > data;
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

#endif // __MATRIX_GRAPH_MAP_H__
