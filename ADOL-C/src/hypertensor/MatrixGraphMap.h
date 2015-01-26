#ifndef __MATRIX_GRAPH_MAP_H__
#define __MATRIX_GRAPH_MAP_H__

#include <map>

#include "MatrixGraph.h"

template <typename T> class VectorGraph;

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


#endif // __MATRIX_GRAPH_MAP_H__
