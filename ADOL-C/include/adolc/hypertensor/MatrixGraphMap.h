#ifndef __MATRIX_GRAPH_MAP_H__
#define __MATRIX_GRAPH_MAP_H__

#include <map>

#include "MatrixGraph.h"

class VectorGraph;

class MatrixGraphMap : public MatrixGraph {
 public:
  MatrixGraphMap();
  MatrixGraphMap(std::map<int, std::map<int, double> >& source);
  MatrixGraphMap(std::map<int, std::map<int, double> >&& source);

  ~MatrixGraphMap();
  void increase(int x, int y, double v);
  VectorGraph* get_and_erase(int x);
  VectorGraph* get(int x);

// private:
  std::map<int, std::map<int, double> > data;
};


#endif // __MATRIX_GRAPH_MAP_H__
