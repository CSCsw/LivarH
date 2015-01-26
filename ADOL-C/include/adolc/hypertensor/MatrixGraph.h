#ifndef __MATRIX_GRAPH_H__
#define __MATRIX_GRAPH_H__

class VectorGraph;

class MatrixGraph {
 public:
  virtual void increase(int x, int y, double v) = 0;
  virtual VectorGraph* get_and_erase(int x) = 0;
  virtual VectorGraph* get(int x) = 0;
};

#endif // __MATRIX_GRAPH_H__
