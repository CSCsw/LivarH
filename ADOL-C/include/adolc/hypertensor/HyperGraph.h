#ifndef __HYPER_GRAPH_H__
#define __HYPER_GRAPH_H__

class MatrixGraph;

class HyperGraph {
 public:
  virtual void increase(int x, int y, int z, double v) = 0;
  virtual MatrixGraph* get_and_erase(int x) = 0;
  virtual MatrixGraph* get(int x) = 0;
};

#endif // __HYPER_GRAPH_H__
