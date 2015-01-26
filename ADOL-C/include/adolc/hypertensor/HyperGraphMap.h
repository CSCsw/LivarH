#ifndef __HYPER_GRAPH_MAP_H__
#define __HYPER_GRAPH_MAP_H__

#include <map>

#include "HyperGraph.h"

class MatrixGraph;

class HyperGraphMap : public HyperGraph {
 public:
  HyperGraphMap();
  ~HyperGraphMap();
  void increase(int x, int y, int z, double v);
  MatrixGraph* get_and_erase(int x);
  MatrixGraph* get(int x);

// private:
  std::map<int, std::map<int, std::map<int, double> > > data;
};


#endif // __HYPER_GRAPH_MAP_H__
