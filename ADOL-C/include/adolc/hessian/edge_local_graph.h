#ifndef _EDGE_LOCAL_GRAPH_H_
#define _EDGE_LOCAL_GRAPH_H_

#include <adolc/hessian/edge_main.h>
#include <adolc/internal/adolc_settings.h>

#define EDGE_LOCAL_SIZE 50

class EdgeLocalGraph {
 public:
  // ctor, dtor
  EdgeLocalGraph();
  ~EdgeLocalGraph();

  size_t AddLiveVar(locint ind);
  size_t find(locint ind, locint* array, size_t size);
  void insert(size_t row, locint ind, double w);

  void reset();
  void erase(size_t ind);

  void Print();
  
  size_t size;
  locint *loc;
  double *adjoints;
  
  size_t *size_array;
  locint **loc_array;
  double **hessian;
};

#endif  // _EDGE_LOCAL_GRAPH_H_
