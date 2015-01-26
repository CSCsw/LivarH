#ifndef __VECTOR_GRAPH_H__
#define __VECTOR_GRAPH_H__

template <typename T>
class VectorGraph {
 public:
  virtual void increase(T x, T v) = 0;
  virtual double get_and_erase(T x) = 0;
  virtual double get(T x) = 0; 
};

#endif // __VECTOR_GRAPH_H__
