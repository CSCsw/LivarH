#ifndef __VECTOR_GRAPH_H__
#define __VECTOR_GRAPH_H__

template <typename T>
class VectorGraph {
 public:
  virtual void increase(T x, double v) = 0;
  virtual double get_and_erase(T x) = 0;
  virtual double get(T x) = 0; 
  virtual bool reset() = 0;
  virtual bool get_next(T& x, double& w) = 0;
  virtual int get_size() = 0;
};

#endif // __VECTOR_GRAPH_H__
