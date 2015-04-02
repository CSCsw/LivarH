#ifndef __HYPER_GRAPH_H__
#define __HYPER_GRAPH_H__

template <typename T> class MatrixGraph;

template <typename T>
class HyperGraph {
 public:
  virtual void increase(T x, T y, T z, double v) = 0;
  virtual MatrixGraph<T>* get_and_erase(T x) = 0;
  virtual MatrixGraph<T>* get(T x) = 0;
  virtual int byte_size() = 0;
  virtual void write_to_byte(char* buf) = 0;
  virtual bool reset() = 0;
  virtual bool get_next(T& x, T& y, T& z, double& w) = 0;
  virtual int get_size() = 0;
  virtual void debug() = 0;
};

#endif // __HYPER_GRAPH_H__
