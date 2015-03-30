#ifndef __MATRIX_GRAPH_H__
#define __MATRIX_GRAPH_H__

template <typename T> class VectorGraph;

template <typename T>
class MatrixGraph {
 public:
  virtual void increase(T x, T y, double v) = 0;
  virtual VectorGraph<T>* get_and_erase(T x) = 0;
  virtual VectorGraph<T>* get(T x) = 0;
  virtual int get_size() const = 0;
  virtual void debug() const = 0;
  virtual int byte_size() const = 0;
  virtual void write_to_byte(char* buf) const = 0;
  class iterator {
   public:
    virtual ~iterator() {};
    virtual bool reset() = 0;
    virtual bool get_next(T& x, T& y, double& w) = 0; 
  };

  virtual typename MatrixGraph<T>::iterator* get_iterator() = 0;
};

#endif // __MATRIX_GRAPH_H__
