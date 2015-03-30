#ifndef __VECTOR_GRAPH_H__
#define __VECTOR_GRAPH_H__

template <typename T>
class VectorGraph {
 public:
  class iterator {
   public:
    virtual ~iterator() {};
    virtual bool reset() = 0;
    virtual bool get_next(T& x, double& w) = 0;
    virtual typename VectorGraph<T>::iterator* copy_iter() = 0;
  };
  virtual bool has_live(T target) const = 0;
  virtual int byte_size() const = 0;
  virtual void write_to_byte(char* buf) const = 0;
  virtual void increase(T x, double v) = 0;
  virtual double get_and_erase(T x) = 0;
  virtual double get(T x) = 0; 
  virtual int get_size() const = 0;
  virtual void debug() const = 0;
  virtual typename VectorGraph<T>::iterator* get_iterator() = 0;
};

#endif // __VECTOR_GRAPH_H__
