#ifndef HYPER_DERIVATIVE_H_
#define HYPER_DERIVATIVE_H_

#include "VectorGraph.h"
#include "VectorGraphMap.h"
#include "MatrixGraph.h"
#include "MatrixGraphMap.h"
#include "HyperGraph.h"
#include "HyperGraphMap.h"


//This is a map implementation

template <typename T>
class HyperDerivative {
 public:
  HyperDerivative();
  HyperDerivative(char* buf, int order);
  ~HyperDerivative();
  
  void init(int order);
  int byte_size();
  void write_to_byte(char* buf);
  void debug();
  void ind_map(T toind, T fromind, int count);

  VectorGraph<T>* adjoints;
  MatrixGraph<T>* hessian;
  HyperGraph<T>* tensor;
};

template <typename T>
HyperDerivative<T>::HyperDerivative() {
  adjoints = NULL;
  hessian = NULL;
  tensor = NULL;
}

template <typename T>
HyperDerivative<T>::~HyperDerivative() {
  delete adjoints;
  delete hessian;
  if (tensor != NULL) {delete tensor;}
}

template <typename T>
int HyperDerivative<T>::byte_size() {
  return adjoints->byte_size() + hessian->byte_size();
}
template <typename T>
void HyperDerivative<T>::write_to_byte(char* buf) {
  char* p = buf;
  adjoints->write_to_byte(p);
  p += adjoints->byte_size();
  hessian->write_to_byte(p);
  if (tensor != NULL) {
    p += hessian->byte_size();
    tensor->write_to_byte(p);
  }
}
template <typename T>
HyperDerivative<T>::HyperDerivative(char* buf, int order) {
  char* p = buf;
  adjoints = new VectorGraphMap<T>(p);

  p += adjoints->byte_size();
  hessian = new MatrixGraphMap<T>(p);

  if (order >= 3) {
    p += hessian->byte_size();
    tensor = new HyperGraphMap<T>(p);
  } else {
    tensor = NULL;
  }
}

template <typename T>
void HyperDerivative<T>::init(int order) {
  if (adjoints == NULL) adjoints = new VectorGraphMap<T>();
  if (hessian == NULL) hessian = new MatrixGraphMap<T>();
  if (order >= 3) {
    if (tensor == NULL) tensor = new HyperGraphMap<T>();
  }
}

template <typename T>
void HyperDerivative<T>::debug() {
  if (adjoints != NULL) adjoints->debug();
  if (hessian != NULL) hessian->debug();
//  if (tensor != NULL) tensor->debug();
}

template <typename T>
void HyperDerivative<T>::ind_map(T toind, T fromind, int count) {
  VectorGraph<T>* t_adjoints = adjoints;
  adjoints = new VectorGraphMap<T>();
  T x,y;
  double w;
  typename VectorGraph<T>::iterator* v_iter = t_adjoints->get_iterator();
  bool has_next = v_iter->reset();
  while(has_next) {
    has_next = v_iter->get_next(x, w);
    if (x > fromind - count && x <= fromind) {
      x = toind + (fromind - x);
    }
    adjoints->increase(x, w);
  }
  delete v_iter;

  MatrixGraph<T>* t_hessian = hessian;
  hessian = new MatrixGraphMap<T>();
  typename MatrixGraph<T>::iterator* h_iter = t_hessian->get_iterator();
  has_next = h_iter->reset();
  while(has_next) {
    has_next = h_iter->get_next(x, y, w);
    if (x > fromind - count && x <= fromind) {
      x = toind + (fromind - x);
    }
    if (y > fromind - count && y <= fromind) {
      y = toind + (fromind - y);
    }
    hessian->increase(x, y, w);
  }
  delete h_iter;
  delete t_adjoints;
  delete t_hessian;
/*
  std::cout << fromind << " --> " << toind << "[" << count << "]" << std::endl;
  t_adjoints->debug();
  t_hessian->debug();
  adjoints->debug();
  hessian->debug();
*/
}


#endif // HYPER_DERIVATIVE_H_
