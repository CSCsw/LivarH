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
  HyperDerivative(char* buf);
  
  void init();
  int byte_size();
  void write_to_byte(char* buf);
  void debug();

  VectorGraph<T>* adjoints;
  MatrixGraph<T>* hessian;
//  HyperGraph<T>* tensor;
};

template <typename T>
HyperDerivative<T>::HyperDerivative() {
  adjoints = NULL;
  hessian = NULL;
//  tensor = NULL;
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
}
template <typename T>
HyperDerivative<T>::HyperDerivative(char* buf) {
  char* p = buf;
  adjoints = new VectorGraphMap<T>(p);
  p += adjoints->byte_size();
  hessian = new MatrixGraphMap<T>(p);

//  p += hessian.get_byte_size();
//  tensor = new HyperGraphMap<T>(p);
}

template <typename T>
void HyperDerivative<T>::init() {
  if (adjoints == NULL) adjoints = new VectorGraphMap<T>();
  if (hessian == NULL) hessian = new MatrixGraphMap<T>();
//  if (tensor == NULL) tensor = new HyperGraphMap<T>();
}

template <typename T>
void HyperDerivative<T>::debug() {
  if (adjoints != NULL) adjoints->debug();
  if (hessian != NULL) hessian->debug();
}
#endif // HYPER_DERIVATIVE_H_
