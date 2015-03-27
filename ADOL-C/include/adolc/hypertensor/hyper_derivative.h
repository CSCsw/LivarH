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
  int get_byte_size();
  void write_to_byte(char* buf);

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
HyperDerivative<T>::HyperDerivative(char* buf) {
  char* p = buf;
  adjoints = new VectorGraphMap<T>(p);
  p += adjoints.get_byte_size();
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

#endif // HYPER_DERIVATIVE_H_
