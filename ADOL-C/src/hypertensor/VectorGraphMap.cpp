#include <map>
#include <iostream>

#include "VectorGraph.h"
#include "VectorGraphMap.h"

template <typename T>
VectorGraphMap<T>::VectorGraphMap() {
  data.clear();
}

// L-value c-tor
template <typename T>
VectorGraphMap<T>::VectorGraphMap(std::map<T, double>& source) {
  std::cout << "VectorGraph (L-ctor)" << std::endl;
  data = source;
}

// R-value c-tor
template <typename T>
VectorGraphMap<T>::VectorGraphMap(std::map<T, double>&& source) {
//  std::cout << "VectorGraph (R-ctor)" << std::endl;
  data = std::move(source);
}

// D-tor
template <typename T>
VectorGraphMap<T>::~VectorGraphMap() {
  data.clear();
}

template <typename T>
void VectorGraphMap<T>::increase(T x, double v) {
  data[x]+=v;
}

template <typename T>
double VectorGraphMap<T>::get_and_erase(T x) {
  double ret = data[x];
  data.erase(x);
  return ret;
}

template <typename T>
double VectorGraphMap<T>::get(T x) {
  return data[x];
}

