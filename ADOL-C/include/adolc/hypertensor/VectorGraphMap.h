#ifndef __VECTOR_GRAPH_MAP_H__
#define __VECTOR_GRAPH_MAP_H__

#include <map>

#include "VectorGraph.h"

template <typename T>
class VectorGraphMap : public VectorGraph<T> {
 public:
  VectorGraphMap();
  VectorGraphMap(std::map<T, double>& source);
  VectorGraphMap(std::map<T, double>&& source);

  ~VectorGraphMap();

  void increase(T x, double v);
  double get_and_erase(T x);
  double get(T x);
  int get_size();

  void debug();

  class iterator : public VectorGraph<T>::iterator {
   public:
    iterator(std::map<T, double>* _data_p):_data(_data_p) {};
    virtual ~iterator();
    bool reset();
    bool get_next(T& x, double& w);
   private:
    const std::map<T, double>* const _data;
    typename std::map<T, double>::const_iterator iter;
  };

  typename VectorGraph<T>::iterator* get_iterator();

 private:
  std::map<T, double> data;
};

template <typename T>
typename VectorGraph<T>::iterator* VectorGraphMap<T>::get_iterator() {
  typename VectorGraph<T>::iterator* ret = new VectorGraphMap<T>::iterator(&data);
  return ret;
}

template <typename T>
VectorGraphMap<T>::iterator::~iterator() {

}

template <typename T>
bool VectorGraphMap<T>::iterator::reset() {
  iter = _data->begin();
  if (iter == _data->end()) {
    return false;
  }
  return true;
}

template <typename T>
bool VectorGraphMap<T>::iterator::get_next(T&x, double& w){
  x = iter->first;
  w = iter->second;
  ++iter;
  if (iter == _data->end()) {
    return false;
  }
  return true;
}

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
  if (v == 0.0) {return;}
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

template <typename T>
int VectorGraphMap<T>::get_size() {
  return data.size();
}

template <typename T>
void VectorGraphMap<T>::debug() {
  typename std::map<T, double>::iterator t_iter;
  t_iter = data.begin();
  while(t_iter != data.end()) {
    std::cout << "A[" << t_iter->first << "]=" << t_iter->second << std::endl;
    ++t_iter;
  }
}

#endif // __VECTOR_GRAPH_MAP_H__
