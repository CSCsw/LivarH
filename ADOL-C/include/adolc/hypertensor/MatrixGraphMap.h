#ifndef __MATRIX_GRAPH_MAP_H__
#define __MATRIX_GRAPH_MAP_H__

#include <map>

#include <adolc/hypertensor/hyper_common.h>

#include "VectorGraph.h"
#include "VectorGraphMap.h"
#include "MatrixGraph.h"

template <typename T>
class MatrixGraphMap : public MatrixGraph<T> {
 public:
  MatrixGraphMap();
  MatrixGraphMap(std::map<T, std::map<T, double> >& source);
  MatrixGraphMap(std::map<T, std::map<T, double> >&& source);

  ~MatrixGraphMap();
  void increase(T x, T y, double v);
  VectorGraph<T>* get_and_erase(T x);
  VectorGraph<T>* get(T x);
  int get_size();
  void debug();
  typename MatrixGraph<T>::iterator* get_iterator();

  class iterator : public MatrixGraph<T>::iterator {
   public:
    iterator(std::map<T, std::map<T, double> >* _data_p): _data(_data_p) {};
    virtual ~iterator();
    bool reset();
    bool get_next(T& x, T& y, double& w);

   private:
    const std::map<T, std::map<T, double> >* const _data;
    typename std::map<T, std::map<T, double> >::const_iterator iter;
    typename std::map<T, double>::const_iterator iter2;
  };

 private:
  std::map<T, std::map<T, double> > data;

};

template <typename T>
MatrixGraphMap<T>::MatrixGraphMap() {
  data.clear();
}

// L-value ctor, copy assignment
template <typename T>
MatrixGraphMap<T>::MatrixGraphMap(std::map<T, std::map<T, double>>& source) {
//  std::cout << "MatrixGraphMap (L-ctor)" << std::endl;
  data = source;
}

// R-value c-tor, move assignment
template <typename T>
MatrixGraphMap<T>::MatrixGraphMap(std::map<T, std::map<T, double>>&& source) {
//  std::cout << "MatrixGraphMap (R-ctor)" << std::endl;
  data = std::move(source);
}

template <typename T>
MatrixGraphMap<T>::~MatrixGraphMap() {
  data.clear();
}

template <typename T>
void MatrixGraphMap<T>::increase(T x, T y, double v) {
  if (v == 0.0) {return;}
  MAX_SWAP(x,y);
  data[x][y]+=v;
}

template <typename T>
VectorGraph<T>* MatrixGraphMap<T>::get_and_erase(T x) {
  VectorGraph<T>* ret = new VectorGraphMap<T>(std::move(data[x]));
  data.erase(x);
  return ret;
}

template <typename T>
VectorGraph<T>* MatrixGraphMap<T>::get(T x) {
  VectorGraph<T>* ret = new VectorGraphMap<T>(data[x]);
  return ret;
}

template <typename T>
typename MatrixGraph<T>::iterator* MatrixGraphMap<T>::get_iterator() {
  typename MatrixGraph<T>::iterator* ret = new MatrixGraphMap<T>::iterator(&data);
  return ret;
}

template <typename T>
MatrixGraphMap<T>::iterator::~iterator() {

}

template <typename T>
bool MatrixGraphMap<T>::iterator::reset() {
  iter = _data->begin();
  while(iter != _data->end()) {
    iter2 = iter->second.begin();
    if (iter2 != iter->second.end()) {
      return true;
    }
    ++iter;
  }
  return false;
}

template <typename T>
bool MatrixGraphMap<T>::iterator::get_next(T& x, T& y, double& w) {
  x = iter->first;
  y = iter2->first;
  w = iter2->second;
  ++iter2;
  while(iter != _data->end()) {
    if (iter2 != iter->second.end()) {
      return true;
    }
    ++iter;
    iter2 = iter->second.begin();
  }
  return false;
}

template <typename T>
int MatrixGraphMap<T>::get_size() {
  int size_count = 0;
  typename std::map<T, std::map<T, double> >::iterator t_iter;
  t_iter = data.begin();
  while (t_iter != data.end()) {
    size_count += t_iter->second.size();
    ++t_iter;
  }
  return size_count;
}

template <typename T>
void MatrixGraphMap<T>::debug() {
  typename std::map<T, std::map<T, double> >::iterator t_iter;
  t_iter = data.begin();
  while (t_iter != data.end()) {
    typename std::map<T, double>::iterator t_iter2;
    t_iter2 = t_iter->second.begin();
    while (t_iter2 != t_iter->second.end()) {
      std::cout << "H[" << t_iter->first << "," <<t_iter2->first << "] = "
                << t_iter2->second << std::endl;
      ++t_iter2;
    }
    ++t_iter;
  }
}
#endif // __MATRIX_GRAPH_MAP_H__
