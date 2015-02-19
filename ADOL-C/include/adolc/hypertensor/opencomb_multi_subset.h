#ifndef OPENCOMB_MULTI_SUBSET_H_
#define OPENCOMB_MULTI_SUBSET_H_

#include <vector>
#include "opencomb.h"

template <typename T>
class OpenCombMultiSubset {
 public:
  OpenCombMultiSubset() {};
  OpenCombMultiSubset(const OpenCombMultiSet<T>& set);
  OpenCombMultiSubset(const OpenCombMultiSet<T>&& set);
  bool init_multi_subset();
  bool has_next_subset();
  bool get_next_subset(OpenCombMultiSet<T>& subset,
                       OpenCombMultiSet<T>& c_subset);
 private:
  std::vector<T> data;
  size_t size;
  size_t p;
  bool has_next;
  std::vector<int> index;
  OpenCombMultiSet<T> multi_set;
  OpenCombMultiSet<T> c_multi_set;
};

template <typename T>
OpenCombMultiSubset<T>::OpenCombMultiSubset(
    const OpenCombMultiSet<T>& set) {
  typename OpenCombMultiSet<T>::iterator iter = set.begin();
  while (iter != set.end()) {
    data.push_back(*iter);
    ++iter;
  }
}

template <typename T>
OpenCombMultiSubset<T>::OpenCombMultiSubset(
    const OpenCombMultiSet<T>&& set) {
  typename OpenCombMultiSet<T>::iterator iter = set.begin();
  while (iter != set.end()) {
    data.push_back(std::move(*iter));
    ++iter;
  }
}

template <typename T>
bool OpenCombMultiSubset<T>::init_multi_subset() {
  size = data.size();
  if (size == 0) {has_next = false;}
  index.clear();
  index.resize(size, 0);
  multi_set.clear();
  c_multi_set.clear();
  for(const T& item: data) {
    c_multi_set.put(item);
  }
  p = 0;
  has_next = true;
  return has_next;
}

template <typename T>
bool OpenCombMultiSubset<T>::has_next_subset() {
  return has_next;
}

template <typename T>
bool OpenCombMultiSubset<T>::get_next_subset(OpenCombMultiSet<T>& subset,
                                             OpenCombMultiSet<T>& c_subset) {
  subset = multi_set;
  c_subset = c_multi_set;
  for (int i = 0; i < p; ++i) {
    multi_set.remove(data[i]);
    c_multi_set.put(data[i]);
    index[i] = 0;
  }
  multi_set.put(data[p]);
  c_multi_set.remove(data[p]);
  index[p] = 1;
  p = 0;
  while (index[p] == 1 && p < size) {++p;}
  if (p < size) {has_next = true;} else {has_next = false;}
  return has_next;
}

#endif // OPENCOMB_MULTI_SUBSET_H_
