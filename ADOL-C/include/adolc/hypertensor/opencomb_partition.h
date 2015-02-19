#ifndef OPENCOMB_PARTITION_H_
#define OPENCOMB_PARTITION_H_

#include "opencomb.h"

template <typename T>
class OpenCombPartition {
 public:
  OpenCombPartition() {};
  OpenCombPartition(const OpenCombMultiSet<T>& multi_set);
  bool init_partition(const OpenCombMultiSet<T>& multi_set);
  bool init_partition();
  bool has_next_partition();
  bool get_next_partition(std::vector<OpenCombMultiSet<T>>& multi_set_vec);
  void generate_partition(std::vector<OpenCombMultiSet<T>>& multi_set_vec);
 private:
  void load_data(const OpenCombMultiSet<T>& multi_set);
  std::vector<T> data;
  std::vector<int> k;
  std::vector<int> m;
  bool has_next;
  bool is_first;
  int size;
};

template <typename T>
void OpenCombPartition<T>::load_data(const OpenCombMultiSet<T>& multi_set) {
  data.clear();
  typename OpenCombMultiSet<T>::iterator iter = multi_set.begin();
  while(iter != multi_set.end()) {
    data.push_back(*iter);
    ++iter;
  }
  size = data.size();

}
template <typename T>
OpenCombPartition<T>::OpenCombPartition(
    const OpenCombMultiSet<T>& multi_set) {
  load_data(multi_set);
}

template <typename T>
bool OpenCombPartition<T>::init_partition(
    const OpenCombMultiSet<T>& multi_set) {
  load_data(multi_set);
  init_partition();
  return has_next;
}

template <typename T>
bool OpenCombPartition<T>::init_partition() {
  k.clear();
  m.clear();
  k.resize(size, 0);
  m.resize(size, 0);
  has_next = false;
  if (size != 0) {has_next = true;}
  return has_next;
}

template <typename T>
bool OpenCombPartition<T>::has_next_partition() {
  return has_next;
}

template <typename T>
void OpenCombPartition<T>::generate_partition(
    std::vector<OpenCombMultiSet<T> >& multi_set_vec) {
  multi_set_vec.clear();
  int i;
  int ssize = m[size - 1] - m[0] + 1;
  multi_set_vec.resize(ssize);
  for(i = 0; i < size; i++) {
    multi_set_vec[k[i]].put(data[i]);
  }
}

template <typename T>
bool OpenCombPartition<T>::get_next_partition(std::vector<OpenCombMultiSet<T> >& multi_set_vec) {
  int i;
  int j;
  generate_partition(multi_set_vec);
  for (i = size - 1; i > 0; i--) {
    if (k[i] <= m[i-1]) {
      k[i] = k[i] + 1;
      m[i] = (k[i]>m[i])?k[i]:m[i];
      for (j = i + 1; j < size; j++) {
        k[j] = k[0];
        m[j] = m[i];
      }
      has_next = true;
      return has_next;
    }
  }
  has_next = false;
  return has_next;
}

#endif // OPENCOMB_PARTITION_H_
