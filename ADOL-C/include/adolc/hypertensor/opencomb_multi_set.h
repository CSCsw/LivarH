#ifndef OPENCOMB_MULTI_SET_H_
#define OPENCOMB_MULTI_SET_H_

#include <set>
#include <iostream>

template <typename T>
class OpenCombMultiSet {
 public:
  OpenCombMultiSet();
  OpenCombMultiSet(const OpenCombMultiSet& rsh);
  OpenCombMultiSet(const OpenCombMultiSet&& rhs);
  OpenCombMultiSet(char* const buf, const int t_size);
  OpenCombMultiSet& operator = (const OpenCombMultiSet& rhs);
  OpenCombMultiSet& operator = (const OpenCombMultiSet&& rhs);
  inline bool operator == (const OpenCombMultiSet<T>& rhs) const;
  inline bool operator != (const OpenCombMultiSet<T>& rhs) const;
  inline bool operator < (const OpenCombMultiSet<T>& rhs) const;

  int put(T element);
  int remove(T element);
  int size() const;
  void debug() const;
  bool find(T target) const;
  int count(T target) const;
  void write_to_byte(char* buf) const;
  void clear();
  class iterator {
   public:
    iterator();
    iterator(const typename OpenCombMultiSet<T>::iterator& rhs);
    iterator(const typename OpenCombMultiSet<T>::iterator&& rhs);

    iterator(const typename std::multiset<T>::iterator& _iter);
    iterator(const typename std::multiset<T>::iterator&& _iter);
    ~iterator();

    OpenCombMultiSet<T>::iterator& operator = 
        (const OpenCombMultiSet<T>::iterator& rhs);
    OpenCombMultiSet<T>::iterator& operator = 
        (const OpenCombMultiSet<T>::iterator&& rhs);

    OpenCombMultiSet<T>::iterator& operator++();
    OpenCombMultiSet<T>::iterator operator++(int);
    const T& operator* ();

    inline bool operator == (const OpenCombMultiSet<T>::iterator& rhs) const;
    inline bool operator != (const OpenCombMultiSet<T>::iterator& rhs) const;
   
   private:
    typename std::multiset<T>::iterator iter;    
  };

  typename OpenCombMultiSet<T>::iterator begin() const;
  typename OpenCombMultiSet<T>::iterator end() const;

 private:
  std::multiset<T> data;
};

template <typename T>
OpenCombMultiSet<T>::iterator::iterator() {
}

template <typename T>
OpenCombMultiSet<T>::iterator::iterator(
    const typename std::multiset<T>::iterator& _iter) {
  this->iter = _iter;
}

template <typename T>
OpenCombMultiSet<T>::iterator::iterator(
    const typename std::multiset<T>::iterator&& _iter) {
  this->iter = std::move(_iter);
}

template <typename T>
OpenCombMultiSet<T>::iterator::iterator(
    const typename OpenCombMultiSet<T>::iterator& rhs) {
  (*this) = rhs;
}

template <typename T>
OpenCombMultiSet<T>::iterator::iterator(
    const typename OpenCombMultiSet<T>::iterator&& rhs) {
  (*this) = std::move(rhs);
}

template <typename T>
OpenCombMultiSet<T>::iterator::~iterator() {
}

template <typename T>
typename OpenCombMultiSet<T>::iterator&
    OpenCombMultiSet<T>::iterator::operator = (
        const OpenCombMultiSet<T>::iterator& rhs) {
  this->iter = rhs.iter;
}

template <typename T>
typename OpenCombMultiSet<T>::iterator&
    OpenCombMultiSet<T>::iterator::operator =(
        const OpenCombMultiSet<T>::iterator&& rhs) {
  if (this != &rhs) {
    this->iter = std::move(rhs.iter);
  }
}

template <typename T>
typename OpenCombMultiSet<T>::iterator&
    OpenCombMultiSet<T>::iterator::operator ++ () {
  ++(this->iter);
  return *this;
}

template <typename T>
typename OpenCombMultiSet<T>::iterator
    OpenCombMultiSet<T>::iterator::operator ++ (int) {
  OpenCombMultiSet<T>::iterator ret(this);
  ++(this->iter);
  return ret;
}

template <typename T>
inline bool OpenCombMultiSet<T>::iterator::operator == (
    const OpenCombMultiSet<T>::iterator& rhs) const {
  return (this->iter == rhs.iter);
}

template <typename T>
inline bool OpenCombMultiSet<T>::iterator::operator != (
    const OpenCombMultiSet<T>::iterator& rhs) const {
  return !((*this) == rhs);
}

template <typename T>
const T& OpenCombMultiSet<T>::iterator::operator * () {
  return (*iter);
}

template <typename T>
OpenCombMultiSet<T>::OpenCombMultiSet() {
  data.clear();
}

template <typename T>
OpenCombMultiSet<T>::OpenCombMultiSet(const OpenCombMultiSet& rhs) {
//  std::cout << "OpenCombMultiSet L-ctor" << std::endl;
  data = rhs.data;
}

template <typename T>
OpenCombMultiSet<T>::OpenCombMultiSet(const OpenCombMultiSet&& rhs) {
//  std::cout << "OpenCombMultiSet R-ctor" << std::endl;
  data = std::move(rhs.data);
}

template <typename T>
OpenCombMultiSet<T>& OpenCombMultiSet<T>::operator = (
    const OpenCombMultiSet<T>& rhs) {
//  std::cout << "OpenCombMultiSet L-assign" << std::endl;
  this->data = rhs.data;
  return *this;
}

template <typename T>
OpenCombMultiSet<T>& OpenCombMultiSet<T>::operator = (
    const OpenCombMultiSet<T>&& rhs) {
//  std::cout << "OpenCombMultiSet R-assign" << std::endl;
  this->data = std::move(rhs.data);
  return *this;
}

template <typename T>
int OpenCombMultiSet<T>::put(T element) {
  data.insert(element);
  return data.count(element);
}

template <typename T>
int OpenCombMultiSet<T>::remove(T element) {
  typename std::multiset<T>::iterator hit(data.find(element));
  if (hit != data.end()) {
    data.erase(hit);
  }
  return data.count(element);
}

template <typename T>
int OpenCombMultiSet<T>::size() const {
  return data.size();
}

template <typename T>
void OpenCombMultiSet<T>::write_to_byte(char* buf) const {
  char* p = buf;
  typename std::multiset<T>::iterator iter;
  iter = data.begin();
  while(iter != data.end()) {
    *((T*)p) = (*iter);
    p += sizeof(T);
    ++iter;
  }
}

template <typename T>
OpenCombMultiSet<T>::OpenCombMultiSet(char* const buf, const int t_size) {
  char* p = buf;
  data.clear();
  for(int i = 0; i < t_size; i++) {
    data.insert(*((T*)p));
    p += sizeof(T);
  }
}

template <typename T>
bool OpenCombMultiSet<T>::find(T target) const {
  if (data.find(target) != data.end()) {
    return true;
  }
  return false;
}

template <typename T>
int OpenCombMultiSet<T>::count(T target) const {
  return data.count(target);
}

template <typename T>
void OpenCombMultiSet<T>::debug() const {
  typename std::multiset<T>::iterator iter;
  std::cout << "{ ";
  for(iter = data.begin(); iter != data.end(); ++iter) {
    std::cout << (*iter) << " ";
  }
  std::cout << "}";
}

template <typename T>
void OpenCombMultiSet<T>::clear() {
  data.clear();
}

template <typename T>
typename OpenCombMultiSet<T>::iterator OpenCombMultiSet<T>::begin() const {
  return OpenCombMultiSet<T>::iterator(data.begin());
}

template <typename T>
typename OpenCombMultiSet<T>::iterator OpenCombMultiSet<T>::end() const {
  return OpenCombMultiSet<T>::iterator(data.end());
}

// kind of redundent because we have < operator.
// But it's more efficient to keep it.
template <typename T>
inline bool OpenCombMultiSet<T>::operator == (
    const OpenCombMultiSet<T>& rhs) const {
  typename std::multiset<T>::iterator iter1;
  typename std::multiset<T>::iterator iter2;
  if (this->data.size() != rhs.data.size()) {
    return false;
  }
  iter1 = this->data.begin();
  iter2 = rhs.data.begin();
  while (iter1 != this->data.end() && iter2 != rhs.data.end()) {
    if (*iter1 != *iter2) {
      return false;
    }
    ++iter1;
    ++iter2;
  }
  return true;
}

template <typename T> 
inline bool OpenCombMultiSet<T>::operator != (
    const OpenCombMultiSet<T>& rhs) const{
  return !((*this) == rhs);
}

template <typename T> 
inline bool OpenCombMultiSet<T>::operator < (
    const OpenCombMultiSet<T>& rhs) const{
  typename std::multiset<T>::iterator iter1;
  typename std::multiset<T>::iterator iter2;
  if (this->data.size() < rhs.data.size()) {
    return true;
  } else if (this->data.size() > rhs.data.size()) {
    return false;
  }
  iter1 = this->data.begin();
  iter2 = rhs.data.begin();
  while (iter1 != this->data.end() && iter2 != rhs.data.end()) {
    if (*iter1 < *iter2) {
      return true;
    } else if (*iter1 > *iter2) {
      return false;
    }
    ++iter1;
    ++iter2;
  }
  return false;
}

#endif // OPENCOMB_MULTI_SET_H_
