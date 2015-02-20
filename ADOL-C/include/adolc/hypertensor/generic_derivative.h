#ifndef GENERIC_DERIVATIVE_H_
#define GENERIC_DERIVATIVE_H_

#include <vector>
#include <map>

#include "opencomb.h"

template <typename T>
class GenericDerivative {
 public:
  GenericDerivative(): d(0) {};
  GenericDerivative(int order): d(order) {info.resize(d);};
  void clear();
  void increase(const OpenCombMultiSet<T>& set, double v);
  double get(const OpenCombMultiSet<T>& set);
  void debug() const;
  void find_and_erase(T target, GenericDerivative& gd);
  class iterator {
   public:
    iterator() {};
    iterator(std::vector<std::map<OpenCombMultiSet<T>, double> >* _p_info);
    bool init_iterator();
    void get_curr_pair(OpenCombMultiSet<T>& multi_set, double& w);
    bool move_to_next();
    typename GenericDerivative<T>::iterator& operator = (
        const GenericDerivative<T>::iterator& rhs);
    
    bool has_next;
   private:
    std::vector<std::map<OpenCombMultiSet<T>, double> >* p_info;
    typename std::vector<std::map<OpenCombMultiSet<T>, double> >::iterator v_iter;
    typename std::map<OpenCombMultiSet<T>, double>::iterator m_iter;
  };

  typename GenericDerivative<T>::iterator get_new_iterator();

 private:
  const int d; // The highest order
  std::vector<std::map<OpenCombMultiSet<T>, double> > info;
};

template <typename T>
void GenericDerivative<T>::clear() {
  info.clear();
  info.resize(d);
}

template <typename T>
void GenericDerivative<T>::increase(const OpenCombMultiSet<T>& set, double v) {
  int size = set.size() - 1;
  if (size >= d) {
    std::cout << "Max order = " << d
              << ", Input size = " << size + 1 << std::endl;
    return;
  }
  info[size][set] += v;
}

template <typename T>
double GenericDerivative<T>::get(const OpenCombMultiSet<T>& set) {
  int size = set.size() - 1;
  if (size >= d) {
    std::cout << "Max order = " << d
              << ", Input size = " << size + 1 << std::endl;
    return 0;
  }
  if (info[size].find(set) == info[size].end()) {
    return 0;
  }
  return info[size][set];
}

template <typename T>
void GenericDerivative<T>::debug() const {
  typename std::vector<std::map<OpenCombMultiSet<T>, double> >::const_iterator v_iter;
  typename std::map<OpenCombMultiSet<T>, double >::const_iterator m_iter;
  v_iter = info.begin();
  while(v_iter != info.end()) {
    m_iter = (*v_iter).begin();
    while(m_iter != (*v_iter).end()) {
      std::cout << "T ";
      m_iter->first.debug();
      std::cout << " = " << m_iter->second << std::endl;
      ++m_iter;
    }
    ++v_iter;
  }
}


template <typename T>
typename GenericDerivative<T>::iterator GenericDerivative<T>::get_new_iterator() {
  typename GenericDerivative<T>::iterator iter(&info);
  return iter; 
}
template <typename T>
GenericDerivative<T>::iterator::iterator(
  std::vector<std::map<OpenCombMultiSet<T>, double> >* _p_info) {
  p_info = _p_info;
}

template <typename T>
typename GenericDerivative<T>::iterator& GenericDerivative<T>::iterator::operator = (
    const GenericDerivative<T>::iterator& rhs) {
  this->p_info = rhs.p_info;
  this->v_iter = rhs.v_iter;
  this->m_iter = rhs.m_iter;
  this->has_next = rhs.has_next;
  return *this;
}

template <typename T>
bool GenericDerivative<T>::iterator::init_iterator() {
  v_iter = p_info->begin();
  while (v_iter != p_info->end()) {
    m_iter = (*v_iter).begin();
    while (m_iter != (*v_iter).end()){
      has_next = true;
      return has_next;
    }
    ++v_iter;
  }
  has_next = false;
  return has_next;
}

template <typename T>
void GenericDerivative<T>::iterator::get_curr_pair(
    OpenCombMultiSet<T>& multi_set, double& w) {
  if (has_next) {
    multi_set = m_iter->first;
    w = m_iter->second;
  }
}

template <typename T>
bool GenericDerivative<T>::iterator::move_to_next() {
  ++m_iter;
  while (v_iter != p_info->end()) {
    if (m_iter != (*v_iter).end()) {
      has_next = true;
      return has_next;
    } else {
      ++v_iter;
      if (v_iter != p_info->end()) {
        m_iter = (*v_iter).begin();
      }
    }
  }
  has_next = false;
  return has_next;
}

template <typename T>
void GenericDerivative<T>::find_and_erase(T target,
                                        GenericDerivative<T>& gd) {
  gd.clear();
  typename std::vector<std::map<OpenCombMultiSet<T>, double> >::iterator v_iter;
  typename std::map<OpenCombMultiSet<T>, double >::iterator m_iter;
  v_iter = info.begin();
  while(v_iter != info.end()) {
    m_iter = (*v_iter).begin();
    while(m_iter != (*v_iter).end()) {
      if (m_iter->first.find(target)) {
        gd.increase(m_iter->first, m_iter->second);
        m_iter->second = 0.0;
        m_iter = (*v_iter).erase(m_iter);
      } else {
        ++m_iter;
      }
    }
    ++v_iter;
  }
}
#endif // GENERIC_DERIVATIVE_H_
