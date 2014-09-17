#ifndef _EDGE_B_TREE_H_
#define _EDGE_B_TREE_H_
#include <adolc/internal/adolc_settings.h>

class EdgeBTreeBlock {
 public:
  // default ctor
  EdgeBTreeBlock();
  
  // d-tor
  ~EdgeBTreeBlock();

  // AppendToArray
  void AppendToArray(size_t *len_array, locint *ind_array, double *weight_array);

  // PrintBlock
  void PrintBlock();

  size_t size;
  bool is_leaf;
  locint *ind;
  double *weight;
  EdgeBTreeBlock **child_block;
};

class EdgeBTree {
 public:
  // default ctor, initialize the root_block;
  EdgeBTree();

  // dtor
  ~EdgeBTree();

  // update an entry, for non-existing entries, create one with 0.0
  void update(locint index, double w);

  // find/create an entry, return the reference to the weight location;
  double* find(locint index);

  void erase(locint index);

  void Print();

  void ToArray(size_t *len_array, locint *ind_array, double *weight_array);
  void ToArrayAlloc(size_t *len_array, locint **ind_array, double **weight_array);

 private:
  EdgeBTreeBlock *root_block;
  size_t size;

};

#endif
