#include <stdlib.h>
#include <stdio.h>
#include "edge_b_tree.h"

#define EdgeBTreeBlockSizeT 20
#define EdgeBTreeBlockSize 41

EdgeBTreeBlock::EdgeBTreeBlock() {
// TODO: Allocate a continguous memory;
  ind = (locint*) malloc(sizeof(locint) * EdgeBTreeBlockSize);
  weight = (double*) malloc(sizeof(double) * EdgeBTreeBlockSize);
  child_block = (EdgeBTreeBlock**) malloc(sizeof(EdgeBTreeBlock*) * (EdgeBTreeBlockSize + 1));
  size = 0;
  is_leaf = false;
}

EdgeBTreeBlock::~EdgeBTreeBlock() {
//  printf("deleting %x:\n", this);
//  fflush(stdout);
  for(size_t i = 0; i <= size; ++i) {
    if (!is_leaf) {
      delete child_block[i];
    }
  }
//  printf("done %x\n", this);
//  fflush(stdout);
  free(ind);
  free(weight);
  free(child_block);
}

void EdgeBTreeBlock::AppendToArray(size_t *len_array, locint *ind_array, double *weight_array) {
  for(size_t i = 0; i < size; ++i) {
    ind_array[*len_array] = ind[i];
    weight_array[*len_array] = weight[i];
    ++(*len_array);
  }
  if (!is_leaf) {
    for(size_t i = 0; i <= size; i++) {
      child_block[i]->AppendToArray(len_array, ind_array, weight_array);
    }
  }
}

void EdgeBTreeBlock::PrintBlock() {
  if (!is_leaf) {
    for(size_t i = 0; i<=size; i++) {
      child_block[i]->PrintBlock();
    }
  }
  printf("B[%x]:",this);
  for(size_t i = 0; i<size; ++i) {
    if (!is_leaf) {
      printf("%x, ", child_block[i]);
    }
    printf("<%d, %.5f>, ", ind[i], weight[i]);
  }
  if (!is_leaf) {
    printf("%x", child_block[size]);
  }
  printf("\n");
}

EdgeBTree::EdgeBTree() {
  root_block = new EdgeBTreeBlock();
  root_block->is_leaf = true;
  size = 0;
}

EdgeBTree::~EdgeBTree() {
  delete root_block;
}

void EdgeBTree::update(locint index, double w) {
  double *weight = find(index);
  *weight+=w;
}

void EdgeBTree::erase(locint index) {
  double *weight = find(index);
  *weight = 0.0;
}

double* EdgeBTree::find(locint index) {
  EdgeBTreeBlock *prev_block = NULL;
  size_t prev_ind = 0;
  EdgeBTreeBlock *curr_block = root_block;

  while (true) {
    if (curr_block->size == EdgeBTreeBlockSize) {
//      printf("splitting block[%x]:\n",curr_block);
//      fflush(stdout);
      //split two blocks;
      EdgeBTreeBlock *new_block = new EdgeBTreeBlock();
      size_t i,j;
      for(i = EdgeBTreeBlockSizeT+1, j=0; i < EdgeBTreeBlockSize; ++i, ++j) {
        new_block->ind[j] = curr_block->ind[i];
        curr_block->ind[i] = 0;
        new_block->weight[j] = curr_block->weight[i];
        curr_block->weight[i] = 0.0;
        new_block->child_block[j] = curr_block->child_block[i];
        curr_block->child_block[i] = NULL;
      }
      new_block->child_block[EdgeBTreeBlockSizeT] = curr_block->child_block[EdgeBTreeBlockSize];
      curr_block->child_block[EdgeBTreeBlockSize] = NULL;
      new_block->is_leaf = curr_block->is_leaf;

      new_block->size = EdgeBTreeBlockSizeT;
      curr_block->size = EdgeBTreeBlockSizeT;
//      curr_block->PrintBlock();
//      new_block->PrintBlock();
//      fflush(stdout);
      // put mid to prev block
      locint mid_ind = curr_block->ind[EdgeBTreeBlockSizeT];
      curr_block->ind[EdgeBTreeBlockSizeT] = 0;
      double mid_weight = curr_block->weight[EdgeBTreeBlockSizeT];
      curr_block->weight[EdgeBTreeBlockSizeT] = 0.0;

//      printf("mid_ind=%d, mid_weight=%.5f\n", mid_ind, mid_weight);
//      fflush(stdout);
      if (prev_block == NULL) {
        //Build a new root;
        root_block = new EdgeBTreeBlock();
        root_block->ind[0] = mid_ind;
        root_block->weight[0] = mid_weight;
        root_block->child_block[0] = curr_block;
        root_block->child_block[1] = new_block;
        root_block->size = 1;
        prev_block = root_block;
//        root_block->PrintBlock();
//        fflush(stdout);
      } else {
        size_t i;
        for(i = prev_block->size; i > prev_ind && i>0; --i) {
          prev_block->ind[i]=prev_block->ind[i-1];
          prev_block->weight[i]=prev_block->weight[i-1];
          prev_block->child_block[i+1] = prev_block->child_block[i];
        }
        prev_block->size++;
        prev_block->ind[prev_ind] = mid_ind;
        prev_block->weight[prev_ind] = mid_weight;
        prev_block->child_block[prev_ind+1] = new_block;
      }
      if (mid_ind == index) {
        return &(prev_block->weight[prev_ind]);
      } else if (mid_ind < index) {
        curr_block = new_block;
      }
    }
    prev_ind=0;
    while(prev_ind < curr_block->size && curr_block->ind[prev_ind] < index) {
      prev_ind++;
    }
    if (prev_ind != curr_block->size && curr_block->ind[prev_ind] == index) {
      return &(curr_block->weight[prev_ind]);
    }
    if (curr_block->is_leaf) {
      size_t i;
      for(i = curr_block->size; i > prev_ind && i > 0; i--) {
//        printf("i=%d\n",i);
//        fflush(stdout);
        curr_block->ind[i] = curr_block->ind[i-1];
        curr_block->weight[i] = curr_block->weight[i-1];
      }
      curr_block->size++;
      curr_block->ind[prev_ind] = index;
      curr_block->weight[prev_ind] = 0.0;
      size++;

      return &(curr_block->weight[prev_ind]);
    } else {
      prev_block = curr_block;
      curr_block = prev_block->child_block[prev_ind];
    }
  }
}

void EdgeBTree::ToArray(size_t *len, locint **ind, double **weight) {
  printf("size=%d\n", size);
  fflush(stdout);
  (*ind) = (locint*) malloc(sizeof(locint) * size);
  (*weight) = (double*) malloc(sizeof(double) * size);
  *len=0;
  root_block->AppendToArray(len, *ind, *weight);
  if (*len != size) {
    fprintf(stderr, "FATAL ERROR: EdgeBTree Broken!\n");
    exit(-1);
  }
}

void EdgeBTree::Print() {
  root_block->PrintBlock();
}

