
#include <stdlib.h>
#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hessian/edge_local_graph.h>

EdgeLocalGraph::EdgeLocalGraph() {
  loc = (locint*) malloc(sizeof(locint) * EDGE_LOCAL_SIZE);
  adjoints = (double*) malloc(sizeof(double) * EDGE_LOCAL_SIZE);
  size_array = (size_t*) malloc(sizeof(size_t) * EDGE_LOCAL_SIZE);
  loc_array = (locint**) malloc(sizeof(locint*) * EDGE_LOCAL_SIZE);
  hessian = (double**) malloc(sizeof(double*) * EDGE_LOCAL_SIZE);
  locint* l_tmp = (locint*) malloc(sizeof(locint) * EDGE_LOCAL_SIZE * EDGE_LOCAL_SIZE);
  double* d_tmp = (double*) malloc(sizeof(double) * EDGE_LOCAL_SIZE * EDGE_LOCAL_SIZE);
  for(size_t i = 0; i < EDGE_LOCAL_SIZE; i++) {
    loc_array[i] = l_tmp;
    hessian[i] = d_tmp;
    l_tmp = l_tmp + EDGE_LOCAL_SIZE;
    d_tmp = d_tmp + EDGE_LOCAL_SIZE;
    size_array[i] = 0;
  }
  size = 0;
  max_size = EDGE_LOCAL_SIZE;
}

EdgeLocalGraph::~EdgeLocalGraph() {
  free(loc);
  free(adjoints);
  free(size_array);
  free(loc_array[0]);
  free(hessian[0]);
  free(loc_array);
  free(hessian);
}

void EdgeLocalGraph::reset() {
  for(size_t i = 0; i < size; ++i) {
    size_array[i] = 0;
  }
  size = 0;
}

void EdgeLocalGraph::ExpandSize() {
// Initialize new arrays
  max_size = max_size * 2;
  locint* n_loc = (locint*) malloc(sizeof(locint) * max_size);
  double* n_adjoints = (double*) malloc(sizeof(double) * max_size);
  size_t* n_size_array = (size_t*) malloc(sizeof(size_t) * max_size);
  locint** n_loc_array = (locint**) malloc(sizeof(locint*) * max_size);
  double** n_hessian = (double**) malloc(sizeof(double*) * max_size);
  locint* n_l_tmp = (locint*) malloc(sizeof(locint) * max_size * max_size);
  double* n_d_tmp = (double*) malloc(sizeof(double) * max_size * max_size);
  for(size_t i = 0; i < max_size; i++) {
    n_loc_array[i] = n_l_tmp;
    n_hessian[i] = n_d_tmp;
    n_l_tmp = n_l_tmp + max_size;
    n_d_tmp = n_d_tmp + max_size;
    n_size_array[i] = 0;
  }
// Copy values from old to new
  for(size_t i = 0; i < size; i++) {
    n_loc[i] = loc[i];
    n_adjoints[i] = adjoints[i];
    n_size_array[i] = size_array[i];
    for(size_t j = 0; j < size_array[i]; j++) {
      n_loc_array[i][j] = loc_array[i][j];
      n_hessian[i][j] = hessian[i][j];
    }
  }
// Deallocate old arrays
  free(loc); loc = n_loc;
  free(adjoints); adjoints = n_adjoints;
  free(size_array); size_array = n_size_array;
  free(loc_array[0]);
  free(loc_array); loc_array = n_loc_array;
  free(hessian[0]);
  free(hessian); hessian = n_hessian;
}
size_t EdgeLocalGraph::AddLiveVar(locint ind) {
  if (ind == NULLLOC) {
    return NULLLOC;
  }
  size_t ret = find(ind, loc, size);
  if (ret == NULLLOC) {
    if (size == max_size) {
      ExpandSize();
    }
    ret = size++;
    loc[ret] = ind;
    adjoints[ret] = 0.0;
    size_array[ret] = 0;
  }
  return ret;
}

void EdgeLocalGraph::erase(size_t ind) {
  if (ind != NULLLOC) {
    adjoints[ind] = 0.0;
    size_array[ind] = 0;
  }
}

size_t EdgeLocalGraph::find(locint ind, locint* array, size_t curr_len) {
  size_t ret = 0;
  while( ret < curr_len && array[ret] != ind) {
    ++ret;
  }
  if (ret == curr_len) {
    return NULLLOC;
  }
  return ret;
}

void EdgeLocalGraph::insert(size_t row, locint ind, double w) {
  size_t ret = find(ind, loc_array[row], size_array[row]);
  if (ret == NULLLOC) {
    ret = size_array[row]++;
    loc_array[row][ret] = ind;
    hessian[row][ret] = 0.0;
  }
  hessian[row][ret] += w;
}

void EdgeLocalGraph::Print() {
  printf("Local size = %d\n", size);
  for(size_t i = 0; i < size; i++) {
    printf("A[%d] (%d) <%.5f> :", loc[i], size_array[i], adjoints[i]);
    for(size_t j = 0; j < size_array[i] ; ++j) {
      printf("[%d]=%.5f, ", loc_array[i][j], hessian[i][j]);
    }
    printf("\n");
  }
  fflush(stdout);
}
