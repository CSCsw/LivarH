
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
}

EdgeLocalGraph::~EdgeLocalGraph() {
  free(loc);
  free(adjoints);
  free(size_array);
  free(loc_array[0]);
  free(hessian[0]);
  free(loc_array);
  free(hessian);
  size = 0;
}

void EdgeLocalGraph::reset() {
  for(size_t i = 0; i < size; ++i) {
    size_array[i] = 0;
  }
  size = 0;
}

size_t EdgeLocalGraph::AddLiveVar(locint ind) {
  if (ind == NULLLOC) {
    return NULLLOC;
  }
  size_t ret = find(ind, loc, size);
  if (ret == NULLLOC) {
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
