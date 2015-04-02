#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/hyper_derivative.h>
#include <adolc/hypertensor/generic_tape.h>
#include <adolc/hypertensor/hyper_mpi_reverse.h>

#include <sys/time.h>

#ifdef ENABLE_GENERIC_MPI
#include "mpi.h"
#endif
// options[0] = 1, first order;
//            = 2, second order;
//            = 3, third order;
int hyper_mpi_tensor(short tag,
                    int ndep,
                    int nindep,
                    const double* basepoint,
                    int order,
                    int** nnz_p,
                    unsigned int**** indices_p,
                    double*** values_p) {
  std::map<locint, HyperDerivative<locint> > global_gd;
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::map<locint, locint> dep_map;
  std::cout << "In hyper_mpi_tensor" << std::endl;
//  hyper_tape(tag, ndep, nindep, basepoint, ind_map, hyper_index, hyper_value);
  generic_tape(tag, ndep, nindep, basepoint,
               ind_map, dep_map, hyper_index, hyper_value);
//  for(const locint& x: hyper_index) {
//    std::cout << x << std::endl;
//  }
  int myid = 0;
#ifdef ENABLE_GENERIC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
  struct timeval tv1, tv2;
  double time_elapsed;
  gettimeofday(&tv1, NULL);
  hyper_mpi_reverse(tag, order, hyper_index, hyper_value, global_gd);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout<<"Proc "<<myid<<" reverse time: "<<time_elapsed<<std::endl;
  typename std::map<locint, HyperDerivative<locint> >::iterator iter;
  locint dep;
/*
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    dep = iter->first;
    std::cout << "Derivatives for Dep " << dep << " : " << std::endl;
    global_gd[dep].debug();
    ++iter;
  }
*/
#ifdef ENABLE_GENERIC_MPI
  gettimeofday(&tv1, NULL);
  hyper_mpi_forward(order, global_gd);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout<<"Proc "<<myid<<" forward time: "<<time_elapsed<<std::endl;
#endif
/*
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    dep = iter->first;
    std::cout << "Derivatives for Dep " << dep << " : " << std::endl;
    global_gd[dep].debug();
    ++iter;
  }
*/
  int num_dep = global_gd.size();
  if (num_dep != ndep) {
    std::cout << "Hyper Tensor: Dependent check error! "
              << "Expected: " << dep << " Actural: " << num_dep << std::endl;
  }
  if (num_dep == 0) {
    return 0;
  }


  int* nnz = (int*)malloc(sizeof(int) * num_dep);
  unsigned int*** indices = (unsigned int***)(malloc(sizeof(unsigned int**)*num_dep));
  double** values = (double**)malloc(sizeof(double*)*num_dep);
  int index = 0;
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    locint dummy_dep = dep_map[iter->first];
    iter->second.debug();
//    std::cout << "index = "<<iter->first << "dep = "<< dummy_dep << std::endl;
    int count = 0;
    if (order == 2) {
      count = iter->second.hessian->get_size();
    } else if (order == 3) {
      count = iter->second.tensor->get_size();
    }
    nnz[dummy_dep] = count;
    indices[dummy_dep] = (unsigned int**)malloc(sizeof(unsigned int*) * count);
    values[dummy_dep] = (double*)malloc(sizeof(double*) * count);
    count = 0;
    if (order == 2) {
      typename MatrixGraph<locint>::iterator* h_iter =
          iter->second.hessian->get_iterator();
      bool has_next = h_iter->reset();
      locint x, y;
      double w;
      while(has_next) {
        has_next = h_iter->get_next(x, y, w);
        indices[dummy_dep][count] = (unsigned int*)malloc(sizeof(unsigned int) * order);
        indices[dummy_dep][count][0] = ind_map[x];
        indices[dummy_dep][count][1] = ind_map[y];
        values[dummy_dep][count] = w;
        ++count;
      }
    } else if (order == 3) {
      bool has_next = iter->second.tensor->reset();
      locint x, y, z;
      double w;
      while(has_next) {
        has_next = iter->second.tensor->get_next(x, y, z, w);
        indices[dummy_dep][count] = (unsigned int*)malloc(sizeof(unsigned int) * order);
        indices[dummy_dep][count][0] = ind_map[x];
        indices[dummy_dep][count][1] = ind_map[y];
        indices[dummy_dep][count][2] = ind_map[z];
        values[dummy_dep][count] = w;
        ++count;
      }
    }
    ++iter;
  }
  *nnz_p = nnz;
  *indices_p = indices;
  *values_p = values;
}
