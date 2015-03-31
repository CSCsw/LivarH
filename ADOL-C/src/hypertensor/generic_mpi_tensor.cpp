#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/generic_tape.h>
#include <adolc/hypertensor/opencomb_multi_set.h>
#include <adolc/hypertensor/generic_derivative_table.h>
#include <adolc/hypertensor/generic_mpi_reverse.h>

#ifdef ENABLE_GENERIC_MPI
#include <sys/time.h>
#include "mpi.h"
#define DEBUG_ID 99

#endif

int generic_mpi_tensor(short tag,
                       int dep,
                       int indep,
                       const double* basepoint,
                       int order,
                       int** nnz_p,
                       unsigned int**** indices_p,
                       double*** values_p) {

#ifdef ENABLE_GENERIC_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif

  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::map<locint, locint> dep_map;
  std::map<locint, std::set<locint> > r_live_set;
  std::map<locint, GenericDerivative<locint> > generic_derivative;
//  std::cout << "In generic_mpi_tensor" << std::endl;

  generic_tape(tag, dep, indep, basepoint,
               ind_map, dep_map, hyper_index, hyper_value);

/*
  if (myid == DEBUG_ID) {
    std::cout << hyper_index.size() << std::endl;
    for(const locint& loc : hyper_index) {
      std::cout << loc << std::endl;
    }
  }
*/
  struct timeval tv1, tv2;
  double time_elapsed;
#ifdef ENABLE_GENERIC_MPI
  gettimeofday(&tv1, NULL);
#endif

  special_derivative_table();
  generic_mpi_reverse(tag, order, hyper_index, hyper_value, r_live_set, generic_derivative);

#ifdef ENABLE_GENERIC_MPI
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout << "Proc "<<myid<<" reverse time: "<<time_elapsed << std::endl;  
#endif
  std::map<locint, GenericDerivative<locint> >::iterator iter;
  iter = generic_derivative.begin();

/*
//  if (myid == DEBUG_ID) {
    while(iter != generic_derivative.end() ){
      std::cout << "Derivative For : " << iter->first << std::endl;
      iter->second.debug();
      iter++;
    }
//  }
*/
#ifdef ENABLE_GENERIC_MPI
  gettimeofday(&tv1, NULL);
  generic_mpi_forward(order, r_live_set, generic_derivative);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout << "Proc "<<myid<<" forward time: "<<time_elapsed << std::endl;  
#endif
/*
  iter = generic_derivative.begin();
  if (myid == 0) {
    while(iter != generic_derivative.end() ){
      std::cout << "Derivative For : " << iter->first << std::endl;
      iter->second.debug();
      iter++;
    }
  }
*/

  int num_dep = generic_derivative.size();
  if (num_dep != dep) {
    std::cout << "Generic Tensor: Dependent check error! "
              << "Expected: " << dep << " Actural: " << num_dep << std::endl; 
  }
  if (num_dep == 0) {
    return 0;
  }

  int* nnz = (int*)malloc(sizeof(int) * num_dep);
  unsigned int*** indices = (unsigned int***)(malloc(sizeof(unsigned int**)*num_dep));
  double** values = (double**)malloc(sizeof(double*)*num_dep);
  int index = 0;
  iter = generic_derivative.begin();
  while(iter != generic_derivative.end()) {
    locint dummy_dep = dep_map[iter->first];
//    std::cout << "index = "<<iter->first << "dep = "<< dummy_dep << std::endl;
    int count = iter->second.get_size_for_order(order);
    nnz[dummy_dep] = count;
    indices[dummy_dep] = (unsigned int**)malloc(sizeof(unsigned int*) * count);
    values[dummy_dep] = (double*)malloc(sizeof(double*) * count);
    count = 0;
    typename GenericDerivative<locint>::iterator s_iter = iter->second.get_new_iterator();
    bool has_next = s_iter.init_iterator();
    OpenCombMultiSet<locint> s_set;
    double w;
    while(has_next) {
      s_iter.get_curr_pair(s_set, w);
      if (s_set.size() == order) {
        indices[dummy_dep][count] = (unsigned int*)malloc(sizeof(unsigned int) * order);
        int lc = 0;
        typename OpenCombMultiSet<locint>::iterator e_iter = s_set.begin();
        while(e_iter != s_set.end()) {
          indices[dummy_dep][count][lc] = ind_map[*e_iter];
          ++e_iter;
          ++lc;
        }
        values[dummy_dep][count] = w;
        ++count;
      }
      has_next = s_iter.move_to_next();
    }
    ++iter;
  }
  *nnz_p = nnz;
  *indices_p = indices;
  *values_p = values; 
}
