#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/hyper_derivative.h>
#include <adolc/hypertensor/hyper_tape.h>
#include <adolc/hypertensor/hyper_mpi_reverse.h>

#include <sys/time.h>
#include "mpi.h"
// options[0] = 1, first order;
//            = 2, second order;
//            = 3, third order;
int hyper_mpi_tensor(short tag,
                    int ndep,
                    int nindep,
                    const double* basepoint,
                    int order) {
  std::map<locint, HyperDerivative<locint> > global_gd;
  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::cout << "In hyper_mpi_tensor" << std::endl;
  hyper_tape(tag, ndep, nindep, basepoint, ind_map, hyper_index, hyper_value);
//  for(const locint& x: hyper_index) {
//    std::cout << x << std::endl;
//  }
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  struct timeval tv1, tv2;
  double time_elapsed;
  gettimeofday(&tv1, NULL);
  hyper_mpi_reverse(tag, hyper_index, hyper_value, global_gd);
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
  MPI_Barrier(MPI_COMM_WORLD);
  gettimeofday(&tv1, NULL);
  hyper_mpi_forward(global_gd);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout<<"Proc "<<myid<<" forward time: "<<time_elapsed<<std::endl;

/*
  iter = global_gd.begin();
  while(iter != global_gd.end()) {
    dep = iter->first;
    std::cout << "Derivatives for Dep " << dep << " : " << std::endl;
    global_gd[dep].debug();
    ++iter;
  }
*/
}
