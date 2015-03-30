#include <vector>
#include <map>
#include <iostream>

#include <stdio.h>

#include <adolc/adolc.h>
#include <adolc/hypertensor/generic_tape.h>
#include <adolc/hypertensor/opencomb.h>
#include <adolc/hypertensor/generic_derivative_table.h>
#include <adolc/hypertensor/generic_mpi_reverse.h>
#include <sys/time.h>
#include "mpi.h"

#define DEBUG_ID 99
int generic_mpi_tensor(short tag,
                       int dep,
                       int indep,
                       const double* basepoint,
                       int d) {
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::vector<locint> hyper_index;
  std::vector<double> hyper_value;
  std::map<locint, locint> ind_map;
  std::map<locint, std::set<locint> > r_live_set;
  std::map<locint, GenericDerivative<locint> > generic_derivative;
  std::cout << "In generic_mpi_tensor" << std::endl;
  generic_tape(tag, dep, indep, basepoint, ind_map, hyper_index, hyper_value);

  if (myid == DEBUG_ID) {
    std::cout << hyper_index.size() << std::endl;
    for(const locint& loc : hyper_index) {
      std::cout << loc << std::endl;
    }
  }

  struct timeval tv1, tv2;
  double time_elapsed;
  gettimeofday(&tv1, NULL);
  special_derivative_table();
  generic_mpi_reverse(tag, d, hyper_index, hyper_value, r_live_set, generic_derivative);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout << "Proc "<<myid<<" reverse time: "<<time_elapsed << std::endl;  
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
  gettimeofday(&tv1, NULL);
  generic_mpi_forward(d, r_live_set, generic_derivative);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                 (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  std::cout << "Proc "<<myid<<" forward time: "<<time_elapsed << std::endl;  

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
}
