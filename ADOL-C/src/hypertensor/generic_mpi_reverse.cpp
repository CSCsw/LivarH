#include <vector>
#include <map>
#include <set>
#include <iostream>

#include "oplate.h"
#include "taping_p.h"
#include <adolc/adolc.h>
#include <adolc/hypertensor/hyper_common.h>
#include <adolc/hypertensor/generic_derivative.h>
#include <adolc/hypertensor/generic_derivative_table.h>
#include <adolc/hypertensor/opencomb.h>
#include <sys/time.h>

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}

void generic_mpi_process_sac(int order,
                             DerivativeInfo<locint>& info,
                             std::map<locint, std::set<locint> >& live_set,
                             GenericDerivative<locint>& local_gd,
                             GenericDerivative<locint>& temp_gd,
                             std::map<locint, GenericDerivative<locint> >& global_gd);

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      std::set<locint>& live_set,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd);

int generic_mpi_reverse(short tag,
                       int order,
                       std::vector<locint>& hyper_index,
                       std::vector<double>& hyper_value,
                       std::map<locint, GenericDerivative<locint> >& global_gd) {
//  std::cout << "In generic reverse" << std::endl;
// static set partition generator
/*
  std::vector<
      std::vector<
          std::vector<OpenCombMultiSet<locint> > > > static_set;
  generate_static_set(order, static_set);
*/
  ADOLC_OPENMP_THREAD_NUMBER;
  ADOLC_OPENMP_GET_THREAD_NUMBER;
  DerivativeInfo<locint> info;
  GenericDerivative<locint> temp_gd(order);
  GenericDerivative<locint> local_gd(order);
  OpenCombMultiSet<locint> d_set;
  std::map<locint, std::set<locint>> live_set;
  unsigned char opcode;
  locint res;
  double coval = 0;
  double* d = NULL;
  int i;
  double r, x, y;
  init_rev_sweep(tag);
  opcode = get_op_r();
  while(opcode != start_of_tape) {
    local_gd.clear();
    info.clear();
    info.opcode = opcode;
    switch (opcode) {
      case end_of_op:
        get_op_block_r();
        opcode = get_op_r();
        break;
      case start_of_tape:
      case end_of_tape:
        break;
      case eq_zero:
      case neq_zero:
      case gt_zero:
      case lt_zero:
      case ge_zero:
      case le_zero:
        break;
      case assign_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case neq_a_a:
      case eq_a_a:
      case le_a_a:
      case ge_a_a:
      case lt_a_a:
      case gt_a_a:
      case assign_d:
      case assign_d_zero:
      case assign_d_one:
        info.r = GET_LAST_INDEX;
        coval = GET_LAST_VALUE;
        break;
      case assign_ind:
        GET_LAST_INDEX;
        GET_LAST_VALUE;
        break;
      case assign_dep:
        res = GET_LAST_INDEX;
        {
          OpenCombMultiSet<locint> dep;
          dep.put(res);
          GenericDerivative<locint> dummy_gd(order);
          dummy_gd.increase(dep, 1.0);
          global_gd.insert(std::make_pair(res, dummy_gd));
          live_set[res].insert(res);
        }
        GET_LAST_VALUE;
        std::cout << "Dummy Dep: " << res << std::endl;
        break;
      case eq_plus_d:
        break;
      case eq_plus_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        populate_derivative_table(order, info, local_gd);
        break;
      case eq_min_d:
        break;
      case eq_min_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        populate_derivative_table(order, info, local_gd);
        break;
      case eq_mult_d:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case eq_mult_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case incr_a:
      case decr_a:
        break;
      case plus_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        populate_derivative_table(order, info, local_gd);
        break;
      case plus_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case min_a_a:
        info.r = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        populate_derivative_table(order, info, local_gd);
        break;
      case min_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case mult_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case mult_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case div_a_a:
        info.r = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;  // r = x / y;
        info.dy = GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case div_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        GET_LAST_VALUE;  // r = coval / x;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case eq_plus_prod:
/*
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        {
          info.opcode = plus_a_a;
          populate_derivative_table(order, info, local_gd);
        }
        local_gd.clear();
        info.opcode = mult_a_a;
        info.r = info.y;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
*/
        break;
      case eq_min_prod:
/*
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        {
          info.opcode = min_a_a;
          populate_derivative_table(order, info, local_gd);
          live_set.erase(info.r);
          global_gd.find_and_erase(info.r, temp_gd);
          generic_d_tuples(order, info, live_set,
                           global_gd, local_gd, temp_gd);
        }
        local_gd.clear();
        info.opcode = mult_a_a;
        info.r = info.y;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
*/
        break;
      case pos_sign_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case neg_sign_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case exp_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case sin_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case cos_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case atan_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case asin_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case acos_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case log_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case pow_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case sqrt_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case min_op:
        info.r = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        if (info.dx <= info.dy) {
          info.x = GET_LAST_INDEX;
          GET_LAST_INDEX;
        } else {
          GET_LAST_INDEX;
          info.x = GET_LAST_INDEX;
        }
        populate_derivative_table(order, info, local_gd);
        break;
      case abs_val:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case cond_assign:
        info.r = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        coval = GET_LAST_VALUE;
        if (coval > 0.0) {
          GET_LAST_INDEX;
          info.x = GET_LAST_INDEX;
          GET_LAST_INDEX;
        } else {
          info.x = GET_LAST_INDEX;
          GET_LAST_INDEX;
          GET_LAST_INDEX;
        }
        populate_derivative_table(order, info, local_gd);
        break; 
      case cond_assign_s:
        info.r = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        coval = GET_LAST_VALUE;
        if (coval > 0.0) {
          info.x = GET_LAST_INDEX;
          GET_LAST_INDEX;
          info.dx = 1.0;
          populate_derivative_table(order, info, local_gd);
        }
        break;
      // Dummy opcodes in reverse
      case death_not:
      case gen_quad:
      case end_of_int:
      case end_of_val:
      case take_stock_op:
      case ceil_op:
      case floor_op:
      case ext_diff:
      case ext_diff_iArr:
        break;
                
#ifdef ATRIG_ERF
      case asinh_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case acosh_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case atanh_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
      case erf_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        populate_derivative_table(order, info, local_gd);
        break;
#endif // ATRIG_ERF
      default:
        fprintf(DIAG_OUT, "HYPER-TENSOR: unimplemented opcode %d\n", opcode);
    }
    opcode = get_op_r();
    // This is where we should do the work
    if (info.r != NULLLOC) {
      generic_mpi_process_sac(order, info, live_set,
                              local_gd, temp_gd, global_gd);
    }
  }
  end_sweep();
}

void generic_mpi_process_sac(int order,
                             DerivativeInfo<locint>& info,
                             std::map<locint, std::set<locint> >& live_set,
                             GenericDerivative<locint>& local_gd,
                             GenericDerivative<locint>& temp_gd,
                             std::map<locint, GenericDerivative<locint> >& global_gd) {
  typename std::map<locint, std::set<locint> >::iterator dep_iter;
  locint dummy_dep;
  dep_iter = live_set.begin();
  while(dep_iter != live_set.end() ) {
    dummy_dep = dep_iter->first;
    if (dep_iter->second.find(info.r) != dep_iter->second.end()) {
      // do the work
      live_set[dummy_dep].erase(info.r);
      global_gd[dummy_dep].find_and_erase(info.r, temp_gd);
      generic_d_tuples(order, info, live_set[dummy_dep],
                       global_gd[dummy_dep], local_gd, temp_gd);
    }
    dep_iter++;
  }
}
