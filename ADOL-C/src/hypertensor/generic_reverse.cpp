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

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}

void generic_derivative_equation(int order,
                                 locint res,
                                 GenericDerivative<locint>& global_gd,
                                 GenericDerivative<locint>& local_gd,
                                 GenericDerivative<locint>& temp_gd,
                                 OpenCombMultiSet<locint>& d_set);

void generate_tuples(int level,
                     typename std::set<locint>::iterator prev_iter,
                     int order,
                     int prec_ind,
                     std::set<locint>& live_set,
                     GenericDerivative<locint>& global_gd,
                     GenericDerivative<locint>& local_gd,
                     GenericDerivative<locint>& temp_gd,
                     OpenCombMultiSet<locint>& d_set);

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      std::set<locint>& live_set,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd);

int generic_reverse(short tag,
                    int order,
                    std::vector<locint>& hyper_index,
                    std::vector<double>& hyper_value,
                    GenericDerivative<locint>& global_gd) {
//  std::cout << "In generic reverse" << std::endl;

  DerivativeInfo<locint> info;
  GenericDerivative<locint> temp_gd(order);
  GenericDerivative<locint> local_gd(order);
  OpenCombMultiSet<locint> d_set;
  std::set<locint> live_set;
  std::vector<std::set<locint>::iterator > iter_stack;
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
// TODO!
//        adjoints->increase(res, 1.0);
        {
          OpenCombMultiSet<locint> dep;
          dep.put(res);
          global_gd.increase(dep, 1.0);
        }
        GET_LAST_VALUE;
        std::cout << "Dep: " << res << std::endl;
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
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        break;
      case eq_min_prod:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
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
//      std::cout << "-------" << std::endl;
//      global_gd.debug();
//      std::cout << "-------" << std::endl;
      live_set.erase(info.r);
      global_gd.find_and_erase(info.r, temp_gd);
      iter_stack.clear();
      iter_stack.resize(order, live_set.begin());
      generic_d_tuples(order, info, live_set,
                       global_gd, local_gd, temp_gd);
    }
  }
  end_sweep();
}

void generate_tuples(int level,
                     typename std::set<locint>::iterator prev_iter,
                     int order,
                     locint prec_ind,
                     locint dep_ind,
                     std::set<locint>& live_set,
                     GenericDerivative<locint>& global_gd,
                     GenericDerivative<locint>& local_gd,
                     GenericDerivative<locint>& temp_gd,
                     OpenCombMultiSet<locint>& d_set) {
  typename std::set<locint>::iterator iter;
  if (level == order) {
    d_set.put(prec_ind);
    generic_derivative_equation(order, dep_ind,
                                global_gd, local_gd, temp_gd, d_set);
    d_set.remove(prec_ind);
  } else {
    iter = prev_iter;
    while (iter != live_set.end()) {
      d_set.put(*iter);
      generate_tuples(level + 1, iter, order, prec_ind, dep_ind,
                      live_set, global_gd, local_gd, temp_gd, d_set);
      d_set.remove(*iter);
      ++iter;
    }
  }
}

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      std::set<locint>& live_set,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd) {
//  std::cout << "global: " << std::endl;
//  global_gd.debug();
//  std::cout << std::endl;

//  std::cout << "local: " << std::endl;
//  local_gd.debug();
//  std::cout << std::endl;

//  std::cout << "temp: " << std::endl;
//  temp_gd.debug();
//  std::cout << std::endl;

  OpenCombMultiSet<locint> d_set;
  // remove x&y from live_set
  if (info.x != NULLLOC) {
    live_set.erase(info.x);
  }
  if (info.y != NULLLOC) {
    live_set.erase(info.y);
  }
  // generate tuples w.r.t x
  if (info.x != NULLLOC) {
    d_set.clear();
    live_set.insert(info.x);
    for (int i = 1; i <= order; i++) {
      generate_tuples(1, live_set.begin(), i, info.x, info.r,
                      live_set, global_gd, local_gd, temp_gd, d_set);
    }
  }
  // generate tuples w.r.t y
  if (info.y != NULLLOC) {
    d_set.clear();
    live_set.insert(info.y);
    for (int i = 1; i <= order; i++) {
      generate_tuples(1, live_set.begin(), i, info.y, info.r,
                      live_set, global_gd, local_gd, temp_gd, d_set);
    }
  }
 
}

void generic_derivative_equation(int order,
                                 locint res,
                                 GenericDerivative<locint>& global_gd,
                                 GenericDerivative<locint>& local_gd,
                                 GenericDerivative<locint>& temp_gd,
                                 OpenCombMultiSet<locint>& d_set) {
//  std::cout << "computeing: ";
//  d_set.debug();
//  std::cout << std::endl;
  OpenCombMultiSubset<locint> d_subset_generator(d_set);
  OpenCombMultiSet<locint> d_subset;
  OpenCombMultiSet<locint> c_subset;
  OpenCombPartition<locint> c_partition;
  std::vector<OpenCombMultiSet<locint> > c_vec;
  bool has_next = d_subset_generator.init_multi_subset();
  bool has_next2;
  int s;
  int t;
  double w;
  double sw = 0;
  int i;
  int l;
  while (has_next) {
    has_next = d_subset_generator.get_next_subset(d_subset, c_subset);
//    std::cout << "WTS";
//    d_subset.debug();
//    c_subset.debug();

    s = d_subset.size();
    has_next2 = c_partition.init_partition(c_subset);
//    std::cout << has_next2 << "d_size = " << s << std::endl;
    while (has_next2) {
      has_next2 = c_partition.get_next_partition(c_vec);
      //compute a value;
      t = c_vec.size();

//      std::cout << "subsets: " << std::endl;
//      for (i=0; i<t; i++) {c_vec[i].debug();}
//      std::cout << std::endl;

      for (i=0; i<t; i++) {d_subset.put(res);}
      w = temp_gd.get(d_subset);
//      d_subset.debug();
//      std::cout << " w = " << w << std::endl;
      if (w != 0.0) {
        l = 0;
        while (w!= 0.0 && l < c_vec.size()) {
          w *= local_gd.get(c_vec[l]);
          l++;
        }
      }
      for (i=0; i<t; i++) {d_subset.remove(res);}
      sw += w;
    }
  }
  if (sw != 0.0) {
    global_gd.increase(d_set, sw);
  }
//  std::cout << "sw = " << sw << " done." << std::endl;
}
