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
#include <adolc/hypertensor/generic_reverse.h>
#include <adolc/hypertensor/opencomb.h>
#include <sys/time.h>

#define DEBUG_ID 99

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}

static const double kFactorial[11] = {1, 1, 2, 6, 24, 120, 720,
                               5040, 40320, 362880, 3628800};

int generic_reverse(short tag,
                    int order,
                    std::vector<locint>& hyper_index,
                    std::vector<double>& hyper_value,
                    GenericDerivative<locint>& global_gd) {
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
        POP_LAST_VALUE(3);
        populate_derivative_table(order, info, local_gd);
        break;
      case min_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        populate_derivative_table(order, info, local_gd);
        break;
      case mult_a_a:
//        std::cout << hyper_index.size() << "=" << hyper_value.size()<<std::endl;
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
        {
          info.opcode = plus_a_a;
          populate_derivative_table(order, info, local_gd);
          global_gd.find_and_erase(info.r, temp_gd);
          generic_d_tuples(order, info,
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
        break;
      case eq_min_prod:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        {
          info.opcode = min_a_a;
          populate_derivative_table(order, info, local_gd);
          global_gd.find_and_erase(info.r, temp_gd);
          generic_d_tuples(order, info,
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
      global_gd.find_and_erase(info.r, temp_gd);
      generic_d_tuples(order, info,
                       global_gd, local_gd, temp_gd);
    }
  }
  end_sweep();
}

double binary_comb_coefficient(int max_level, int* dx, int* dy) {
  double ret = 1;
  int n = 0;
  int m = 0;
/*
  std::cout << "D: ";
  for (int i=0; i<= max_level; i++){
    std::cout << " " << dx[i] << ","<< dy[i];
  }
  std::cout << std::endl;
*/
  for (int i=0; i<= max_level; i++) {n+=dx[i];m+=dy[i];}
  ret = kFactorial[n];
  for (int i=0; i<= max_level; i++) {
    ret = ret / kFactorial[dx[i]];
  }
  ret = ret * kFactorial[m];
  for (int i=0; i<= max_level; i++) {
    ret = ret / kFactorial[dy[i]];
  }
  int count = 1;
  int l = 2;
  while (l <= max_level) {
    if (dx[l] == dx[l-1] && dy[l] == dy[l-1]) {
      count++;
    } else {
      ret = ret / kFactorial[count];
      count = 1;
    }
    l++;
  }
  ret = ret / kFactorial[count];
//  std::cout << "binary coeff = " << ret << std::endl; 
  return ret;
}

double unary_comb_coefficient(int max_level, int* dx) {
  double ret = 1;
  int n = 0;
/*
  std::cout << "D: ";
  for (int i=0; i<= max_level; i++){
    std::cout << " " << dx[i];
  }
  std::cout << std::endl;
*/
  for (int i=0; i<= max_level; i++) {n+=dx[i];}
  ret = kFactorial[n];
  for (int i=0; i<= max_level; i++) {ret = ret / kFactorial[dx[i]];}
  int count = 1;
  int l = 2;
  while (l <= max_level) {
    if (dx[l] == dx[l-1]) {
      count++;
    } else {
      ret = ret / kFactorial[count];
      count = 1;
    }
    l++;
  }
  ret = ret / kFactorial[count];
  
  return ret;
}

void generate_binary_tuples(int curr_level,
                            int max_level,
                            int curr_order,
                            int max_order,
                            double cw,
                            typename GenericDerivative<locint>::iterator& prev,
                            int* dx,
                            int* dy,
                            double* ssw, int order, locint x, locint y) {

  if (curr_level > max_level) {
    //compute the coefficient and increase
    double coeff = binary_comb_coefficient(max_level, dx, dy);
    int cx = 0;
    int cy = 0;
    for (int i=1; i<=max_level; i++) {cx+=dx[i];cy+=dy[i];}
    ssw[cx*(order+1)+cy] += coeff* cw;
/*
    std::cout << "coeff = " << coeff
              << " cw = " << cw
              << " curr_order = " << curr_order
              << " ssw[ "<<cx<<","<<cy<<"] = "
              << ssw[cx*(order+1)+cy] << std::endl; 
*/
    return;
  }

  typename GenericDerivative<locint>::iterator iter = prev;
  bool has_next = true;
  OpenCombMultiSet<locint> dc;
  double w;
  int this_order;
  while (has_next) {
    iter.get_curr_pair(dc, w);
//    dc.debug();
//    std::cout << "w = " << w << std::endl;
    if (w != 0) {
      this_order = dc.size();
/*
      std::cout << " this_order " << this_order
                << " curr_order " << curr_order
                << " max_order" << max_order << std::endl;
      std::cout << " curr_level " << curr_level
                << " max_level " << max_level << std::endl;
*/
      if (this_order * (max_level - curr_level + 1) + curr_order <= max_order) {
        dx[curr_level] = dc.count(x);
        dy[curr_level] = dc.count(y);
        generate_binary_tuples(curr_level + 1,
                              max_level,
                              curr_order + this_order,
                              max_order,
                              cw * w,
                              iter,
                              dx,
                              dy,
                              ssw, order, x, y);
        dx[curr_level] = 0;
      }
    }
    has_next = iter.move_to_next();
  }
}

void generate_unary_tuples(int curr_level,
                           int max_level,
                           int curr_order,
                           int max_order,
                           double cw,
                           typename GenericDerivative<locint>::iterator& prev,
                           int* dx,
                           double* sw) {

  if (curr_level > max_level) {
    //compute the coefficient and increase
    double coeff = unary_comb_coefficient(max_level, dx);
    sw[curr_order] += coeff* cw;
/*
    std::cout << "coeff = " << coeff
              << " cw = " << cw
              << " sw[ " << curr_order << "] = " << sw[curr_order] << std::endl;
*/
    return;
  }
  typename GenericDerivative<locint>::iterator iter = prev;
  bool has_next = true;
  OpenCombMultiSet<locint> dc;
  double w;
  int this_order;
  while (has_next) {
    iter.get_curr_pair(dc, w);
//    dc.debug();
//    std::cout << "w = " << w << std::endl;
    if (w != 0) {
      this_order = dc.size();
/*
      std::cout << " this_order " << this_order
                << " curr_order " << curr_order
                << " max_order" << max_order << std::endl;
      std::cout << " curr_level " << curr_level
                << " max_level" << max_level << std::endl;
*/
      if (this_order * (max_level - curr_level + 1) + curr_order <= max_order) {
        dx[curr_level] = dc.size();
        generate_unary_tuples(curr_level + 1,
                              max_level,
                              curr_order + this_order,
                              max_order,
                              cw * w,
                              iter,
                              dx,
                              sw);
        dx[curr_level] = 0;
      }
    }
    has_next = iter.move_to_next();
  }
}

void generic_d_tuples(int order,
                      DerivativeInfo<locint>& info,
                      GenericDerivative<locint>& global_gd,
                      GenericDerivative<locint>& local_gd,
                      GenericDerivative<locint>& temp_gd) {
/*
//  int myid;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//if (myid == DEBUG_ID) {
  std::cout << "global: " << std::endl;
  global_gd.debug();
  std::cout << std::endl;

  std::cout << "local: " << std::endl;
  local_gd.debug();
  std::cout << std::endl;

  std::cout << "temp: " << std::endl;
  temp_gd.debug();
  std::cout << std::endl;
//}
*/
  typename GenericDerivative<locint>::iterator local_iter;
  local_iter = local_gd.get_new_iterator();
  local_iter.init_iterator();

  typename GenericDerivative<locint>::iterator temp_iter;
  temp_iter = temp_gd.get_new_iterator();
  bool temp_has_next = temp_iter.init_iterator();
  OpenCombMultiSet<locint> s_set;
  double w;
  int dx[MAX_ORDER+1];
  int dy[MAX_ORDER+1];
  double sw[MAX_ORDER+1];
  double ssw[(MAX_ORDER+1)*(MAX_ORDER+1)];
  if (info.y != NULLLOC) {
    //binary case
    while (temp_has_next) {
      temp_iter.get_curr_pair(s_set, w);
      int t_size = s_set.count(info.r);
      int s_size = s_set.size() - t_size;
//      s_set.debug();
//      std::cout << "s_size = " << s_size << ", t_size = "<< t_size<< std::endl;
      while (s_set.count(info.r) != 0) {
        s_set.remove(info.r);
      }
      dx[0] = s_set.count(info.x);
      dy[0] = s_set.count(info.y);
      for (int i=0; i<=(order+1)*(order+1); i++) {ssw[i] = 0;}
      generate_binary_tuples(1,
                             t_size,
                             0,
                             order - s_size,
                             w,
                             local_iter,
                             dx,
                             dy,
                             ssw, order, info.x, info.y);

      for(int i=0; i <= order; i++) {
        OpenCombMultiSet<locint> ss_set(s_set);
        for(int j=0; j<= order; j++) {
//          ss_set.debug();
//          std::cout << "s_weight = " << ssw[i*(order+1)+j] << std::endl;
          if (ssw[i*(order+1) + j] != 0) {
            global_gd.increase(ss_set, ssw[i*(order+1)+j]);
          }
          ss_set.put(info.y);
        }
        s_set.put(info.x);
      }
      temp_has_next = temp_iter.move_to_next();
    }
  } else if (info.x != NULLLOC) {
    //unary case
    while (temp_has_next) {
      temp_iter.get_curr_pair(s_set, w);
      int t_size = s_set.count(info.r);
      int s_size = s_set.size() - t_size;
//      s_set.debug();
//      std::cout << "s_size = " << s_size << ", t_size = "<< t_size<< std::endl;
      while (s_set.count(info.r) != 0) {
        s_set.remove(info.r);
      }
      dx[0] = s_set.count(info.x);
      for (int i=0; i<= order; i++) {sw[i] = 0;}
      generate_unary_tuples(1,
                            t_size,
                            0,
                            order - s_size,
                            w,
                            local_iter,
                            dx,
                            sw);

      for(int i=1; i <= order; i++) {
        s_set.put(info.x);
//        s_set.debug();
//        std::cout << "s_weight = " << sw[i] << std::endl;
        if (sw[i] != 0) {
          global_gd.increase(s_set, sw[i]);
        }
      }
      temp_has_next = temp_iter.move_to_next();
    }
  }
/*
//if (myid == DEBUG_ID) {
  std::cout << "new global:" << std::endl;
  global_gd.debug();
  std::cout << std::endl; 
//}
*/
}
