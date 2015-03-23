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
#include <adolc/hypertensor/generic_mpi_trace.h>
#include <adolc/hypertensor/generic_reverse.h>
#include <adolc/hypertensor/opencomb.h>
#include <sys/time.h>
#include "mpi.h"

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}

#define DEBUG_ID 0
extern std::vector<SRinfo> sr_stack;

void generic_mpi_process_sac(int order,
                             DerivativeInfo<locint>& info,
                             std::map<locint, std::set<locint> >& live_set,
                             GenericDerivative<locint>& local_gd,
                             GenericDerivative<locint>& temp_gd,
                             std::map<locint, GenericDerivative<locint> >& global_gd);

int generic_mpi_reverse(short tag,
                       int order,
                       std::vector<locint>& hyper_index,
                       std::vector<double>& hyper_value,
                       std::map<locint, std::set<locint> >& live_set,
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
  unsigned char opcode;
  locint res;
  double coval = 0;
  double* d = NULL;
  int i;
  double r, x, y;
  std::vector<SRinfo>::reverse_iterator sr_riter;
  sr_riter = sr_stack.rbegin();
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
//        std::cout << info.r << " = " << info.x << std::endl; 
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
//        std::cout << "Dummy Dep: " << res << std::endl;
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
      case ampi_send:
        if (sr_riter->SR_tag == RMPI_SEND_TAG) {
        } else {
          std::cout << "Generic MPI reverse trace Send ERROR!" << std::endl;
        }
        ++sr_riter;        
        break;
      case ampi_recv:
        if (sr_riter->SR_tag == RMPI_RECV_TAG) {
          for(int i = 0; i < sr_riter->count; i++) {
            opcode = get_op_r();
            if (opcode == assign_ind) {
              GET_LAST_INDEX;
              GET_LAST_VALUE;
            } else {
              std::cout << "Generic MPI reverse Recv opcode error!" << std::endl;
            }
          }
          for(int i = 0; i < sr_riter->count; i++) {
            opcode = get_op_r();
            if (opcode == assign_d_zero) {
              GET_LAST_INDEX;
              GET_LAST_VALUE;
            } else {
              std::cout << "Generic MPI reverse Recv opcode error!" << std::endl;
            }
          }
        } else {
          std::cout << "Generic MPI reverse trace RecvERROR!" << std::endl;
        }
        ++sr_riter;
        break;
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
  if (info.r == NULLLOC) {
    return;
  }
/*
  std::cout << " opcode = " << (int)info.opcode
            << " r = " << info.r
            << " x = " << info.x
            << " y = " << info.y << std::endl;
*/
  typename std::map<locint, std::set<locint> >::iterator dep_iter;
  locint dummy_dep;
  dep_iter = live_set.begin();
  while(dep_iter != live_set.end() ) {
    dummy_dep = dep_iter->first;
/*
    std::cout << "check dep: " << dummy_dep << "size = " << live_set[dummy_dep].size() << std::endl;
    typename std::set<locint>::iterator set_iter;
    set_iter = live_set[dummy_dep].begin();
    while(set_iter != live_set[dummy_dep].end()) {
      std::cout << *set_iter << std::endl;
      ++set_iter;
    }
*/
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

void print_live_set(std::set<locint>& live_set) {

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
if (myid == DEBUG_ID) {
  typename std::set<locint>::iterator iter;
  iter = live_set.begin();
  std::cout << "SET: <";
  while(iter != live_set.end()) {
    std::cout << " " << *iter;
    ++iter;
  }
  std::cout << ">" << std::endl;
}
}

double symmetric_coeff(
    double w,
    OpenCombMultiSet<locint>& s_set,
    std::vector<OpenCombMultiSet<locint> >& set_vec,
    std::vector<int>& set_count,
    GenericDerivative<locint>& global_gd) {
  std::map<locint, int> total_count;
  std::map<locint, int> local_count;
  double coeff = 1.0;
  if (set_vec.size() != set_count.size()) {
    std::cout << "SIZE mismatch in symmetric coeff" << std::endl;
    return 0;
  }
  OpenCombMultiSet<locint> ss_set(s_set);
  locint loc;
  typename OpenCombMultiSet<locint>::iterator iter;
  iter = s_set.begin();
  double tg = 1.0;
  double tl = 1.0;
  while(iter != s_set.end()) {
    loc = *iter;
    total_count[loc]++;
    local_count[loc]++;
    tg *= total_count[loc];
    tl *= local_count[loc];
    ++iter;
  }
  coeff = coeff * tg / tl;

  for(int i = 0; i < set_vec.size(); i++) {
    for(int j = 0; j < set_count[i]; j++) {
      local_count.clear();
      iter = set_vec[i].begin();
      tg = 1.0;
      tl = 1.0;
      while(iter != set_vec[i].end()) {
        loc = *iter;
        ss_set.put(loc);
        total_count[loc]++;
        local_count[loc]++;
        tg *= total_count[loc];
        tl *= local_count[loc];
        ++iter;
      }
      coeff = coeff * tg / (j+1) / tl;
    }
  }
// Debug
/*
  std::cout << " coeff = " << coeff
            << " w = " << w << std::endl;
  for (int i=0; i< set_vec.size(); i++) {
    set_vec[i].debug();
    std::cout << std::endl;
  }
*/
  global_gd.increase(ss_set, w * coeff);
}

void compute_multi_set(
    int curr_level,
    int max_level,
    int curr_order,
    int max_order,
    double w,
    typename GenericDerivative<locint>::iterator& prev,
    std::vector<OpenCombMultiSet<locint> >& set_vec,
    std::vector<int>& set_count,
    OpenCombMultiSet<locint>& s_set,
    GenericDerivative<locint>& global_gd) {

  if (curr_level > max_level) {
// compute the symmetric coefficient
    symmetric_coeff(w, s_set, set_vec, set_count, global_gd);    
    return;
  }

  int this_order;
  typename GenericDerivative<locint>::iterator iter = prev;
  OpenCombMultiSet<locint> dc;
  double cw;
  iter.get_curr_pair(dc, cw);
//  dc.debug();
//  std::cout << " cw = " << cw << std::endl;
//  this_order = dc.size();
  if (curr_level == 1) {
    set_vec.push_back(dc);
    set_count.push_back(1);
  } else  {
    set_count[set_count.size() - 1]++;
  }
  compute_multi_set(curr_level + 1, max_level,
                    curr_order + this_order, max_order,
                    cw * w, iter, set_vec, set_count, s_set, global_gd);
  if (curr_level == 1) {
    set_vec.pop_back();
    set_count.pop_back();
  } else {
    set_count[set_count.size() - 1]--;
  }
  bool has_next = iter.move_to_next();
  while(has_next) {
    iter.get_curr_pair(dc, cw);
//    dc.debug();
//    std::cout << " cw = " << cw << std::endl;
    if (w != 0) {
      this_order = dc.size();
      set_vec.push_back(std::move(dc));
      set_count.push_back(1);
      if (this_order * (max_level - curr_level + 1) + curr_order <= max_order) {
        compute_multi_set(curr_level + 1, max_level,
                          curr_order + this_order, max_order,
                          cw * w, iter,                          
                          set_vec, set_count, s_set, global_gd);
      }
      set_vec.pop_back();
      set_count.pop_back();
    }
    has_next = iter.move_to_next();
  }
}

void generic_tuples(int order,
                    locint res,
                    std::set<locint>& live_set,
                    GenericDerivative<locint>& global_gd,
                    GenericDerivative<locint>& local_gd,
                    GenericDerivative<locint>& temp_gd) {
//populare live set
//  print_live_set(live_set);
  typename GenericDerivative<locint>::iterator local_iter;
  local_iter = local_gd.get_new_iterator();
  bool local_has_next = local_iter.init_iterator();
  OpenCombMultiSet<locint> s_set;
  typename OpenCombMultiSet<locint>::iterator s_iter;
  double w;
  while(local_has_next) {
    local_iter.get_curr_pair(s_set, w);
    s_iter = s_set.begin();
    while(s_iter != s_set.end()) {
      live_set.insert(*s_iter);
      ++s_iter;
    }
    local_has_next = local_iter.move_to_next();
  }
//  print_live_set(live_set);

// compute derivatives
  typename GenericDerivative<locint>::iterator temp_iter;
  temp_iter = temp_gd.get_new_iterator();
  bool temp_has_next = temp_iter.init_iterator();
  while(temp_has_next) {
    temp_iter.get_curr_pair(s_set, w);
    int t_size = s_set.count(res);
    int s_size = s_set.size() - t_size;
//    s_set.debug();
//    std::cout << " w = " << w << std::endl;

    while (s_set.count(res) != 0) {
      s_set.remove(res);
    }

    local_iter = local_gd.get_new_iterator();
    std::vector<OpenCombMultiSet<locint> > set_vec;
    std::vector<int> set_count;
    local_has_next = local_iter.init_iterator();
    if (local_has_next) {
      compute_multi_set(1, t_size, 0, order - s_size, w,
                        local_iter, set_vec, set_count, s_set, global_gd); 
    }

    temp_has_next = temp_iter.move_to_next();
  }

}

bool generic_is_linear(GenericDerivative<locint>& local_gd) {
//  std::cout << "In checkint : ";
//  local_gd.debug();
//  std::cout << std::endl;

  int total_num_elements = 0;
  OpenCombMultiSet<locint> s_set;
  double w;

  typename GenericDerivative<locint>::iterator local_iter;
  local_iter = local_gd.get_new_iterator();
  bool local_has_next = local_iter.init_iterator();
  if (local_has_next) {
    local_iter.get_curr_pair(s_set, w);
    local_has_next = local_iter.move_to_next();
    if ((!local_has_next) && (s_set.size() == 1)) {
      return true;
    }
  }
  return false;
}

void generic_linear_tuples(int order,
                           locint res,
                           std::set<locint>& live_set,
                           GenericDerivative<locint>& global_gd,
                           GenericDerivative<locint>& local_gd,
                           GenericDerivative<locint>& temp_gd) {
//populare live set
  typename GenericDerivative<locint>::iterator local_iter;
  local_iter = local_gd.get_new_iterator();
  bool local_has_next = local_iter.init_iterator();
  OpenCombMultiSet<locint> s_set;
  locint x;
  typename OpenCombMultiSet<locint>::iterator s_iter;
  double w;
  local_iter.get_curr_pair(s_set, w);
  s_iter = s_set.begin();
  x = *s_iter;
  live_set.insert(x);

// compute derivatives
  typename GenericDerivative<locint>::iterator temp_iter;
  temp_iter = temp_gd.get_new_iterator();
  bool temp_has_next = temp_iter.init_iterator();
  double lw = 0;
  double cw = 0;
  while(temp_has_next) {
    temp_iter.get_curr_pair(s_set, lw);
    int t_size = s_set.count(res);
    cw = 1;
    for(int i=0; i < t_size; i++) {
      s_set.remove(res);
      s_set.put(x);
      cw = cw * w;
    }
    global_gd.increase(s_set, lw * cw);
    temp_has_next = temp_iter.move_to_next();
  }
}

void generic_mpi_process_recv_gd(
    int order,
    locint res,
    std::map<locint, std::set<locint> >& live_set,
    GenericDerivative<locint>& recv_gd,
    std::map<locint, GenericDerivative<locint> >& global_gd) {
  typename std::map<locint, std::set<locint> >::iterator dep_iter;
  locint dummy_dep;
  GenericDerivative<locint> temp_gd(order);

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (generic_is_linear(recv_gd)) {
// For linear derivatives, we have a more efficient method
//    std::cout << " res = " << res << " is linear : ";
//    recv_gd.debug();
    std::cout << std::endl;
    dep_iter = live_set.begin();
    while(dep_iter != live_set.end() ) {
      dummy_dep = dep_iter->first;
      if (dep_iter->second.find(res) != dep_iter->second.end()) {
        // do the work
        live_set[dummy_dep].erase(res);
        global_gd[dummy_dep].find_and_erase(res, temp_gd);
        generic_linear_tuples(order, res, live_set[dummy_dep],
                              global_gd[dummy_dep], recv_gd, temp_gd);
      }
      dep_iter++;
    }
    return;
  }




  dep_iter = live_set.begin();
  while(dep_iter != live_set.end() ) {
    dummy_dep = dep_iter->first;
    if (dep_iter->second.find(res) != dep_iter->second.end()) {
      // do the work
      live_set[dummy_dep].erase(res);
      global_gd[dummy_dep].find_and_erase(res, temp_gd);
      generic_tuples(order, res, live_set[dummy_dep],
                       global_gd[dummy_dep], recv_gd, temp_gd);
    }
    dep_iter++;
  }
}


void generic_mpi_forward(int order,
                         std::map<locint, std::set<locint> >& live_set,
                         std::map<locint, GenericDerivative<locint> >& global_gd) {
  int myid;
  locint loc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  for(const SRinfo& sr_info: sr_stack) {
    if (sr_info.SR_tag == RMPI_SEND_TAG) {
      int total_buf_size = 0;
      for(int i = 0; i < sr_info.count; i++) {
        loc = sr_info.loc + i;
        total_buf_size += global_gd[loc].byte_size();
      }
      char* buf = (char*)malloc(sizeof(char) * total_buf_size);
      total_buf_size = 0;
      for(int i = 0; i < sr_info.count; i++) {
        loc = sr_info.loc + i;
        total_buf_size += global_gd[loc].to_byte(&(buf[total_buf_size]));
        global_gd.erase(loc);
      }
      MPI_Send(&total_buf_size, 1, MPI_INT, sr_info.peer,
               sr_info.tag, sr_info.comm);
      MPI_Send((void*)buf, total_buf_size, MPI_CHAR, sr_info.peer,
               sr_info.tag, sr_info.comm);
//      std::cout << myid << " send size " << total_buf_size << std::endl;
      free(buf);
    } else {
      int total_buf_size = 0;
      int buf_size;
      MPI_Recv(&total_buf_size, 1, MPI_INT, sr_info.peer, sr_info.tag,
               sr_info.comm, MPI_STATUS_IGNORE);
//      std::cout << myid << " recv size " << total_buf_size << std::endl;
      char* buf = (char*)malloc(sizeof(char) * total_buf_size);
      MPI_Recv((void*)buf, total_buf_size, MPI_CHAR, sr_info.peer,
               sr_info.tag, sr_info.comm, MPI_STATUS_IGNORE);
      total_buf_size = 0;
      for(int i = 0; i < sr_info.count; i++) {
        loc = sr_info.loc + i;
        GenericDerivative<locint> recv_gd(&(buf[total_buf_size]), buf_size);
        total_buf_size += buf_size;
//        recv_gd.debug();      
        generic_mpi_process_recv_gd(order, loc,
                                    live_set, recv_gd, global_gd);
      }
      free(buf);
    }
  }
}
