#include <vector>
#include <map>
#include <set>
#include <iostream>

#include "oplate.h"
#include "taping_p.h"
#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>
#include <adolc/hypertensor/hyper_pushing.h>
#include <adolc/hypertensor/hyper_derivative.h>

#ifdef ENABLE_GENERIC_MPI

#include "mpi.h"
#include <adolc/hypertensor/generic_mpi_trace.h>
extern std::vector<SRinfo> sr_stack;

#endif

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}

#define COMBINE_D_1 info.dx += info.dy; info.dy = 0.0
#define COMBINE_D_2 info.pxx += 2.0 * info.pxy + info.pyy;\
                    info.pxy = 0.0; info.pyy = 0.0;\
                    COMBINE_D_1;
#define COMBINE_D_3 info.pxxx += 3.0 * info.pxyy + 3.0 * info.pxxy + info.pyy;\
                    info.pxyy = 0; info.pxxy = 0.0; info.pyyy = 0.0;\
                    COMBINE_D_2

#define PSEUDO_BINARY_1 if (info.y == info.x) {info.y = NULLLOC; COMBINE_D_1;}
#define PSEUDO_BINARY_2 if (info.y == info.x) {info.y = NULLLOC; COMBINE_D_2;}
#define PSEUDO_BINARY_3 if (info.y == info.x) {info.y = NULLLOC; COMBINE_D_3;}

int hyper_mpi_reverse(short tag,
                      int order,
                      std::vector<locint>& hyper_index,
                      std::vector<double>& hyper_value,
                      std::map<locint, HyperDerivative<locint> >& global_gd) {
  ADOLC_OPENMP_THREAD_NUMBER;
  ADOLC_OPENMP_GET_THREAD_NUMBER;
  std::cout << "In hyper_third_reverse " << std::endl;
#ifdef ENABLE_GENERIC_MPI
  std::vector<SRinfo>::reverse_iterator sr_riter;
  sr_riter = sr_stack.rbegin();
#endif
  std::map<locint, std::set<locint> > reverse_live;
  unsigned char opcode;
  locint res;
  double coval = 0;
  double* d = NULL;
  int i;

  double r, x, y;

  DerivativeInfo<locint> info;
  
  init_rev_sweep(tag);
  opcode = get_op_r();
  while(opcode != start_of_tape) {
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
        info.dx = 1.0;
        POP_LAST_VALUE(2);
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
        GET_LAST_VALUE;
        global_gd[res].init(order);
        global_gd[res].adjoints->increase(res, 1.0);
        reverse_live[res].insert(res); 
//        std::cout << "Dep: " << res << std::endl;
        break;
      case eq_plus_d:
        break;
      case eq_plus_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = 1.0; info.dy = 1.0;
        POP_LAST_VALUE(3);
        PSEUDO_BINARY_1;
        break;
      case eq_min_d:
        break;
      case eq_min_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = 1.0; info.dy = -1.0;
        POP_LAST_VALUE(3);
        PSEUDO_BINARY_1;
        break;
      case eq_mult_d:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        info.dx = GET_LAST_VALUE;
        break;
      case eq_mult_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.pxy = 1.0;
        PSEUDO_BINARY_2;
        break;
      case incr_a:
      case decr_a:
        break;
      case plus_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = 1.0;
        info.dy = 1.0;
        POP_LAST_VALUE(3);
        PSEUDO_BINARY_1;
        break;
      case plus_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = 1.0;
        POP_LAST_VALUE(2);
        break;
      case min_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = -1.0;
        info.dy = 1.0;
        POP_LAST_VALUE(3);
        PSEUDO_BINARY_1;
        break;
      case min_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = -1.0;
        POP_LAST_VALUE(2);
        break;
      case mult_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        GET_LAST_VALUE;
        info.dy = GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.pxy = 1.0;
        PSEUDO_BINARY_2;
        break;
      case mult_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        info.dx = GET_LAST_VALUE;
        break;
      case div_a_a:
        info.r = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;  // r = x / y;
        y = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 1.0 / y;
        info.dy = -r / y;
        info.pyy = 2.0 * r / (y * y);
        info.pxy = -1.0 / (y * y);
        // TODO: third order
        info.pxyy = 2.0 / (y * y * y);
        info.pyyy = -6.0 * r / (y * y * y);
        PSEUDO_BINARY_3;
        break;
      case div_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        r = GET_LAST_VALUE;  // r = coval / x;
        x = GET_LAST_VALUE;
        coval = GET_LAST_VALUE;
        info.dx = -r / x;
        info.pxx = 2.0 * r / (x * x);
        // TODO: third order
        info.pxxx = -6.0 * r / (x * x);
        PSEUDO_BINARY_3;
        break;
      case eq_plus_prod:

        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        info.dx = 1.0; info.dy = 1.0;
        {
          locint dep;
          std::set<locint> dep_set = std::move(reverse_live[info.r]);
          reverse_live.erase(info.r);
          typename std::set<locint>::iterator iter;
          iter = dep_set.begin();
          while(iter != dep_set.end()) {
            dep = *iter;
            hyper_process_sac(info, order, global_gd[dep]);
            if (info.x != NULLLOC) {
              reverse_live[info.x].insert(dep);
            }
            if (info.y != NULLLOC) {
              reverse_live[info.y].insert(dep);
            }
            ++iter;
          }
        }
        
        res = info.y;
        info.clear();
        info.r = res;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dy = GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.pxy = 1.0;
        PSEUDO_BINARY_2;

        break;
      case eq_min_prod:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        POP_LAST_VALUE(3);
        info.dx = 1.0; info.dy = -1.0;
        {
          locint dep;
          std::set<locint> dep_set = std::move(reverse_live[info.r]);
          reverse_live.erase(info.r);
          typename std::set<locint>::iterator iter;
          iter = dep_set.begin();
          while(iter != dep_set.end()) {
            dep = *iter;
            hyper_process_sac(info, order, global_gd[dep]);
            if (info.x != NULLLOC) {
              reverse_live[info.x].insert(dep);
            }
            if (info.y != NULLLOC) {
              reverse_live[info.y].insert(dep);
            }
            ++iter;
          }
        }
        
        res = info.y;
        info.clear();
        info.r = res;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dy = GET_LAST_VALUE;
        info.dx = GET_LAST_VALUE;
        info.pxy = 1.0;
        PSEUDO_BINARY_2;
        break;
      case pos_sign_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = 1.0;
        POP_LAST_VALUE(2);
        break;
      case neg_sign_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = -1.0;
        POP_LAST_VALUE(2);
        break;
      case exp_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = r;
        info.pxx = r;
        info.pxxx = r;
        break;
      case sin_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = cos(x);
        info.pxx = -r;
        info.pxxx = -info.dx;
        break;
      case cos_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = -sin(x);
        info.pxx = -r;
        info.pxxx = - info.dx;
        break;
      case atan_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = 1.0 + x * x;
        info.dx = 1.0 / coval;
        info.pxx = -2.0 * x / (coval * coval);
        // TODO: third order
        info.pxxx = (6.0 * x * x - 2.0) / (coval * coval * coval);
        break;
      case asin_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = sqrt(1.0 - x * x);
        info.dx = 1.0 / coval;
        info.pxx = x / (coval * coval * coval);
        // TODO: third order
        info.pxxx = (2.0 * x * x + 1.0) / (coval * coval * coval * coval * coval);
        break;
      case acos_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = - sqrt(1.0 - x * x);
        info.dx = 1.0 / coval;
        info.pxx = x / (coval * coval * coval);
        // TODO: third order
        info.pxxx = (2.0 * x * x + 1.0) / (coval * coval * coval * coval * coval);
        break;
      case log_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 1.0 / x;
        info.pxx = -info.dx / x;
        // TODO: third order
        info.pxxx = -2.0 * info.pxx / x;
        break;
      case pow_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = GET_LAST_VALUE;
        if (x == 0.0) {
          info.dx = 0.0;
          info.pxx = 0.0;
          info.pxxx = 0.0;
        } else {
          info.dx = coval * (r / x);
          info.pxx = (coval - 1) * (info.dx / x);
          info.pxxx = (coval - 2) * (info.pxx / x);
        }
        // TODO: third order
        break;
      case sqrt_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        if (x == 0.0) {
          info.dx = 0.0;
          info.pxx = 0.0;
          info.pxxx = 0.0;
        } else {
          info.dx = 0.5 * r / x;
          info.pxx = -0.5 * info.dx / x;
          info.pxxx = -1.5 * info.pxx / x;
        }
        // TODO: third order;
        break;
      case min_op:
        info.r = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        y = GET_LAST_VALUE;
        if (x < y) {
          info.x = GET_LAST_INDEX;
          info.y = GET_LAST_INDEX;
        } else if (x > y) {
          info.y = GET_LAST_INDEX;
          info.x = GET_LAST_INDEX;
        } else {
          info.x = GET_LAST_INDEX;
          info.y = GET_LAST_INDEX;
          fprintf(DIAG_OUT, "WARNING: Tie in min(a,b)!\n");
        }
        info.dx = 1.0; info.y = NULLLOC;
        break;
      case abs_val:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        if (x != 0.0) {
          info.dx = x / abs(x);
        } else {
          fprintf(DIAG_OUT, "WARNING: abs(x), x==0!\n");
        }
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
        info.dx = 1.0;
        break; 
      case cond_assign_s:
        info.r = GET_LAST_INDEX;
        POP_LAST_VALUE(2);
        coval = GET_LAST_VALUE;
        if (coval > 0.0) {
          info.x = GET_LAST_INDEX;
          GET_LAST_INDEX;
          info.dx = 1.0;
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
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = sqrt(1.0 + x * x);
        info.dx = 1.0 / coval;
        info.pxx = -x / (coval * coval * coval);
        // TODO: third order
        break;
      case acosh_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = sqrt(x * x - 1.0);
        info.dx = 1.0 / coval;
        info.pxx = -x / (coval * coval * coval);
        // TODO: third order
        break;
      case atanh_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 1.0 / (1.0 - x * x);
        info.pxx = 2.0 * x * info.dx * info.dx;
        // TODO: third order
        break;
      case erf_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 2.0 / sqrt(acos(-1.0)) * exp(-x * x);
        info.pxx = -2.0 * x * info.dx;
        // TODO: third order
        break;
#endif // ATRIG_ERF
#ifdef ENABLE_GENERIC_MPI
      case ampi_send:
        if (sr_riter->SR_tag == RMPI_SEND_TAG) {
        } else {
          std::cout << "HYPER MPI reverse trace Send ERROR!" << std::endl;
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
              std::cout << "HYPER MPI reverse Recv opcode error!" << std::endl;
            }
          }
        } else {
          std::cout << "HYPER MPI reverse trace RecvERROR!" << std::endl;
        }
        ++sr_riter;
        break;
#endif
      default:
        fprintf(DIAG_OUT, "HYPER-TENSOR: unimplemented opcode %d\n", opcode);
    }
    // This is when we should do the work
/*
    std::cout << (int)info.opcode << " : "
              << info.r << "<---" << info.x << ", " << info.y << std::endl
              << info.dx << "," << info.dy << std::endl
              << info.pxx << "," << info.pxy << "," << info.pyy << std::endl;
*/
    if (info.r != NULLLOC) {
      locint dep;
      std::set<locint> dep_set = std::move(reverse_live[info.r]);
      reverse_live.erase(info.r);
      typename std::set<locint>::iterator iter;
      iter = dep_set.begin();
      while(iter != dep_set.end()) {
        dep = *iter;
        hyper_process_sac(info, order, global_gd[dep]);
//        global_gd[dep].debug();
        if (info.x != NULLLOC) {
          reverse_live[info.x].insert(dep);
        }
        if (info.y != NULLLOC) {
          reverse_live[info.y].insert(dep);
        }
        ++iter;
      }
    }
    opcode = get_op_r(); 
  }
  end_sweep();
}

void hyper_mpi_forward(int order,
                       std::map<locint, HyperDerivative<locint> >& global_gd) {
#ifdef ENABLE_GENERIC_MPI
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
      char* buf = (char*)malloc(sizeof(char*) * total_buf_size);
      total_buf_size = 0;
      for(int i = 0; i < sr_info.count; i++) {
        loc = sr_info.loc + i;
        global_gd[loc].write_to_byte(&buf[total_buf_size]);
        total_buf_size += global_gd[loc].byte_size();
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
        loc = sr_info.loc - i;
        HyperDerivative<locint> recv_gd(&(buf[total_buf_size]), order);
        total_buf_size += recv_gd.byte_size();
//        recv_gd.debug();      
        hyper_process_recv_gd(loc, order, recv_gd, global_gd);
      }
      free(buf);
    }
  }
#endif
}
