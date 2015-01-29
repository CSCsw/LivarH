#include <vector>
#include <map>
#include <iostream>

#include "oplate.h"
#include "taping_p.h"
#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>
#include <adolc/hypertensor/hyper_third_reverse.h>
#include <adolc/hypertensor/hyper_reverse.h>

#define GET_LAST_INDEX hyper_index.back(); hyper_index.pop_back();
#define GET_LAST_VALUE hyper_value.back(); hyper_value.pop_back();
#define POP_LAST_VALUE(n) for(int i = 0; i < n; ++i) {hyper_value.pop_back();}
int hyper_third_reverse(short tag,
                        std::vector<locint>& hyper_index,
                        std::vector<double>& hyper_value,
                        VectorGraph<locint>* adjoints,
                        MatrixGraph<locint>* hessian,
                        HyperGraph<locint>* tensor) {
  std::cout << "In hyper_third_reverse " << std::endl;
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
        adjoints->increase(res, 1.0);
        GET_LAST_VALUE;
        std::cout << "Dep: " << res << std::endl;
        break;
      case eq_plus_d:
        break;
      case eq_plus_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = 1.0; info.dy = 1.0;
        POP_LAST_VALUE(3);
        break;
      case eq_min_d:
        break;
      case eq_min_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = 1.0; info.dy = -1.0;
        POP_LAST_VALUE(3);
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
        break;
      case plus_d_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.dx = 1.0;
        POP_LAST_VALUE(3);
        break;
      case min_a_a:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        info.y = GET_LAST_INDEX;
        info.dx = -1.0;
        info.dy = 1.0;
        POP_LAST_VALUE(3);
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
        break;
      case acos_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        coval = - sqrt(1.0 - x * x);
        info.dx = 1.0 / coval;
        info.pxx = 1.0 / (coval * coval * coval);
        // TODO: third order
        break;
      case log_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 1.0 / x;
        info.pxx = -info.dx / x;
        // TODO: third order
        break;
      case sqrt_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        if (x = 0.0) {
          info.dx = 0.0;
          info.pxx = 0.0;
        } else {
          info.dx = 0.5 * r / x;
          info.pxx = -0.5 * info.dx / x;
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
        info.dx = 0;
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
      case cond_assign_s:
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

#ifdef ARTIG_ERF
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
        info.dx = 2.0 * x * info.dx * info.dx;
        // TODO: third order
        break;
      case erf_op:
        info.r = GET_LAST_INDEX;
        info.x = GET_LAST_INDEX;
        r = GET_LAST_VALUE;
        x = GET_LAST_VALUE;
        info.dx = 2.0 / sqrt(acos(-1.0)) * exp(-x * x);
        info.pxx = -2.0 * x * info->dx;
        // TODO: third order
        break;
#endif // ATRIG_ERF

      default:
        fprintf(DIAG_OUT, "HYPER-TENSOR: unimplemented opcode %d\n", opcode);
    }
    // This is when we should do the work
    if (info.r != NULLLOC) {
      // TODO: third order
      // Hessian
      hyper_hessian(info, adjoints, hessian);
      // Adjoints
      hyper_adjoints(info, adjoints);
    }
    opcode = get_op_r(); 
  }
  end_sweep();
}
