#include <iostream>
#include <cmath>

#include "adolc/adolc.h"
#include "adolc/hypertensor/hyper_common.h"
#include "adolc/hypertensor/opencomb.h"
#include "adolc/hypertensor/generic_derivative.h"
#include "oplate.h"


void populate_derivative_table(int order,
                               DerivativeInfo<locint>& info,
                               GenericDerivative<locint>& local_gd) {
//  std::cout << "in populate_derivative_table: " << std::endl;
  double coval;
  OpenCombMultiSet<locint> term;
  unsigned char opcode = info.opcode;
  switch(opcode) {
    case cond_assign:
    case cond_assign_s:
    case min_op:
    case pos_sign_a:
    case plus_d_a:
    case assign_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0);
      break;
    case plus_a_a:
    case eq_plus_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0);
      term.clear();
      term.put(info.y);
      local_gd.increase(term, 1.0);
      break;
    case min_a_a:
    case eq_min_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0);
      term.clear();
      term.put(info.y);
      local_gd.increase(term, -1.0);
      break;
    case mult_d_a:
    case eq_mult_d:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, info.dy);
      break;
    case mult_a_a:
    case eq_mult_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, info.dy);
      term.clear();
      term.put(info.y);
      local_gd.increase(term, info.dx);
      if (order > 1) {
        term.put(info.x);
        if (info.x == info.y) { // r = x * x
          info.y = NULLLOC;
          local_gd.increase(term, 2.0);
        } else {
          local_gd.increase(term, 1.0);
        }
      }
      break;
    case neg_sign_a:
    case min_d_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -1.0);
      break;
    case div_a_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0/info.dy);
      term.clear();
      term.put(info.y);
      local_gd.increase(term, -info.dx/(info.dy * info.dy));
      if (order > 1) {
        if (info.x != info.y) {
          term.put(info.x);
          local_gd.increase(term, -1.0 / (info.dy * info.dy));
          term.remove(info.x);
          term.put(info.y);
          local_gd.increase(term, -2.0*info.dx/(info.dy*info.dy*info.dy));
        }
      }
      break;
    case div_d_a:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -info.dy/(info.dx * info.dx));
      if (order > 1){
        term.put(info.x);
        local_gd.increase(term, 2.0*info.dy/(info.dx*info.dx*info.dx));
      }
      break;
    case eq_plus_prod:
    case eq_min_prod:
      std::cout << "unimplementd part" << std::endl;
      break;
    case exp_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, exp(info.dx));
      if (order > 1) {
        for(int i = 2; i < order; i++){
          term.put(info.x);
          local_gd.increase(term, exp(info.dx));
        }
      }
      break;
    case sin_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, cos(info.dx));
      if (order > 1) {
        double dsin[4] = {sin(info.dx), cos(info.dx),
                          -sin(info.dx), -cos(info.dx)};
        for (int i = 2; i <= order; i++) {
          term.put(info.x);
          local_gd.increase(term, dsin[i % 4]);
        }
      }
      break;
    case cos_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -sin(info.dx));
      if (order > 1) {
        double dcos[4] = {cos(info.dx), -sin(info.dx),
                          -cos(info.dx), sin(info.dx)};
        for (int i = 2; i <= order; i++) {
          term.put(info.x);
          local_gd.increase(term, dcos[i % 4]);
        }
      }
      break;
    case atan_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / (1.0 + info.dx * info.dx));
      break;
    case asin_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / sqrt(1.0 - info.dx * info.dx));
      break;
    case acos_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -1.0 / sqrt(1.0 - info.dx * info.dx));
      break;
    case log_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / info.dx);
      if (order > 1) {
        double dlog = 1.0 / info.dx;
        for(int i = 2; i <= order; i++) {
          term.put(info.dx);
          dlog = dlog * ((1-i) / info.dx);
          local_gd.increase(term, dlog);
        }
      }
      break;
    case pow_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, info.dy * pow(info.dx, info.dy - 1));
      if (order > 1) {
        double dpow = info.dy * pow(info.dx, info.dy - 1);
        for (int i = 2; i <= order; i++) {
          term.put(info.x);
          dpow = (dpow * (info.dy + 1 - i)) / info.dy;
          local_gd.increase(term, dpow);
        }
      }
      break;
    case sqrt_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 0.5 / sqrt(info.dx));
      if (order > 1) {
        double dsqrt = 0.5 /sqrt(info.dx);
        for(int i = 2; i <= order; i++) {
          term.put(info.x);
          dsqrt = dsqrt * ((1.5-i) / info.dx);
          local_gd.increase(term, dsqrt);
        }
      }
      break;
    case abs_val:
      // 1:
      term.clear();
      term.put(info.x);
      coval = 1.0;
      if (info.dx < 0) { coval = -1.0;}
      local_gd.increase(term, coval); 
      break;
    case asinh_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / sqrt(1.0 + info.dx * info.dx));
      break;
    case acosh_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / sqrt(info.dx * info.dx - 1.0));
      break;
    case atanh_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / (1.0 + info.dx * info.dx));
      break;
    case erf_op:
      // 1:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 2.0 / sqrt(acos(-1.0)) * exp (-info.x * info.x));
      break;
    default:
      fprintf(DIAG_OUT, "Generic Derivative Error: Unknown opcode = %d\n", opcode);
  }
//  local_gd.debug();
}
