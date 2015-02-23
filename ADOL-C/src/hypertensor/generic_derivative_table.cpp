#include <iostream>
#include <cmath>

#include "adolc/adolc.h"
#include "adolc/hypertensor/hyper_common.h"
#include "adolc/hypertensor/opencomb.h"
#include "adolc/hypertensor/generic_derivative.h"
#include "oplate.h"

double c_atan[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];
double c_asin[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];
double c_acos[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];
double c_atanh[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];
double c_asinh[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];
double c_acosh[MAX_ORDER+1][MAX_ORDER+1][MAX_ORDER+1];



void special_derivative_table() {
  for (int i=0; i<=MAX_ORDER; i++) {
    for (int j=0; j<=MAX_ORDER; j++) {
      for (int k=0; k<=MAX_ORDER; k++) {
        c_atan[i][j][k] = 0.0;
        c_asin[i][j][k] = 0.0;
      }
    }  
  }
  c_atan[1][1][0] = 1.0;
  c_asin[1][1][0] = 1.0;
  for (int i=2; i<=MAX_ORDER; i++) {
    for (int j = MAX_ORDER; j>=2; j--) {
      for (int k=0; k < MAX_ORDER; k++) {
        c_atan[i][j][k] = c_atan[i-1][j][k+1] * (k+1);
        c_asin[i][j][k] = c_asin[i-1][j][k+1] * (k+1);
      }
      for (int k=1; k <= MAX_ORDER; k++) {
        c_atan[i][j][k] += c_atan[i-1][j-1][k-1] * 2 * (1-j);
        c_asin[i][j][k] += c_asin[i-1][j-1][k-1] * 2 * (j - 1.5);
      }
    }
  }
  for (int i=0; i<=MAX_ORDER; i++) {
    for (int j=0; j<=MAX_ORDER; j++) {
      for (int k=0; k<=MAX_ORDER; k++) {
        std::cout << "c_asin[" << i
                  << ", " << j
                  << ", " << k
                  << "] = " << c_asin[i][j][k] << std::endl;
      }
    }  
  }
}

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
      // 2:
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
      // n:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0/info.dy);
      term.clear();
      term.put(info.y);
      local_gd.increase(term, -info.dx/(info.dy * info.dy));
      if (order > 1) {
        if (info.x != info.y) {
          coval = 1.0 / info.dy;
          OpenCombMultiSet<locint> t_term;
          t_term.put(info.x);
          for (int i=2; i<=order; i++) {
            coval = coval * (1-i) / info.dy;
            t_term.put(info.y);
            local_gd.increase(t_term, coval);
          }
          coval = -info.dx / (info.dy * info.dy);
          for (int i=2; i<=order; i++) {
            coval = coval * (-i) / info.dy;
            term.put(info.y);
            local_gd.increase(term, coval);
          }
        }
      }
      break;
    case div_d_a:
      // n:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -info.dy/(info.dx * info.dx));
      if (order > 1){
        coval = -info.dy / (info.dx * info.dx);
        for (int i = 2; i <= order; i++){
          coval = coval * (-i) / info.dx;
          term.put(info.x);
          local_gd.increase(term, coval);
        }
      }
      break;
    case eq_plus_prod:
    case eq_min_prod:
      std::cout << "should never see this." << std::endl;
      break;
    case exp_op:
      // n:
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
      // n:
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
      // n:
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
      // n:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / (1.0 + info.dx * info.dx));
      if (order > 1) {
        coval = 1.0 / (1.0 + info.dx * info.dx);
        double sw = 0;
        double w = 0;
        double s = 0;
        for (int i=2; i <= order; i++) {
          term.put(info.x);
          sw = 0;
          s = 1;
          for (int j=0; j <= i; j++) {
            w = 1;
            for (int k = 0; k <= j; k++) {
              sw += c_atan[i][j][k] * w * s;
              w = w * info.dx;
            }
            s = s * coval;
          }  
          local_gd.increase(term, sw);
        }
      }
      break;
    case asin_op:
      // n:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, 1.0 / sqrt(1.0 - info.dx * info.dx));
      if (order > 1) {
        coval = 1.0 / (1.0 - info.dx * info.dx);
        double sw = 0;
        double w = 0;
        double s = 0;
        for (int i=2; i <= order; i++) {
          term.put(info.x);
          sw = 0;
          s = sqrt(1.0 - info.dx * info.dx);
          for (int j=0; j <= i; j++) {
            w = 1;
            for (int k = 0; k <= j; k++) {
              sw += c_asin[i][j][k] * w * s;
              w = w * info.dx;
            }
            s = s * coval;
          }
          local_gd.increase(term, sw);
        }
      }
      break;
    case acos_op:
      // n:
      term.clear();
      term.put(info.x);
      local_gd.increase(term, -1.0 / sqrt(1.0 - info.dx * info.dx));
      if (order > 1) {
        coval = 1.0 / (1.0 - info.dx * info.dx);
        double sw = 0;
        double w = 0;
        double s = 0;
        for (int i=2; i <= order; i++) {
          term.put(info.x);
          sw = 0;
          s = sqrt(1.0 - info.dx * info.dx);
          for (int j=0; j <= i; j++) {
            w = 1;
            for (int k = 0; k <= j; k++) {
              sw += c_asin[i][j][k] * w * s;
              w = w * info.dx;
            }
            s = s * coval;
          }
          local_gd.increase(term, -sw);
        }
      }
      break;
    case log_op:
      // n:
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
      // n:
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
      // n:
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
