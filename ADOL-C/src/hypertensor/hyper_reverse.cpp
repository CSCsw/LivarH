#include <iostream>

#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>

int hyper_hessian(DerivativeInfo<locint>& info,
                  VectorGraph<locint>* adjoints,
                  MatrixGraph<locint>* hessian) {
  std::cout << "In hessian" << std::endl;
  double w = adjoints->get(info.r);
  VectorGraph<locint>* r = hessian->get_and_erase(info.r);
  bool has_next = r->reset();
  locint p;
  double pw;
  while (has_next) {
    has_next = r->get_next(p, pw);
    if (w != 0.0) {
      if (info.y != NULLLOC){
        if (p != info.r){
          if (p == info.x){
            hessian->increase(p, p, 2 * info.dx * pw);
          } else {
            hessian->increase(info.x, p, info.dx * pw);
          }
          if(p == info.y){
            hessian->increase(p, p, 2 * info.dy * pw);
          } else {
            hessian->increase(info.y, p, info.dy * pw);
          }
        } else {
          if (info.x != info.y){
            hessian->increase(info.x, info.x, info.dx * info.dx * pw);
            hessian->increase(info.x, info.y, info.dx * info.dy * pw);
            hessian->increase(info.y, info.y, info.dy * info.dy * pw);
          } else {
            hessian->increase(info.x, info.x, (info.dx * info.dx + 
                                               2 * info.dx * info.dy +
                                               info.dy * info.dy) * pw);
          }
        }
      } else {
        if (p != info.r){
          if (p == info.x){
            hessian->increase(p, p, 2 * info.dx * pw);
          } else {
            hessian->increase(info.x, p, info.dx * pw);
          }
        } else {
          hessian->increase(info.x,info.x,info.dx * info.dx * pw);
        }
      }
    } // if (pw != 0)
  } // while (has_next)

  if (w != 0.0) {
    if (info.pxx != 0.0) {
      hessian->increase(info.x, info.x, w * info.pxx);
    }
    if (info.pyy != 0.0) {
      hessian->increase(info.y, info.y, w * info.pyy);
    }
    if (info.pxy != 0.0) {
      if (info.x != info.y) {
        hessian->increase(info.x, info.y, w * info.pxy);
      } else {
        hessian->increase(info.x, info.y, 2.0 * w * info.pxy);
      }
    }
  } // if (w != 0)
}

int hyper_adjoints(DerivativeInfo<locint>& info,
                   VectorGraph<locint>* adjoints) {
  std::cout << "In adjoints" << std::endl;
  double w = adjoints->get_and_erase(info.r);
  if (w == 0.0) {
    return 0;
  }
  if (info.dx != 0.0) {
    adjoints->increase(info.x, w * info.dx);
  }
  if (info.dy != 0.0) {
    adjoints->increase(info.y, w * info.dy);
  }
}

