#include <iostream>

#include <adolc/adolc.h>
#include <adolc/hypertensor/VectorGraph.h>
#include <adolc/hypertensor/MatrixGraph.h>
#include <adolc/hypertensor/HyperGraph.h>
#include <adolc/hypertensor/hyper_common.h>
#include <adolc/hypertensor/hyper_derivative.h>


#ifndef ENABLE_GENERIC_MPI

#define IF_MY_DEBUG if(1==0){
#define END_MY_DEBUG }

#else // ENABLE_GENERIC_MPI
#include "mpi.h"
#define DEBUG_ID 99

#define IF_MY_DEBUG if(myid==DEBUG_ID){

#define END_MY_DEBUG }

#endif // ENABLE_GENERIC_MPI

void hyper_third(DerivativeInfo<locint>& info,
                 VectorGraph<locint>* adjoints,
                 MatrixGraph<locint>* hessian,
                 HyperGraph<locint>* tensor,
                 double w,
                 VectorGraph<locint>* r,
                 MatrixGraph<locint>* e) {
//  std::cout << "In tensor " << std::endl;

  locint p;
  locint q;
  double pw;
  double coeff = 0.0;
  double xcoeff = 0.0;
  double ycoeff = 0.0;

  MatrixGraph<locint>::iterator* e_iter = e->get_iterator();
  VectorGraph<locint>::iterator* r_iter = r->get_iterator();

  bool has_next = e_iter->reset();
  while(has_next) {
    // p >= q
    has_next = e_iter->get_next(p, q, pw);
    if (p != info.r and q!= info.r) {
      if (info.x != NULLLOC && info.dx != 0.0) {
        xcoeff = 1.0;
        if (info.x == p) {xcoeff += 1.0;}
        if (info.x == q) {xcoeff += 1.0;}
        tensor->increase(p, q, info.x, xcoeff * pw * info.dx);
      }
      if (info.y != NULLLOC && info.dy != 0.0) {
        ycoeff = 1.0;
        if (info.y == p) {ycoeff += 1.0;}
        if (info.y == q) {ycoeff += 1.0;}
        tensor->increase(p, q, info.y, ycoeff * pw * info.dy);
      }
    } else { // p == r or q == r;
      if (p != q) {
        if (q == info.r) {T_SWAP(p,q);}  // p == r != q
        if (info.x != NULLLOC && info.y != NULLLOC) {
          coeff = 1.0; xcoeff = 1.0; ycoeff = 1.0;
          if (info.x == q) {coeff += 1.0; xcoeff += 2.0;}
          if (info.y == q) {coeff += 1.0; ycoeff += 2.0;}
          tensor->increase(q, info.x, info.x, xcoeff * pw * info.dx * info.dx);
          tensor->increase(q, info.x, info.y, coeff * pw * info.dx * info.dy);
          tensor->increase(q, info.y, info.y, ycoeff * pw * info.dy * info.dy);
        } else if (info.x != NULLLOC) {
          xcoeff = 1.0;
          if (info.x == q) {xcoeff += 2.0;}
          tensor->increase(q, info.x, info.x, xcoeff * pw * info.dx * info.dx);
        }
      } else {  // p == q == r
        if (info.x != NULLLOC && info.dx != 0.0) {
          tensor->increase(info.x, info.x, info.x,
                           pw * info.dx * info.dx * info.dx);
          if (info.y != NULLLOC && info.dy != 0.0) {
            tensor->increase(info.x, info.x, info.y,
                             pw * info.dx * info.dx * info.dy);
            tensor->increase(info.x, info.y, info.y,
                             pw * info.dx * info.dy * info.dy);
            tensor->increase(info.y, info.y, info.y,
                             pw * info.dy * info.dy * info.dy);
  
          }
        }
      }
    }
  }

  has_next = r_iter->reset();
  while (has_next) {
    has_next = r_iter->get_next(p, pw);
    if (p != info.r) {
      if (info.x != NULLLOC && info.y != NULLLOC) {
        coeff = 1.0; xcoeff = 1.0; ycoeff = 1.0;
        if (info.x == p) {coeff += 1.0; xcoeff += 2.0;}
        if (info.y == p) {coeff += 1.0; ycoeff += 2.0;}
        tensor->increase(p, info.x, info.x, xcoeff * pw * info.pxx);
        tensor->increase(p, info.x, info.y, coeff * pw * info.pxy);
        tensor->increase(p, info.y, info.y, ycoeff * pw * info.pyy);
      } else if (info.x != NULLLOC) {
        xcoeff = 1.0;
        if (info.x == p) {xcoeff += 2.0;}
        tensor->increase(p, info.x, info.x, xcoeff * pw * info.pxx);
      }
    } else { // p == r;
      if (info.x != NULLLOC) {
        if (info.y != NULLLOC) {
          tensor->increase(info.x, info.x, info.x,
                           3.0 * pw * info.dx * info.pxx);
          tensor->increase(info.x, info.x, info.y, 
                           2.0 * pw * info.dx * info.pxy);
          tensor->increase(info.x, info.y, info.y,
                           1.0 * pw * info.dx * info.pyy);
          tensor->increase(info.y, info.x, info.x,
                           1.0 * pw * info.dy * info.pxx);
          tensor->increase(info.y, info.x, info.y,
                           2.0 * pw * info.dy * info.pxy);
          tensor->increase(info.y, info.y, info.y,
                           3.0 * pw * info.dy * info.pyy);
        } else {
          tensor->increase(info.x, info.x, info.x,
                           3.0 * pw * info.dx * info.pxx);
        }
      }
    }
  }
  if (w != 0.0) {
    if (info.x != NULLLOC && info.y != NULLLOC) {
      tensor->increase(info.x, info.x, info.x, w * info.pxxx);
      tensor->increase(info.x, info.x, info.y, w * info.pxxy);
      tensor->increase(info.x, info.y, info.y, w * info.pxyy);
      tensor->increase(info.y, info.y, info.y, w * info.pyyy);
    } else if (info.x != NULLLOC) {
      tensor->increase(info.x, info.x, info.x, w * info.pxxx);
    }
  }
  delete e_iter;
  delete r_iter;
//  tensor->debug();
}

void hyper_hessian(DerivativeInfo<locint>& info,
                  VectorGraph<locint>* adjoints,
                  MatrixGraph<locint>* hessian,
                  double w,
                  VectorGraph<locint>* r) {
//  std::cout << "In hessian" << std::endl;
  typename VectorGraph<locint>::iterator* r_iter = r->get_iterator();
  bool has_next = r_iter->reset();
  locint p;
  double pw;
  while (has_next) {
    has_next = r_iter->get_next(p, pw);
//    std::cout << p << "," << pw << std::endl;
    if (pw != 0.0) {
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
          hessian->increase(info.x, info.x, info.dx * info.dx * pw);
          hessian->increase(info.x, info.y, info.dx * info.dy * pw);
          hessian->increase(info.y, info.y, info.dy * info.dy * pw);
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
//  hessian->debug();
  delete r_iter;
}

void hyper_adjoints(DerivativeInfo<locint>& info,
                   VectorGraph<locint>* adjoints,
                   double w) {
//  std::cout << "In adjoints" << std::endl;
  if (info.x != NULLLOC) {
    adjoints->increase(info.x, w * info.dx);
  }
  if (info.y != NULLLOC) {
    adjoints->increase(info.y, w * info.dy);
  }
//  adjoints->debug();
}

void hyper_process_sac(DerivativeInfo<locint>& info,
                       int order,
                       HyperDerivative<locint>& global_gd) {
  double w = global_gd.adjoints->get_and_erase(info.r);
  VectorGraph<locint>* r = global_gd.hessian->get_and_erase(info.r);
  if (order >= 3) {
    MatrixGraph<locint>* e = global_gd.tensor->get_and_erase(info.r);
    hyper_third(info, global_gd.adjoints, global_gd.hessian, global_gd.tensor,
                w, r, e);
    delete e;
  }
  // Hessian
  hyper_hessian(info, global_gd.adjoints, global_gd.hessian, w, r);
  // Adjoints
  hyper_adjoints(info, global_gd.adjoints, w);
  delete r;
}

void hyper_hessian_gd(locint dep,
                      VectorGraph<locint>* l_a,
                      MatrixGraph<locint>* l_h,
                      double w,
                      VectorGraph<locint>* r,
                      MatrixGraph<locint>* g_h) {
  typename VectorGraph<locint>::iterator* r_iter = r->get_iterator();
  bool has_next = r_iter->reset();
  locint p;
  double pw;
  locint v;
  double vw;
  locint v2;
  double vw2;
  while (has_next) {
    has_next = r_iter->get_next(p, pw);
//    std::cout << p << "," << pw << std::endl;
    if (p != dep) {
      typename VectorGraph<locint>::iterator* a_iter = l_a->get_iterator();
      bool a_has_next = a_iter->reset();
      while (a_has_next) {
        a_has_next = a_iter->get_next(v, vw);
        if (p != v) {
          g_h->increase(p, v, pw * vw);
        } else {
          g_h->increase(p, v, 2.0 * pw * vw);
        }
      }
    } else { // p == dep
      typename VectorGraph<locint>::iterator* a_iter = l_a->get_iterator();
      bool a_has_next = a_iter->reset();
      while(a_has_next) {
        VectorGraph<locint>::iterator* a2_iter = a_iter->copy_iter();
        bool a2_has_next = a_has_next;
        a_has_next = a_iter->get_next(v, vw);
        while(a2_has_next) {
          a2_has_next = a2_iter->get_next(v2, vw2);
          g_h->increase(v, v2, pw * vw * vw2);
        }
        delete a2_iter;
      }
    }
  } // while
  if (w != 0.0) {
    typename MatrixGraph<locint>::iterator* h_iter = l_h->get_iterator();
    bool h_has_next = h_iter->reset();
    while(h_has_next) {
      h_has_next = h_iter->get_next(v, v2, vw);
      g_h->increase(v, v2, w * vw);
    }
  }
//  hessian->debug();
  delete r_iter;
}

void hyper_adjoints_gd(VectorGraph<locint>* l_a,
                       VectorGraph<locint>* g_a,
                       double w) {
  if (w != 0.0) {
    locint v;
    double vw;
    typename VectorGraph<locint>::iterator* a_iter = l_a->get_iterator();
    bool a_has_next = a_iter->reset();
    while(a_has_next) {
      a_has_next = a_iter->get_next(v, vw);
      g_a->increase(v, w * vw);
    }
  }
}

void hyper_third_gd(locint dep,
                    VectorGraph<locint>* l_a,
                    MatrixGraph<locint>* l_h,
                    HyperGraph<locint>* l_t,
                    double w,
                    VectorGraph<locint>* r,
                    MatrixGraph<locint>* e,
                    HyperGraph<locint>* g_t) {

  typename MatrixGraph<locint>::iterator* e_iter = e->get_iterator();
  bool e_has_next = e_iter->reset();
  locint p, q;
  double ww;
  locint x, y, z;
  double wx, wy, wz;
  double coeff;
  while(e_has_next) {
    e_has_next = e_iter->get_next(p, q, ww);
    if ((p != dep) && (q != dep)) {
      typename VectorGraph<locint>::iterator* a1_iter = l_a->get_iterator();
      bool a1_has_next = a1_iter->reset();
      while(a1_has_next) {
        a1_has_next = a1_iter->get_next(x, wx);
        coeff = 1.0;
        if (p == x) {coeff += 1.0;}
        if (q == x) {coeff += 1.0;}
        g_t->increase(x, p, q, coeff * ww * wx);
      }
      delete a1_iter;
    } else if ((p != dep) || (q != dep)) { // p == dep != q
      if (q == dep) {T_SWAP(p, q);}
      typename VectorGraph<locint>::iterator* a1_iter = l_a->get_iterator();
      typename VectorGraph<locint>::iterator* a2_iter;
      bool a1_has_next = a1_iter->reset();
      while(a1_has_next) {
        a2_iter = a1_iter->copy_iter();
        a1_has_next = a1_iter->get_next(x, wx);
        bool a2_has_next = true;
        while(a2_has_next) {
          a2_has_next = a2_iter->get_next(y, wy);
          coeff = 1.0;
          if (x == y) {
            if (q == x) {coeff += 2.0;}
          } else {
            if ((q == x)||(q == y)) {coeff += 1.0;}
          }
          g_t->increase(q, x, y, coeff * ww * wx * wy);
        }
        delete a2_iter;
      }
      delete a1_iter;
    } else { // p == q == dep
      typename VectorGraph<locint>::iterator* a1_iter = l_a->get_iterator();
      typename VectorGraph<locint>::iterator* a2_iter;
      typename VectorGraph<locint>::iterator* a3_iter;
      bool a1_has_next = a1_iter->reset();
      while(a1_has_next) {
        a2_iter = a1_iter->copy_iter();
        a1_has_next = a1_iter->get_next(x, wx);
        bool a2_has_next = true;
        while(a2_has_next) {
          a3_iter = a2_iter->copy_iter();
          a2_has_next = a2_iter->get_next(y, wy);
          bool a3_has_next = true;
          while(a3_has_next) {
            a3_has_next = a3_iter->get_next(z, wz);
            g_t->increase(x, y, z, ww * wx * wy * wz);
          }
          delete a3_iter;
        }
        delete a2_iter;
      }
      delete a1_iter;
    }
  }
  delete e_iter;
  
  typename VectorGraph<locint>::iterator* r_iter = r->get_iterator();
  bool r_has_next = r_iter->reset();
  while(r_has_next) {
    r_has_next = r_iter->get_next(p, ww);
    if (p != dep) {
      typename MatrixGraph<locint>::iterator* h_iter = l_h->get_iterator();
      bool h_has_next = h_iter->reset();
      while(h_has_next) {
        h_has_next = h_iter->get_next(x, y, wx);
        coeff = 1.0;
        if (p == x) {coeff += 1.0;}
        if (p == y) {coeff += 1.0;}
        g_t->increase(p, x, y, coeff * ww * wx);
      }
      delete h_iter;
    } else {
      typename MatrixGraph<locint>::iterator* h_iter = l_h->get_iterator();
      typename VectorGraph<locint>::iterator* a_iter;
      bool h_has_next = h_iter->reset();
      while(h_has_next) {
        h_has_next = h_iter->get_next(x, y, wx);
        a_iter = l_a->get_iterator();
        bool a_has_next = a_iter->reset();
        while(a_has_next) {
          a_has_next = a_iter->get_next(z, wz);
          coeff = 1.0;
          if (x == z) {coeff += 1.0;}
          if (y == z) {coeff += 1.0;}
          g_t->increase(x, y, z, coeff * ww * wx * wz);
        }
        delete a_iter;
      }
      delete h_iter;
    }
  }
  delete r_iter;
  
  bool t_has_next = l_t->reset();
  while(t_has_next) {
    t_has_next = l_t->get_next(x, y, z, ww);
    g_t->increase(x, y, z, w * ww);
  }
}

void hyper_process_gd(locint dep,
                      int order,
                      HyperDerivative<locint>& local_gd,
                      HyperDerivative<locint>& global_gd) {
#ifdef ENABLE_GENERIC_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (myid == DEBUG_ID) {
    global_gd.debug();
  }
#endif
  IF_MY_DEBUG
  std::cout << "Processing dep: " << dep << std::endl;
  END_MY_DEBUG
  double w = global_gd.adjoints->get_and_erase(dep);
  IF_MY_DEBUG
  std::cout << " w = " << w << std::endl;
  END_MY_DEBUG
  VectorGraph<locint>* r = global_gd.hessian->get_and_erase(dep);
  IF_MY_DEBUG
  r->debug();
  global_gd.debug();
  END_MY_DEBUG
  // third order tensor
  if (order >= 3) {
    MatrixGraph<locint>* e = global_gd.tensor->get_and_erase(dep);
    hyper_third_gd(dep, local_gd.adjoints, local_gd.hessian, local_gd.tensor,
                   w, r, e, global_gd.tensor);
    delete e;
  }
  // Hessian
  hyper_hessian_gd(dep, local_gd.adjoints, local_gd.hessian,
                   w, r, global_gd.hessian);
  // Adjoints
  hyper_adjoints_gd(local_gd.adjoints, global_gd.adjoints, w);
  delete r;
}

void hyper_process_recv_gd(
    locint dep,
    int order,
    HyperDerivative<locint>& local_gd,
    std::map<locint, HyperDerivative<locint> >& global_gd) {
  locint loc;
#ifdef ENABLE_GENERIC_MPI
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
  IF_MY_DEBUG
  std::cout << "dep = " << dep << std::endl;
  std::cout << "local_gd :" << std::endl;
  local_gd.debug();
  std::cout << "---------" << std::endl;
  END_MY_DEBUG
  std::map<locint, HyperDerivative<locint> >::iterator t_iter;
  t_iter = global_gd.begin();
  while(t_iter != global_gd.end()) {
    loc = t_iter->first;
    if (global_gd[loc].adjoints->has_live(dep)) {
      hyper_process_gd(dep, order, local_gd, global_gd[loc]);
    }
    ++t_iter;
  }
}

void hyper_process_recv_ind(locint toind,
                            locint fromind,
                            int count,
                            int order,
                            std::map<locint, HyperDerivative<locint> >& global_gd) {
  std::map<locint, HyperDerivative<locint> >::iterator t_iter;
  t_iter = global_gd.begin();
  while(t_iter != global_gd.end()) {
    t_iter->second.ind_map(toind, fromind, count);
    ++t_iter;
  }
}
