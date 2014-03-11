#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_graph.h>
#include <adolc/hessian/edge_check.h>
#include <adolc/hessian/edge_tape.h>
#include <cmath>
#include "oplate.h"
#include <adolc/adalloc.h>
#include <adolc/interfaces.h>
#include "taping_p.h"

extern int inc_edge_count;
extern int del_edge_count;
extern void (*increase_edge)(int ,int ,double , map<int, map<int,double> >*);

#ifdef NO_ASSIGN_BYPASS

#define ASYSMMETRIC_MATRIX 1
#define PRE_ACC 1
#include "edge_uni5_push.cpp"
#undef ASYSMMETRIC_MATRIX
#undef PRE_ACC

#endif
