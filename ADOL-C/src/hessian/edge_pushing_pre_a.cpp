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
#include <adolc/adolc.h>
#include "taping_p.h"


#ifdef PREACC

#define PRE_ACC
#define ASYSMMETRIC_MATRIX 1
#include "edge_uni5_push.cpp"
#undef ASYSMMETRIC_MATRIX
#undef PRE_ACC

#endif
