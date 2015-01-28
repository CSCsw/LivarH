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

int hyper_third_reverse(short tag,
                        std::vector<locint>& hyper_index,
                        std::vector<double>& hyper_value,
                        VectorGraph<locint>* adjoints,
                        MatrixGraph<locint>* hessian,
                        HyperGraph<locint>* tensor) {
  std::cout << "In hyper_third_reverse " << std::endl;
  unsigned char opcode;
  locint size = 0;
  locint res = 0;
  locint arg = 0;
  locint arg1 = 0;
  locint arg2 = 0;
  double coval = 0;
  double* d = NULL;
  int i;

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
      default:
        fprintf(DIAG_OUT, "HYPER-TENSOR: unimplemented opcode %d\n", opcode);
    }
    opcode = get_op_r();    
  }
  end_sweep();
}
