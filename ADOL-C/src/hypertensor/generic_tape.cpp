#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include "oplate.h"
#include <adolc/adolc.h>
#include <adolc/adalloc.h>
#include <adolc/interfaces.h>
#include "taping_p.h"
#include <adolc/hypertensor/hyper_common.h>
#include <adolc/hypertensor/generic_tape.h>
#include <adolc/hypertensor/generic_mpi_trace.h>
#include "mpi.h"
/*
#define TRANSLATE_ARG(arg) (index_translate[arg])
#define TRANSLATE_RES(res) (translate_result(index_translate, max_ind, res))

#define TRANSLATE_IND(ind) TRANSLATE_RES(ind)
#define TRANSLATE_DEP(dep) TRANSLATE_ARG(dep)
*/

#define TRANSLATE_ARG(arg) (arg + max_ind)
#define TRANSLATE_RES(res) (res + max_ind)
#define TRANSLATE_IND(ind) (ind + max_ind)
#define TRANSLATE_DEP(dep) (dep + max_ind)

#define TEMP_INDEX (NULLLOC - 1)

extern std::vector<SRinfo> sr_stack;
extern std::vector<double> dummy_ind_vec;

/*
static locint translate_result(std::map<locint, locint>& index_translate,
                               locint& max_ind,
                               locint res) {
  locint ret = max_ind;
  ++max_ind;
  index_translate[res] = ret;
  return ret;
}
*/
int generic_tape(short tag,
                 int depcheck,
                 int indcheck,
                 const double* basepoint,
                 std::map<locint, locint>& ind_map,
                 std::vector<locint>& hyper_index,
                 std::vector<double>& hyper_value) {
  std::cout << "In generic_tape" << std::endl;
  unsigned char opcode;
  locint size = 0;
  locint res = 0;
  locint arg = 0;
  locint arg1 = 0;
  locint arg2 = 0;
  double coval = 0;
  double* d = NULL;
  int i;

  ADOLC_OPENMP_THREAD_NUMBER;
  ADOLC_OPENMP_GET_THREAD_NUMBER;
  locint index_ind = 0;
  locint max_ind = 0;
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  max_ind = myid * LOCINT_PER_PROC;

  std::map<locint, locint> index_translate; 
  std::vector<SRinfo>::iterator sr_iter;
  sr_iter = sr_stack.begin();
  init_for_sweep(tag);

  if (indcheck+dummy_ind_vec.size() != ADOLC_CURRENT_TAPE_INFOS.stats[NUM_INDEPENDENTS]) {
    fprintf(DIAG_OUT,"ADOL-C error: Tape_doc on tape %d  aborted!\n",tag);
    fprintf(DIAG_OUT,"Number of dependent (%d) and/or independent (%d) "
            "variables passed to Tape_doc is\ninconsistent with "
            "number recorded on tape %d (%d:%d)\n", depcheck,
            indcheck, tag, (int)ADOLC_CURRENT_TAPE_INFOS.stats[NUM_DEPENDENTS],
            (int)ADOLC_CURRENT_TAPE_INFOS.stats[NUM_INDEPENDENTS]);
    exit (-1);
  }

  double *dp_T0 = NULL;
  dp_T0 = new double[ADOLC_CURRENT_TAPE_INFOS.stats[NUM_MAX_LIVES]];
  opcode = get_op_f();
  while(opcode != end_of_tape) {
    switch (opcode) {
      // CONTROL OPCODES
      case end_of_op:
        get_op_block_f();
        opcode = get_op_f();
        break;
      case end_of_int:
        get_loc_block_f();
        break;
      case end_of_val:
        get_val_block_f();
        break;
      case start_of_tape:
      case end_of_tape:
        break;

      // COMPARISON
      case eq_zero:
      case neq_zero:
      case le_zero:
      case gt_zero:
      case ge_zero:
      case lt_zero:
        arg = get_locint_f();
        break;
      
      // ASSIGNMENTS
      case assign_a:
        arg = get_locint_f();
        res = get_locint_f();
        dp_T0[res] = dp_T0[arg];
//        std::cout << res << " <--- " << arg << std::endl;
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[arg]);
        hyper_value.push_back(dp_T0[res]);
//        std::cout << index_translate[res] << " <--- " << index_translate[arg] << std::endl;
        break;
      case assign_d:
        res = get_locint_f();
        coval = get_val_f();
        dp_T0[res] = coval;
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case assign_d_one:
        res = get_locint_f();
        dp_T0[res] = 1.0;
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case assign_d_zero:
        res = get_locint_f();
        dp_T0[res] = 0.0;
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case assign_ind:
        res = get_locint_f();
        if (index_ind < indcheck) {
          dp_T0[res] = basepoint[index_ind];
        } else if (index_ind < indcheck + dummy_ind_vec.size()){
          dp_T0[res] = dummy_ind_vec[index_ind - indcheck];
        } else {
          std::cout << "FATAL error in dummy independent!" << std::endl;
        }
        index_ind++;
        hyper_index.push_back(TRANSLATE_IND(res));
        hyper_value.push_back(dp_T0[res]);
//        std::cout << res << " I--> " << index_translate[res] << " = " << dp_T0[res] << std::endl;
        break;
      case assign_dep:
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_DEP(res));
        hyper_value.push_back(dp_T0[res]);
//        std::cout << res << " D--> " << hyper_index.back() << std::endl;
        break;
      case eq_plus_d:
        res = get_locint_f();
        coval = get_val_f();
        dp_T0[res] += coval;
        break;
      case eq_plus_a:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] += dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case eq_min_d:
        res = get_locint_f();
        coval = get_val_f();
        dp_T0[res] -= coval;
        break;
      case eq_min_a:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] -= dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case eq_mult_d:
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(NULLLOC);
        hyper_value.push_back(coval);
        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] *= coval;
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case eq_mult_a:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] *= dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case incr_a:
        res = get_locint_f();
//        std::cout << "D[" << res << "]" << dp_T0[res] << std::endl;
        dp_T0[res]++;
//        std::cout << "D[" << res << "]" << dp_T0[res] << std::endl;
        break;
      case decr_a:
        res = get_locint_f();
        dp_T0[res]--;
        break;
      case plus_a_a:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        dp_T0[res] = dp_T0[arg1] + dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
//        std::cout << arg1 << " ---> " << index_translate[arg1] << std::endl;
//        std::cout << arg2 << " ---> " << index_translate[arg2] << std::endl;
//        std::cout << res << " ---> " << index_translate[res] << std::endl;
        break;
      case plus_d_a:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = dp_T0[arg] + coval;
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case min_a_a:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        dp_T0[res] = dp_T0[arg1] - dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case min_d_a:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = coval - dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case mult_a_a:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
//        std::cout << res << " <--- " << arg1 << " * " << arg2 << std::endl;
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        dp_T0[res] = dp_T0[arg1] * dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case mult_d_a:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
//        std::cout << res << " <--- " << arg << " * " << coval << std::endl;
        hyper_index.push_back(NULLLOC);
        hyper_value.push_back(coval);
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = coval * dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case div_a_a:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        dp_T0[res] = dp_T0[arg1] / dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case div_d_a:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(NULLLOC);
        hyper_value.push_back(coval);
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = coval / dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case eq_plus_prod:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        hyper_index.push_back(TEMP_INDEX);
        hyper_value.push_back(dp_T0[arg1] * dp_T0[arg2]);

        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] += dp_T0[arg1] * dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case eq_min_prod:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        hyper_index.push_back(TEMP_INDEX);
        hyper_value.push_back(dp_T0[arg1] * dp_T0[arg2]);

        hyper_index.push_back(TRANSLATE_ARG(res));
        hyper_value.push_back(dp_T0[res]);
        dp_T0[res] -= dp_T0[arg1] * dp_T0[arg2];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case pos_sign_a:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case neg_sign_a:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = -dp_T0[arg];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case exp_op:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = exp(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case log_op:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = log(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case pow_op:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(NULLLOC);
        hyper_value.push_back(coval);
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = pow(dp_T0[arg], coval);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case sqrt_op:
        arg = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = sqrt(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case sin_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = sin(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case cos_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
//        std::cout << res << " = cos " << arg1 << std::endl;
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = cos(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case atan_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = atan(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case asin_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = asin(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case acos_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = acos(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case min_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        if (dp_T0[arg1] > dp_T0[arg2])
          dp_T0[res] = dp_T0[arg2];
        else
          dp_T0[res] = dp_T0[arg1];
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case abs_val:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        dp_T0[res] = fabs(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case ceil_op:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        dp_T0[res] = ceil(dp_T0[arg]);
        break;
      case floor_op:
        arg = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        dp_T0[res] = floor(dp_T0[arg]);
        break;
      case cond_assign:
        arg = get_locint_f();
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_ARG(arg2));
        hyper_value.push_back(dp_T0[arg2]);
        if (dp_T0[arg]>0) {
          if (coval<=0){
            fprintf(DIAG_OUT,"Inconsistency in cond_assign. Retape?\n");
          }
          dp_T0[res]=dp_T0[arg1];
        } else {
          if (coval>0){
            fprintf(DIAG_OUT,"Inconsistency in cond_assign. Retape?\n");
          }
          dp_T0[res]=dp_T0[arg2];
        }
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case cond_assign_s:
        arg = get_locint_f();
        arg1 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_value.push_back(dp_T0[arg]);
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        if (dp_T0[arg] > 0 ) {
          if (coval <= 0) {
            fprintf(DIAG_OUT,"Inconsistency in cond_assign_s. Retape?\n");
          }
          dp_T0[res] = dp_T0[arg1];
        }
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;

#ifdef ATRIG_ERF
      case asinh_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = asinh(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case acosh_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = acosh(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case atanh_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = atanh(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
      case erf_op:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        hyper_index.push_back(TRANSLATE_ARG(arg1));
        hyper_value.push_back(dp_T0[arg1]);
        dp_T0[res] = erf(dp_T0[arg1]);
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[res]);
        break;
        break;
#endif // ATRIG_ERF

#ifdef ADOLC_ADVANCED_BRANCHING
      case neq_a_a:
      case eq_a_a:
      case le_a_a:
      case ge_a_a:
      case lt_a_a:
      case gt_a_a:
        break;
#endif // ADOLC_ADVANCED_BRANCHING
      // Control opcodes
      case take_stock_op:
        size = get_locint_f();
        res = get_locint_f();
        d = get_val_v_f(size);
        for(i = 0; i < size; i++) {
          dp_T0[res+i] = d[i];
        }
        break;
      case death_not:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        break;
      case gen_quad:
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res = get_locint_f();
        coval = get_val_f();
        if (coval != dp_T0[arg1]) {
          fprintf(DIAG_OUT, "HYPER-TENSOR: forward sweep aborted; tape invalid!\n");
        }
        coval = get_val_f();
      case ext_diff:
        get_locint_f();
        get_locint_f();
        get_locint_f();
        get_locint_f();
        get_locint_f();
        get_locint_f();
      case ignore_me:
        break;
      case ampi_send:
        if (sr_iter->SR_tag == RMPI_SEND_TAG) {
          res = sr_iter->loc;
          sr_iter->loc = TRANSLATE_ARG(res); 
        } else {
          std::cout << "GENERIC MPI: Send trace error." << std::endl;
        }
        sr_iter++;  
        break;
      case ampi_recv:
        if (sr_iter->SR_tag == RMPI_RECV_TAG) {
          res = sr_iter->loc;
          sr_iter->loc = TRANSLATE_ARG(res); 
        } else {
          std::cout << "GENERIC MPI: Recv trace error." << std::endl;
        }
        sr_iter++;
        break;
      default:
        fprintf(DIAG_OUT, "GENERIC-TENSOR: unimplemented opcode %d\n", opcode);
    }
    opcode = get_op_f();
  }
  end_sweep();
}
