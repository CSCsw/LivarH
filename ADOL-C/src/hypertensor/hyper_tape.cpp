#include <vector>
#include <map>


int hyper_tape(short tag,
               int dep,
               int indep,
               const double* basepoint,
               std::map<int, int> ind_map,
               std::vector<int>& hyper_index,
               std::vector<double>& hyper_value) {

  unsigned char opcode;
  locint size = 0;
  locint res = 0;
  locint arg = 0;
  locint arg1 = 0;
  locint arg2 = 0;
  double coval = 0;
  
  init_for_sweep(tnum);

  if ((depcheck != ADOLC_CURRENT_TAPE_INFOS.stats[NUM_DEPENDENTS]) ||
     (indcheck != ADOLC_CURRENT_TAPE_INFOS.stats[NUM_INDEPENDENTS]) ) {
    fprintf(DIAG_OUT,"ADOL-C error: Tape_doc on tape %d  aborted!\n",tnum);
    fprintf(DIAG_OUT,"Number of dependent (%d) and/or independent (%d) "
            "variables passed to Tape_doc is\ninconsistent with "
            "number recorded on tape %d (%d:%d)\n", depcheck,
            indcheck, tnum, (int)ADOLC_CURRENT_TAPE_INFOS.stats[NUM_DEPENDENTS],
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
        arg = get_locinf_f();
        res = get_locint_f();
        dp_T0[res] = dp_T0[arg];
        hyper_index.push_back(TRANSLATE_ARG(arg));
        hyper_index.push_back(TRANSLATE_RES(res));
        hyper_value.push_back(dp_T0[arg]);
        hyper_value.push_back(dp_T0[res]);
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
        hyper_index.push_back(TRANSLATE_ 
      default:
        fprintf(DIAG_OUT, "HYPER-TENSOR: unimplemented opcode %d\n", opcode);
    }
  }
}
