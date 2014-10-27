#include <iostream>
#include <vector>
#include <cmath>
#include "oplate.h"
#include <adolc/adalloc.h>
#include <adolc/interfaces.h>
#include "taping_p.h"
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_check.h>

//#define BUG_TEST

using namespace std;


#define EDGE_TRANSLATE_INDEX    \
{                                                                \
  edge_t_new_len = edge_value_len - 1;                           \
  if (edge_t_new_len >=edge_t_old_len) {                        \
    unsigned int i;                                              \
    for(i = edge_t_old_len; i < edge_t_new_len; i++){            \
      if (edge_index[i] != NULLLOC)                              \
        edge_index[i] = edge_t_index[edge_index[i]];             \
    }                                                            \
    edge_t_index[edge_index[edge_t_new_len]] = edge_t_cur_index; \
    edge_index[edge_t_new_len] = edge_t_cur_index;               \
    edge_t_cur_index++;                                          \
  }                                                              \
  edge_t_old_len=edge_value_len;                                 \
}


#define EDGE_TRANSLATE_INDEPENDENT  \
{                                                                \
  edge_t_new_len = edge_value_len - 1;                           \
  edge_t_index[edge_index[edge_t_new_len]] = edge_t_cur_index;   \
  edge_index[edge_t_new_len] = edge_t_cur_index;                 \
  indmap[edge_t_cur_index] = indexi++;                           \
  edge_t_cur_index++;                                            \
  edge_t_old_len = edge_value_len;                               \
}


//For hessian, there is a unique dependent variable
#define EDGE_TRANSLATE_DEPENDENT  \
{                                                                         \
  edge_t_new_len = edge_value_len - 1;                                    \
  edge_index[edge_t_new_len] = edge_t_index[edge_index[edge_t_new_len]];  \
}



int edge_tape(short tnum,                         /* tape id */
              int depcheck,                       /* consistency chk on # of dependents */
              int indcheck,                       /* consistency chk on # of independents */
              const double*   basepoint,          /* independent variable values */
              unsigned int**  indmap_p,           /* Mapping from location to independent index */
              locint**        edge_index_p,       /* The translated index (increasing) */
              double**        edge_value_p,       /* The value of variables (for reverse) */
              unsigned int*   edge_index_len_p,   /* The length of edge_index[] */
              unsigned int*   edge_value_len_p,   /* The length of edge_value[] */
              locint*         max_index_p)        /* The max index (after translation)    */
{
    unsigned int edge_t_new_len;
    unsigned int edge_t_old_len;
    unsigned int *edge_t_index;
    locint edge_t_cur_index;
    int ret_val = 1;
    int max_tmp;
    unsigned char operation;
    unsigned int* indmap;
    locint *edge_index;
    double *edge_value;
    unsigned int edge_index_len = 0;
    unsigned int edge_value_len = 0;
    locint max_index = 0;
    locint tmp_index = 0;

    locint size = 0;
    locint res  = 0;
    locint arg  = 0;
    locint arg1 = 0;
    locint arg2 = 0;

    double coval = 0, *d = 0;

    int indexi = 0;

    /****************************************************************************/
    /*                                                                    INITs */

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

//Sweep the tape once, get the edge_tape_size
//Is there a better way to do this?
    int edge_op_cnt=0;
    unsigned int edge_tape_size=0;

    operation=get_op_f();
    while (operation !=end_of_tape) {
      switch (operation){
            case end_of_op:                                          /* end_of_op */
                get_op_block_f();
                operation=get_op_f();
                break;
            case assign_d:            /* assign an adouble variable a    assign_d */
            case assign_d_one:    /* assign an adouble variable a    assign_d_one */
            case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
            case assign_ind:       /* assign an adouble variable an    assign_ind */
            case assign_dep:           /* assign a float variable a    assign_dep */
            case neq_a_a:
            case eq_a_a:
            case le_a_a:
            case ge_a_a:
            case lt_a_a:
            case gt_a_a:
                edge_tape_size+=1;
                edge_op_cnt++;
                break;
            case assign_a:           /* assign an adouble variable an    assign_a */
            case plus_d_a:             /* Add an adouble and a double    plus_d_a */
            case min_d_a:                /* Subtract an adouble from a    min_d_a */
            case pos_sign_a:                                        /* pos_sign_a */
            case neg_sign_a:                                        /* neg_sign_a */
            case exp_op:                          /* exponent operation    exp_op */
            case log_op:                                                /* log_op */
            case sqrt_op:                                              /* sqrt_op */
            case sin_op:                              /* sine operation    sin_op */
            case cos_op:                            /* cosine operation    cos_op */
            case atan_op:                                              /* atan_op */
            case asin_op:                                              /* asin_op */
            case acos_op:                                              /* acos_op */
            case abs_val:                                              /* abs_val */
            case asinh_op:                                            /* asinh_op */	
            case acosh_op:                                            /* acosh_op */
            case atanh_op:                                            /* atanh_op */
            case erf_op:                                                /* erf_op */
                edge_tape_size+=2;
                edge_op_cnt++;
                break;
            case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
            case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
            case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
            case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
            case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
            case min_a_a:              /* Subtraction of two adoubles     min_a_a */
            case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
            case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
            case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
            case div_d_a:             /* Division double - adouble (/)    div_d_a */
            case pow_op:                                                /* pow_op */
            case min_op:                                                /* min_op */
            case cond_assign_s:                                  /* cond_assign_s */
                edge_tape_size+=3;
                edge_op_cnt++;
                break;
            case cond_assign:                                      /* cond_assign */
                edge_tape_size+=4;
                edge_op_cnt++;
                break;
            case eq_plus_prod:    /* Add an product to an            eq_plus_prod */
            case eq_min_prod:     /* Subtract an product from an      eq_min_prod */
                edge_tape_size+=5;
                edge_op_cnt+=2;
                break;
            case subscript_ref:
            case ref_incr_a:
            case ref_decr_a:
            case ref_eq_plus_d:
            case ref_eq_min_d:
                edge_op_cnt++;
                break;
            case ref_assign_ind:
                edge_tape_size++;
                edge_op_cnt++;
                break;
            case subscript:
            case ref_assign_d:
            case ref_assign_d_zero:
            case ref_assign_d_one:
            case ref_assign_a:
                edge_tape_size+=2;
                edge_op_cnt++;
                break;
            case ref_eq_plus_a:
            case ref_eq_min_a:
            case ref_eq_mult_d:
            case ref_eq_mult_a:
                edge_tape_size+=3;
                edge_op_cnt++;
                break;
            case ref_copyout:
                edge_tape_size+=2;
                edge_op_cnt++;
                break;
            case ref_cond_assign:
                edge_tape_size+=4;
                edge_op_cnt++;
                break;
            case ref_cond_assign_s:
                edge_tape_size+=3;
                edge_op_cnt++;
                break;
            default:
                edge_op_cnt++;
                break;
        }
        operation=get_op_f();
    }
    end_sweep();
#ifdef EDGE_DEBUG
    printf("edge_tape_size = %d\n",edge_tape_size);
    printf("edge_op_size = %d\n", edge_op_cnt);
    fflush(stdout);
#endif  // EDGE_DEBUG

    tmp_index=edge_tape_size+1;
    edge_tape_size+=10;
    max_index=edge_tape_size;
    int i;
    indmap=new unsigned int[edge_tape_size];
    for(i=0;i<edge_tape_size;i++){indmap[i]=0;}

//For translating the ADOLC locint into monotonic indexing
    edge_t_old_len=0;
    edge_t_new_len=0;
    edge_t_cur_index=0;
    edge_t_index=new locint[edge_tape_size];
    for(i=0;i<edge_tape_size;i++){edge_t_index[i]=0;}

    /* globals */
    edge_value=new double[edge_tape_size];
    edge_index=new locint[edge_tape_size];
    edge_value_len=0;
    edge_index_len=0;


//Because dp_T0 works on ADOL-C locint indexint
    double *dp_T0=NULL;
//    dp_T0 = myalloc1(ADOLC_CURRENT_TAPE_INFOS.stats[NUM_MAX_LIVES]);
    dp_T0=new double[edge_tape_size];

//This sweep stores index&value
    init_for_sweep(tnum);
    operation=get_op_f();
    while (operation !=end_of_tape) {
        switch (operation) {
                /****************************************************************************/
                /*                                                                  MARKERS */
                /*--------------------------------------------------------------------------*/
            case end_of_op:                                          /* end_of_op */
                get_op_block_f();
                operation=get_op_f();
                /* Skip next operation, it's another end_of_op */
                break;

                /*--------------------------------------------------------------------------*/
            case end_of_int:                                        /* end_of_int */
                get_loc_block_f();
                break;

                /*--------------------------------------------------------------------------*/
            case end_of_val:                                        /* end_of_val */
                get_val_block_f();
                break;

                /*--------------------------------------------------------------------------*/
            case start_of_tape:                                  /* start_of_tape */
                break;

                /*--------------------------------------------------------------------------*/
            case end_of_tape:                                      /* end_of_tape */
                break;

                /****************************************************************************/
                /*                                                               COMPARISON */
                /*--------------------------------------------------------------------------*/
//Also No Check Here. May need Retape
//TO DO
            case eq_zero  :                                            /* eq_zero */
            case neq_zero :                                           /* neq_zero */
            case le_zero  :                                            /* le_zero */
            case gt_zero  :                                            /* gt_zero */
            case ge_zero  :                                            /* ge_zero */
            case lt_zero  :                                            /* lt_zero */
                arg  = get_locint_f();
                break;

                /****************************************************************************/
                /*                                                              ASSIGNMENTS */
                /*--------------------------------------------------------------------------*/
            case assign_a:           /* assign an adouble variable an    assign_a */
                /* adouble value. (=) */
                arg = get_locint_f();
                res = get_locint_f();
                dp_T0[res]= dp_T0[arg];
                edge_index[edge_index_len++]=arg;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case assign_d:            /* assign an adouble variable a    assign_d */
                /* double value. (=) */
                res  = get_locint_f();
                coval=get_val_f();
                dp_T0[res]= coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case assign_d_one:    /* assign an adouble variable a    assign_d_one */
                /* double value. (1) (=) */
                res  = get_locint_f();
                dp_T0[res]= 1.0;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
                /* double value. (0) (=) */
                res  = get_locint_f();
                dp_T0[res]= 0.0;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case assign_ind:       /* assign an adouble variable an    assign_ind */
                /* independent double value (<<=) */
                res  = get_locint_f();
                dp_T0[res]= basepoint[indexi];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEPENDENT
                break;

                /*--------------------------------------------------------------------------*/
            case assign_dep:           /* assign a float variable a    assign_dep */
                /* dependent adouble value. (>>=) */
                res = get_locint_f();
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_DEPENDENT
                break;


                /****************************************************************************/
                /*                                                   OPERATION + ASSIGNMENT */

                /*--------------------------------------------------------------------------*/
            case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
                /* adouble. (+=) */
                res   = get_locint_f();
                coval = get_val_f();
                dp_T0[res] += coval;
                //EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
                /* adouble. (+=) */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]+= dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_plus_prod:    /* Add an product to an            eq_plus_prod */
                /* adouble. (+= x1*x2) */
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                edge_index[edge_index_len++]=tmp_index;
                edge_value[edge_value_len++]=dp_T0[arg1]*dp_T0[arg2];
                EDGE_TRANSLATE_INDEX

                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res] += dp_T0[arg1]*dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                /* adouble. (-=) */
                res   = get_locint_f();
                coval = get_val_f();
                dp_T0[res] -= coval;
                //EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
                /* adouble. (-=) */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]-= dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_min_prod:     /* Subtract an product from an      eq_min_prod */
                /* adouble. (-= x1*x2) */
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                edge_index[edge_index_len++]=tmp_index;
                edge_value[edge_value_len++]=dp_T0[arg1]*dp_T0[arg2];
                EDGE_TRANSLATE_INDEX

                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res] -= dp_T0[arg1]*dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
                /* flaoting point. (*=) */
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=NULLLOC;
                edge_value[edge_value_len++]=coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res] *= coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
                /* (*=) */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]*= dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
                
                /*--------------------------------------------------------------------------*/
            case incr_a:                        /* Increment an adouble    incr_a */
                res = get_locint_f();
//Should do nothing... Is that right?
                dp_T0[res]++;
                break;

                /*--------------------------------------------------------------------------*/
            case decr_a:                        /* Increment an adouble    decr_a */
                res = get_locint_f();
//Should do nothing... Is that right?
                dp_T0[res]--;
                break;

                /****************************************************************************/
                /*                                                        BINARY OPERATIONS */
                /*--------------------------------------------------------------------------*/
            case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                dp_T0[res]=dp_T0[arg1]+dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case plus_d_a:             /* Add an adouble and a double    plus_d_a */
                /* (+) */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= dp_T0[arg] + coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case min_a_a:              /* Subtraction of two adoubles     min_a_a */
                /* (-) */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                dp_T0[res]=dp_T0[arg1]-dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case min_d_a:                /* Subtract an adouble from a    min_d_a */
                /* double (-) */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]  = coval - dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                dp_T0[res]=dp_T0[arg1]*dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                /* (*) */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=NULLLOC;
                edge_value[edge_value_len++]=coval;
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]  = coval * dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                /* (/) */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                dp_T0[res]=dp_T0[arg1]/dp_T0[arg2];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case div_d_a:             /* Division double - adouble (/)    div_d_a */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=NULLLOC;
                edge_value[edge_value_len++]=coval;
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]  = coval / dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;


                /****************************************************************************/
                /*                                                         SIGN  OPERATIONS */
                /*--------------------------------------------------------------------------*/
            case pos_sign_a:                                        /* pos_sign_a */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case neg_sign_a:                                        /* neg_sign_a */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= -dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /****************************************************************************/
                /*                                                         UNARY OPERATIONS */

                /*--------------------------------------------------------------------------*/
            case exp_op:                          /* exponent operation    exp_op */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= exp(dp_T0[arg]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case log_op:                                                /* log_op */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= log(dp_T0[arg]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case pow_op:                                                /* pow_op */
                arg  = get_locint_f();
                res  = get_locint_f();
                coval   = get_val_f();
                edge_index[edge_index_len++]=NULLLOC;
                edge_value[edge_value_len++]=coval;
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res] = pow(dp_T0[arg],coval);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case sqrt_op:                                              /* sqrt_op */
                arg  = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]= sqrt(dp_T0[arg]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case sin_op:                              /* sine operation    sin_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                dp_T0[arg2]= cos(dp_T0[arg1]);
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = sin(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case cos_op:                            /* cosine operation    cos_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                dp_T0[arg2]= sin(dp_T0[arg1]);
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = cos(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case atan_op:                                              /* atan_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = atan(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case asin_op:                                              /* asin_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = asin(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case acos_op:                                              /* acos_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = acos(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

#ifdef ATRIG_ERF

                /*--------------------------------------------------------------------------*/
            case asinh_op:                                            /* asinh_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = asinh(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case acosh_op:                                           /* acosh_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = acosh(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case atanh_op:                                            /* atanh_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = atanh(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case erf_op:                                                /* erf_op */
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                res  = get_locint_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res] = erf(dp_T0[arg1]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

#endif

#ifdef ADOLC_ADVANCED_BRANCHING
//We do not check the consistancy here
//Different inputs would destroy the consistancy! Retape is needed!
//TO DO: 
            case neq_a_a:
            case eq_a_a:
            case le_a_a:
            case ge_a_a:
            case lt_a_a:
            case gt_a_a:
                coval =	get_val_f();
                arg   = get_locint_f();
                arg1  = get_locint_f();
                res   = get_locint_f();
                dp_T0[res]=coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
#endif

                /*--------------------------------------------------------------------------*/
            case min_op:                                                /* min_op */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];
                if (dp_T0[arg1] > dp_T0[arg2])
                    dp_T0[res] = dp_T0[arg2];
                else
                    dp_T0[res] = dp_T0[arg1];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
		        EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case abs_val:                                              /* abs_val */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]  = fabs(dp_T0[arg]);
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /*--------------------------------------------------------------------------*/
            case ceil_op:                                              /* ceil_op */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                dp_T0[res]  = ceil(dp_T0[arg]);
                break;

                /*--------------------------------------------------------------------------*/
            case floor_op:                 /* Compute ceil of adouble    floor_op */
                arg   = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                dp_T0[res]  = floor(dp_T0[arg]);
                break;


                /****************************************************************************/
                /*                                                             CONDITIONALS */

                /*--------------------------------------------------------------------------*/
            case cond_assign:                                      /* cond_assign */
                arg   = get_locint_f();
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                edge_index[edge_index_len++]=arg2;
                edge_value[edge_value_len++]=dp_T0[arg2];

                if (dp_T0[arg]>0){
                    if (coval<=0){
                        fprintf(DIAG_OUT,"Inconsistency in cond_assign. Retape?\n");
                    }
                    dp_T0[res]=dp_T0[arg1];
                }
                else{
                    if (coval>0){
                        fprintf(DIAG_OUT,"Inconsistency in cond_assign. Retape?\n");
                    }
                    dp_T0[res]=dp_T0[arg2];
                }
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;


                /*--------------------------------------------------------------------------*/
            case cond_assign_s:                                  /* cond_assign_s */
                arg   = get_locint_f();
                arg1  = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                if (dp_T0[arg]>0) {
                    if (coval<=0){
                        fprintf(DIAG_OUT,"Inconsistency in cond_assign_s. Retape?\n");
                    }
                    dp_T0[res]=dp_T0[arg1];
                }
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
/*
advector ops
Probably buggy :(
*/
            // advector[adouble] as r-value
            // In the constructor of advector(n), the locations are ensured to be contiguous
            // The implementation is copying it out via something like assign_a
            case subscript:                                             /*   advector[adouble]  */
                coval = get_val_f();
                arg   = get_locint_f();
                res   = get_locint_f();
                arg1  = get_locint_f();
                {
                    locint idx,numvar=(locint)trunc(fabs(coval));
                    idx=dp_T0[arg];
                    if (idx>=numvar){
                        fprintf(DIAG_OUT, "ADOL-C warning: index out of bounds while subscripting n=%z, idx=%z\n", numvar, idx); 
                    }
                    arg1=arg1+idx;
                    edge_index[edge_index_len++]=arg1;
                    edge_value[edge_value_len++]=dp_T0[arg1];
                    dp_T0[res]=dp_T0[arg1];
                    edge_index[edge_index_len++]=res;
                    edge_value[edge_value_len++]=dp_T0[res];
                    EDGE_TRANSLATE_INDEX
                }
		        break;

            // advector[adouble] as l-value?
            // return adubref
            // afterwards, dp_T0[res] stores the location it references to
            // then other ref_* op read it
            case subscript_ref:						/* &advector[adouble]  */
                coval = get_val_f();
                arg   = get_locint_f();
                res   = get_locint_f();
                arg1  = get_locint_f();
                {
                    locint idx,numvar=(locint)trunc(fabs(coval));
                    idx=dp_T0[arg];
                    if (idx>=numvar){
                        fprintf(DIAG_OUT, "ADOL-C warning: index out of bounds while subscripting n=%z, idx=%z\n", numvar, idx); 
                    }
                    arg1=arg1+idx;
                    dp_T0[res]=arg1;
                }
                break;

            case ref_assign_d_zero:
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]=0.0;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=0.0;
                EDGE_TRANSLATE_INDEX;
                break;

            case ref_assign_d_one:
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]=1.0;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=1.0;
                EDGE_TRANSLATE_INDEX
                break;
		
            case ref_assign_d:
                coval = get_val_f();
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]=coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=coval;
                EDGE_TRANSLATE_INDEX
                break;
                
            case ref_assign_a:
                arg  = get_locint_f();
                arg1 = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg1]));
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                dp_T0[res]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
                
            case ref_assign_ind:
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]=basepoint[indexi];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEPENDENT
                break;
                
            case ref_incr_a:
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]++;
                break;
            
            case ref_decr_a:
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]--;
                break;
                
            case ref_eq_plus_d:
                coval=get_val_f();
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]+=coval;
                break;
                
            case ref_eq_min_d:
                coval=get_val_f();
                arg = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                dp_T0[res]-=coval;
                break;
                
            case ref_eq_plus_a:
                arg = get_locint_f();
                arg1 = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg1]));
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]+=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
                
            case ref_eq_min_a:
                arg = get_locint_f();
                arg1 = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg1]));
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]-=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
            
            case ref_eq_mult_d:
                coval = get_val_f();
                arg  = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg]));
                edge_index[edge_index_len++]=NULLLOC;
                edge_value[edge_value_len++]=coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]*=coval;
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;
		
            case ref_eq_mult_a:
                arg = get_locint_f();
                arg1 = get_locint_f();
                res=(locint)trunc(fabs(dp_T0[arg1]));
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                dp_T0[res]*=dp_T0[arg];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;		

//works like assign_a
//(adub)advector[adouble]?
            case ref_copyout:
                arg = get_locint_f();
                res = get_locint_f();
                arg1=(locint)trunc(fabs(dp_T0[arg]));
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                dp_T0[res]=dp_T0[arg1];
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;		
//Store both the cond.loc and the value of that?
//What if they are inconsistant? Retape?
            case ref_cond_assign:
                arg = get_locint_f();
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                coval= get_val_f();
                {
                    locint ref = get_locint_f();
                    res=(locint)trunc(fabs(dp_T0[ref]));
                    edge_index[edge_index_len++]=arg;
                    edge_value[edge_value_len++]=dp_T0[arg];
                    edge_index[edge_index_len++]=arg1;
                    edge_value[edge_value_len++]=dp_T0[arg1];
                    edge_index[edge_index_len++]=arg2;
                    edge_value[edge_value_len++]=dp_T0[arg2];
                    if (dp_T0[arg]>0.0){
                        if (coval<=0.0){
                            fprintf(DIAG_OUT,"Inconsistancy in ref_cond_assign, Retape?\n");
                        }
                        dp_T0[res]=dp_T0[arg1];
                    }
                    else{
                        if (coval>0.0){
                            fprintf(DIAG_OUT,"Inconsistancy in ref_cond_assign, Retape?\n");
                        }
                        dp_T0[res]=dp_T0[arg2];
                    }
                    edge_index[edge_index_len++]=res;
                    edge_value[edge_value_len++]=dp_T0[res];
                }
                EDGE_TRANSLATE_INDEX
                break;

            case ref_cond_assign_s:
                arg = get_locint_f();
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                coval = get_val_f();
                res=(locint)trunc(fabs(dp_T0[arg2]));
                edge_index[edge_index_len++]=arg;
                edge_value[edge_value_len++]=dp_T0[arg];
                edge_index[edge_index_len++]=arg1;
                edge_value[edge_value_len++]=dp_T0[arg1];
                if (dp_T0[arg]>0.0){
                    if (coval<=0.0){
                        fprintf(DIAG_OUT, "Inconsistancy in ref_cond_assign_s, Retape?\n");
                    }
                    dp_T0[res]=dp_T0[arg1];
                }
                edge_index[edge_index_len++]=res;
                edge_value[edge_value_len++]=dp_T0[res];
                EDGE_TRANSLATE_INDEX
                break;

                /****************************************************************************/
                /*                                                          REMAINING STUFF */
                /*--------------------------------------------------------------------------*/
            case take_stock_op:                                  /* take_stock_op */
                size = get_locint_f();
                res  = get_locint_f();
                d    = get_val_v_f(size);
                for (i=0; i<size; i++)
                    dp_T0[res+i] = d[i];
                break;

                /*--------------------------------------------------------------------------*/
            case death_not:                                          /* death_not */
                arg1 = get_locint_f();
                arg2 = get_locint_f();
                break;

                /*--------------------------------------------------------------------------*/
            case gen_quad:                                            /* gen_quad */
                arg1  = get_locint_f();
                arg2  = get_locint_f();
                res   = get_locint_f();
                coval = get_val_f();
                if (coval!=dp_T0[arg1]){
                    fprintf(DIAG_OUT,"ADOL-C Warning: forward sweep aborted; tape invalid!\n");
                }
                coval = get_val_f();
                NOT_IMPLEMENTED_YET
                break;

                /****************************************************************************/
            case ext_diff:
                get_locint_f(); 
                get_locint_f(); 
                get_locint_f(); 
                get_locint_f(); 
                get_locint_f(); 
                get_locint_f();
                NOT_IMPLEMENTED_YET;
                break;

            case ignore_me:
                break;

                /*--------------------------------------------------------------------------*/
            default:                                                   /* default */
                /* Die here, we screwed up */
                fprintf(DIAG_OUT,"Edge_Hess: Fatal error in tape_doc for op %d\n",operation);
                break;
        } /* endswitch */
        /* Read the next operation */
        operation=get_op_f();
    }  /* endwhile */
    delete[] dp_T0;
    dp_T0 = NULL;
    end_sweep();
#ifdef EDGE_DEBUG
    printf("edge_value_len=%d\n",edge_value_len);
#endif  // EDGE_EDBUG
    delete[] edge_t_index;
    *edge_index_p=edge_index;
    *edge_value_p=edge_value;
    *edge_index_len_p=edge_index_len;
    *edge_value_len_p=edge_value_len;
    *max_index_p=max_index;
    *indmap_p=indmap;
    return ret_val;
}
