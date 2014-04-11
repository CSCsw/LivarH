#ifndef __EDGE_TAPE_H__
#define __EDGE_TAPE_H__
#include <adolc/adolc_settings.h>
/*  edge_tape.cpp */
int edge_tape(  short tnum,         /* tape id */
                int depcheck,       /* consistency chk on # of dependents */
                int indcheck,       /* consistency chk on # of independents */
                const double*       basepoint, /* independent variable values */
                int                 edge_translate_flag,/* 1=Translate the index or not */
                unsigned int**      indmap_p,           /* Mapping from location to independent index */
                locint**            edge_index_p,       /* The translated index (increasing) */
                double**            edge_value_p,       /* The value of variables (for reverse) */                
                unsigned int*       edge_index_len_p,   /* The length of edge_index[] */
                unsigned int*       edge_value_len_p,   /* The length of edge_value[] */
                locint*             max_index_p                /* The max index (after translation */
);

//void translate_tape(vector<derivative_info*> *tape_info);
 
#endif
