#ifndef __EDGE_TAPE_H__
#define __EDGE_TAPE_H__
#include <vector>
/*  edge_tape.cpp */
int edge_tape(short tnum,         /* tape id */
              int depcheck,       /* consistency chk on # of dependents */
              int indcheck,       /* consistency chk on # of independents */
              const double *basepoint, /* independent variable values */
              vector<derivative_info*> *tape_info_p,  /* derivative info on the tape */
              unsigned int **indmap); /* Mapping from location to independent index */

void translate_tape(vector<derivative_info*> *tape_info);
 
#endif
