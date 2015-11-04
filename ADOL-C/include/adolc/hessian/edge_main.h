#ifndef __EDGE_MAIN_H__
#define __EDGE_MAIN_H__
#include <map>
#include <limits.h>
#include <adolc/adolc.h>

#define NOT_IMPLEMENTED_YET fprintf(stderr,"edge_hess: Not implemented yet\n");fflush(stderr);
#define NULLLOC UINT_MAX
#define MAX_TEMP_ARRAY_SIZE 1000

#define EDGE_DEBUG

class derivative_info{
  public:
    unsigned char opcode;
    locint r,x,y;
    double dx,dy;
    double px,py,pxy;

    derivative_info(){
      opcode=-1;
      r=NULLLOC;x=NULLLOC;y=NULLLOC;
      dx=0.0;dy=0.0;
      px=0.0;py=0.0;pxy=0.0;
    };
    ~derivative_info(){};
};

/*  edge_main.cpp */
int edge_hess(
    short          tag,        /* tape identification                     */
    int            dep,        /* consistency chk on # of dependents      */
    int            indep,      /* number of independent variables         */
    const double  *basepoint,  /* independant variable values             */
    int           *nnz,        /* number of nonzeros                      */
    unsigned int **rind,       /* row index                               */
    unsigned int **cind,       /* column index                            */
    double       **values,     /* non-zero values                         */
    int           *options     /* control options                         */
);

void edge_retrive(std::map<locint, std::map<locint, double> > *graph, unsigned int *indmap, int *nnz, unsigned int **rind, unsigned int **cind, double **values);

#endif
