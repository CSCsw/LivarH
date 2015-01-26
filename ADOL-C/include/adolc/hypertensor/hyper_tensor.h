#ifndef __HYPER_MAIN_H__
#define __HYPER_MAIN_H__

int hyper_tensor(short tag,
                 int dep,
                 int indep,
                 const double* basepoint,
                 int* nnz,
                 unsigned int** rind,
                 unsigned int** cind,
                 unsigned int** xind,
                 double** values,
                 int* optinos);




#endif // __HYPER_MAIN_H__
