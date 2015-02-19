#ifndef HYPER_TENSOR_H_
#define HYPER_TENSOR_H_

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




#endif // HYPER_TENSOR_H_
