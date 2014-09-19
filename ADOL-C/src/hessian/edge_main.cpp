#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <adolc/adolc.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_check.h>
#include <adolc/hessian/edge_tape.h>
#include <adolc/hessian/edge_uni5_push.h>


#define ADD_INC_EDGE_COUNT	if (edge_is_local==0) { local_inc_edge_count++;} else {global_inc_edge_count++;}

using namespace std;

int __edge_is_symmetric__=1;

//void (*increase_edge)(int ,int ,double , map<int, map<int,double> >*);

int edge_count_global = 0;
int edge_count_local = 0;
int edge_hess(
    short           tag,        /* tape identification                     */
    int             dep,        /* consistency chk on # of dependents      */
    int             indep,      /* number of independent variables         */
    const double*   basepoint,  /* independant variable values             */
    int*            nnz,        /* number of nonzeros                      */
    unsigned int**  rind,       /* row index                               */
    unsigned int**  cind,       /* column index                            */
    double**        values,     /* non-zero values                         */
    int*            options     /* control options                         */
                                /* options[0]=0,    disable statement preaccumulation    */
                                /*           =1,    enable  statement preaccumulation    */
                                /* options[1]=0,    ADOL-C indexing        */
                                /*           =1,    Monotonic indesing     */

)
{
    if ((options[0]!=0) && (options[0]!=1)){
        fprintf(stderr,"Edge_Hess requires options[0]=0 or 1 (Default is 0)\n");
        options[0]=0;
    }
    if ((options[1]!=0) && (options[1]!=1)){
        fprintf(stderr,"Edge_Hess requires options[1]=0 or 1 (Default is 1)\n");
        options[0]=1;
    }
    map<locint,map<locint,double> > *graph=new map<locint, map<locint, double> >;
    unsigned int *indmap;
    locint *edge_index;
    double *edge_value;
    unsigned int edge_index_len;
    unsigned int edge_value_len;
    unsigned int max_index;
//Step 1: translate index, forward
    edge_tape(tag,dep,indep,basepoint,options[1],&indmap,&edge_index,&edge_value,&edge_index_len,&edge_value_len,&max_index);

    edge_count_global = 0;
    edge_count_local = 0;

//Step 2: edge_pushing, reverse
    __edge_is_symmetric__=options[1];
    if (options[0]==0){
        if (options[1]==0){
            edge_pushing_a(tag,graph,edge_index,edge_value,edge_index_len,edge_value_len,max_index);
        }
        else{
            edge_pushing_s(tag,graph,edge_index,edge_value,edge_index_len,edge_value_len,max_index);
        }
    }
    else{
#ifdef PREACC
        if (options[1]==0){
            edge_pushing_pre_a(tag,graph,edge_index,edge_value,edge_index_len,edge_value_len,max_index);
        }
        else{
            edge_pushing_pre_s(tag,graph,edge_index,edge_value,edge_index_len,edge_value_len,max_index);
        }
#endif
#ifndef PREACC
        fprintf(stderr, "Preaccumulation in Hessian must be enabled WITH --enable-preacc when configure\n");
#endif
    }
    printf("global edge count = %d\n", edge_count_global);
    printf("local edge count = %d\n", edge_count_local);
//Step 3: retrive results
    edge_retrive(graph,indmap,nnz,rind,cind,values);
    delete[] indmap;
    delete[] edge_index;
    delete[] edge_value;
    delete graph;
    return 1;
}

void edge_retrive(map<locint, map<locint, double> > *graph, unsigned int* indmap, int *nnz, unsigned int **rind, unsigned int **cind, double **values){
    unsigned int n=0;
    map<locint, double> *edge;
    for(map<locint, map<locint, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
    {
        edge=&(ii->second);
        for(map<locint, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
            if (indmap[ii->first]>=indmap[mi->first]){
                n++;
            }
        }
    }
//  cout<<"nnz="<<n<<endl;
    *nnz=n;
    if (*rind==NULL){
        *rind=(unsigned int*)malloc(sizeof(unsigned int)*n);
    }
    if (*cind==NULL){
        *cind=(unsigned int*)malloc(sizeof(unsigned int)*n);
    }
    if (*values==NULL){
        *values=(double*)malloc(sizeof(double)*n);
    }
    n=0;
    for(map<locint, map<locint, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
    {
        edge=&(ii->second);
        for(map<locint, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
            if (indmap[ii->first]>=indmap[mi->first]){
                (*rind)[n]=indmap[ii->first];
                (*cind)[n]=indmap[mi->first];
                (*values)[n]=mi->second;
                n++;
            }
        }
    }
}


/*
void increase_edge_a(int i,int j,double w,map<int, map<int,double> > *graph){
//  printf("increase_edge <%d,%d><%10.10f>\n",i,j,w);
  if (i!=j){
    ADD_INC_EDGE_COUNT;
    (*graph)[i][j]+=w;
    ADD_INC_EDGE_COUNT;
    (*graph)[j][i]+=w;
  }
  else{
    ADD_INC_EDGE_COUNT;
    (*graph)[i][j]+=w;
  }
}
void increase_edge_s(int i,int j,double w,map<int, map<int,double> > *graph){
//  printf("increase_edge <%d,%d><%10.10f>\n",i,j,w);
  ADD_INC_EDGE_COUNT;
  if (i>=j){
    (*graph)[i][j]+=w;
  }
  else{
    (*graph)[j][i]+=w;
  }
}
*/

