#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <adolc/adalloc.h>
#include <adolc/interfaces.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_check.h>
#include <adolc/hessian/edge_tape.h>
#include <adolc/hessian/edge_uni5_push.h>


#define ADD_INC_EDGE_COUNT	if (edge_is_local==0) { local_inc_edge_count++;} else {global_inc_edge_count++;}

using namespace std;
int local_inc_edge_count;
int global_inc_edge_count;
int del_edge_count;
int edge_is_local=0;
int is_symmetric=0;

int edge_op_cnt;
int edge_tape_size;
double *edge_value;
int *edge_index;
int edge_value_len;
int edge_index_len;
int max_active;
int edge_translate_flag;

void (*increase_edge)(int ,int ,double , map<int, map<int,double> >*);

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
                               /* options[0]=0,    no pre-accumulation    */
                               /*           =1,    pre-accumulation       */
                               /* options[1]=0,    ADOL-C indexing        */
                               /*           =1,    Computational Graph    */

)
{
  if ((options[0]!=0) && (options[0]!=1)){
    fprintf(stderr,"edge_pushing requires options[0]=0 or 1\n");
    return 0;
  }
  if ((options[1]!=0) && (options[1]!=1)){
    fprintf(stderr,"edge_pushing requires options[1]=0 or 1\n");
    return 0;
  }
  vector<derivative_info*> *tape_info=new vector<derivative_info*>;
  map<int,map<int,double> > *graph=new map<int, map<int, double> >;
  unsigned int *indmap;
//Step 1: reverse sweep
  edge_translate_flag=options[1];
  edge_tape(tag,dep,indep,basepoint,tape_info,&indmap);
//Step 2: edge pushing
  local_inc_edge_count=0;
  global_inc_edge_count=0;
  edge_is_local=0;
  del_edge_count=0;
  if (options[0]==0){
    if (options[1]==0){
      is_symmetric=0;
      edge_pushing_a(tag,graph);
    }
    else{
      is_symmetric=1;
      edge_pushing_s(tag,graph);
    }
  }
  else{
#ifdef NO_ASSIGN_BYPASS
    if (options[1]==0){
      is_symmetric=0;
      edge_pushing_pre_a(tag,graph);
    }
    else{
      is_symmetric=1;
      edge_pushing_pre_s(tag,graph);
    }
#endif
#ifndef NO_ASSIGN_BYPASS
    fprintf(stderr, "Pre-accumulation in Hessian must be enabled without assign_a shortpath.\n");
#endif
  }
//Step 3: retrive results
//  int i;
//  for(i=0;i<20;i++){
//    printf("indmap[%d]=%d\n",i,indmap[i]);
//  }
//  edge_check_graph(graph);
  edge_retrive(graph,indmap,nnz,rind,cind,values);
  cout<<"Number of local updates: "<<local_inc_edge_count<<endl;
  cout<<"Number of global updates: "<<global_inc_edge_count<<endl;
  cout<<"Number of deletions: "<<del_edge_count<<endl;

//  cout<<*nnz<<endl;
  delete[] indmap;
  delete tape_info;
  delete graph;
  return 1;
}

void edge_retrive(map<int, map<int, double> > *graph, unsigned int* indmap, int *nnz, unsigned int **rind, unsigned int **cind, double **values){
  int n=0;
  map<int, double> *edge;
  for(map<int, map<int, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
  {
    edge=&(ii->second);
    for(map<int, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
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
  for(map<int, map<int, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
  {
    edge=&(ii->second);
    for(map<int, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
      if (indmap[ii->first]>=indmap[mi->first]){
        (*rind)[n]=indmap[ii->first];
        (*cind)[n]=indmap[mi->first];
        (*values)[n]=mi->second;
        n++;
      }
    }
  }
}

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


