#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <adolc/hessian/edge_check.h>
#include <adolc/hessian/edge_b_tree.h>
//#define PRINT_CHECK_EDGES
//#define PRINT_CHECK_GRAPH
//#define PRINT_CHECK_INFO
//#define PRINT_CHECK_TAPE
//#define PRINT_CHECK_ADJOINTS

void edge_check_info(derivative_info* ri){
#ifdef PRINT_CHECK_INFO
     cout << "op="<<(int)ri->opcode<<": res="<<ri->r<<" x="<<ri->x<<" y="<<ri->y<<std::endl;
     cout<<"dx="<<ri->dx<<"  dy="<<ri->dy<<endl;
     cout<<"px="<<ri->px<<"  py="<<ri->py<<"  pxy="<<ri->pxy<<endl;
     fflush(stdout);
#endif
}
void edge_check_index(locint *edge_index, unsigned int edge_index_len){
#ifdef PRINT_CHECK_INDEX
  unsigned int i;
  for(i=0;i<edge_index_len;i++){
    cout<<"ind["<<i<<"]="<<edge_index[i]<<std::endl;
  }
#endif
}
void edge_check_graph(map<locint, EdgeBTree* > *graph){
#ifdef PRINT_CHECK_GRAPH
  EdgeBTree *edge;
  for(map<locint, EdgeBTree*>::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
  {
    edge=ii->second;
    printf("Row[%d] ",ii->first);
    fflush(stdout);
    edge_check_edges(edge);
  }
#endif
}
void edge_check_edges(EdgeBTree *edge){
#ifdef PRINT_CHECK_EDGES
  size_t size;
  locint *tp;
  double *tw;
  edge->ToArrayAlloc(&size, &tp, &tw);
  printf(" (%d): ", size);
  for(size_t i=0; i<size; i++) {
    printf("[%d]=%.5f, ", tp[i], tw[i]);
  }
  printf("\n");
  fflush(stdout);
  free(tp);
  free(tw);
#endif
}
void edge_check_adjoints(map<locint, double> *Adjoints, locint max_index){
#ifdef PRINT_CHECK_ADJOINTS
  unsigned int i;
  cout<<"--------------------------ADJOINT-----------------------------------------"<<endl;
  for(i=0;i<max_index;i++){
    cout<<"["<<i<<"]="<<(*Adjoints)[i]<<"  ";
  }
  cout<<endl;
#endif
}
