#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <adolc/hessian/edge_check.h>
//#define PRINT_CHECK_EDGES
//#define PRINT_CHECK_GRAPH
#define PRINT_CHECK_INFO
//#define PRINT_CHECK_TAPE
//#define PRINT_CHECK_ADJOINTS

void edge_check_info(derivative_info* ri){
#ifdef PRINT_CHECK_INFO
     cout << "op="<<(int)ri->opcode<<": res="<<ri->r<<" x="<<ri->x<<" y="<<ri->y<<std::endl;
     cout<<"dx="<<ri->dx<<"  dy="<<ri->dy<<endl;
     cout<<"px="<<ri->px<<"  py="<<ri->py<<"  pxy="<<ri->pxy<<endl;
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
void edge_check_graph(map<locint, map<locint, double> > *graph){
#ifdef PRINT_CHECK_GRAPH
  map<locint, double> *edge;
  for(map<locint, map<locint, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
  {
    edge=&(ii->second);
    cout<<ii->first<<": "<<" size="<<edge->size()<<"  :";
    for(map<locint, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
      cout<<"["<<mi->first<<"]="<<mi->second<<"  ";
    }
    cout<<endl;
  }
#endif
}
void edge_check_edges(map<locint, double> *edge){
#ifdef PRINT_CHECK_EDGES
    cout<<"A: "<<" size="<<edge->size()<<"  :";
    for(map<locint, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
      cout<<"["<<mi->first<<"]="<<mi->second<<"  ";
    }
    cout<<endl;
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
