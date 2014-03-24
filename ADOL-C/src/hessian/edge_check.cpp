#include <iostream>
#include <vector>
#include <cmath>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_check.h>

//#define PRINT_CHECK_PRE
//#define PRINT_CHECK_EDGES
//#define PRINT_CHECK_GRAPH
#define PRINT_CHECK_INFO
//#define PRINT_CHECK_TAPE
//#define PRINT_CHECK_ADJOINTS
//#define PRINT_CHECK_INDEX

void edge_check_pre(vector<global_trace*> *gTrace){
#ifdef PRINT_CHECK_PRE
  vector<global_trace*>::reverse_iterator ri;
  global_trace* rt;
  cout<<"pre-accumulation size="<<gTrace->size()<<endl;
  for(ri=gTrace->rbegin();ri!=gTrace->rend();ri++){
    rt=(*ri);
    cout<<"target="<<rt->r<<endl;
    edge_check_edges(rt->first);
    edge_check_graph(rt->second);
  }
#endif
}

void edge_check_tape(vector<derivative_info*> *tape_info){
#ifdef PRINT_CHECK_TAPE
   cout<<"Checking the informations on the tape..."<<std::endl;
   cout << std::endl << "Reverse Iterator:" <<tape_info->size()<< std::endl;	 
   vector<derivative_info*>::reverse_iterator ri;
   derivative_info* info;
   for(ri=tape_info->rbegin(); ri!=tape_info->rend(); ri++)
   {
     info=(*ri);
     cout << "op="<<(int)info->opcode<<": res="<<info->r<<" x="<<info->x<<" y="<<info->y<<std::endl;
     cout<<"dx="<<info->dx<<"  dy="<<info->dy<<endl;
     cout<<"px="<<info->px<<"  py="<<info->py<<"  pxy="<<info->pxy<<endl;
   }
#endif
}
void edge_check_info(derivative_info* ri){
#ifdef PRINT_CHECK_INFO
     cout << "op="<<(int)ri->opcode<<": res="<<ri->r<<" x="<<ri->x<<" y="<<ri->y<<std::endl;
     cout<<"dx="<<ri->dx<<"  dy="<<ri->dy<<endl;
     cout<<"px="<<ri->px<<"  py="<<ri->py<<"  pxy="<<ri->pxy<<endl;
#endif
}
void edge_check_index(int *edge_index, int edge_index_len){
#ifdef PRINT_CHECK_INDEX
  int i;
  for(i=0;i<edge_index_len;i++){
    cout<<"ind["<<i<<"]="<<edge_index[i]<<std::endl;
  }
#endif
}
void edge_check_graph(map<int, map<int, double> > *graph){
#ifdef PRINT_CHECK_GRAPH
  map<int, double> *edge;
  for(map<int, map<int, double> >::iterator ii=graph->begin(); ii!=graph->end(); ++ii)
  {
    edge=&(ii->second);
    cout<<ii->first<<": "<<" size="<<edge->size()<<"  :";
    for(map<int, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
      cout<<"["<<mi->first<<"]="<<mi->second<<"  ";
    }
    cout<<endl;
  }
#endif
}
void edge_check_edges(map<int, double> *edge){
#ifdef PRINT_CHECK_EDGES
    cout<<"A: "<<" size="<<edge->size()<<"  :";
    for(map<int, double>::iterator mi=edge->begin(); mi!=edge->end(); ++mi){
      cout<<"["<<mi->first<<"]="<<mi->second<<"  ";
    }
    cout<<endl;
#endif
}
void edge_check_adjoints(map<int, double> *Adjoints, int max_active){
#ifdef PRINT_CHECK_ADJOINTS
  int i;
  cout<<"--------------------------ADJOINT-----------------------------------------"<<endl;
  for(i=0;i<max_active;i++){
    cout<<"["<<i<<"]="<<(*Adjoints)[i]<<"  ";
  }
  cout<<endl;
#endif
}
