#define isCreating (ri->px!=0.0)||(ri->py!=0.0)||(ri->pxy!=0.0)

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_graph.h>
#include <adolc/hessian/edge_b_tree.h>
#include <adolc/adolc.h>

extern int __edge_is_symmetric__;

void increase_edge(locint i, locint j, double w, map<locint, EdgeBTree*> *graph){
//  printf("insert <%d,%d, %.5f>\n", i,j,w);
//  fflush(stdout);
    if (__edge_is_symmetric__==0){
        if (i!=j){
            (*graph)[i]->update(j, w);
            (*graph)[j]->update(i, w);
        }
        else{
            (*graph)[i]->update(j, w);
        }
    }
    else{
        if (i>=j){
            (*graph)[i]->update(j, w);
        }
        else{
            (*graph)[j]->update(i, w);
        }
    }
}



void compute_pushing(locint *tp,
                     double *tw,
                     derivative_info* ri,
                     map<locint, EdgeBTree* > *graph){
    size_t i, size;
    locint p;
    double w;

//    printf("in pushing...\n");
//    fflush(stdout);
    EdgeBTree *edge;
    if ( graph->find(ri->r) == graph->end() ){
       graph->insert(std::pair<locint, EdgeBTree*>(ri->x, new EdgeBTree()));
    }
    edge = graph->find(ri->r)->second;
    
//    printf("edge = %x\n", edge);
//    fflush(stdout);
    edge->ToArray(&size, tp, tw);
//    printf("1");
//    fflush(stdout);
//remove the row from the graph and symmetric term if needed
    if (__edge_is_symmetric__==0){
        for(i=0;i<size;i++){
          if (tw[i] != 0.0) {
	    (*graph)[tp[i]]->erase(ri->r);
          }
        }
    }
//    printf("2");
//    fflush(stdout);
    graph->erase(ri->r);
    delete edge;

//    printf("retrieve done!\n");
//    fflush(stdout);
    if ( ri->x != NULLLOC && graph->find(ri->x) == graph->end() ) {
      graph->insert(std::pair<locint, EdgeBTree*>(ri->x, new EdgeBTree()));
    }
    if ( ri->y != NULLLOC && graph->find(ri->y) == graph->end() ) {
      graph->insert(std::pair<locint, EdgeBTree*>(ri->y, new EdgeBTree()));
    }

    if (ri->x!=NULLLOC){
        for(i=0;i<size;i++){
            p=tp[i];w=tw[i];
            if (w == 0.0) {continue;}
            if (ri->y!=-1){
                if (p!=ri->r){
                    if (ri->x==p){
                        increase_edge(p,p,2*ri->dx*w,graph);
                    }
                    else{
                        increase_edge(ri->x,p,ri->dx*w,graph);
                    }
                    if(ri->y==p){
                        increase_edge(p,p,2*ri->dy*w,graph);
                    }
                    else{
                        increase_edge(ri->y,p,ri->dy*w,graph);
                    }
                }
                else{
                    if (ri->x!=ri->y){
                        increase_edge(ri->x,ri->x,ri->dx*ri->dx*w,graph);
                        increase_edge(ri->x,ri->y,ri->dx*ri->dy*w,graph);
                        increase_edge(ri->y,ri->y,ri->dy*ri->dy*w,graph);
                    }
                    else{
                        increase_edge(ri->x,ri->x,(ri->dx*ri->dx+2*ri->dx*ri->dy+ri->dy*ri->dy)*w,graph);
                    }
                }
            }
            else{
                if (p!=ri->r){
                    if (ri->x==p){
                        increase_edge(p,p,2*ri->dx*w,graph);
                    }
                    else{
                        increase_edge(ri->x,p,ri->dx*w,graph);
                    }
                }
                else{
                    increase_edge(ri->x,ri->x,ri->dx*ri->dx*w,graph);
                }
            }
        }//for
    }
    else{
    }//nothing to be pushed
}

void compute_creating(derivative_info* ri,
                      map<locint, double > *Adjoints,
                      map<locint, EdgeBTree* > *graph){
    if ((isCreating) && ((*Adjoints)[ri->r]!=0.0)){
        if (ri->px!=0.0){
            increase_edge(ri->x,ri->x,(*Adjoints)[ri->r]*ri->px,graph);
        }
        if (ri->py!=0.0){
            increase_edge(ri->y,ri->y,(*Adjoints)[ri->r]*ri->py,graph);
        }
        if (ri->pxy!=0.0){
            if (ri->x!=ri->y){
                increase_edge(ri->x,ri->y,(*Adjoints)[ri->r]*ri->pxy,graph);
            }
            else{
            increase_edge(ri->x,ri->y,(*Adjoints)[ri->r]*ri->pxy*2.0,graph);
            }
        }	    
    }//creating
}

void compute_adjoints(derivative_info* ri,
                      map<locint, double> *Adjoints){
    double w;
    w=(*Adjoints)[ri->r];
    Adjoints->erase(ri->r);
    if (ri->dx!=0.0){
        (*Adjoints)[ri->x]+=w*ri->dx;
    }
    if (ri->dy!=0.0){
        (*Adjoints)[ri->y]+=w*ri->dy;
    }
}


#ifdef PREACC
void compute_global_pushing(locint *tp,
                            double *tw,
                            locint r,
                            map<locint, double > *first,
                            map<locint, EdgeBTree > *gGraph){
    size_t i;
    size_t tl;
    locint p;
    double w;
//remove the row from graph and symmetric part if needed
    (*gGraph)[r].ToArray(&tl, &tp, &tw);
    if (__edge_is_symmetric__==0){
        for(i=0;i<tl;i++){
            (*gGraph)[tp[i]].erase(r);
        }
    }
    gGraph->erase(r);

    for(i=0;i<tl;i++){
        p=tp[i];w=tw[i];
        if (p!=r){
//loop all adjoints
            for(std::map<locint, double>::iterator it=first->begin();it!=first->end();it++){
                if (it->first!=p){
//px*w
                    increase_edge(it->first,p,it->second*w,gGraph);
                }
                else{
//2*px*w
                    increase_edge(it->first,p,2.0*it->second*w,gGraph);
                }
            }
        }
        else{
//for all unordered set
            for(std::map<locint, double>::iterator it1=first->begin();it1!=first->end();it1++){
                for(std::map<locint, double>::iterator it2=it1;it2!=first->end();it2++){
                    increase_edge(it1->first,it2->first,it1->second*it2->second*w,gGraph);
                }
            }
        }
    }
}

void compute_global_creating(locint *tp,
                             double *tw,
                             locint r,
                             map<locint, EdgeBTree > *second,
                             map<locint, double> *Adjoints,
                             map<locint, Edge > *gGraph){
    locint x,y;
    size_t tl;
    double w;
    double a=(*Adjoints)[r];
    EdgeBTree *edges;
    for(map<locint, EdgeBTree >::iterator it=second->begin();it!=second->end();it++){
        x=it->first;
        edges=&(it->second);
        edges->ToArray(&tl, tp, tw);
        for(size_t i = 0; i < tl; ++i) {
          if (x > tl[i]) {
            increase_edge(x, y tw[i] * a, gGraph);
          }
        }
    }
}

void compute_global_adjoints(locint r,
                             map<locint, double> *first,
                             map<locint, double> *Adjoints){
    double w;
    w=(*Adjoints)[r];
    Adjoints->erase(r);
    for(map<locint, double>::iterator it=first->begin();it!=first->end();it++){
        (*Adjoints)[it->first]+=it->second*w;
    }
}
#endif
