#define isCreating (ri->px!=0.0)||(ri->py!=0.0)||(ri->pxy!=0.0)

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/hessian/edge_local_graph.h>
#include <adolc/internal/adolc_settings.h>
#include <adolc/adolc.h>

extern int __edge_is_symmetric__;

extern int edge_count_global;
extern int edge_count_local;

void increase_edge(locint i, locint j, double w, map<locint, map<locint, double> > *graph){
    edge_count_global++;
    if (__edge_is_symmetric__==0){
        if (i!=j){
            (*graph)[i][j]+=w;
            (*graph)[j][i]+=w;
        }
        else{
            (*graph)[i][j]+=w;
        }
    }
    else{
        if (i>=j){
            (*graph)[i][j]+=w;
        }
        else{
            (*graph)[j][i]+=w;
        }
    }
}



void compute_pushing(unsigned int tl, locint *tp, double *tw, derivative_info* ri, map<locint,map<locint, double> > *graph){
    unsigned int i;
    map<locint, double> *edges;
    locint p;
    double w;
    
    edges = &(*graph)[ri->r];
    tl=0;
    for (std::map<locint,double>::iterator it=edges->begin(); it!=edges->end(); it++){
        tp[tl]=it->first;tw[tl]=it->second;tl++;
    }
//get rid of symmetric term
	graph->erase(ri->r);
    if (__edge_is_symmetric__==0){
        for(i=0;i<tl;i++){
	        (*graph)[tp[i]].erase(ri->r);
        }
    }
	if (ri->x!=NULLLOC){
        for(i=0;i<tl;i++){
            p=tp[i];w=tw[i];
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

void compute_creating(derivative_info* ri, map<locint, double> *Adjoints, map<locint, map<locint, double> > *graph){
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

void compute_adjoints(derivative_info* ri, map<locint, double> *Adjoints){
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

void increase_local_edge(locint i, locint j, double w, EdgeLocalGraph *local_graph){
    edge_count_local++;
    size_t ind;
    if (i>=j){
      ind = local_graph->AddLiveVar(i);
      local_graph->insert(ind, j, w);
    }
    else{
      ind = local_graph->AddLiveVar(j);
      local_graph->insert(ind, i, w);
    }
}

void compute_local(locint *tp,
                   double *tw,
                   derivative_info* ri,
                   EdgeLocalGraph *local_graph) {

    size_t i;
    locint p;
    double w;
    double aw;
    size_t ind_r = local_graph->AddLiveVar(ri->r);
    size_t tl = 0;
    for(tl = 0; tl < local_graph->size_array[ind_r]; ++tl) {
      tp[tl] = local_graph->loc_array[ind_r][tl];
      tw[tl] = local_graph->hessian[ind_r][tl];
    }
    aw = local_graph->adjoints[ind_r];
    local_graph->erase(ind_r);

    

    size_t ind_x = local_graph->AddLiveVar(ri->x);
    size_t ind_y = local_graph->AddLiveVar(ri->y);

    // pushing
    if (ri->x!=NULLLOC){
        for(i=0;i<tl;i++){
            p=tp[i];w=tw[i];
            if (ri->y!=-1){
                if (p!=ri->r){
                    if (ri->x==p){
                        increase_local_edge(p,p,2*ri->dx*w,local_graph);
                    }
                    else{
                        increase_local_edge(ri->x,p,ri->dx*w,local_graph);
                    }
                    if(ri->y==p){
                        increase_local_edge(p,p,2*ri->dy*w,local_graph);
                    }
                    else{
                        increase_local_edge(ri->y,p,ri->dy*w,local_graph);
                    }
                }
                else{
                    if (ri->x!=ri->y){
                        increase_local_edge(ri->x,ri->x,ri->dx*ri->dx*w,local_graph);
                        increase_local_edge(ri->x,ri->y,ri->dx*ri->dy*w,local_graph);
                        increase_local_edge(ri->y,ri->y,ri->dy*ri->dy*w,local_graph);
                    }
                    else{
                        increase_local_edge(ri->x,ri->x,(ri->dx*ri->dx+2*ri->dx*ri->dy+ri->dy*ri->dy)*w,local_graph);
                    }
                }
            }
            else{
                if (p!=ri->r){
                    if (ri->x==p){
                        increase_local_edge(p,p,2*ri->dx*w,local_graph);
                    }
                    else{
                        increase_local_edge(ri->x,p,ri->dx*w,local_graph);
                    }
                }
                else{
                    increase_local_edge(ri->x,ri->x,ri->dx*ri->dx*w,local_graph);
                }
            }
        }//for
    }
    else{
    }//nothing to be pushed

    // creating
    if ( (isCreating) && ( aw != 0.0 ) ){
        if (ri->px!=0.0){
            increase_local_edge(ri->x,ri->x, aw*ri->px, local_graph);
        }
        if (ri->py!=0.0){
            increase_local_edge(ri->y,ri->y, aw*ri->py, local_graph);
        }
        if (ri->pxy!=0.0){
            if (ri->x!=ri->y){
                increase_local_edge(ri->x,ri->y, aw*ri->pxy, local_graph);
            }
            else{
                increase_local_edge(ri->x,ri->y, aw*ri->pxy*2.0, local_graph);
            }
        }	    
    } //creating

    // adjoints
    if (ri->dx!=0.0){
        local_graph->adjoints[ind_x]+=aw*ri->dx;
    }
    if (ri->dy!=0.0){
        local_graph->adjoints[ind_y]+=aw*ri->dy;
    } // adjoints
}


void compute_global(locint *tp,
                    double *tw,
                    EdgeLocalGraph *local_graph,
                    locint r,
                    map<locint, double> *Adjoints,
                    map<locint, map<locint, double> > *graph) {

    unsigned int i,j;
    map<locint, double> *edges;
    locint p;
    double w;
    edges = &(*graph)[r];
    size_t tl = 0;
    for (std::map<locint,double>::iterator it=edges->begin(); it!=edges->end(); it++){
        tp[tl]=it->first;tw[tl]=it->second;tl++;
    }
//get rid of symmetric term
    graph->erase(r);

    for(i=0;i<tl;i++){
      p=tp[i];w=tw[i];
      if (w != 0.0) {
        if (p!=r){
//loop all adjoints
            for(size_t ii = 0; ii< local_graph->size; ++ii) {
                if (local_graph->adjoints[ii] != 0.0) {
                  if (local_graph->loc[ii]!=p){
//px*w
                      increase_edge(local_graph->loc[ii], p, local_graph->adjoints[ii]*w, graph);
                  }
                  else{
//2*px*w
                      increase_edge(local_graph->loc[ii], p, 2.0*local_graph->adjoints[ii]*w,graph);
                  }
                }
            }
        }
        else{
//for all unordered set
            for(size_t ii = 0; ii < local_graph->size; ++ii) {
              if (local_graph->adjoints[ii] != 0.0) {
                for(size_t jj =0; jj < local_graph->size; ++jj) {
                  if (local_graph->adjoints[jj] != 0.0) {
                      increase_edge(local_graph->loc[ii], local_graph->loc[jj],
                                    local_graph->adjoints[ii]*local_graph->adjoints[jj], graph);
                  }
                }
              }
            }
        }
      }
    }

    //creating
    w = (*Adjoints)[r];
    Adjoints->erase(r);
    if (w != 0.0) {
      for(size_t ii = 0; ii < local_graph->size; ++ii) {
        for(size_t jj =0; jj < local_graph->size_array[ii]; ++jj) {
            increase_edge(local_graph->loc[ii], local_graph->loc_array[ii][jj],
                          w * local_graph->hessian[ii][jj], graph);
        }
      }
      for(size_t ii = 0; ii< local_graph->size; ++ii) {
        if(local_graph->adjoints[ii] != 0.0) {
          (*Adjoints)[local_graph->loc[ii]] += w * local_graph->adjoints[ii];
        }
      }
    }
}
#endif
