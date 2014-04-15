#define isCreating (ri->px!=0.0)||(ri->py!=0.0)||(ri->pxy!=0.0)

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <adolc/hessian/edge_main.h>
#include <adolc/adolc.h>

extern int __edge_is_symmetric__;

void increase_edge(locint i, locint j, double w, map<locint, map<locint, double> > *graph){
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
void compute_global_pushing(unsigned int tl, locint *tp, double *tw, locint r, map<locint, double> *first, map<locint, map<locint, double> > *gGraph){
    unsigned int i,j;
    map<locint, double> *edges;
    locint p;
    double w;
    edges = &(*gGraph)[r];
    tl=0;
    for (std::map<locint,double>::iterator it=edges->begin(); it!=edges->end(); it++){
        tp[tl]=it->first;tw[tl]=it->second;tl++;
    }
//get rid of symmetric term
    gGraph->erase(r);
    if (__edge_is_symmetric__==0){
        for(i=0;i<tl;i++){
            (*gGraph)[tp[i]].erase(r);
        }
    }

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

void compute_global_creating(locint r, map<locint, map<locint, double> > *second, map<locint, double> *Adjoints, map<locint, map<locint, double> > *gGraph){
    locint x,y;
    double w;
    double a=(*Adjoints)[r];
    map<locint, double> *edges;
    for(map<locint, map<locint, double> >::iterator it=second->begin();it!=second->end();it++){
        x=it->first;
        edges=&(it->second);
        for(map<locint, double>::iterator it2=edges->begin();it2!=edges->end();it2++){
            y=it2->first;
            if (x>=y){
                w=it2->second;
                increase_edge(x,y,w*a,gGraph);      
            }
        }
    }
}

void compute_global_adjoints(locint r, map<locint, double> *first, map<locint, double> *Adjoints){
    double w;
    w=(*Adjoints)[r];
//  printf("Adjoints[%d]=%10.5f\n",r,w);
//  (*Adjoints)[r]=0.0;
    Adjoints->erase(r);
    for(map<locint, double>::iterator it=first->begin();it!=first->end();it++){
        (*Adjoints)[it->first]+=it->second*w;
//    printf("Adjoints[%d]=%10.5f\n",it->first,(*Adjoints)[it->first]);
    }
}
#endif
