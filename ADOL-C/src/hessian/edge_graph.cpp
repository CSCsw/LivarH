#define isCreating (ri->px!=0.0)||(ri->py!=0.0)||(ri->pxy!=0.0)
#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <adolc/hessian/edge_main.h>

extern int inc_edge_count;
extern int del_edge_count;
extern int is_symmetric;
extern void (*increase_edge)(int ,int ,double , map<int, map<int,double> >*);

void compute_pushing(size_t tl, int *tp, double *tw, derivative_info* ri, map<int,map<int, double> > *graph){
  size_t i;
  map<int, double> *edges;
  int p;
  double w;

	    edges = &(*graph)[ri->r];
            tl=0;
	    for (std::map<int,double>::iterator it=edges->begin(); it!=edges->end(); it++){
              tp[tl]=it->first;tw[tl]=it->second;tl++;
	    }
//get rid of symmetric term
	    graph->erase(ri->r);
            if (is_symmetric==0){
              for(i=0;i<tl;i++){
                del_edge_count++; 
	        (*graph)[tp[i]].erase(ri->r);
              }
            }
	if (ri->x!=NULLLOC){
            for(i=0;i<tl;i++){
//a interesting bug
//	      cout<<ri->r;edge_check_edges(edges);
//	      p=it->first;w=it->second;
//end bug
              p=tp[i];w=tw[i];
//	      printf("p=%d,w=%10.10f\n",p,w);
		/*   insected code from original version */
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
		/*             copied code ends                  */

	    }//for
	}
	else{
	}//nothing to be pushed
}

void compute_creating(derivative_info* ri, map<int, double> *Adjoints, map<int,map<int, double> > *graph){
	if ((isCreating) && ((*Adjoints)[ri->r]!=0.0)){
/* old code */
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
/* copied code end */
	    
	}//creating
}

void compute_adjoints(derivative_info* ri, map<int, double> *Adjoints){
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


#ifdef NO_ASSIGN_BYPASS
void compute_global_pushing(size_t tl, int *tp, double *tw, int r, map<int, double> *first, map<int, map<int, double> > *gGraph){
  size_t i,j;
  map<int, double> *edges;
  int p;
  double w;
  edges = &(*gGraph)[r];
  tl=0;
  for (std::map<int,double>::iterator it=edges->begin(); it!=edges->end(); it++){
    tp[tl]=it->first;tw[tl]=it->second;tl++;
  }
//get rid of symmetric term
  gGraph->erase(r);
  if (is_symmetric==0){
    for(i=0;i<tl;i++){
      del_edge_count++; 
      (*gGraph)[tp[i]].erase(r);
    }
  }

  for(i=0;i<tl;i++){
    p=tp[i];w=tw[i];
    if (p!=r){
//loop all adjoints
      for(std::map<int, double>::iterator it=first->begin();it!=first->end();it++){
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
//set cardinal
      for(std::map<int, double>::iterator it1=first->begin();it1!=first->end();it1++){
        for(std::map<int, double>::iterator it2=it1;it2!=first->end();it2++){
          increase_edge(it1->first,it2->first,it1->second*it2->second*w,gGraph);
        }
      }
    }
  }
}

void compute_global_creating(int r, map<int, map<int, double> > *second, map<int, double> *Adjoints, map<int, map<int, double> > *gGraph){
  int x,y;
  double w;
  double a=(*Adjoints)[r];
//  printf("Adjoints[%d]=%10.5f\n",r,a);
  map<int, double> *edges;
  for(map<int, map<int, double> >::iterator it=second->begin();it!=second->end();it++){
    x=it->first;
    edges=&(it->second);
    for(map<int, double>::iterator it2=edges->begin();it2!=edges->end();it2++){
      y=it2->first;
      if (x>=y){
        w=it2->second;
        increase_edge(x,y,w*a,gGraph);      
      }
    }
  }
}

void compute_global_adjoints(int r, map<int, double> *first, map<int, double> *Adjoints){
  double w;
  w=(*Adjoints)[r];
//  printf("Adjoints[%d]=%10.5f\n",r,w);
//  (*Adjoints)[r]=0.0;
  Adjoints->erase(r);
  for(map<int, double>::iterator it=first->begin();it!=first->end();it++){
    (*Adjoints)[it->first]+=it->second*w;
//    printf("Adjoints[%d]=%10.5f\n",it->first,(*Adjoints)[it->first]);
  }
}
#endif
