#ifndef __EDGE_GRAPH_H__
#define __EDGE_GRAPH_H__
#include <iostream>
#include <vector>
#include <map>
#include <adolc/adolc.h>

class EdgeLocalGraph;

/*  edge_graph.cpp */
void compute_pushing(unsigned int tl, locint *tp, double *tw, derivative_info* ri, map<locint,map<locint, double> > *graph);
void compute_creating(derivative_info* ri, map<locint, double> *Adjoints, map<locint,map<locint, double> > *graph);
void compute_adjoints(derivative_info* ri, map<locint, double> *Adjoints);

#ifdef PREACC
void compute_local(locint *tp,
                   double *tw,
                   derivative_info* ri,
                   EdgeLocalGraph *local_graph);

void compute_global(locint *tp,
                    double *tw,
                    EdgeLocalGraph *local_graph,
                    locint r,
                    map<locint, double> *Adjoints,
                    map<locint, map<locint, double> > *graph);
//void compute_global_pushing(unsigned tl, locint *tp, double *tw, locint r, map<locint, double> *first, map<locint, map<locint, double> > *gGraph);
//void compute_global_creating(locint r, map<locint, map<locint, double> > *second, map<locint, double> *Adjoints, map<locint, map<locint, double> > *gGraph);
//void compute_global_adjoints(locint r, map<locint, double> *first, map<locint, double> *Adjoints);
#endif
#endif
