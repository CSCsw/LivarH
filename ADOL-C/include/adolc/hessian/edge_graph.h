#ifndef __EDGE_GRAPH_H__
#define __EDGE_GRAPH_H__
#include <iostream>
#include <vector>
#include <map>
#include <adolc/adolc.h>

class EdgeLocalGraph;

/*  edge_graph.cpp */
void compute_pushing(unsigned int tl,
                     locint *tp,
                     double *tw,
                     derivative_info* ri,
                     std::map<locint,std::map<locint, double> > *graph);

void compute_creating(derivative_info* ri,
                      std::map<locint, double> *Adjoints,
                      std::map<locint, std::map<locint, double> > *graph);

void compute_adjoints(derivative_info* ri,
                      std::map<locint, double> *Adjoints);

#ifdef PREACC
void compute_local(locint *tp,
                   double *tw,
                   derivative_info* ri,
                   EdgeLocalGraph *local_graph);

void compute_global(locint *tp,
                    double *tw,
                    EdgeLocalGraph *local_graph,
                    locint r,
                    std::map<locint, double> *Adjoints,
                    std::map<locint, std::map<locint, double> > *graph);
#endif  // PREACC

#endif  // __EDGE_GRAPH_H__
