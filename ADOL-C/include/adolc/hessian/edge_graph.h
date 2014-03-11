#ifndef __EDGE_GRAPH_H__
#define __EDGE_GRAPH_H__
#include <iostream>
#include <vector>
#include <map>
/*  edge_graph.cpp */
void compute_pushing(size_t tl, int *tp, double *tw, derivative_info* ri, map<int,map<int, double> > *graph);
void compute_creating(derivative_info* ri, map<int, double> *Adjoints, map<int,map<int, double> > *graph);
void compute_adjoints(derivative_info* ri, map<int, double> *Adjoints);
#ifdef NO_ASSIGN_BYPASS
void compute_global_pushing(size_t tl, int *tp, double *tw, int r, map<int, double> *first, map<int, map<int, double> > *gGraph);
void compute_global_creating(int r, map<int, map<int, double> > *second, map<int, double> *Adjoints, map<int, map<int, double> > *gGraph);
void compute_global_adjoints(int r, map<int, double> *first, map<int, double> *Adjoints);
#endif
#endif
