#ifndef __EDGE_CHECK_H__
#define __EDGE_CHECK_H__
#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <adolc/hessian/edge_main.h>
/*  edge_check.cpp */
void edge_check_graph(map<locint, map<locint, double > > *graph);
void edge_check_adjoints(map<locint, double > *Adjoints, locint max_active);
void edge_check_info(derivative_info* ri);
void edge_check_edges(map<locint, double > *edges);
void edge_check_index(locint* edge_index, unsigned int edge_index_len);
#endif
