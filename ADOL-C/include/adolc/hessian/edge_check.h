#ifndef __EDGE_CHECK_H__
#define __EDGE_CHECK_H__
#include <map>
#include <adolc/hessian/edge_main.h>
#include <adolc/internal/adolc_settings.h>
using namespace std;

class EdgeBTree;

/*  edge_check.cpp */
void edge_check_graph(map<locint, EdgeBTree* > *graph);

void edge_check_adjoints(map<locint, double > *Adjoints, locint max_active);

void edge_check_info(derivative_info* ri);

void edge_check_edges(EdgeBTree *edges);

void edge_check_index(locint* edge_index, unsigned int edge_index_len);
#endif
