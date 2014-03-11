#ifndef __EDGE_CHECK_H__
#define __EDGE_CHECK_H__
#include <iostream>
#include <vector>
#include <map>
/*  edge_check.cpp */
void edge_check_pre(vector<global_trace*> *gTrace);
void edge_check_tape(vector<derivative_info*> *tape_info);
void edge_check_graph(map<int,map<int, double> > *graph);
void edge_check_adjoints(map<int, double> *Adjoints, int max_active);
void edge_check_info(derivative_info* ri);
void edge_check_edges(map<int, double> *edges);
void edge_check_index(int* edge_index, int edge_index_len);
#endif
