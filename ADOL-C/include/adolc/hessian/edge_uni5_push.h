#ifndef __EDGE_UNI5_PUSH_H__
#define __EDGE_UNI5_PUSH_H__
#include <iostream>
#include <vector>
#include <map>
#include <adolc/hessian/edge_main.h>
/*  edge_uni5_push.cpp */

void edge_pushing_a(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);

void edge_pushing_s(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);

#ifdef NO_ASSIGN_BYPASS
void edge_pushing_pre_a(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);
void edge_pushing_pre_s(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);
#endif
#endif

