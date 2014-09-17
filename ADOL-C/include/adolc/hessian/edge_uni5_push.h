#ifndef __EDGE_UNI5_PUSH_H__
#define __EDGE_UNI5_PUSH_H__
#include <map>
/*  edge_uni5_push.cpp */

class EdgeBTree;

void edge_pushing_a(short           tnum,
                    map<locint, EdgeBTree* > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);

void edge_pushing_s(short           tnum,
                    map<locint, EdgeBTree* > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);

#ifdef PREACC
void edge_pushing_pre_a(short           tnum,
                    map<locint, EdgeBTree > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);
void edge_pushing_pre_s(short           tnum,
                    map<locint, EdgeBTree > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
);
#endif
#endif

