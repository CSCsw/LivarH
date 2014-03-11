#ifndef __EDGE_UNI5_PUSH_H__
#define __EDGE_UNI5_PUSH_H__
#include <iostream>
#include <vector>
#include <map>
#include <adolc/hessian/edge_main.h>
/*  edge_uni5_push.cpp */
void edge_pushing_a(short tnum,map<int,map<int,double> > *graph);
void edge_pushing_s(short tnum,map<int,map<int,double> > *graph);
#ifdef NO_ASSIGN_BYPASS
void edge_pushing_pre_a(short tnum,map<int,map<int,double> > *graph);
void edge_pushing_pre_s(short tnum,map<int,map<int,double> > *graph);
#endif
#endif

