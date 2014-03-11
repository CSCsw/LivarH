##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am 471 2014-02-19 18:11:34Z kulshres $
##
## Copyright (C)  Andrea Walther, Andreas Kowarz
##
## This file is part of ADOL-C. This software is provided as open source.
## Any use, reproduction, or distribution of the software constitutes 
## recipient's acceptance of the terms of the accompanying license file.
## 
##############################################################################

AM_CFLAGS                = @ac_adolc_cflags@ -DADOLC_INTERNAL=1
AM_CXXFLAGS              = @ac_adolc_cxxflags@ -DADOLC_INTERNAL=1
AM_CPPFLAGS               = -I$(top_srcdir)/ADOL-C/include -I$(top_srcdir)/ADOL-C/src

MAINTAINERCLEANFILES     = Makefile.in *~ *.orig

noinst_LTLIBRARIES       = libhessian.la

libhessian_la_SOURCES     = edge_main.cpp edge_graph.cpp edge_tape.cpp edge_check.cpp edge_pushing_a.cpp edge_pushing_s.cpp edge_pushing_pre_a.cpp edge_pushing_pre_s.cpp

