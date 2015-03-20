#include <vector>
#include <iostream>
#include <mpi.h>
#include <adolc/adolc.h>
#include <adolc/hypertensor/generic_mpi_trace.h>

// This is a toy implementation of the SendRecv Stack
std::vector<SRinfo> sr_stack;

void RMPI_trace_init() {
  sr_stack.clear();
}

void put_mpi_trace(int sr_tag,
                   locint loc,
                   int count,
                   int peer,
                   int tag,
                   MPI_Comm comm) {
  SRinfo sr_info(sr_tag, loc, count, peer, tag, comm);
  sr_stack.push_back(sr_info);
}

/*
void get_mpi_trace(SRInfo& sr_info) {
  sr_info = sr_stack.back();
  sr_stack.pop_back();
}
*/
