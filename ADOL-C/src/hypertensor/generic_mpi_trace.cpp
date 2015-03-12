#include <vector>
#include <iostream>
#include "mpi.h"


static std::vector<SRinfo> srinfo;


void RMPI_trace_init() {
  srinfo.clear();
}

void put_mpi_trace(int sr_tag,
                   locint loc,
                   int count,
                   int peer,
                   int tag,
                   MPI_Comm comm) {
  SRinfo sr_info(sr_tag, loc, count, peer, tag, comm);
  srinfo.push_back(sr_info);
}
