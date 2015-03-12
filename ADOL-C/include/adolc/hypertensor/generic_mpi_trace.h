#ifndef GENERIC_MPI_TRACE_H_
#define GENERIC_MPI_TRACE_H_

#include <adolc/adolc.h>
#include "mpi.h"
#define RMPI_SEND_TAG 0
#define RMPI_RECV_TAG 1

class SRinfo {
 public:
  SRinfo(int sr_tag, locint loc, int count, int peer, int tag, MPI_Comm comm) {
    this->SR_tag = sr_tag;
    this->loc = loc;
    this->count = count;
    this->peer = peer;
    this->tag = tag;
    this->comm = comm;
  };
  int SR_tag;
  locint loc;
  int count;
  int peer;
  int tag;
  MPI_Comm comm;
};


#endif // GENERIC_MPI_TRACE_H_
