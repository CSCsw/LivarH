#include "taping_p.h"
#include "oplate.h"
#include <adolc/adolc.h>


// For the current implementation, we only send 1 variable!
void* ADOLC_rawData_p(void* buf, int count) {
  void* ret = 0;
  if (cout > 0){
    adouble* adouble_p = (adouble*)buf;
    ret = (void*)(&(ADOLC_GLOBAL_TAPE_VARS.store[adouble_p->loc()]));
  }
  return ret;
}
int RMPI_Send(void* buf,
              int count,
              MPI_Datatype datatype,
              int dest,
              int tag,
              MPI_Comm comm) {
  if (count != 1) {
    std::cout << "Only 1 variable can be sent in RMPI_Send()" << std::endl;
  }
// 1 Generate a dummy dependent variable
  double dummy_dep;
  adouble* adouble_p = (adouble*)buf;
  adouble[0] >> dummy_dep;
// 2 Do the send
  void* rawData_p = ADOLC_rawData_p(buf, count);
  int rc = MPI_Send(rawData_p, count, MPI_DOUBLE, dest, tag, comm);
// 3 Record on MPI trace
  put_op(ampi_send);
  put_mpi_trace(RMPI_SEND_TAG, adouble_p->loc(), count, dest, tag, comm);
  return rc;
}

int RMPI_Recv(void* buf,
              int count,
              MPI_Datatype datatype,
              MPI_Comm comm,
              MPI_Status* status) {
  if (count != 1) {
    std::cout << "Only 1 variable can be recv in RMPI_Recv()" << std::endl;
  }
// 1 Do the recv
  void* rawData_p = ADOLC_rawData_p(buf, count);
  int rc = MPI_Recv(rawData_p, count, MPI_DOUBLE, src, tag, comm, status);
// 2 Generate a dummy independent variable
  adouble* adouble_p = (adouble*)buf;
  adouble dummy_ind;
  dummy_ind <<= ((double*)rawData_p)[0];
  adouble_p[0] = dummy_ind;
// 3 Record on MPI trace
  put_op(ampi_recv);
  put_mpi_trace(RMPI_RECV_TAG, adouble_p->loc(), count, src, tag, comm);
  return rc;
}
