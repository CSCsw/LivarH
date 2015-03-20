#include "taping_p.h"
#include "oplate.h"
#include <adolc/adolc.h>
#include <mpi.h>
#include <adolc/hypertensor/generic_mpi_trace.h>

MPI_Datatype RMPI_ADOUBLE;

std::vector<double> dummy_ind_vec;

void RMPI_Init(int* argc, char** argv[]) {
  MPI_Init(argc, argv);
  MPI_Type_contiguous(1, MPI_DOUBLE, &RMPI_ADOUBLE);
  MPI_Type_commit(&RMPI_ADOUBLE);
  RMPI_trace_init();  
  dummy_ind_vec.clear();
}

void RMPI_Finalize() {
  if (RMPI_ADOUBLE != MPI_DATATYPE_NULL) {MPI_Type_free(&RMPI_ADOUBLE);}
  MPI_Finalize();
}

// For the current implementation, we only send 1 variable!
static void* ADOLC_rawData_p(void* buf, int count) {
  void* ret = 0;
  if (count > 0){
    adouble* adouble_p = (adouble*)buf;
    std::cout << "adouble_p->loc = " << adouble_p->loc() << std::endl;
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
  std::cout << "In RMPI_Send:" << std::endl;
  if (count != 1) {
    std::cout << "Only 1 variable can be sent in RMPI_Send()" << std::endl;
  }
  void* send_buf;
  MPI_Datatype send_type;
  if (datatype == RMPI_ADOUBLE) {
    std::cout << "Active type" << std::endl;
// 1 Generate a dummy dependent variable
    double dummy_dep_v;
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_dep = new adouble[count];
    dummy_dep[0] = adouble_p[0];
    dummy_dep[0] >>= dummy_dep_v;
    std::cout << "dummy_dep = " << dummy_dep << std::endl;
// 2 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_send);
      put_mpi_trace(RMPI_SEND_TAG, dummy_dep[0].loc(), count, dest, tag, comm);
    }
// 3 Redirect the send buffer
    send_buf = ADOLC_rawData_p(buf, count);
    send_type = MPI_DOUBLE;
    std::cout << adouble_p->loc() << " ---> [" << dest << "]" << std::endl;
  } else {
    send_buf = buf;
    send_type = datatype;
  }
  int rc = MPI_Send(send_buf, count, send_type, dest, tag, comm);
  return rc;
}

int RMPI_Recv(void* buf,
              int count,
              MPI_Datatype datatype,
              int src,
              int tag,
              MPI_Comm comm,
              MPI_Status* status) {
  if (count != 1) {
    std::cout << "Only 1 variable can be recv in RMPI_Recv()" << std::endl;
  }
  void* recv_buf;
  MPI_Datatype recv_type;
  if (datatype == RMPI_ADOUBLE) {
    recv_buf = ADOLC_rawData_p(buf, count);
    recv_type = MPI_DOUBLE;
  } else {
    recv_buf = buf;
    recv_type = datatype;
  }
// 1 Do the recv
  int rc = MPI_Recv(recv_buf, count, recv_type, src, tag, comm, status);
  if (datatype == RMPI_ADOUBLE) {
// 2 Generate a dummy independent variable
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_ind = new adouble[count];
    dummy_ind[0] <<= ((double*)recv_buf)[0];

// 3 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_recv);
      put_mpi_trace(RMPI_RECV_TAG, dummy_ind[0].loc(), count, src, tag, comm);
    }

    adouble_p[0] = dummy_ind[0];
    dummy_ind_vec.push_back(((double*)recv_buf)[0]);

    std::cout << adouble_p->loc() << " <--- [" << src << "]" << std::endl; 
  }
  return rc;
}
