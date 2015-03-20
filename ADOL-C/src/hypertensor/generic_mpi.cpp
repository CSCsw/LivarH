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
  int rc;
  std::cout << "In RMPI_Send:" << std::endl;
  if (datatype == RMPI_ADOUBLE) {
    void* send_buf;
    std::cout << "Active type" << std::endl;
// 1 Generate a dummy dependent variable
    double dummy_dep_v;
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_dep = new adouble[count];
    for(int i = 0; i < count; i++) {
      dummy_dep[i] = adouble_p[i];
      dummy_dep[i] >>= dummy_dep_v;
    }
// 2 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_send);
      put_mpi_trace(RMPI_SEND_TAG, dummy_dep[0].loc(), count, dest, tag, comm);
    }
// 3 Redirect the send buffer
    send_buf = ADOLC_rawData_p(buf, count);
    rc = MPI_Send(send_buf, count, MPI_DOUBLE, dest, tag, comm);
//    std::cout << adouble_p->loc() << " ---> [" << dest << "]" << std::endl;
  } else {
    rc = MPI_Send(buf, count, datatype, dest, tag, comm);
  }

  return rc;
}

int RMPI_Recv(void* buf,
              int count,
              MPI_Datatype datatype,
              int src,
              int tag,
              MPI_Comm comm,
              MPI_Status* status) {
  int rc;


  if (datatype == RMPI_ADOUBLE) {
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_ind = new adouble[count];
    void* recv_buf;
    recv_buf = ADOLC_rawData_p(buf, count);
// 1 Do the recv
    rc = MPI_Recv(recv_buf, count, MPI_DOUBLE, src, tag, comm, status);
// 2 Generate a dummy independent variable
    for(int i = 0; i < count; i++) {
      std::cout << "recved: " << ((double*)recv_buf)[i] << std::endl;
      dummy_ind[i] <<= ((double*)recv_buf)[i];
    }
// 3 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_recv);
      put_mpi_trace(RMPI_RECV_TAG, dummy_ind[0].loc(), count, src, tag, comm);
    }
    for(int i = 0; i < count; i++) {
      adouble_p[i] = dummy_ind[i];
      dummy_ind_vec.push_back(((double*)recv_buf)[i]);
    }

//    std::cout << adouble_p->loc() << " <--- [" << src << "]" << std::endl; 
  } else {
    rc = MPI_Recv(buf, count, datatype, src, tag, comm, status);
  }
  return rc;
}
