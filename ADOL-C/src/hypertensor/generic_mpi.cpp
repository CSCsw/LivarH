#include "taping_p.h"
#include "oplate.h"
#include <adolc/adolc.h>

#ifdef ENABLE_GENERIC_MPI

#include <mpi.h>
#include <adolc/hypertensor/generic_mpi_trace.h>
MPI_Datatype RMPI_ADOUBLE;

std::vector<double> dummy_ind_vec;

#define MAX_DUMMY_SIZE 100000

static adouble* dummy_reverse;
static int dummy_reverse_size = 0;

static adouble* get_dummy(int count) {
//  std::cout << "In get dummy, count = " << count << std::endl;
  if (count + dummy_reverse_size >= MAX_DUMMY_SIZE) {
    std::cout << "Not enough dummy adoubles reversed! "
              << "Current size = " << MAX_DUMMY_SIZE << std::endl;
    exit(-1);
  }
  adouble* ret = &(dummy_reverse[dummy_reverse_size]);
  dummy_reverse_size += count;
  return ret;
}

void RMPI_Init(int* argc, char** argv[]) {
  MPI_Init(argc, argv);
  MPI_Type_contiguous(1, MPI_DOUBLE, &RMPI_ADOUBLE);
  MPI_Type_commit(&RMPI_ADOUBLE);
  RMPI_trace_init();  
  dummy_ind_vec.clear();
  dummy_reverse = new adouble[MAX_DUMMY_SIZE];
  dummy_reverse_size = 0;
}

void RMPI_Finalize() {
  if (RMPI_ADOUBLE != MPI_DATATYPE_NULL) {MPI_Type_free(&RMPI_ADOUBLE);}
  MPI_Finalize();
  delete[] dummy_reverse;
}

// For the current implementation, we only send 1 variable!
static void* ADOLC_rawData_p(void* buf, int count) {
  void* ret = 0;
  if (count > 0){
    adouble* adouble_p = (adouble*)buf;
//    std::cout << "adouble_p->loc = " << adouble_p->loc() << std::endl;
    ret = (void*)(&(ADOLC_GLOBAL_TAPE_VARS.store[adouble_p->loc()]));
  }
  return ret;
}

// For Send_ind we do not need the dummy_dep counter part
// Because, we assume they'll always stay in the fixed location in the memory.
int RMPI_Send_ind(void* buf,
                  int count,
                  MPI_Datatype datatype,
                  int dest,
                  int tag,
                  MPI_Comm comm) {
  int rc = 0;
//  std::cout << "In RMPI_Send:" << std::endl;
  if (datatype == RMPI_ADOUBLE) {
    void* send_buf;
//    std::cout << "Active type" << std::endl;
// 1 Generate a dummy dependent variable
    double dummy_dep_v;
    adouble* dummy_dep = (adouble*)buf;
// 2 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_send);
      put_mpi_trace(RMPI_SEND_IND_TAG, dummy_dep[0].loc(), count, dest, tag, comm);
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


int RMPI_Send(void* buf,
              int count,
              MPI_Datatype datatype,
              int dest,
              int tag,
              MPI_Comm comm) {
  int rc = 0;
//  std::cout << "In RMPI_Send:" << std::endl;
  if (datatype == RMPI_ADOUBLE) {
    void* send_buf;
//    std::cout << "Active type" << std::endl;
// 1 Generate a dummy dependent variable
    double dummy_dep_v;
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_dep = get_dummy(count);
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

int RMPI_Recv_ind(void* buf,
                  int count,
                  MPI_Datatype datatype,
                  int src,
                  int tag,
                  MPI_Comm comm,
                  MPI_Status* status) {
  int rc = 0;
//  std::cout << "In Recv:" << std::endl;
  if (datatype == RMPI_ADOUBLE) {
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_ind = get_dummy(count);
    void* recv_buf;
    recv_buf = ADOLC_rawData_p(buf, count);
// 1 Do the recv
    rc = MPI_Recv(recv_buf, count, MPI_DOUBLE, src, tag, comm, status);
// 2 Generate a dummy independent variable
    for(int i = 0; i < count; i++) {
//      std::cout << "recved: " << ((double*)recv_buf)[i] << std::endl;
      dummy_ind[i] <<= ((double*)recv_buf)[i];
    }
// 3 Record on MPI trace
    if (ADOLC_CURRENT_TAPE_INFOS.traceFlag) {
      put_op(ampi_recv);
      put_mpi_trace(RMPI_RECV_IND_TAG, dummy_ind[0].loc(), count, src, tag, comm);
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


int RMPI_Recv(void* buf,
              int count,
              MPI_Datatype datatype,
              int src,
              int tag,
              MPI_Comm comm,
              MPI_Status* status) {
  int rc = 0;
//  std::cout << "In Recv:" << std::endl;
  if (datatype == RMPI_ADOUBLE) {
    adouble* adouble_p = (adouble*)buf;
    adouble* dummy_ind = get_dummy(count);
    void* recv_buf;
    recv_buf = ADOLC_rawData_p(buf, count);
// 1 Do the recv
    rc = MPI_Recv(recv_buf, count, MPI_DOUBLE, src, tag, comm, status);
// 2 Generate a dummy independent variable
    for(int i = 0; i < count; i++) {
//      std::cout << "recved: " << ((double*)recv_buf)[i] << std::endl;
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

//for simplicity, we assume root = 0;
int RMPI_Reduce(void* sendbuf,
                void* recvbuf,
                int count,
                MPI_Datatype datatype,
                MPI_Op op,
                int root,
                MPI_Comm comm) {
  int rc = 0;
  if (datatype == RMPI_ADOUBLE) {
    adouble* t_buf = new adouble[count];
    adouble* s_buf = (adouble*)sendbuf;
    int size;
    int myid;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &myid);
//    std::cout << myid << "in reduction" << std::endl;
    int d = 0;
    int p = 1;
    int mask;
    int peer;
    while(p < size) {
      mask = p - 1;
//      std::cout << myid << " mask : " << mask << std::endl;
      if ((myid & mask) == 0) {
        peer = myid ^ p;
//        std::cout << myid << " peer : " << peer << std::endl;
        if (peer < size) {
          if ((myid & p) == 0) {
            // recv from peer
//            std::cout << myid << "recv from : " << peer << std::endl;
            RMPI_Recv((void*)t_buf, count, RMPI_ADOUBLE, peer, 0, comm, MPI_STATUS_IGNORE);
            for(int i = 0; i < count; i++) {
              if (op == MPI_SUM) {
                s_buf[i] += t_buf[i];
              } else if(op == MPI_PROD) {
                s_buf[i] *= t_buf[i];
              } else {
                std::cout<<"Unsupported operator in RMPI_Reduce!" << std::endl;
              }
            }
          } else {
            // send to peer
//            std::cout << myid << " send to : " << peer << std::endl;
            RMPI_Send((void*)sendbuf, count, RMPI_ADOUBLE, peer, 0, comm);
          }
        }
      }
      d = d + 1;
      p = p * 2;
    }
    if (myid == 0) {
      adouble* r_buf = (adouble*)recvbuf;
      for(int i=0; i < count; i++) {
        r_buf[i] = s_buf[i];
      }
    }
//    delete t_buf;  
  } else {
    rc = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
  }
  return rc;
}
#endif // ENABLE_GENERIC_MPI
