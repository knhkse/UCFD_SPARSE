#pragma once

#include "sparsemat.h"
#include "mpicontext.h"


typedef struct {
    MPI_Comm        comm;
    UCFDInt         blocksize;                          // Only used in BSR
    UCFDInt         nrecv, tag;
    UCFDInt         *recv_nei, *recv_count, *recv_off;
    
    UCFDInt         nsend;
    UCFDInt         *send_nei, *send_count, *send_off;
    UCFDInt         send_total, *send_idx;

    UCFDReal        *sbuf;
    UCFDReal        *lvec;
    UCFDInt         nghost;
    MPI_Request     *reqs;
    UCFDInt         nreq;
} UCFDSpMVContext;

typedef struct {
    BaseCSR          A;
    BaseCSR          B;
    UCFDInt          n_local, n_ghost, n_boundary;
    UCFDInt          *garray;
    UCFDInt          *boundary_rows;
    UCFDSpMVContext  spmvctx;
} MPICSR;


typedef struct {
    BaseBSR          A;
    BaseBSR          B;
    UCFDInt          n_local, n_ghost, n_boundary;
    UCFDInt          *garray;
    UCFDInt          *boundary_rows;
    UCFDSpMVContext  spmvctx;
} MPIBSR;


#if defined(__cplusplus)
extern "C" {
#endif

/* Internal functions */
UCFD_INTERN ucfd_status_t UCFDSpMVContextDestroy(UCFDSpMVContext *c);

#if defined(__cplusplus)
}
#endif
