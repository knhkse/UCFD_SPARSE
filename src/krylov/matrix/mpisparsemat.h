#pragma once

#include <mpi.h>
#include "sparsemat.h"
#include "mpicontext.h"


typedef struct {
    MPI_Comm        comm;
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




