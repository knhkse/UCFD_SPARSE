/**
 * @file        mpicontext.h
 * @brief       MPI context for splitted CSR matrix
 * 
 */
#pragma once

#include "ucfdtypes.h"

typedef struct {
    MPI_Comm        comm;
    UCFDInt         nrecv;
    UCFDInt         *recv_nei, *recv_count, *recv_off;
    
    UCFDInt         nsend;
    UCFDInt         *send_nei, *send_count, *send_off;
    UCFDInt         send_total, *send_idx;

    UCFDReal        *sbuf;
    UCFDReal        *lvec;
    UCFDInt         nghost;
    MPI_Request     *reqs;
    UCFDInt         nreq;
} UCFDMPIContext;

