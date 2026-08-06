/**
 * @file        mpicontext.h
 * @brief       MPI context for splitted CSR matrix
 * 
 */
#pragma once

#include <mpi.h>
#include "ucfdmpi.h"


struct _Ctx {
    MPI_Comm    comm;
    UCFDInt     next_tag;
    UCFDInt     tag_bound;
};

UCFD_INTERN void ContextNextTag(Ctx ctx, UCFDInt *tag_out);

