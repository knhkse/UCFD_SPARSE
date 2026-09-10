/**
 * @file        mpicontext.h
 * @brief       MPI context for splitted CSR matrix
 * 
 */
#pragma once

#include "ucfdmpi.h"


struct _Ctx {
    MPI_Comm    comm;
    UCFDInt     rank;
    UCFDInt     size;
    UCFDInt     next_tag;
    UCFDInt     tag_bound;
};

UCFD_INTERN void ContextNextTag(Ctx ctx, UCFDInt *tag_out);

