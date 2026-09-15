/**
 * @file        mpicontext.h
 * @brief       MPI context for splitted CSR matrix
 * 
 */
#pragma once

#include "ucfdmpi.h"


struct _Ctx {
    MPI_Comm    comm;
    int         rank, size;
    UCFDInt     next_tag;
    UCFDInt     tag_bound;
};

static inline void ContextNextTag(Ctx ctx, UCFDInt *tag_out)
{
    if (ctx->next_tag > ctx->tag_bound)
        UCFDRaiseError("Tag space exhausted\n");
    *tag_out = ctx->next_tag++;
}
