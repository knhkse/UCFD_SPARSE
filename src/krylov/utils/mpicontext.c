#include "mpicontext.h"


static int mpi_is_active(void)
{
    int initialized = 0;
    int finalized = 0;

    MPI_Initialized(&initialized);
    if (initialized)
        MPI_Finalized(&finalized);
    return initialized && !finalized;
}

ucfd_mpi_t UCFDMPIContextCreate(MPI_Fint fcomm, Ctx *ctx)
{
    if (!mpi_is_active()) return UCFD_MPI_NOT_ACTIVE;

    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    if (comm == MPI_COMM_NULL) return UCFD_MPI_INVALID_COMM;

    Ctx c = (Ctx)calloc(1, sizeof(*c));
    UCFDCheckNull(c, "MPI Context allocation failed\n");

    /* Duplicate MPI context */
    if (MPI_Comm_dup(comm, &c->comm) != MPI_SUCCESS) {
        free(c); return UCFD_MPI_ERROR;
    }

    /* Query the */
    /* Query the largest legal tag rather than assuming 32767.        */
    void *attr_val = NULL;
    int   flag = 0;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr_val, &flag);
    c->tag_bound = flag ? *(int *)attr_val : 32767;
    c->next_tag  = 0;

    *ctx = c;

    UCFDFunctionReturn(UCFD_MPI_SUCCESS);
}

ucfd_mpi_t UCFDMPIContextDestroy(Ctx *ctx)
{
    if (!ctx || !*ctx) UCFDFunctionReturn(UCFD_MPI_SUCCESS);

    /* Guard against calling this after MPI_Finalize */
    if (!mpi_is_active()) UCFDFunctionReturn(UCFD_MPI_NOT_ACTIVE);
    if ((*ctx)->comm != MPI_COMM_NULL)
        MPI_Comm_free(&(*ctx)->comm);
    free(*ctx);
    UCFDFunctionReturn(UCFD_MPI_SUCCESS);
}


inline void ContextNextTag(Ctx ctx, UCFDInt *tag_out)
{
    if (ctx->next_tag > ctx->tag_bound)
        UCFDRaiseError(UCFD_MPI_ERROR, "Tag space exhausted\n");
    *tag_out = ctx->next_tag++;
}
