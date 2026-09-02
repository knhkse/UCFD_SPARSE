#include <mpi.h>
#include "config.h"


static UCFDInt cmp_idx(const void *a, const void *b)
{
    const UCFDInt x = *(const UCFDInt *)a;
    const UCFDInt y = *(const UCFDInt *)b;
    return (x > y) - (x < y);
}

static UCFDInt bsearch_idx(const UCFDInt *a, UCFDInt m, UCFDInt g)
{
    UCFDInt lo = 0;
    UCFDInt hi = m;

    while (lo < hi)
    {
        const UCFDInt md = lo + (hi - lo) / 2;
        if (a[md] < g)
            lo = md + 1;
        else if (a[md] > g)
            hi = md;
        else
            return md;
    }
    return -1;
}

static size_t owner(const UCFDInt *range, size_t size, UCFDInt g)
{
    size_t lo = 0;
    size_t hi = size;

    while (hi - lo > 1)
    {
        const int md = lo + (hi - lo) / 2;
        if (g < range[md])
            hi = md;
        else
            lo = md;
    }
    return lo;
}

static void build_range(MPI_Comm comm, UCFDInt n_local, UCFDInt *range)
{
    UCFDInt size;

    MPI_Comm_size(comm, &size);
    range[0] = 0;
    MPI_Allgather(&n_local, 1, MPI_INT, range + 1, 1, MPI_INT, comm);
    for (int p = 1; p <= size; ++p)
        range[p] += range[p - 1];
}

static inline void halo_wait(UCFDSpMVContext *c)
{
    if (c->nreq)
        MPI_Waitall(c->nreq, c->reqs, MPI_STATUSES_IGNORE);
}
