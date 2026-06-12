#ifndef SPARSE_H
#define SPARSE_H

#include <ucfd_types.h>


typedef struct _spmat {
    UCFD_INT bn;
    UCFD_INT block;
    const UCFD_INT *row_ptr;
    const UCFD_INT *col_ind;
    const UCFD_INT *diag_ind;
    UCFD_FLOAT *nnz_data;
} ucfd_spmat;


typedef enum {
    UCFD_SPARSE_SUCCESS   = 0,
    UCFD_SPARSE_ERR_ARG   = -1,
    UCFD_SPARSE_ERR_ALLOC = -2,
    UCFD_SPARSE_ERR_SIZE = -3
} ucfd_sparse_status_t;

typedef void (*ucfd_precon_solve)(ucfd_spmat*, UCFD_FLOAT*);

#endif /* SPARSE_H */
