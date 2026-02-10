/**
 * @file        ucfd_types.h
 * @brief       Header file for UCFD_SPRASE Library types
 */
#ifndef _UCFD_TYPES_H
#define _UCFD_TYPES_H
#include "config.h"

typedef enum
{
    UCFD_STATUS_RHO_BREAKDOWN = -1,
    UCFD_STATUS_SUCCESS = 0,
    UCFD_STATUS_FAILED = 1,
    UCFD_MKL_FAILED = 2,
    UCFD_STATUS_CONVERGED = 3,
    UCFD_STATUS_DIVERGED = 4,
    UCFD_MAX_ITER = 5,
    UCFD_STATUS_NOT_SUPPORTED = 6
} ucfd_status_t;

typedef enum
{
    NONE = 0,
    BILU = 1,
    LUSGS = 2
} ucfd_precon_type_t;

typedef void (*ucfd_precon_solve)(int, int*, int*, int*, double*, double*);



// TODO
typedef struct ucfd_bsrmat {
    UCFD_INT    *row_ptr;
    UCFD_INT    *col_ind;
    UCFD_INT    *diag_ind;
    UCFD_FLOAT  *nnz_data;
} ucfd_bsrmat_t;

#endif  // _UCFD_TYPES_H