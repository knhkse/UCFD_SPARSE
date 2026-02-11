#ifndef PRECON_H
#define PRECON_H
#include "ucfd_types.h"
#include "config.h"


/**
 * @file        precon.h
 * @brief       Header file for preconditioners for Krylov subspace methods
 * @details     Declaration of each prepare/psolve function.
 *              `prepare` must be executed before Krylov method routine starts.
 *              `psolve` function is executed in Krylov method to solve Px = b.
 */

/**
 * @brief       Block fill-in Incomplete LU preconditioner for BSR matrix format.
 * @param       bn          Number of element cells
 * @param       blk         Block size of BSR matrix
 * @param       iw          Working array
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 */
ucfd_status_t bilu_prepare(int bn, int *iw,
                           int *row_ptr, int *col_ind, int *diag_ind, double *nnz_data);

/**
 * @brief       Solver function for BILU preconditioner.
 * @param       bn          Number of element cells
 * @param       blk         Block size of BSR matrix
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 * @param       b           Right-hand-side
 */
void bilu_psolve(int bn, int *row_ptr,
                 int *col_ind, int *diag_ind, double *nnz_data, double *b);

/**
 * @brief       LU-SGS preconditioner for BSR matrix format.
 * @param       bn          Number of element cells
 * @param       blk         Block size of BSR matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 */
ucfd_status_t lusgs_prepare(int bn, int *diag_ind, double *nnz_data);

/**
 * @brief       Solver function for LU-SGS preconditioner.
 * @param       bn          Number of element cells
 * @param       blk         Block size of BSR matrix
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 * @param       b           Right-hand-side
 */
void lusgs_psolve(int bn, int *row_ptr,
                  int *col_ind, int *diag_ind, double *nnz_data, double *b);


/**
 * @brief       Unpreconditioned solver
 */
void none_psolve(int bn, int *row_ptr,
                 int *col_ind, int *diag_ind, double *nnz_data, double *b);

#endif // PRECON.H