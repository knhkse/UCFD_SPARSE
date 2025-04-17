#ifndef PRECON_H
#define PRECON_H
#include "types.h"

/**
 * @file        precon.h
 * @brief       Header file for preconditioners for Krylov subspace methods
 * @details     Declaration of each prepare/psolve function.
 *              `prepare` must be executed before Krylov method routine starts.
 *              `psolve` function is executed in Krylov method to solve Px = b.
 */

/**
 * @brief       Block fill-in Incomplete LU preconditioner for BSR matrix format.
 * @param       neles       Number of element cells
 * @param       nvars       Block size of BSR matrix
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 */
ucfd_status_t bilu_prepare(const int neles, const int nvars,
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *nnz_data);

/**
 * @brief       Solver function for BILU preconditioner.
 * @param       neles       Number of element cells
 * @param       nvars       Block size of BSR matrix
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 * @param       b           Right-hand-side
 */
ucfd_status_t bilu_psolve(const int neles, const int nvars, const int *row_ptr,
                          const int *col_ind, const int *diag_ind, double *nnz_data, double *b);

/**
 * @brief       LU-SGS preconditioner for BSR matrix format.
 * @param       neles       Number of element cells
 * @param       nvars       Block size of BSR matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 */
ucfd_status_t lusgs_prepare(const int neles, const int nvars, const int *diag_ind, double *nnz_data);

/**
 * @brief       Solver function for LU-SGS preconditioner.
 * @param       neles       Number of element cells
 * @param       nvars       Block size of BSR matrix
 * @param       row_ptr     Row-directional index pointer array of the matrix
 * @param       col_ind     Column index array of the matrix
 * @param       diag_ind    Diagonal matrix index array based on `row_ptr`
 * @param       nnz_data    Non-zero value array of the matrix
 * @param       b           Right-hand-side
 */
ucfd_status_t lusgs_psolve(const int neles, const int nvars, const int *row_ptr,
                           const int *col_ind, const int *diag_ind, double *nnz_data, double *b);

#endif // PRECON.H