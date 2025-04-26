#ifndef KRYLOV_H
#define KRYLOV_H
#include "mkl.h"
#include "types.h"

/**
 * @file        krylov.h
 * @brief       Header file for Krylov subspace methods
*/

/**
 * @brief       Serial GMRES routine
 * @param       op              `Intel MKL` sparse matrix handler (System matrix)
 * @param       precon_type     Preconditioner type
 * @param       neles           Number of element cells (= Number of block rows)
 * @param       nvars           Block size of BSR matrix (= Size of blocks)
 * @param       m               Restart number
 * @param       row_ptr         Row-directional index pointer array of the matrix
 * @param       col_ind         Column index array of the matrix
 * @param       diag_ind        Diagonal matrix index array based on `row_ptr`
 * @param       pre_nnz_data    Non-zero value array of the **preconditioner** matrix
 * @param       tol             Tolerence of convergence
 * @param       itmax           Maximum iteration number
 * @param       x               Solution array
 * @param       b               Right-hand-side array
 * @param       H               Upper Hessenberg matrix
 * @param       V               Orthogonal matrix
 * @param       g               Givens rotation prameters
 * @param       y               Least squares array for argmin(||beta - Hy||)
 * @param       w               Arnoldi iteration array
 * @param       r               Residual array
 */
ucfd_status_t serial_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m,
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data,
                           const double tol, const double itmax, double *x, double *b, double *H, double *V, double *g, double *y, double *w, double *r);


/**
* @brief       Single GMRES iteration routine
* @param       op              `Intel MKL` sparse matrix handler (System matrix)
* @param       precon_type     Preconditioner type
* @param       neles           Number of element cells (= Number of block rows)
* @param       nvars           Block size of BSR matrix (= Size of blocks)
* @param       m               Restart number
* @param       row_ptr         Row-directional index pointer array of the matrix
* @param       col_ind         Column index array of the matrix
* @param       diag_ind        Diagonal matrix index array based on `row_ptr`
* @param       pre_nnz_data    Non-zero value array of the **preconditioner** matrix
* @param       x               Solution array
* @param       b               Right-hand-side array
* @param       H               Upper Hessenberg matrix
* @param       V               Orthogonal matrix
* @param       g               Givens rotation prameters
* @param       y               Least squares array for argmin(||beta - Hy||)
* @param       w               Arnoldi iteration array
* @param       r               Residual array
*/
ucfd_status_t step_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m,
                         const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data,
                         double *x, double *b, double *H, double *V, double *g, double *y, double *w, double *r);


/**
 * @brief       Serial BiCGstab routine
 * @param       op              `Intel MKL` sparse matrix handler (System matrix)
 * @param       precon_type     Preconditioner type
 * @param       neles           Number of element cells (= Number of block rows)
 * @param       nvars           Block size of BSR matrix (= Size of blocks)
 * @param       row_ptr         Row-directional index pointer array of the matrix
 * @param       col_ind         Column index array of the matrix
 * @param       diag_ind        Diagonal matrix index array based on `row_ptr`
 * @param       pre_nnz_data    Non-zero value array of the **preconditioner** matrix
 * @param       tol             Tolerence of convergence
 * @param       itmax           Maximum iteration number
 * @param       x               Solution array
 * @param       b               Right-hand-side array
 * @param       r               Residual array
 * @param       p               
 */
ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, \
                              const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                              const double tol, const double itmax, double *x, double *b, double *r, double *p, double *v, double *s, double *t);


#endif // KRYLOV_H