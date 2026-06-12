/**
 * @file        krylov.h
 * @brief       Header file for Krylov subspace methods
*/
#ifndef KRYLOV_H
#define KRYLOV_H
#include <mkl.h>
#include "ucfd_types.h"
#include "sparse.h"
#include "config.h"

#ifndef EPS
    #define eps 2.22e-16 // Machine epsilon of double precision
#endif


/**
 * @brief       Serial GMRES routine
 * @param       op                  `Intel MKL` sparse matrix handler (System matrix)
 * @param       precon_type         Preconditioner type
 * @param       bn                  Number of element cells (= Number of block rows)
 * @param       block               Block size of BSR matrix (= Size of blocks)
 * @param       m                   Restart number
 * @param       iter                Maximum iteration number
 * @param       tol                 Tolerence of convergence
 * @param       row_ptr             Row-directional index pointer array of the matrix
 * @param       col_ind             Column index array of the matrix
 * @param       diag_ind            Diagonal matrix index array based on `row_ptr`
 * @param       precon_nnz_data     Non-zero value array of the **preconditioner** matrix
 * @param       x                   Solution array
 * @param       b                   Right-hand-side array
 * @param       H                   Upper Hessenberg matrix
 * @param       V                   Orthogonal matrix
 * @param       g                   Givens rotation prameters
 * @param       y                   Least squares array for argmin(||beta - Hy||)
 * @param       w                   Arnoldi iteration array
 * @param       r                   Residual array
 */
ucfd_status_t serial_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, UCFD_INT bn, UCFD_INT block, UCFD_INT m, UCFD_INT *iter, UCFD_FLOAT tol,
                           ucfd_spmat *precon, UCFD_FLOAT *x, UCFD_FLOAT *b, UCFD_FLOAT *H, UCFD_FLOAT *V, UCFD_FLOAT *g, UCFD_FLOAT *y, UCFD_FLOAT *w, UCFD_FLOAT *r);


/**
* @brief       Single GMRES iteration routine
* @param       op                   `Intel MKL` sparse matrix handler (System matrix)
* @param       precon_type          Preconditioner type
* @param       bn                   Number of element cells (= Number of block rows)
* @param       block                Block size of BSR matrix (= Size of blocks)
* @param       m                    Restart number
* @param       row_ptr              Row-directional index pointer array of the matrix
* @param       col_ind              Column index array of the matrix
* @param       diag_ind             Diagonal matrix index array based on `row_ptr`
* @param       precon_nnz_data      Non-zero value array of the **preconditioner** matrix
* @param       x                    Solution array
* @param       b                    Right-hand-side array
* @param       H                    Upper Hessenberg matrix
* @param       V                    Orthogonal matrix
* @param       g                    Givens rotation prameters
* @param       y                    Least squares array for argmin(||beta - Hy||)
* @param       w                    Arnoldi iteration array
* @param       r                    Residual array
*/
ucfd_status_t step_gmres(sparse_matrix_t op, ucfd_precon_solve psolve, const struct matrix_descr descr,
                         UCFD_INT bn, UCFD_INT m, UCFD_INT *flag, ucfd_spmat *precon,
                         UCFD_FLOAT *x, UCFD_FLOAT *b, UCFD_FLOAT *H, UCFD_FLOAT *V, UCFD_FLOAT *g, UCFD_FLOAT *y, UCFD_FLOAT *w, UCFD_FLOAT *r);


/**
 * @brief       Serial BiCGstab routine
 * @param       op                  `Intel MKL` sparse matrix handler (System matrix)
 * @param       precon_type         Preconditioner type
 * @param       bn                  Number of element cells (= Number of block rows)
 * @param       block               Block size of BSR matrix (= Size of blocks)
 * @param       iter                Maximum iteration number
 * @param       tol                 Tolerence of convergence
 * @param       row_ptr             Row-directional index pointer array of the matrix
 * @param       col_ind             Column index array of the matrix
 * @param       diag_ind            Diagonal matrix index array based on `row_ptr`
 * @param       precon_nnz_data     Non-zero value array of the **preconditioner** matrix
 * @param       x                   Solution array
 * @param       b                   Right-hand-side array
 * @param       r                   Residual array
 * @param       p                   
 */
ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, UCFD_INT bn, UCFD_INT *iter, UCFD_FLOAT tol,
                              ucfd_spmat *precon, UCFD_FLOAT *x, UCFD_FLOAT *b, UCFD_FLOAT *r, UCFD_FLOAT *p, UCFD_FLOAT *v, UCFD_FLOAT *s, UCFD_FLOAT *t);

#endif // KRYLOV_H