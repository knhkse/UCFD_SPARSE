

#include <stdio.h>
#include <math.h>
#include "mkl.h"
#include "bicgstab.h"
#include "precon.h"

#define eps 2.22e-16    // Machine epsilon


/**
 * @note        
 */
ucfd_status_t serial_bicgstab_mkl(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, \
                                  const double tol, const double itmax, double *dub, double *rhsb, double *r, double *p, double *v, double *s, double *t)
{
    int it, i, j, err;
    double rho, rhoprev, alpha, beta, omega, resid;
    double rv, ts, tt;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    // Computes residual : r := rhsb - A @ dub
    cblas_dcopy(n, rhsb, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");

    // Choose r\tilde as r
    /* r[:n] = r
       r[n:] = r\tilde */
    cblas_dcopy(n, r, 1, &r[n], 1);

    // Outer iteration
    for (it=0; it<itmax; it++) {
        resid = cblas_dnrm2(n, r, 1);
        printf("%d : %.10f\n", it, resid);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        rho = cblas_ddot(n, r, 1, &r[n], 1);
        
        // Rho breakdown
        if (fabs(rho) < eps) return -5;

        if (it == 0) cblas_dcopy(n, r, 1, p, 1);
        else {
            beta = (rho/rhoprev) * (alpha/omega);

            // p = r + beta*(p - omega*v)
            cblas_daxpy(n, -omega, v, 1, p, 1);
            cblas_dscal(n, beta, p, 1);
            cblas_daxpy(n, 1.0, r, 1, p, 1);
        }

        //phat = inv(M) @ p
        cblas_dcopy(n, p, 1, &p[n], 1);

        // Apply preconditioner to p[n]~
        if (precon_type == ILU0) {
            ilu0_sweep_bsr_mkl(precon, &p[n]);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
        }

        // v = A @ phat
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &p[0], 0.0, v);
        if (status != SPARSE_STATUS_SUCCESS) return -1;

        rv = cblas_ddot(n, v, 1, &r[n], 1);
        alpha = rho/rv;

        // s = r - alpha*v
        // 1) r := r - alpha*v (r = s)
        cblas_daxpy(n, -alpha, v, 1, r, 1);
        // 2) s := r
        cblas_dcopy(n, r, 1, s, 1);

        // x := x + alpha*phat
        cblas_daxpy(n, alpha, &p[n], 1, dub, 1);
        resid = cblas_dnrm2(n, s, 1);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        // shat = inv(M) @ s
        cblas_dcopy(n, s, 1, &s[n], 1);

        if (precon_type == ILU0) {
            ilu0_sweep_bsr_mkl(precon, &s[n]);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
        }

        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &s[0], 0.0, t);
        if (status != SPARSE_STATUS_SUCCESS) return -1;

        ts = cblas_ddot(n, t, 1, s, 1);
        tt = cblas_ddot(n, t, 1, t, 1);
        omega = ts/tt;

        // Update solution
        cblas_daxpy(n, omega, &s[0], 1, dub, 1);

        // r = s - omega*t
        cblas_daxpy(n, -omega, t, 1, r, 1);

        // Update rho
        rhoprev = rho;
    }

    if (it == itmax){
        printf("Not converged in %d iteration, residual=%.5f", it, beta);
        return it;
    }
    return UCFD_STATUS_ERROR;
}


ucfd_status_t serial_bicgstab(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, \
                              const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                              const double tol, const double itmax, double *dub, double *rhsb, double *r, double *p, double *v, double *s, double *t)
{
    int it, i, j, err;
    const int n = neles * nvars;
    double rho, rhoprev, alpha, beta, omega, resid;
    double rv, ts, tt;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    // Computes residual : r := rhsb - A @ dub
    cblas_dcopy(n, rhsb, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");

    // Choose r\tilde as r
    /* r[:n] = r
       r[n:] = r\tilde */
    cblas_dcopy(n, r, 1, &r[n], 1);

    // Outer iteration
    for (it=0; it<itmax; it++) {
        resid = cblas_dnrm2(n, r, 1);
        printf("%d : %.10f\n", it, resid);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        rho = cblas_ddot(n, r, 1, &r[n], 1);
        
        // Rho breakdown
        if (fabs(rho) < eps) return -5;

        if (it == 0) cblas_dcopy(n, r, 1, p, 1);
        else {
            beta = (rho/rhoprev) * (alpha/omega);

            // p = r + beta*(p - omega*v)
            cblas_daxpy(n, -omega, v, 1, p, 1);
            cblas_dscal(n, beta, p, 1);
            cblas_daxpy(n, 1.0, r, 1, p, 1);
        }

        //phat = inv(M) @ p
        cblas_dcopy(n, p, 1, &p[n], 1);

        // Apply preconditioner to p[n]~
        if (precon_type == ILU0) {
            ilu0_sweep_bsr(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &p[n]);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
        }

        // v = A @ phat
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &p[0], 0.0, v);
        if (status != SPARSE_STATUS_SUCCESS) return -1;

        rv = cblas_ddot(n, v, 1, &r[n], 1);
        alpha = rho/rv;

        // s = r - alpha*v
        // 1) r := r - alpha*v (r = s)
        cblas_daxpy(n, -alpha, v, 1, r, 1);
        // 2) s := r
        cblas_dcopy(n, r, 1, s, 1);

        // x := x + alpha*phat
        cblas_daxpy(n, alpha, &p[n], 1, dub, 1);
        resid = cblas_dnrm2(n, s, 1);
        if (resid < tol) return UCFD_STATUS_CONVERGED;

        // shat = inv(M) @ s
        cblas_dcopy(n, s, 1, &s[n], 1);

        if (precon_type == ILU0) {
            ilu0_sweep_bsr(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, &s[n]);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
        }

        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &s[0], 0.0, t);
        if (status != SPARSE_STATUS_SUCCESS) return -1;

        ts = cblas_ddot(n, t, 1, s, 1);
        tt = cblas_ddot(n, t, 1, t, 1);
        omega = ts/tt;

        // Update solution
        cblas_daxpy(n, omega, &s[0], 1, dub, 1);

        // r = s - omega*t
        cblas_daxpy(n, -omega, t, 1, r, 1);

        // Update rho
        rhoprev = rho;
    }

    if (it == itmax){
        printf("Not converged in %d iteration, residual=%.5f", it, beta);
        return it;
    }
    return UCFD_STATUS_ERROR;
}

