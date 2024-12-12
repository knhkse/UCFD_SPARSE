/** ======================================================================================================================
 * @file        gmres.c
 * @brief       Generalized Minimal RESidual method
 * @details     GMRES method is introduced by Saad and Schultz, which is the one of the Krylov
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        Dec 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "gmres.h"
#include "precon.h"


ucfd_status_t serial_gmres_mkl(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, const int m, \
                               const double tol, const double itmax, double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r)
{   
    int it, i, j, err;
    double beta, tmp, c, s, h1, h2, rr;
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

    for (it=0; it<itmax; it++) {
        beta = cblas_dnrm2(n, r, 1);
        printf("%d iteration, %.10f residual\n", it, beta);
        if (beta < tol) return UCFD_STATUS_CONVERGED;

        // Apply Preconditioner
        if (precon_type == ILU0) {
            ilu0_sweep_bsr_mkl(precon, r);
            beta = cblas_dnrm2(n, r, 1);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
            beta = cblas_dnrm2(n, r, 1);
        }

        y[0] = beta;

        // V = r/beta
        cblas_dcopy(n, r, 1, V, 1);
        cblas_dscal(n, 1/beta, V, 1);

        for (j=0; j<m; j++) {
            status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j*n], 0.0, w);
            if (status != SPARSE_STATUS_SUCCESS) return -1;
            
            if (precon_type == ILU0) {
                ilu0_sweep_bsr_mkl(precon, w);
            }
            else if (precon_type == LUSGS) {
                // TODO : LU-SGS preconditioner function
            }
            
            // Arnoldi iteration
            for (i=0; i<(j+1); i++) {
                tmp = cblas_ddot(n, w, 1, &V[i*n], 1);
                H[j+m*i] = tmp;
                cblas_daxpy(n, -tmp, &V[i*n], 1, w, 1);
            }

            tmp = cblas_dnrm2(n, w, 1);
            H[j+(j+1)*m] = tmp;
            cblas_dcopy(n, w, 1, &V[(j+1)*n], 1);
            cblas_dscal(n, 1/tmp, &V[(j+1)*n], 1);

            // Givens Rotation
            for (i=0; i<j; i++) {
                c = g[i*2];         // g[i, 0]
                s = g[i*2+1];       // g[i, 1]
                h1 = H[j+i*m];
                h2 = H[j+(i+1)*m];
                H[j+i*m] = c*h1 - s*h2;
                H[j+(i+1)*m] = s*h1 + c*h2;
            }

            h1 = H[j*(m+1)];
            h2 = H[j+(j+1)*m];
            rr = sqrt(h1*h1 + h2*h2);
            c = h1/rr;
            s = -h2/rr;
            H[j*(m+1)] = c*h1 - s*h2;
            g[j*2] = c;
            g[j*2+1] = s;

            // Modify e1 vector
            y[j+1] = y[j];
            y[j] *= c;
            y[j+1] *= s;
        }

        // Back substitution
        cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

        // Update
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, dub, 1);

        // Computes next iteration residual
        cblas_dcopy(n, rhsb, 1, r, 1);
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
        if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");
    }

    if (it == itmax){
        printf("Not converged in %d iteration, residual=%.5f", it, beta);
        return it;
    }

    return UCFD_STATUS_ERROR;
}


ucfd_status_t serial_gmres(sparse_matrix_t op, ucfd_precon_type_t precon_type, const int neles, const int nvars, const int m, \
                           const int *row_ptr, const int *col_ind, const int *diag_ind, double *pre_nnz_data, \
                           const double tol, const double itmax, double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r)
{   
    int it, i, j, err;
    const int n = neles * nvars;
    double beta, tmp, c, s, h1, h2, rr;
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

    for (it=0; it<itmax; it++) {
        beta = cblas_dnrm2(n, r, 1);
        // printf("%d iteration, %.10f residual\n", it, beta);
        if (beta < tol) return UCFD_STATUS_CONVERGED;

        if (precon_type == ILU0) {
            ilu0_sweep_bsr_raw(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);
            beta = cblas_dnrm2(n, r, 1);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
            beta = cblas_dnrm2(n, r, 1);
        }

        y[0] = beta;

        // V = r/beta
        cblas_dcopy(n, r, 1, V, 1);
        cblas_dscal(n, 1/beta, V, 1);

        for (j=0; j<m; j++) {
            status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j*n], 0.0, w);
            if (status != SPARSE_STATUS_SUCCESS) return -1;
            
            if (precon_type == ILU0) {
                ilu0_sweep_bsr_raw(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);
            }
            else if (precon_type == LUSGS) {
                // TODO : LU-SGS preconditioner function
            }
            
            // Arnoldi iteration
            for (i=0; i<(j+1); i++) {
                tmp = cblas_ddot(n, w, 1, &V[i*n], 1);
                H[j+m*i] = tmp;
                cblas_daxpy(n, -tmp, &V[i*n], 1, w, 1);
            }

            tmp = cblas_dnrm2(n, w, 1);
            H[j+(j+1)*m] = tmp;
            cblas_dcopy(n, w, 1, &V[(j+1)*n], 1);
            cblas_dscal(n, 1/tmp, &V[(j+1)*n], 1);

            // Givens Rotation
            for (i=0; i<j; i++) {
                c = g[i*2];         // g[i, 0]
                s = g[i*2+1];       // g[i, 1]
                h1 = H[j+i*m];
                h2 = H[j+(i+1)*m];
                H[j+i*m] = c*h1 - s*h2;
                H[j+(i+1)*m] = s*h1 + c*h2;
            }

            h1 = H[j*(m+1)];
            h2 = H[j+(j+1)*m];
            rr = sqrt(h1*h1 + h2*h2);
            c = h1/rr;
            s = -h2/rr;
            H[j*(m+1)] = c*h1 - s*h2;
            g[j*2] = c;
            g[j*2+1] = s;

            // Modify e1 vector
            y[j+1] = y[j];
            y[j] *= c;
            y[j+1] *= s;
        }

        // Back substitution
        cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

        // Update
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, dub, 1);

        // Computes next iteration residual
        cblas_dcopy(n, rhsb, 1, r, 1);
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
        if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");
    }

    if (it == itmax){
        printf("Not converged in %d iteration, residual=%.5f", it, beta);
        return it;
    }

    return UCFD_STATUS_ERROR;
}


// * Temporal arrays embedded ver.
// * Slightly slower than serial_gmres
ucfd_status_t serial_gmres2(const int neles, const int nvars, const int m, sparse_matrix_t op, const int *row_ptr, const int *col_ind, const int *diag_ind, \
                            ucfd_precon_type_t precon_type, double *pre_nnz_data, const double tol, const double itmax, double *dub, double *rhsb)
{
    int it, i, j, err;
    const int n = neles * nvars;
    double beta, tmp, c, s, h1, h2, rr;
    sparse_status_t status;

    // Sparse matrix description for system matrix A(op)
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = 0;
    descr.diag = 0;

    // Temporal arrays
    double *H, *V, *g, *y, *w, *r;
    H = (double *)malloc(sizeof(double)*m*(m+1));
    V = (double *)malloc(sizeof(double)*n*(m+1));
    g = (double *)malloc(sizeof(double)*2*m);
    y = (double *)malloc(sizeof(double)*(m+1));
    w = (double *)malloc(sizeof(double)*n);
    r = (double *)malloc(sizeof(double)*n);
    

    // Computes residual : r := rhsb - A @ dub
    cblas_dcopy(n, rhsb, 1, r, 1);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
    if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");

    for (it=0; it<itmax; it++) {
        beta = cblas_dnrm2(n, r, 1);
        if (beta < tol) return UCFD_STATUS_CONVERGED;

        if (precon_type == ILU0) {
            ilu0_sweep_bsr_raw(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, r);
            beta = cblas_dnrm2(n, r, 1);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
            beta = cblas_dnrm2(n, r, 1);
        }

        y[0] = beta;

        // V = r/beta
        cblas_dcopy(n, r, 1, V, 1);         // ! Segmentation Fault 왜?????????????
        cblas_dscal(n, 1/beta, V, 1);

        for (j=0; j<m; j++) {
            status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j*n], 0.0, w);
            if (status != SPARSE_STATUS_SUCCESS) return -1;
            
            if (precon_type == ILU0) {
                ilu0_sweep_bsr_raw(neles, nvars, row_ptr, col_ind, diag_ind, pre_nnz_data, w);
            }
            else if (precon_type == LUSGS) {
                // TODO : LU-SGS preconditioner function
            }
            
            // Arnoldi iteration
            for (i=0; i<(j+1); i++) {
                tmp = cblas_ddot(n, w, 1, &V[i*n], 1);
                H[j+m*i] = tmp;
                cblas_daxpy(n, -tmp, &V[i*n], 1, w, 1);
            }

            tmp = cblas_dnrm2(n, w, 1);
            H[j+(j+1)*m] = tmp;
            cblas_dcopy(n, w, 1, &V[(j+1)*n], 1);
            cblas_dscal(n, 1/tmp, &V[(j+1)*n], 1);

            // Givens Rotation
            for (i=0; i<j; i++) {
                c = g[i*2];         // g[i, 0]
                s = g[i*2+1];       // g[i, 1]
                h1 = H[j+i*m];
                h2 = H[j+(i+1)*m];
                H[j+i*m] = c*h1 - s*h2;
                H[j+(i+1)*m] = s*h1 + c*h2;
            }

            h1 = H[j*(m+1)];
            h2 = H[j+(j+1)*m];
            rr = sqrt(h1*h1 + h2*h2);
            c = h1/rr;
            s = -h2/rr;
            H[j*(m+1)] = c*h1 - s*h2;
            g[j*2] = c;
            g[j*2+1] = s;

            // Modify e1 vector
            y[j+1] = y[j];
            y[j] *= c;
            y[j+1] *= s;
        }

        // Back substitution
        cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

        // Update
        cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, dub, 1);

        // Computes next iteration residual
        cblas_dcopy(n, rhsb, 1, r, 1);
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, -1.0, op, descr, dub, 1.0, r);
        if (status != SPARSE_STATUS_SUCCESS) printf("Residual computation error\n");
    }

    if (it == itmax){
        printf("Not converged in %d iteration, residual=%.5f", it, beta);
        return it;
    }

    free(H);
    free(V);
    free(g);
    free(y);
    free(w);
    free(r);

    return UCFD_STATUS_ERROR;
}                            







ucfd_status_t single_gmres_mkl(sparse_matrix_t op, ucfd_precon_type_t precon_type, sparse_matrix_t precon, const int n, const int m, \
                               double *dub, double *rhsb, double *H, double *V, double *g, double *y, double *w, double *r)
{
    int it, i, j, err;
    double beta, tmp, c, s, h1, h2, rr;
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

    beta = cblas_dnrm2(n, r, 1);
    
    if (precon_type == ILU0) {
        ilu0_sweep_bsr_mkl(precon, r);
        beta = cblas_dnrm2(n, r, 1);
    }
    else if (precon_type == LUSGS) {
        // TODO : LU-SGS preconditioner function
        beta = cblas_dnrm2(n, r, 1);
    }
    
    // beta = cblas_dnrm2(n, r, 1);
    y[0] = beta;

    // V = r/beta
    cblas_dcopy(n, r, 1, V, 1);
    cblas_dscal(n, 1/beta, V, 1);

    for (j=0; j<m; j++) {
        status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, op, descr, &V[j*n], 0.0, w);
        if (status != SPARSE_STATUS_SUCCESS)
            return -1;
        
        if (precon_type == ILU0) {
            ilu0_sweep_bsr_mkl(precon, w);
        }
        else if (precon_type == LUSGS) {
            // TODO : LU-SGS preconditioner function
        }
        
        // Arnoldi iteration
        for (i=0; i<j+1; i++) {
            tmp = cblas_ddot(n, w, 1, &V[i*n], 1);
            H[j+m*i] = tmp;
            cblas_daxpy(n, -tmp, &V[i*n], 1, w, 1);
        }

        tmp = cblas_dnrm2(n, w, 1);
        H[j+(j+1)*m] = tmp;
        cblas_dcopy(n, w, 1, &V[(j+1)*n], 1);
        cblas_dscal(n, 1/tmp, &V[(j+1)*n], 1);

        // Givens Rotation
        for (i=0; i<j; i++) {
            c = g[i*2];
            s = g[i*2+1];
            h1 = H[j+i*m];
            h2 = H[j+(i+1)*m];
            H[j+i*m] = c*h1 - s*h2;
            H[j+(i+1)*m] = s*h1 + c*h2;
        }

        h1 = H[j*(m+1)];
        h2 = H[j+(j+1)*m];
        rr = sqrt(h1*h1 + h2*h2);
        c = h1/rr;
        s = -h2/rr;
        H[j*(m+1)] = c*h1 - s*h2;
        g[j*2] = c;
        g[j*2+1] = s;

        // Modify e1 vector
        y[j+1] = y[j];
        y[j] *= c;
        y[j+1] *= s;
    }

    // Back substitution
    cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H, m, y, 1);

    // Update
    cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, V, n, y, 1, 1.0, dub, 1);

    return UCFD_STATUS_SUCCESS;
}