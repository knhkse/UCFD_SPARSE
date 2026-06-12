/**
 * @file        precon.h
 * @brief       Header file for preconditioners for Krylov subspace methods
 * @details     Declaration of each prepare/psolve function.
 *              `prepare` must be executed before Krylov method routine starts.
 *              `psolve` function is executed in Krylov method to solve Px = b.
 */
#ifndef PRECON_H
#define PRECON_H
#include "ucfd_types.h"
#include "sparse.h"
#include "config.h"


/**
 * @brief       Block fill-in Incomplete LU preconditioner for BSR matrix format.
 * @param       precon      Preconditioner structure
 * @param       iw          Working array to compute BILU preparation
 */
ucfd_sparse_status_t bilu_prepare(ucfd_spmat *precon, UCFD_INT *iw);

/**
 * @brief       Solver function for BILU preconditioner.
 * @param       precon      Preconditioner structure
 * @param       b           Right-hand-side
 */
void bilu_psolve(ucfd_spmat *precon, UCFD_FLOAT *b);


/**
 * @brief       LU-SGS preconditioner for BSR matrix format.
 * @param       precon      Preconditioner structure
 */
ucfd_status_t lusgs_prepare(ucfd_spmat *precon);

/**
 * @brief       Solver function for LU-SGS preconditioner.
 * @param       precon      Preconditioner structure
 * @param       b           Right-hand-side
 */
void lusgs_psolve(ucfd_spmat *precon, UCFD_FLOAT *b);


/**
 * @brief       Unpreconditioned solver
 * @param       precon      Preconditioner structure
 * @param       b           Right-hand-side
 */
void none_psolve(ucfd_spmat *precon, UCFD_FLOAT *b);

#endif // PRECON.H