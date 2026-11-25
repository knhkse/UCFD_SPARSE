#pragma once

#include <assert.h>
#include <math.h>
#include "config.h"

static inline void block_matvec(UCFDInt nvars, const UCFDReal *A, UCFDReal *b)
{
    UCFDReal work[nvars];

    /* Preserve the original b until all output components are computed. */
    for (UCFDInt row = 0; row < nvars; ++row) {
        UCFDReal value = A[row*nvars]*b[0];
        for (UCFDInt col = 1; col < nvars; ++col) {
            value += A[row*nvars + col]*b[col];
        }
        work[row] = value;
    }
    for (UCFDInt row = 0; row < nvars; ++row) {
        b[row] = work[row];
    }
}

static inline void pivoted_inverse(UCFDInt nvars, UCFDReal *A, UCFDReal pivot_tol)
{
    UCFDReal B[nvars * nvars];

    assert(nvars >= 1 && nvars <= nvars);

    for (UCFDInt row = 0; row < nvars; ++row) {
        for (UCFDInt col = 0; col < nvars; ++col) {
            B[row*nvars + col] = row == col ? 1.0 : 0.0;
        }
    }

    /* Compute PA = LU, applying the same row swaps to the identity. */
    for (UCFDInt k = 0; k < nvars; ++k) {
        UCFDInt pivot = k;
        UCFDReal pivot_abs = fabs(A[k*nvars + k]);
        for (UCFDInt row = k + 1; row < nvars; ++row) {
            UCFDReal value_abs = fabs(A[row*nvars + k]);
            if (value_abs > pivot_abs) {
                pivot = row;
                pivot_abs = value_abs;
            }
        }

        if (pivot != k) {
            for (UCFDInt col = 0; col < nvars; ++col) {
                UCFDReal tmp = A[k*nvars + col];
                A[k*nvars + col] = A[pivot*nvars + col];
                A[pivot*nvars + col] = tmp;
                tmp = B[k*nvars + col];
                B[k*nvars + col] = B[pivot*nvars + col];
                B[pivot*nvars + col] = tmp;
            }
        }

        if (fabs(A[k*nvars + k]) < pivot_tol) {
            A[k*nvars + k] = A[k*nvars + k] >= 0.0
                          ? pivot_tol : -pivot_tol;
        }

        for (UCFDInt row = k + 1; row < nvars; ++row) {
            A[row*nvars + k] /= A[k*nvars + k];
            for (UCFDInt col = k + 1; col < nvars; ++col) {
                A[row*nvars + col] -= A[row*nvars + k]*A[k*nvars + col];
            }
        }
    }

    /* Solve LY = P (the diagonal of L is implicitly one). */
    for (UCFDInt row = 1; row < nvars; ++row) {
        for (UCFDInt col = 0; col < nvars; ++col) {
            UCFDReal sum = A[row*nvars]*B[col];
            for (UCFDInt k = 1; k < row; ++k) {
                sum += A[row*nvars + k]*B[k*nvars + col];
            }
            B[row*nvars + col] -= sum;
        }
    }

    /* Solve UX = Y; this reverse loop avoids unsigned index underflow. */
    for (UCFDInt remaining = nvars; remaining > 0; --remaining) {
        UCFDInt row = remaining - 1;
        for (UCFDInt col = 0; col < nvars; ++col) {
            UCFDReal value = B[row*nvars + col];
            for (UCFDInt k = row + 1; k < nvars; ++k) {
                value -= A[row*nvars + k]*B[k*nvars + col];
            }
            B[row*nvars + col] = value/A[row*nvars + row];
        }
    }

    for (UCFDInt row = 0; row < nvars; ++row) {
        for (UCFDInt col = 0; col < nvars; ++col) {
            A[row*nvars + col] = B[row*nvars + col];
        }
    }
}
