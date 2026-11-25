#pragma once

#include <float.h>
#include <stddef.h>
#include <math.h>

#include "config.h"

template <UCFDInt nvars>
__device__ __forceinline__
void block_matvec(const UCFDReal *A, UCFDReal *b)
{
    UCFDReal work[nvars];
    /* Preserve the original b until all output components are computed. */
    for (size_t row = 0; row < nvars; ++row) {
        UCFDReal value = A[row*nvars]*b[0];
        for (size_t col = 1; col < nvars; ++col) {
            value += A[row*nvars + col]*b[col];
        }
        work[row] = value;
    }
    for (size_t row = 0; row < nvars; ++row) {
        b[row] = work[row];
    }
}


template <UCFDInt nvars>
__device__ __forceinline__
void pivoted_inverse(UCFDReal *A, UCFDReal pivot_tol = DBL_EPSILON)
{
    UCFDReal B[nvars * nvars];

    static_assert(nvars >= 1 && nvars <= 7, "nvars must be between 1 and 7");

    for (size_t row = 0; row < nvars; ++row) {
        for (size_t col = 0; col < nvars; ++col) {
            B[row*nvars + col] = row == col ? 1.0 : 0.0;
        }
    }

    /* Compute PA = LU, applying the same row swaps to the identity. */
    for (size_t k = 0; k < nvars; ++k) {
        size_t pivot = k;
        UCFDReal pivot_abs = fabs(A[k*nvars + k]);
        for (size_t row = k + 1; row < nvars; ++row) {
            UCFDReal value_abs = fabs(A[row*nvars + k]);
            if (value_abs > pivot_abs) {
                pivot = row;
                pivot_abs = value_abs;
            }
        }

        if (pivot != k) {
            for (size_t col = 0; col < nvars; ++col) {
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

        for (size_t row = k + 1; row < nvars; ++row) {
            A[row*nvars + k] /= A[k*nvars + k];
            for (size_t col = k + 1; col < nvars; ++col) {
                A[row*nvars + col] -= A[row*nvars + k]*A[k*nvars + col];
            }
        }
    }

    /* Solve LY = P (the diagonal of L is implicitly one). */
    for (size_t row = 1; row < nvars; ++row) {
        for (size_t col = 0; col < nvars; ++col) {
            UCFDReal sum = A[row*nvars]*B[col];
            for (size_t k = 1; k < row; ++k) {
                sum += A[row*nvars + k]*B[k*nvars + col];
            }
            B[row*nvars + col] -= sum;
        }
    }

    /* Solve UX = Y; this reverse loop avoids unsigned index underflow. */
    for (size_t remaining = nvars; remaining > 0; --remaining) {
        size_t row = remaining - 1;
        for (size_t col = 0; col < nvars; ++col) {
            UCFDReal value = B[row*nvars + col];
            for (size_t k = row + 1; k < nvars; ++k) {
                value -= A[row*nvars + k]*B[k*nvars + col];
            }
            B[row*nvars + col] = value/A[row*nvars + row];
        }
    }

    for (size_t row = 0; row < nvars; ++row) {
        for (size_t col = 0; col < nvars; ++col) {
            A[row*nvars + col] = B[row*nvars + col];
        }
    }
}