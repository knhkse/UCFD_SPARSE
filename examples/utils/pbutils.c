/** ======================================================================================================================
 * @file        pbutils.c
 * @brief       
 * @details     This file contains utility functions to construct arrays used in example.
 *              Shape and size of each array are identical to the pyBaram formulation,
 *              and each function was modified for the unstructured hexahedral mesh format.
 * 
 * @author
 *              - Namhyoung Kim (knhkse@inha.edu), Department of Aerospace Engineering, Inha University
 *              - Jin Seok Park (jinseok.park@inha.ac.kr), Department of Aerospace Engineering, Inha University
 * 
 * @date        July 2024
 * @version     1.0
 * @par         Copyright
 *              Copyright (c) 2024, Namhyoung Kim and Jin Seok Park, Inha University, All rights reserved.
 * @par         License
 *              This project is release under the terms of the MIT License (see LICENSE file).
 * 
 * =======================================================================================================================
 */
#include "pbutils.h"
#include <stdio.h> 
#include <stdlib.h>


enum {
    FACE_XM = 0,  // -x
    FACE_XP = 1,  // +x
    FACE_YM = 2,  // -y
    FACE_YP = 3,  // +y
    FACE_ZM = 4,  // -z
    FACE_ZP = 5   // +z
};

/**
 * @brief       Computes neighbor cell elements array.
 * 
 * @param       m           number of elements in x-direction
 * @param       n           number of elements in y-direction
 * @param       l           number of elements in z-direction
 * @param       nei_ele     array storing neighbor elements
 */
int make_nei_ele(const int m, const int n, const int l, int **adj)
{
    if (!adj || m <= 0 || n <= 0 || l <= 0) return -1;

    // Overflow-safe products for N = m*n*l
    int mn = m * n;
    if (m != 0 && mn / m != n) return -2;

    int N = mn * l;
    if (mn != 0 && N / mn != l) return -2;

    // Optional OpenMP: parallelize over (z,y) slabs if desired.
    // Compile with -fopenmp and uncomment pragmas.
    // #pragma omp parallel for collapse(2)
    for (int z = 0; z < l; ++z) {
        const int base_z = z * mn;
        for (int y = 0; y < n; ++y) {
            const int base_zy = base_z + y * m;
            for (int x = 0; x < m; ++x) {
                const int idx = base_zy + x;

                const int xm = (x > 0)     ? (idx - 1)  : idx;
                const int xp = (x + 1 < m) ? (idx + 1)  : idx;

                const int ym = (y > 0)     ? (idx - m)  : idx;
                const int yp = (y + 1 < n) ? (idx + m)  : idx;

                const int zm = (z > 0)     ? (idx - mn) : idx;
                const int zp = (z + 1 < l) ? (idx + mn) : idx;

                adj[(long)FACE_XM][idx] = xm;
                adj[(long)FACE_XP][idx] = xp;
                adj[(long)FACE_YM][idx] = ym;
                adj[(long)FACE_YP][idx] = yp;
                adj[(long)FACE_ZM][idx] = zm;
                adj[(long)FACE_ZP][idx] = zp;
            }
        }
    }

    return 0;
}


/**
 * @brief       Coloring algorithm for unstructured grid
 * @note        In this function, original multi-coloring algorithm
 *              is simplified into Red-Black coloring (checkerboard) scheme
 *              due to the properties of hexahedral mesh.
 * 
 * @param       m           number of elements in x-direction
 * @param       n           number of elements in y-direction
 * @param       l           number of elements in z-direction
 * @param       icolor      Color Index of each cell
 */
void make_coloring(int m, int n, int l, int *icolor)
{
    int N = m*n*l;
    int mn = m*n;
    int half = N / 2;
    int ka = 0;         /* write cursor for A */
    int kb = 0;         /* write cursor for B (offset by half) */

    for (int z = 0; z < l; ++z) {
        int base_z = z * mn;
        for (int y = 0; y < n; ++y) {
            int idx = base_z + y * m;          /* x=0 */
            int parity = (int)((y + z) & 1L);   /* parity at x=0 */

            for (int x = 0; x < m; ++x) {
                if (parity == 0) icolor[ka++] = idx;
                else             icolor[half + kb++] = idx;

                parity ^= 1; /* toggles when x increments */
                idx += 1;
            }
        }
    }
}
