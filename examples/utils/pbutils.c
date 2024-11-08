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

#include <stdio.h> 
#include <stdlib.h>
#include "pbutils.h"
#include "queue.h"


/**
 * @brief       Computes neighbor cell elements array.
 * 
 * @param       m           number of elements in x-direction
 * @param       n           number of elements in y-direction
 * @param       l           number of elements in z-direction
 * @param       nei_ele     array storing neighbor elements
 */
void make_nei_ele(const int m, const int n, const int l, int **nei_ele)
{

    int idx = 0;
    int ml = m*l;
    int nml = n*m*l;
    int ny, py, nx, px, nz, pz;

    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
            for (int k=0; k<l; k++) {
                // Next or previous neighbor cell count
                ny = idx-ml;
                py = idx+ml;
                nx = idx-l;
                px = idx+l;
                nz = idx-1;
                pz = idx+1;

                // Handling grid boundary
                if (ny<0) nei_ele[0][idx] = idx; else nei_ele[0][idx] = ny;
                if (nx<j*ml) nei_ele[1][idx] = idx; else nei_ele[1][idx] = nx;
                if (pz >= j*ml+(i+1)*l) nei_ele[2][idx] = idx; else nei_ele[2][idx] = pz;
                if (px >= (j+1)*ml) nei_ele[3][idx] = idx; else nei_ele[3][idx] = px;
                if (nz < j*ml+i*l) nei_ele[4][idx] = idx; else nei_ele[4][idx] = nz;
                if (py >= nml) nei_ele[5][idx] = idx; else nei_ele[5][idx] = py;

                idx++;
            }
        }
    }
}


/**
 * @brief       Search element in the given array.
 *              If element exists, return 1. 
 * 
 * @param       val         Value for searching
 * @param       arr         Object array
 * @param       size        Size of array
 */
int searcharr(int val, int *arr, int size)
{
    for (int i=0; i<size; i++) {
        // If value is found, return 1
        if (arr[i] == val)
            return 1;
    }
    // If not found, return 0
    return 0;
}


/**
 * @brief       Reverse Cuthill-McKee algorithm using neighbor elements.
 * 
 * @param       nele        Number of element cells
 * @param       nface       Number of faces depends on element type
 * @param       nei_ele     Indices for neighbor cells [nface, neles]
 * @param       mapping     Reordered index by Reverse Cuthill-McKee algorithm [neles]
 * @param       unmapping   Original index before Reverse Cuthill-McKee algorithm [neles]
 */
void make_reordering(const int nele, const int nface, int **nei_ele, int *mapping, int *unmapping)
{
    struct Queue q;             /** Queue data structure */
    int fe = nele - 1;          /** Last element index */
    int idx, jdx, ele, neib;    /** Indexing variables */

    // Initialize queue structure
    initQueue(&q);
    // Starts from the last element index
    enqueue(&q, fe);
    mapping[0] = fe;
    idx = 1;

    // Compute `mapping` array
    while (!isEmpty(&q)) {
        ele = dequeue(&q);
        for (jdx=0; jdx<nface; jdx++) {
            neib = nei_ele[jdx][ele];

            // If neighbor element index is not stored in the mapping array,
            // append index in mapping and queue
            if (!searcharr(neib, mapping, nele)) {
                mapping[idx] = neib;
                enqueue(&q, neib);
                idx++;
            }
        }
    }
    
    // Compute `unmapping` array
    for (int i=0; i<nele; i++) {
        idx = 0;
        
        // Find each mapping element index
        while (mapping[idx] != i)
            idx ++;
        
        unmapping[i] = idx;
    }
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
 * @param       lcolor      Color level of each cell     
 */
void make_coloring(const int m, const int n, const int l, \
                    int *icolor, int *lcolor)
{
    int neles = m*n*l;
    int cb = neles/2;
    int denom = cb;
    int ml = m*l;
    int ca = 0;
    int ele = 0;
    int ind, j, ist, iend, tmp;

    // Starts from zero-index element
    while (ele < neles) {
        for (j=0; j<n; j++) {
            ist = (j*ml)/2;
            iend = (j+1)*ml/2;
            for (ind=ist; ind<iend; ind++) {
                icolor[ca+ind] = ele;
                lcolor[ele] = ca/denom;
                ele++;
                icolor[cb+ind] = ele;
                lcolor[ele] = cb/denom;
            }
            tmp = ca;
            ca = cb;
            cb = tmp;
        }
    }
}
