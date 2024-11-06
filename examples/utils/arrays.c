/** ======================================================================================================================
 * @file        arrays.c
 * @brief       Array allocation functions for C example files
 * @details     This file defines some functions for allocating multi-dimensional array.
 *              Each function needs shape of array and size of data type.
 *              One-dimensional array can be allocated simply with `malloc` or `calloc` function in <stdlib.h>. 
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
#include "arrays.h"


/**
 * @brief       Allocate 2D array
 * @param       rows        number of rows
 * @param       cols        number of columns
 * @param       T           size of data
 */
void **malloc_2d(const size_t rows, const size_t cols, const size_t T)
{
	void *data= (void *) malloc(rows*cols*T);
    void **ar = (void **)malloc(rows*sizeof(void*));
	int i;
	char *p = (char*)data;
    for (i=0; i<rows; i++)
        ar[i] = &(p[cols*i*T]);
    return ar;
}


/**
 * @brief       Allocate 3D array
 * @param       rows        number of rows
 * @param       cols        number of columns
 * @param       depth       depth length, equals with z-direction elements
 * @param       T           size of data
 */
void ***malloc_3d(const size_t rows, const size_t cols, const size_t depth, const size_t T)
{
	void *data = (void *) malloc(rows*cols*depth*T);  // using calloc
	void **ar2 = (void **)malloc(rows*cols*sizeof(void*));
    void ***ar = (void ***)malloc(rows*sizeof(void**));
	int i,j;
	char *p = (char*)data;
    for (i=0; i<rows; i++){
		for (j=0; j<cols; j++){
			ar2[cols*i+j] = &(p[(cols*i+j)*depth*T]);
		}
		ar[i] = &(ar2[cols*i]);
	}
    return ar;
}


/**
 * @brief       Deallocate 2D array
 * @param       mat         2D array
 */
void dealloc_2d(void **mat)
{
    free(*mat);
    free(mat);
}


/**
 * @brief       Deallocate 3D array
 * @param       mat         3D array
 */
void dealloc_3d(void ***mat)
{
    free(**mat);
    free(*mat);
    free(mat);
}