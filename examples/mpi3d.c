/** ======================================================================================================================
 * @file        mpi3d.c
 * @brief       Example for LU-SGS method using MPI parallelism.
 * @details     This file is an executable example file for LU-SGS method with MPI parallelism.
 *              Basic hexahedral gird is used, and 
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

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "pbutils.h"
#include "readinput.h"
#include "arrays.h"
#include "lusgs.h"

#define nvars 5
#define ndims 3
#define ncellface 6
#define root 0

double run(int rank, int *params);
void write_output(int nprocs, int *params, double *times);

int main(int argc, char **argv)
{
    int nprocs, rank;
    int dm, dn, dl;
    int err = 0;
    int params[7];
    double avtimei;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *times = (double *)calloc(nprocs, sizeof(double));

    if (rank == 0) {
        if (argc < 2) {
            printf("[Error] Input file argument is needed\n");
            err = 1;
        }
        else {
            FILE *inpfile = fopen(argv[1], "r");
            if (inpfile == NULL) {
                printf("[Error] Input file not found\n");
                err = 1;
            }
            else {
                input_params(inpfile, params);
                dm = params[3];
                dn = params[4];
                dl = params[5];
                fclose(inpfile);
                
                if (nprocs != dm*dn*dl) {
                    printf("[Error] Partitioning number mismatch\n");
                    printf("%d processors, but total %d partitioning\n", nprocs, dm*dn*dl);
                    err = 1;
                }
                else {
                    printf("[Main] 3D hexahedral example starts...\n");
                }
            }
        }
    }

    MPI_Bcast(&err, 1, MPI_INT, root, MPI_COMM_WORLD);
    if (err == 1) return -1;

    MPI_Bcast(&params, 7, MPI_INT, root, MPI_COMM_WORLD);

    avtimei = run(rank, params);
    
    MPI_Gather(&avtimei, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("[Main] Run finished. Writing result...\n");
        write_output(nprocs, params, times);
    }

    MPI_Finalize();
    free(times);
    return 0;
}


double run(int rank, int *params)
{
    const double gamma = 1.4;
    const double rho = gamma;
    const double u = 1.5;
    const double v = 0.0;
    const double w = 0.0;
    const double p = 1.0;
    const double et = p / (gamma - 1) + 0.5*rho*u*u;
    const int bm = params[0]/params[3];
    const int bn = params[1]/params[4];
    const int bl = params[2]/params[5];
    const int neles = bm*bn*bl;
    const int nstep = params[6];
    double start;
    double avtime = 0.0;

    if (rank == 0) printf("[Run] Allocating arrays...\n");
    double **upts = (double **) malloc_2d(nvars, neles, sizeof(double));
    double **rhs = (double **) malloc_2d(nvars, neles, sizeof(double));
    double **dub = (double **) malloc_2d(nvars, neles, sizeof(double));
    double *diag = (double *) calloc(neles, sizeof(double));
    double **fspr = (double **) malloc_2d(ncellface, neles, sizeof(double));
    double *dt = (double *) calloc(neles, sizeof(double));
    double **fnorm_vol = (double **) malloc_2d(ncellface, neles, sizeof(double));
    double ***vec_fnorm = (double ***) malloc_3d(ncellface, ndims, neles, sizeof(double));
    int **nei_ele = (int **) malloc_2d(ncellface, neles, sizeof(int));
    int *mapping = (int *)calloc(neles, sizeof(int));
    int *unmapping = (int *)calloc(neles, sizeof(int));
    double *tarr = (double *)malloc(sizeof(double)*nstep);

    /* Pointer of each array */
    double *uptsp = &upts[0][0];
    double *rhsp = &rhs[0][0];
    double *dubp = &dub[0][0];
    double *fsprp = &fspr[0][0];
    double *fnp = &fnorm_vol[0][0];
    double *vfp = &vec_fnorm[0][0][0];
    int *nep = &nei_ele[0][0];

    /* Initialize */
    if (rank == 0) printf("[Run] Initializing arrays...\n");
    for (int i=0; i < neles; i++) {
        upts[0][i] = rho;
        upts[1][i] = rho*u;
        upts[2][i] = rho*v;
        upts[3][i] = rho*w;
        upts[4][i] = et;

        rhs[0][i] = rho;
        rhs[1][i] = rho*u;
        rhs[2][i] = rho*v;
        rhs[3][i] = rho*w;
        rhs[4][i] = et;

        dub[0][i] = rho;
        dub[1][i] = rho*u;
        dub[2][i] = rho*v;
        dub[3][i] = rho*w;
        dub[4][i] = et;

        dt[i] = 0.1;

        for (int j=0; j<ncellface; j++) {
            fnorm_vol[j][i] = 1.0;
            fspr[j][i] = 1.0;
            for (int k=0; k<ndims; k++) {
                vec_fnorm[j][k][i] = 0.0;
            }
        }

        vec_fnorm[0][1][i] = -1.0;
        vec_fnorm[1][0][i] = -1.0;
        vec_fnorm[2][2][i] = 1.0;
        vec_fnorm[3][0][i] = 1.0;
        vec_fnorm[4][2][i] = -1.0;
        vec_fnorm[5][1][i] = 1.0;
    }

    if (rank == 0) printf("[Run] Constructing nei_ele array...\n");
    make_nei_ele(bm, bn, bl, nei_ele);

    if (rank == 0) printf("[Run] Processing Reverse Cuthill-McKee...\n");
    make_reordering(neles, ncellface, nei_ele, mapping, unmapping);

    if (rank == 0) printf("[Run] Starting iteration...\n");
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0; i<nstep; i++) {
        start = MPI_Wtime();
        serial_pre_lusgs(neles, ncellface, 1.0, fnp, dt, diag, fsprp);
        ns_serial_lower_sweep(neles, nvars, ncellface, ndims,
                    nep, mapping, unmapping, fnp, vfp, \
                    uptsp, rhsp, dubp, diag, fsprp);
        ns_serial_upper_sweep(neles, nvars, ncellface, ndims,
                    nep, mapping, unmapping, fnp, vfp, \
                    uptsp, rhsp, dubp, diag, fsprp);
        serial_update(neles, nvars, uptsp, rhsp);
        MPI_Barrier(MPI_COMM_WORLD);
        tarr[i] = (double) MPI_Wtime() - start;
    }

    for (int i=0; i<nstep; i++)
        avtime += tarr[i];
    avtime = avtime*1000/nstep;

    //deallocate array
    dealloc_2d((void **) upts);
    dealloc_2d((void **) rhs);
    dealloc_2d((void **) dub);
    dealloc_2d((void **) fspr);
    dealloc_2d((void **) fnorm_vol);
    dealloc_3d((void ***) vec_fnorm);
    dealloc_2d((void **) nei_ele);
    free(mapping);
    free(unmapping);
    free(diag);
    free(dt);
    free(tarr);

    return avtime;
}

/**
 * @brief       Write output file as txt format
 */
void write_output(int nprocs, int *params, double *times)
{
    int m = params[0];
    int n = params[1];
    int l = params[2];
    int dm = params[3];
    int dn = params[4];
    int dl = params[5];
    int nstep = params[6];
    int neles = m*n*l;
    double avtime = 0.0;
    char filename[100];

    // Make filename
    sprintf(filename, "MPIOUTPUT_c_%d_%d.txt", nprocs, neles);
    FILE *outfile = fopen(filename, "w");
    fprintf(outfile, "========= LU-SGS example output =========\n\n");
    fprintf(outfile, "*** Problem setup ***\n");
    fprintf(outfile, "Number of cores: %d\n", nprocs);
    fprintf(outfile, "Number of iteration step: %d\n", nstep);
    fprintf(outfile, "Number of elements: %d = %d x %d x %d\n", neles, m, n, l);
    fprintf(outfile, "Partitioning with %d: %d x %d x %d\n", dm*dn*dl, dm, dn, dl);
    fprintf(outfile, "Elements per core: %d = %d x %d x %d\n\n", neles/nprocs, m/dm, n/dn, l/dl);
    fprintf(outfile, "*** Average runtime for each processor [ms] ***\n");
    for (int i=0; i<nprocs; i++) {
        avtime += times[i];
        fprintf(outfile, "%lf\n", times[i]);
    }
    fprintf(outfile, "\n*** Average runtime for entire processors [ms] ***\n");
    fprintf(outfile, "%lf\n", avtime/nprocs);

    fclose(outfile);
}

