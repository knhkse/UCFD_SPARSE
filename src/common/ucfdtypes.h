/**
 * @file        ucfdtypes.h
 * @brief       Header file for UCFD_SPRASE Library types
 */
#pragma once

#include <stdint.h>
#include "config.h"
#include "macros.h"

#define UCFD_FALSE false
#define UCFD_TRUE true

typedef enum
{
    UCFD_SUCCESS = 0,
    UCFD_FAILED  = 1
} ucfd_status_t;

typedef enum
{
    INITIALIZED         = -2,
    ITERATING           = -1,
    CONVERGED           = 0,
    HAPPYBREAKDOWN      = 1,
    REACH_ITERMAX       = 2,
    DIVERGED_BREAKDOWN  = 3,
    DIVERGED_DTOL       = 4,
    RHO_BREAKDOWN        = 5,
    PI_BREAKDOWN         = 6
} ucfd_solver_t;

typedef enum
{
    UCFD_MPI_SUCCESS            = 0,
    UCFD_MPI_INVALID_ARGUMENT   = 1,
    UCFD_MPI_NOT_ACTIVE         = 2,
    UCFD_MPI_INVALID_COMM       = 3,
    UCFD_MPI_ALLOCATION_FAILED  = 4,
    UCFD_MPI_ERROR              = 5
} ucfd_mpi_t;
