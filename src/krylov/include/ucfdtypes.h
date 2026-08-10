/**
 * @file        ucfdtypes.h
 * @brief       Header file for UCFD_SPRASE Library types
 */
#pragma once

#include <stdbool.h>
#include <stdint.h>
#include "config.h"
#include "macros.h"

#define UCFD_FALSE false
#define UCFD_TRUE true

typedef int8_t UCFDInt8;

typedef enum
{
    UCFD_SUCCESS = 0,
    UCFD_FAILED  = 1
} ucfd_status_t;

typedef enum
{
    INITIALIZED      = -1,
    CONVERGED        = 0,
    HAPPYBREAKDOWN   = 1,
    REACH_ITERMAX    = 2,
    RHOBREAKDOWN     = 3,
    PIBREAKDOWN      = 4
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


typedef bool UCFDBool;
