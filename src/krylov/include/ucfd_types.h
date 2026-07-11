/**
 * @file        ucfd_types.h
 * @brief       Header file for UCFD_SPRASE Library types
 */
#ifndef _UCFD_TYPES_H
#define _UCFD_TYPES_H
#include <stdbool.h>
#if defined(USE_MKL)
    #include <mkl.h>
#endif

#include "config.h"
#include "macros.h"

#define UCFD_FALSE false
#define UCFD_TRUE true

#define INITFLAG 1234


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
    NONE = 0,
    BILU = 1,
    LUSGS = 2
} ucfd_precon_type_t;

typedef bool UCFDBool;


#endif  // _UCFD_TYPES_H