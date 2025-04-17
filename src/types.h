#ifndef TYPES_H
#define TYPES_H

/**
 * @file        types.h
 * @brief       Header file for UCFD_SPRASE Library types
 */

typedef enum
{
    UCFD_STATUS_RHO_BREAKDOWN = -2,
    UCFD_STATUS_NOT_CONVERGED = -1,
    UCFD_STATUS_CONVERGED = 0,
    UCFD_STATUS_SUCCESS = 1,
    UCFD_STATUS_ERROR = 2
} ucfd_status_t;

typedef enum
{
    NONE = 0,
    BILU = 1,
    LUSGS = 2
} ucfd_precon_type_t;

#endif // TYPES_H