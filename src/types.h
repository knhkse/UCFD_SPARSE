#ifndef TYPES_H
#define TYPES_H

typedef enum
{
    UCFD_STATUS_NOT_CONVERGED = -1,
    UCFD_STATUS_CONVERGED = 0,
    UCFD_STATUS_SUCCESS = 1,
    UCFD_STATUS_ERROR = 2
} ucfd_status_t;


typedef enum
{
    NONE = 0,
    ILU0 = 1,
    LUSGS = 2
} ucfd_precon_type_t;

typedef enum
{
    UCFD_PRECON_SPARSE_FAILED = -1,
    UCFD_PRECON_SUCCESS = 0,
    UCFD_PRECON_FAILED = 1
} ucfd_precon_status_t;

#endif // TYPES_H