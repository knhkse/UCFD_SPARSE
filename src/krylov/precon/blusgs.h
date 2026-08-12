#pragma once

#include "ucfdpcimpl.h"

typedef struct {
    UCFDInt bn;
    UCFDInt block;
    UCFDReal *diagvalues;
} Precon_BLUSGS;

typedef struct {
    Precon_BLUSGS base;
    UCFDInt ncolors;
    UCFDInt *icolors;
} Precon_PBLUSGS;

UCFD_INTERN ucfd_status_t BLUSGSPreconPrepare(Precon);
UCFD_INTERN ucfd_status_t BLUSGSPreconDestroy(Precon);
