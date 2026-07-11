#pragma once

#include "ucfdpcimpl.h"

typedef struct {
    UCFDInt bn;
    UCFDInt block;
    UCFDInt *iw;
} Precon_BILU;


typedef struct {
    Precon_BILU base;
    UCFDInt ncolors;
    UCFDInt *icolors;
} Precon_PBILU;


UCFD_INTERN ucfd_status_t BILUPreconDestroy(Precon precon);