#pragma once

#include "ucfdpcimpl.h"

typedef struct {
    UCFDInt *iw;
    UCFDInt n;
} Precon_ILU;

typedef struct {
    UCFDInt *iw;
    UCFDInt bn;
    UCFDInt block;
} Precon_BILU;


typedef struct {
    Precon_BILU base;
    UCFDInt ncolors;
    UCFDInt *icolors;
} Precon_PBILU;

/* Destroy inner structure */
UCFD_INTERN ucfd_status_t ILUPreconDestroy(Precon precon);