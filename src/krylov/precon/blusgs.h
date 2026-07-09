#pragma once

#include "ucfdpcimpl.h"

typedef struct {
    UCFDInt bn;
    UCFDInt block;
    UCFDReal *diagvalues;
    UCFDReal *values;    /* Same with system matrix */
} Precon_BLUSGS;

