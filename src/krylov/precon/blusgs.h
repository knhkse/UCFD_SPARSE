#pragma once

#include "ucfdpcimpl.h"

typedef struct {
    UCFDInt bn;
    UCFDInt block;
    UCFDReal *diagvalues;
} Precon_BLUSGS;

