#pragma once

#include "ucfdkrylovimpl.h"

typedef struct {
    UCFDInt restart;    // Restart number
    UCFDReal *H, *V, *y, *w, *sn, *cs, *htmp, *r;      // Working arrays
} Solver_GMRES;
