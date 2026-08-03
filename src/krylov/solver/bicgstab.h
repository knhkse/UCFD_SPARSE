#pragma once

#include "ucfdsolverimpl.h"


typedef struct {
    UCFDInt     n;
    UCFDReal    *r, *rt, *p, *pt, *v, *s, *shat, *t;
} Solver_BICGSTAB;