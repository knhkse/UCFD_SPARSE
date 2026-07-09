#pragma once

#include "ucfdkrylovimpl.h"


typedef struct {
    UCFDReal *r, *rt, *p, *pt, *v, *s, *shat, *t;
} Solver_BICGSTAB;