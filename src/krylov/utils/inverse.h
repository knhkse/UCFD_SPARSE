/**
 * @file        inverse.h
 * @brief       Header file for LU decomposition/substitution
 */
#pragma once

#include "config.h"
#include "macros.h"


MAYBE_UNUSED void ludcmp(const UCFDInt, UCFDReal*);
MAYBE_UNUSED void lusub(const UCFDInt, const UCFDReal*, UCFDReal*);
MAYBE_UNUSED void lusubmattrans(const UCFDInt, const UCFDReal*, UCFDReal*);
