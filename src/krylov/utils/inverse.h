/**
 * @file        inverse.h
 * @brief       Header file for LU decomposition/substitution
 */
#pragma once

#include "config.h"
#include "macros.h"


MAYBE_UNUSED void ludcmp(UCFDInt, UCFDReal*);
MAYBE_UNUSED void lusub(UCFDInt, UCFDReal*, UCFDReal*);
MAYBE_UNUSED void lusubmattrans(UCFDInt, UCFDReal*, UCFDReal*);
