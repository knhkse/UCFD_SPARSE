#pragma once

#include "ucfdtypes.h"


#if defined(UCFD_FLOAT32)
    #ifndef BETAST
        #define BETAST 0.09f
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4f
    #endif
    #ifndef PMIN
        #define PMIN 1e-13f
    #endif
#else
    #ifndef BETAST
        #define BETAST 0.09
    #endif
    #ifndef GAMMA
        #define GAMMA 1.4
    #endif
    #ifndef PMIN
        #define PMIN 1e-13
    #endif
#endif

typedef void (*fluxfunc)(UCFDInt, UCFDInt, UCFDInt, UCFDReal*, UCFDReal*, UCFDReal*);

void ns_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f);
void rans_flux_container(UCFDInt nfvars, UCFDInt nturbvars, UCFDInt ndims, UCFDReal *u, UCFDReal *nf, UCFDReal *f);




