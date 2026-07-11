#pragma once

#include "ucfdprecon.h"


typedef struct _PreconOps *PreconOps;
struct _PreconOps {
    ucfd_status_t (*prepare)(Precon);
    ucfd_status_t (*apply)(Precon, UCFDReal*);
    ucfd_status_t (*destroy)(Precon);
};


struct _Precon {
    PreconType type_name;
    UCFDInt *rowptr;
    UCFDInt *colidx;
    UCFDInt *diagslots;
    UCFDReal *values;
    void *data;
    struct _PreconOps ops[1];
};
