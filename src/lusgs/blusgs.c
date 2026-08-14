#include "flowsys.h"
#include "flux.h"


/**
 * Pack & update kernels => per-element execution
 */
ucfd_status_t UCFDBLUSGS_Pack(FlowSys sys, UCFDInt eidx,
                              UCFDReal a0, UCFDReal turb_factor)
{
    FlowElem *e  = &sys->eles[eidx];
    const UCFDInt nvars = sys->nvars, nfvars = sys->nfvars;
    const UCFDInt *cell_ids = e->cell_ids;
    // UCFDReal *rank_rhs = sys->du, 




}














