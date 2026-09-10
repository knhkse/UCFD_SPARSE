#include "flowsys.h"
#include "flux.h"
#include "inverse_cuda.cuh"


/**
 * ! Notes !
 * In colored (Block) LU-SGS scheme,
 * computation kernels are created in element-unit.
 */

__global__ static void
_rank_blusgs_diag_pack(const UCFDInt neles, const UCFDInt nfvars,
                       const UCFDReal a0,
                       const UCFDInt *__restrict__ cell_ids,
                       const UCFDReal *__restrict__ dt,
                       UCFDReal *__restrict__ diag)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;
    
    const UCFDInt ridx = cell_ids[idx];
    const UCFDInt dim2 = nfvars*nfvars;
    UCFDInt row, col;

    for (row=0; row<nfvars; ++row) {
        for (col=0; col<nfvars; ++col)
            diag[col + row*nfvars + ridx*dim2] = 0.0;
        diag[row + row*nfvars + ridx*dim2] = 1/dt[idx] + a0;
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Pack(PFlowSys sys, UCFDInt eidx, UCFDReal a0)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;
    _rank_blusgs_diag_pack<<<bpg, TPB>>>(e->base.neles, sys->nfvars, a0,
                                         e->base.cell_ids, e->base.dt, blusgs->diag);
    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt nvars, rans_t rt>
__global__ static void
_rank_tblusgs_pack(const UCFDInt neles, const UCFDInt nturbvars,
                   const UCFDReal factor, const UCFDReal a0,
                   const UCFDInt *__restrict__ cell_ids,
                   const UCFDReal *__restrict__ upts,
                   const UCFDReal *__restrict__ dsrc,
                   const UCFDReal *__restrict__ dt,
                   UCFDReal *__restrict__ tdiag)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    const UCFDInt ridx = cell_ids[idx];
    const UCFDInt dim2 = nturbvars*nturbvars;
    UCFDInt kdx, row, col;
    UCFDReal u[nvars], d[nvars];

    for (row=0; row<nturbvars; ++row) {
        for (col=0; col<nturbvars; ++col)
            tdiag[col + row*nturbvars + ridx*dim2] = 0.0;
    }

    // Prepare cell uf and dsrc
    for (kdx=0; kdx<nvars; ++kdx) {
        u[kdx] = upts[idx + kdx*neles];
        d[kdx] = dsrc[idx + kdx*neles];
    }

    if (rt == RANS_KWSST)
        cuda_kwsst_src_jacobian(nvars, nturbvars, u, &tdiag[ridx*dim2], d);
    else if (rt == RANS_SA)
        cuda_sa_src_jacobian(nvars, nturbvars, u, &tdiag[ridx*dim2], d);

    for (row=0; row<nturbvars; ++row)
        tdiag[row + row*nturbvars + ridx*dim2] += 1/(dt[idx]*factor) + a0;
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_KWSST_Pack(PFlowSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;

    switch (sys->nvars) {
        case 4:
            _rank_tblusgs_pack<4, RANS_KWSST><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 5:
            _rank_tblusgs_pack<5, RANS_KWSST><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 6:
            _rank_tblusgs_pack<6, RANS_KWSST><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 7:
            _rank_tblusgs_pack<7, RANS_KWSST><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        default:
            fprintf(stderr, "Unsupported `nvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_SA_Pack(PFlowSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;

    switch (sys->nvars) {
        case 4:
            _rank_tblusgs_pack<4, RANS_SA><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 5:
            _rank_tblusgs_pack<5, RANS_SA><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 6:
            _rank_tblusgs_pack<6, RANS_SA><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        case 7:
            _rank_tblusgs_pack<7, RANS_SA><<<bpg, TPB>>>(
            e->base.neles, sys->nturbvars, turb_factor, a0,
            e->base.cell_ids, e->base.uptsb, e->base.dsrc, e->base.dt, blusgs->tdiag
        ); break;
        default: fprintf(stderr, "Unsupported `nvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


__global__ static void
_rank_blusgs_update(const UCFDInt neles, const UCFDInt nvars,
                    const UCFDInt *__restrict__ cell_ids,
                    const UCFDReal *__restrict__ rank_du,
                    UCFDReal *__restrict__ upts)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    UCFDInt kdx;
    const UCFDInt ridx = cell_ids[idx];
    for (kdx=0; kdx<nvars; ++kdx)
        upts[idx + kdx*neles] += rank_du[kdx + ridx*nvars];
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Update(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;

    _rank_blusgs_update<<<bpg, TPB>>>(
        e->base.neles, sys->nvars, e->base.cell_ids,
        blusgs->du, e->base.uptsb
    );
    UCFDFunctionReturn(UCFD_SUCCESS);
}

__global__ static void
_rank_sub_residual(const UCFDInt neles, const UCFDInt nvars,
                   const UCFDInt *__restrict__ cell_ids,
                   const UCFDReal *__restrict__ vol,
                   const UCFDReal *__restrict__ du,
                   UCFDReal *__restrict__ dup,
                   UCFDReal *__restrict__ res)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    UCFDInt kdx, offset;
    UCFDReal diff;
    const UCFDInt ridx = cell_ids[idx];
    
    for (kdx=0; kdx<nvars; ++kdx) {
        offset = kdx + ridx*nvars;
        diff = du[offset] - dup[offset];
        res[idx + kdx*neles] = diff*diff*vol[idx];
        dup[offset] = du[offset];
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_SubResidual(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;

    _rank_sub_residual<<<bpg, TPB>>>(
        e->base.neles, sys->nvars, e->base.cell_ids,
        e->base.vol, blusgs->du, blusgs->dup, e->base.resid_out
    );
    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt nvars>
__global__ static void
_ele_pre_blusgs(const UCFDInt neles, const UCFDInt nface, const UCFDInt nfaces,
                const UCFDInt *__restrict__ cell_ids,
                const UCFDInt *__restrict__ face_refs,
                const UCFDReal *__restrict__ face_factors,
                UCFDReal *__restrict__ diag,
                UCFDReal *__restrict__ offdiag,
                const UCFDReal *__restrict__ jmat
                )
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    const UCFDInt ridx = cell_ids[idx];
    const UCFDInt dim2 = nvars*nvars;
    UCFDInt fidx, row, col, face, slot, ind0, ind1;
    UCFDReal dmat[dim2], fv, dval, oval;

    for (row=0; row<nvars; ++row) {
        for (col=0; col<nvars; ++col)
            dmat[col + row*nvars] = diag[col + row*nvars + ridx*dim2];
    }

    for (fidx=0; fidx<nface; ++fidx)
    {
        face = face_refs[idx + fidx*neles];
        if (face > 0) slot = face - 1;
        else slot = -face - 1;
        fv = face_factors[idx + fidx*neles];

        for (row=0; row<nvars; ++row) {
            for (col=0; col<nvars; ++col) {
                ind0 = slot + col*nfaces + row*nvars*nfaces;
                ind1 = ind0 + dim2*nfaces;
                if (face > 0) { dval = jmat[ind0]; oval = jmat[ind1]; }
                else { dval = -jmat[ind1]; oval = -jmat[ind0]; }

                dmat[col + row*nvars] += dval*fv;
                offdiag[col + row*nvars + fidx*dim2 + idx*dim2*nface] = oval*fv;
            }
        }
    }

    ludcmp(nvars, dmat);

    for (row=0; row<nvars; ++row) {
        for (col=0; col<nvars; ++col)
            diag[col + row*nvars + ridx*dim2] = dmat[col + row*nvars];
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSPrepare(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;
    
    switch (sys->nfvars) {
        case 4:
            _ele_pre_blusgs<4><<<bpg, TPB>>>(
                e->base.neles, e->base.nface, sys->nfaces,
                e->base.cell_ids, e->face_refs, e->face_factors,
                blusgs->diag, e->flow_offdiag, blusgs->jmat
            ); break;
        case 5:
            _ele_pre_blusgs<5><<<bpg, TPB>>>(
                e->base.neles, e->base.nface, sys->nfaces,
                e->base.cell_ids, e->face_refs, e->face_factors,
                blusgs->diag, e->flow_offdiag, blusgs->jmat
            ); break;
        case 6:
            _ele_pre_blusgs<6><<<bpg, TPB>>>(
                e->base.neles, e->base.nface, sys->nfaces,
                e->base.cell_ids, e->face_refs, e->face_factors,
                blusgs->diag, e->flow_offdiag, blusgs->jmat
            ); break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSPrepare(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt bpg = (e->base.neles + TPB - 1)/TPB;
    
    switch (sys->nturbvars) {
        case 1:
            _ele_pre_blusgs<1><<<bpg, TPB>>>(
                e->base.neles, e->base.nface, sys->nfaces,
                e->base.cell_ids, e->face_refs, e->face_factors,
                blusgs->tdiag, e->turb_offdiag, blusgs->tjmat
            ); break;
        case 2:
            _ele_pre_blusgs<2><<<bpg, TPB>>>(
                e->base.neles, e->base.nface, sys->nfaces,
                e->base.cell_ids, e->face_refs, e->face_factors,
                blusgs->tdiag, e->turb_offdiag, blusgs->tjmat
            ); break;
        default:
            fprintf(stderr, "Unsupported `nturbvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt nv>
__global__ static void
_ele_colored_blusgs_sweep(const UCFDInt nstart, const UCFDInt nend, const UCFDInt var0,
                          const UCFDInt neles, const UCFDInt nvars, const UCFDInt nface,
                          const UCFDInt *__restrict__ color_order,
                          const UCFDInt *__restrict__ cell_ids,
                          const UCFDInt *__restrict__ face_neighbors,
                          const UCFDReal *__restrict__ diag ,
                          const UCFDReal *__restrict__ offdiag,
                          const UCFDReal *__restrict__ rhsb,
                          UCFDReal *__restrict__ dub)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= (nend-nstart)) return;

    const UCFDInt idx = color_order[_idx + nstart];
    const UCFDInt ridx = cell_ids[idx];
    const UCFDInt dim2 = nv*nv;
    
    UCFDInt row, col, fidx, neib;
    UCFDReal rhs[nv], dmat[dim2], val;

    for (row=0; row<nv; ++row) {
        rhs[row] = rhsb[idx + (var0+row)*neles];
        for (col=0; col<nv; ++col)
            dmat[col + row*nv] = diag[col + row*nv + ridx*dim2];
    }

    for (fidx=0; fidx<nface; ++fidx)
    {
        neib = face_neighbors[idx + fidx*neles];
        if (neib >= 0) {
            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col)
                    val += (
                        // offdiag[idx, fidx, row, col]
                        offdiag[col + row*nv + fidx*dim2 + idx*dim2*nface]
                        * dub[var0+col + neib*nvars]
                    );
                rhs[row] -= val;
            }
        }
    }
    lusub(nv, dmat, rhs);
    for (row=0; row<nv; ++row)
        dub[var0+row + ridx*nvars] = rhs[row];
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSLowerSweep(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;

    switch (sys->nfvars) {
        case 4:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<4><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 5:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<5><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 6:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<6><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSUpperSweep(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;

    switch (sys->nfvars) {
        case 4:
            for (i=ncolors-1; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<4><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 5:
            for (i=ncolors-1; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<5><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 6:
            for (i=ncolors-1; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<6><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], 0, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->diag, e->flow_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSLowerSweep(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e  = &sys->eles[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;

    switch (sys->nturbvars) {
        case 1:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<1><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], sys->nfvars, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->tdiag, e->turb_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 2:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<2><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], sys->nfvars, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->tdiag, e->turb_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSUpperSweep(PFlowSys sys, UCFDInt eidx)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    PFlowElem *e = &sys->eles[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;

    switch (sys->nturbvars) {
        case 1:
            for (i=ncolors-1; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<1><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], sys->nfvars, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->tdiag, e->turb_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        case 2:
            for (i=ncolors-1; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<2><<<bpg, TPB>>>(
                    e->color_offsets[i], e->color_offsets[i+1], sys->nfvars, e->base.neles,
                    sys->nvars, e->base.nface, e->color_order, e->base.cell_ids,
                    e->neighbors, blusgs->tdiag, e->turb_offdiag, e->base.rhs, blusgs->du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n");
            UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Reset(PFlowSys sys)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    const UCFDInt nlocal = sys->nlocal, nvars = sys->nvars;

    CUDACall(cudaMemset(blusgs->du, 0, nlocal*nvars*sizeof(UCFDReal)));
    CUDACall(cudaMemset(blusgs->dup, 0, nlocal*nvars*sizeof(UCFDReal)));

    UCFDFunctionReturn(UCFD_SUCCESS);
}

static ucfd_status_t CUBLUSGSDestroy(PFlowSys sys)
{
    CUBLUSGSSys *blusgs = (CUBLUSGSSys *)sys->data;
    CUDACall(cudaFree(blusgs->du));
    CUDACall(cudaFree(blusgs->dup));
    CUDACall(cudaFree(blusgs->diag));
    CUDACall(cudaFree(blusgs->tdiag));

    PFlowElem *eles = sys->eles;
    for (UCFDInt i=0; i<sys->nelem; ++i) {
        free(eles[i].nblocks);
        CUDACall(cudaFree(eles[i].flow_offdiag));
        CUDACall(cudaFree(eles[i].turb_offdiag));
    }

    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDPFlowSysSetCUBLUSGS(PFlowSys *sys, UCFDReal *jmat, UCFDReal *tjmat)
{
    PFlowSys s          = *sys;
    CUBLUSGSSys *blu    = (CUBLUSGSSys *)calloc(1, sizeof(*blu));
    UCFDInt nlocal      = s->nlocal;
    UCFDInt nvars       = s->nvars;
    UCFDInt nfvars      = s->nfvars;
    UCFDInt ntvars      = s->nturbvars;

    CUDACall(cudaMalloc(&blu->du, nlocal*nvars*sizeof(UCFDReal)));
    CUDACall(cudaMalloc(&blu->dup, nlocal*nvars*sizeof(UCFDReal)));
    CUDACall(cudaMalloc(&blu->diag, nlocal*nfvars*nfvars*sizeof(UCFDReal)));
    if (ntvars != 0)
        CUDACall(cudaMalloc(&blu->tdiag, nlocal*ntvars*ntvars*sizeof(UCFDReal)));
    else blu->tdiag     = NULL;

#if defined(DEBUG)
    CheckCUDAPointer(jmat);
    if (ntvars != 0)
        CheckCUDAPointer(tjmat);
#endif

    blu->jmat           = jmat;
    blu->tjmat          = tjmat;

    s->data             = blu;
    s->destroy          = CUBLUSGSDestroy;

    UCFDFunctionReturn(UCFD_SUCCESS);
}

ucfd_status_t UCFDPFlowSysSetBLUSGSElement(PFlowSys sys, UCFDInt eidx,
                                           UCFDInt neles, UCFDInt nface, UCFDInt *cell_ids,
                                           UCFDReal *uptsb, UCFDReal *rhs,
                                           UCFDReal *dt, UCFDReal *dsrc,
                                           UCFDReal *vol, UCFDReal *resid_out,
                                           UCFDInt rank_ncolors, UCFDInt *face_refs, UCFDReal *face_factors,
                                           UCFDInt *color_offsets, UCFDInt *color_order,
                                           UCFDInt *neighbors)
{
    PFlowElem *eles = sys->eles;
    PFlowElem *e = &eles[eidx];
    const UCFDInt nfvars = sys->nfvars, ntvars = sys->nturbvars;

    ((FlowElem *)e)->neles          = neles;
    ((FlowElem *)e)->nface          = nface;
    ((FlowElem *)e)->cell_ids       = cell_ids;
    ((FlowElem *)e)->uptsb          = uptsb;
    ((FlowElem *)e)->rhs            = rhs;
    ((FlowElem *)e)->dt             = dt;
    ((FlowElem *)e)->dsrc           = dsrc;
    ((FlowElem *)e)->vol            = vol;
    ((FlowElem *)e)->resid_out      = resid_out;

    /* Get number of private colors for each element */
    UCFDInt ncolors = 0;
    for (UCFDInt i=0; i<rank_ncolors; ++i)
        if (color_offsets[i] != color_offsets[i+1]) ncolors++;
    
    /* Construct nblocks */
    e->nblocks = (UCFDInt *)malloc((size_t)ncolors * sizeof(*e->nblocks));
    if (ncolors != 0)
        UCFDCheckNull(e->nblocks, "Color-block allocation failed\n");
    for (UCFDInt i=0; i<ncolors; ++i)
        e->nblocks[i] = color_offsets[i+1] - color_offsets[i];

    e->ncolors                      = ncolors;
    e->face_refs                    = face_refs;
    e->face_factors                 = face_factors;
    e->color_offsets                = color_offsets;
    e->color_order                  = color_order;

    /* Colored BLU-SGS attributes */
    e->neighbors                    = neighbors;
    CUDACall(cudaMalloc(&e->flow_offdiag, nfvars*nfvars*nface*neles*sizeof(UCFDReal)));
    if (ntvars != 0)
        CUDACall(cudaMalloc(&e->turb_offdiag, ntvars*ntvars*nface*neles*sizeof(UCFDReal)));
    else e->turb_offdiag            = NULL;

    /* Colored LU-SGS attributes */
    e->lower_neighbors              = NULL;
    e->upper_neighbors              = NULL;
    e->face_normal                  = NULL;
    e->wave_factors                 = NULL;

    UCFDFunctionReturn(UCFD_SUCCESS);
}
