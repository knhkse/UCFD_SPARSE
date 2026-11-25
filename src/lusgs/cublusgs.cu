#include "flowsys.h"
#include "flux.h"
#include "smallmat.cuh"
#include "source_jacobian.cuh"


template<UCFDInt nfvars>
__global__ static void
_rank_blusgs_diag_pack(const UCFDInt neles, const UCFDInt nlocal,
                       const UCFDReal a0,
                       const UCFDInt *__restrict__ cell_ids,
                       const UCFDReal *__restrict__ dt,
                       UCFDReal *__restrict__ diag)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;
    
    const UCFDInt ridx = cell_ids[idx];
    UCFDInt row, col;

    #pragma unroll
    for (row=0; row<nfvars; ++row) {
        #pragma unroll
        for (col=0; col<nfvars; ++col)
            // diag[row, col, ridx] = 0.0;
            diag[(size_t)(row*nfvars + col)*nlocal + ridx] = 0.0;
        diag[(size_t)(row*nfvars + row)*nlocal + ridx] = 1/dt[idx] + a0;
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nfvars = sys->nfvars;

    switch (nfvars) {
        case 4:
            _rank_blusgs_diag_pack<4><<<bpg, TPB>>>(
                e->neles, sys->nlocal, a0,
                e->cell_ids, e->dt, blusgs->flow_diag
            ); break;
        case 5:
            _rank_blusgs_diag_pack<5><<<bpg, TPB>>>(
                e->neles, sys->nlocal, a0,
                e->cell_ids, e->dt, blusgs->flow_diag
            ); break;
        case 6:
            _rank_blusgs_diag_pack<6><<<bpg, TPB>>>(
                e->neles, sys->nlocal, a0,
                e->cell_ids, e->dt, blusgs->flow_diag
            ); break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt nvars, UCFDInt nturbvars, rans_t rans_type>
__global__ static void
_rank_tblusgs_pack(const UCFDInt neles, const UCFDInt nlocal,
                   const UCFDReal factor, const UCFDReal a0,
                   const UCFDInt *__restrict__ cell_ids,
                   const UCFDReal *__restrict__ upts,
                   const UCFDReal *__restrict__ dt,
                   const UCFDReal *__restrict__ dsrc,
                   UCFDReal *__restrict__ tdiag)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    const UCFDInt ridx = cell_ids[idx];
    UCFDInt kdx, row, col;
    UCFDReal u[nvars], ds[nvars], dmat[nturbvars*nturbvars];

    #pragma unroll
    for (kdx=0; kdx<nvars; ++kdx) {
        u[kdx] = upts[idx + kdx*neles];
        ds[kdx] = dsrc[idx + kdx*neles];
    }

    /* Source Jacobian contributes to the turbulence diagonal block. */
    #pragma unroll
    for (row=0; row<nturbvars; ++row) {
        #pragma unroll
        for (col=0; col<nturbvars; ++col)
            dmat[col + nturbvars*row] = 0.0;
    }

    if (rans_type == RANS_KWSST)
        cuda_kwsst_src_jacobian<nvars>(u, dmat, ds);
    else if (rans_type == RANS_SA)
        cuda_sa_src_jacobian<nvars>(dmat, ds);

    #pragma unroll
    for (row=0; row<nturbvars; ++row)
        #pragma unroll
        for (col=0; col<nturbvars; ++col)
            tdiag[(size_t)(row*nturbvars + col)*nlocal + ridx] = dmat[col + nturbvars*row];
        tdiag[(size_t)(row*nturbvars + row)*nlocal + ridx] += 1.0/(dt[idx]*factor) + a0;
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_KWSST_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nvars = sys->nvars;

    switch (nvars) {
        case 6:
            _rank_tblusgs_pack<6, 2, RANS_KWSST><<<bpg, TPB>>>(
            e->neles, sys->nlocal, turb_factor, a0,
            e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
        ); break;
        case 7:
            _rank_tblusgs_pack<7, 2, RANS_KWSST><<<bpg, TPB>>>(
            e->neles, sys->nlocal, turb_factor, a0,
            e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
        ); break;
        default:
            fprintf(stderr, "Unsupported `nvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_SA_Pack(UCFDSys sys, UCFDInt eidx, UCFDReal turb_factor, UCFDReal a0)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nvars = sys->nvars;

    switch (nvars) {
        case 5:
            _rank_tblusgs_pack<5, 1, RANS_SA><<<bpg, TPB>>>(
            e->neles, sys->nlocal, turb_factor, a0,
            e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
        ); break;
        case 6:
            _rank_tblusgs_pack<6, 1, RANS_SA><<<bpg, TPB>>>(
            e->neles, sys->nlocal, turb_factor, a0,
            e->cell_ids, e->uptsb, e->dsrc, e->dt, blusgs->turb_diag
        ); break;
        default: fprintf(stderr, "Unsupported `nvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

template<UCFDInt nvars>
__global__ static void
_rank_blusgs_update(const UCFDInt neles, const UCFDInt nlocal,
                    const UCFDInt *__restrict__ cell_ids,
                    const UCFDReal *__restrict__ du,
                    UCFDReal *__restrict__ upts)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    UCFDInt kdx;
    const UCFDInt ridx = cell_ids[idx];

    #pragma unroll
    for (kdx=0; kdx<nvars; ++kdx)
        upts[idx + kdx*neles] += du[kdx*nlocal + ridx];
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Update(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nvars = sys->nvars;

    switch (nvars) {
        case 4:
            _rank_blusgs_update<4><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids, blusgs->rank_du, e->uptsb
            ); break;
        case 5:
            _rank_blusgs_update<5><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids, blusgs->rank_du, e->uptsb
            ); break;
        case 6:
            _rank_blusgs_update<6><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids, blusgs->rank_du, e->uptsb
            ); break;
        case 7:
            _rank_blusgs_update<7><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids, blusgs->rank_du, e->uptsb
            ); break;
        default: fprintf(stderr, "Unsupported `nvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

template<UCFDInt nvars>
__global__ static void
_rank_sub_residual(const UCFDInt neles, const UCFDInt nlocal,
                   const UCFDInt *__restrict__ cell_ids,
                   const UCFDReal *__restrict__ vol,
                   const UCFDReal *__restrict__ du,
                   UCFDReal *__restrict__ dup,
                   UCFDReal *__restrict__ res)
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    UCFDInt kdx;
    UCFDReal diff;
    const UCFDInt ridx = cell_ids[idx];

    #pragma unroll
    for (kdx=0; kdx<nvars; ++kdx) {
        size_t offset = ridx + kdx*nlocal;
        diff = du[offset] - dup[offset];
        res[(size_t)kdx*nlocal + idx] = diff*diff*vol[idx];
        dup[offset] = du[offset];
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_SubResidual(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nvars = sys->nvars;

    switch (nvars) {
        case 4:
            _rank_sub_residual<4><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids,
                e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
            ); break;
        case 5:
            _rank_sub_residual<5><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids,
                e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
            ); break;
        case 6:
            _rank_sub_residual<6><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids,
                e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
            ); break;
        case 7:
            _rank_sub_residual<7><<<bpg, TPB>>>(
                e->neles, sys->nlocal, e->cell_ids,
                e->vol, blusgs->rank_du, blusgs->rank_dup, e->resid_out
            ); break;
        default: fprintf(stderr, "Unsupported `nvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


/**
 * Prepare & sweep functions for colored BLU-SGS
 * nv : nfvars or nturbvars
 * var0 : 0 or nfvars
 */

template<UCFDInt nv>
__global__ static void
_ele_pre_blusgs(const UCFDInt neles, const UCFDInt nface,
                const UCFDInt nlocal, const UCFDInt rank_nfaces,
                const UCFDInt *__restrict__ cell_ids,
                const UCFDInt *__restrict__ face_refs,
                const UCFDInt *__restrict__ neighbors,
                const UCFDReal *__restrict__ face_factors,
                UCFDReal *__restrict__ diag,
                UCFDReal *__restrict__ offdiag,
                const UCFDReal *__restrict__ jmat
                )
{
    const UCFDInt idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= neles) return;

    const UCFDInt ridx = cell_ids[idx];
    const UCFDInt dim2 = nv*nv;
    UCFDInt fidx, row, col, face, slot;
    size_t ind0, ind1;
    UCFDReal dmat[dim2], fv, dval, oval;

    #pragma unroll
    for (row=0; row<nv; ++row) {
        #pragma unroll
        for (col=0; col<nv; ++col)
            dmat[col + row*nv] = diag[(size_t)(row*nv + col)*nlocal + ridx];
    }

    for (fidx=0; fidx<nface; ++fidx) {
        face = face_refs[idx + fidx*neles];
        UCFDBool has_neighbor = neighbors[idx + fidx*neles] >= 0;
        if (face > 0) slot = face - 1;
        else slot = -face - 1;
        fv = face_factors[idx + fidx*neles];

        for (row=0; row<nv; ++row) {
            for (col=0; col<nv; ++col) {
                ind0 = slot + (size_t)col*rank_nfaces + (size_t)row*nv*rank_nfaces;
                ind1 = ind0 + (size_t)dim2*rank_nfaces;
                if (face > 0) dval = jmat[ind0];
                else dval = -jmat[ind1];
                dmat[col + row*nv] += dval*fv;

                /* Colored sweeps use off-diagonal blocks */
                if (has_neighbor) {
                    if (face > 0) oval = jmat[ind1];
                    else oval = -jmat[ind0];
                    offdiag[idx + fidx*neles + (size_t)col*nface*neles + (size_t)row*nv*nface*neles] = oval*fv;
                }
            }
        }
    }

    pivoted_inverse<nv>(dmat);

    #pragma unroll
    for (row=0; row<nv; ++row) {
        #pragma unroll
        for (col=0; col<nv; ++col)
            diag[(size_t)(row*nv + col)*nlocal + ridx] = dmat[col + row*nv];
    }
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSPrepare(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nfvars = sys->nfvars;
    
    switch (nfvars) {
        case 4:
            _ele_pre_blusgs<4><<<bpg, TPB>>>(
                e->neles, e->nface, sys->nlocal, sys->nfaces,
                e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
                blusgs->flow_diag, e->flow_offdiag, blusgs->jmat
            ); break;
        case 5:
            _ele_pre_blusgs<5><<<bpg, TPB>>>(
                e->neles, e->nface, sys->nlocal, sys->nfaces,
                e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
                blusgs->flow_diag, e->flow_offdiag, blusgs->jmat
            ); break;
        case 6:
            _ele_pre_blusgs<6><<<bpg, TPB>>>(
                e->neles, e->nface, sys->nlocal, sys->nfaces,
                e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
                blusgs->flow_diag, e->flow_offdiag, blusgs->jmat
            ); break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSPrepare(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    const UCFDInt bpg = (e->neles + TPB - 1)/TPB;
    const UCFDInt nturbvars = sys->nturbvars;
    
    switch (nturbvars) {
        case 1:
            _ele_pre_blusgs<1><<<bpg, TPB>>>(
                e->neles, e->nface, sys->nlocal, sys->nfaces,
                e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
                blusgs->turb_diag, e->turb_offdiag, blusgs->tjmat
            ); break;
        case 2:
            _ele_pre_blusgs<2><<<bpg, TPB>>>(
                e->neles, e->nface, sys->nlocal, sys->nfaces,
                e->cell_ids, e->face_refs, e->neighbors, e->face_factors,
                blusgs->turb_diag, e->turb_offdiag, blusgs->tjmat
            ); break;
        default:
            fprintf(stderr, "Unsupported `nturbvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}


template<UCFDInt nv>
__global__ static void
_ele_colored_blusgs_sweep(const UCFDInt i_begin, const UCFDInt i_end,
                          const UCFDInt neles, const UCFDInt nface, const UCFDInt var0,
                          const UCFDInt nlocal, const UCFDInt cell_offset,
                          const UCFDInt *__restrict__ neighbors,
                          const UCFDReal *__restrict__ diag ,
                          const UCFDReal *__restrict__ offdiag,
                          const UCFDReal *__restrict__ rhsb,
                          UCFDReal *__restrict__ dub)
{
    const UCFDInt _idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (_idx >= (i_end-i_begin)) return;

    const UCFDInt idx = _idx + i_begin;
    const UCFDInt ridx = cell_offset + idx;
    const UCFDInt dim2 = nv*nv;
    
    UCFDInt row, col, fidx, neib;
    UCFDReal rhs[nv], dmat[dim2], val;

    #pragma unroll
    for (row=0; row<nv; ++row) {
        rhs[row] = rhsb[idx + (size_t)(var0+row)*neles];
        #pragma unroll
        for (col=0; col<nv; ++col)
            dmat[col + row*nv] = diag[(size_t)(row*nv + col)*nlocal + ridx];
    }

    for (fidx=0; fidx<nface; ++fidx) {
        neib = neighbors[idx + fidx*neles];
        if (neib >= 0) {
            for (row=0; row<nv; ++row) {
                val = 0.0;
                for (col=0; col<nv; ++col)
                    val += (
                        // offdiag[row, col, fidx, idx]
                        offdiag[idx + fidx*neles + (size_t)col*nface*neles + (size_t)row*nv*nface*neles]
                        * dub[neib + (var0+col)*nlocal]
                    );
                rhs[row] -= val;
            }
        }
    }
    block_matvec<nv>(dmat, rhs);

    #pragma unroll
    for (row=0; row<nv; ++row)
        dub[ridx + (var0+row)*nlocal] = rhs[row];
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSLowerSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;
    const UCFDInt nfvars = sys->nfvars;

    switch (nfvars) {
        case 4:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<4><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 5:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<5><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 6:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<6><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_NSUpperSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;
    const UCFDInt nfvars = sys->nfvars;

    switch (nfvars) {
        case 4:
            for (i=ncolors-2; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<4><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 5:
            for (i=ncolors-2; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<5><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 6:
            for (i=ncolors-2; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<6><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, 0, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->flow_diag, e->flow_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSLowerSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;
    const UCFDInt nturbvars = sys->nturbvars;

    switch (nturbvars) {
        case 1:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<1><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, sys->nfvars, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 2:
            for (i=0; i<ncolors; ++i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<2><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, sys->nfvars, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_RANSUpperSweep(UCFDSys sys, UCFDInt eidx)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    PFlowEle *e = &((PFlowEle *)sys->eles)[eidx];
    UCFDInt i, bpg;
    const UCFDInt ncolors = e->ncolors;
    const UCFDInt *nblocks = e->nblocks;
    const UCFDInt nturbvars = sys->nturbvars;

    switch (nturbvars) {
        case 1:
            for (i=ncolors-2; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<1><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, sys->nfvars, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        case 2:
            for (i=ncolors-2; i>=0; --i) {
                bpg = (nblocks[i] + TPB - 1)/TPB;
                _ele_colored_blusgs_sweep<2><<<bpg, TPB>>>(
                    e->color_offset[i], e->color_offset[i+1], e->neles,
                    e->nface, sys->nfvars, sys->nlocal, e->cell_offset,
                    e->neighbors, blusgs->turb_diag, e->turb_offdiag, e->rhs, blusgs->rank_du
                );
            }
            break;
        default:
            fprintf(stderr, "Unsupported `nfvars` size\n"); UCFDFunctionReturn(UCFD_FAILED);
    }
    UCFDFunctionReturn(UCFD_SUCCESS);
}

extern "C" ucfd_status_t
UCFDCUBLUSGS_Reset(UCFDSys sys)
{
    PFlowSys *blusgs = (PFlowSys *)sys->flowsys;
    const UCFDInt nlocal = sys->nlocal, nvars = sys->nvars;

    CUDACall(cudaMemset(blusgs->rank_du, 0, nlocal*nvars*sizeof(UCFDReal)));
    CUDACall(cudaMemset(blusgs->rank_dup, 0, nlocal*nvars*sizeof(UCFDReal)));

    UCFDFunctionReturn(UCFD_SUCCESS);
}
