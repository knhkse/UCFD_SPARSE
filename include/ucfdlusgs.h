#pragma once

#include "ucfdtypes.h"

typedef struct _FlowSys *FlowSys;
typedef struct _PFlowSys *PFlowSys;

#if defined(__cplusplus)
extern "C" {
#endif

/* Common functions */
UCFD_EXTERN ucfd_status_t UCFDFlowSysCreate(FlowSys*, UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt*, UCFDInt*, UCFDInt*, UCFDInt8*,
                                            UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetElement(FlowSys, UCFDInt, UCFDInt, UCFDInt, UCFDInt*, UCFDReal*,
                                                UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysDestroy(FlowSys*);

UCFD_EXTERN ucfd_status_t UCFDPFlowSysCreate(PFlowSys*,
                                             UCFDInt, UCFDInt, UCFDInt,
                                             UCFDInt, UCFDInt, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDPFlowSysSetBLUSGSElement(PFlowSys, UCFDInt,
                                                       UCFDInt, UCFDInt, UCFDInt*,
                                                       UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*,
                                                       UCFDReal*, UCFDReal*, UCFDInt, UCFDInt*, UCFDReal*,
                                                       UCFDInt*, UCFDInt*, UCFDInt*);
UCFD_EXTERN ucfd_status_t UCFDPFlowSysDestroy(PFlowSys*);


/* LU-SGS functions */
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetLUSGS(FlowSys*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_Pack(FlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_Update(FlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSPrepare(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSPrepare(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSLowerSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSUpperSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSLowerSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSUpperSweep(FlowSys, UCFDReal);

/* BLU-SGS functions */
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetBLUSGS(FlowSys*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Pack(FlowSys, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_KWSST_Pack(FlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_SA_Pack(FlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Update(FlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_SubResidual(FlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSPrepare(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSPrepare(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSLowerSweep(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSLowerSweep(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSUpperSweep(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSUpperSweep(FlowSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Reset(FlowSys);

/* CUDA functions */
#if defined(__CUDACC__)
UCFD_EXTERN ucfd_status_t UCFDPFlowSysSetCUBLUSGS(PFlowSys*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Pack(PFlowSys, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_KWSST_Pack(PFlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_SA_Pack(PFlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Update(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_SubResidual(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSPrepare(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSPrepare(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSLowerSweep(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSUpperSweep(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSLowerSweep(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSUpperSweep(PFlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Reset(PFlowSys);
#endif

#if defined(__cplusplus)
}
#endif
