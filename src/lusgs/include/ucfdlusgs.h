#pragma once

#include "ucfdtypes.h"

typedef struct _FlowSys *FlowSys;

/* Common functions */
UCFD_EXTERN ucfd_status_t UCFDFlowSysCreate(FlowSys*, UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt*, UCFDInt*, UCFDInt*, UCFDInt8*,
                                            UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetElement(FlowSys*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*, UCFDReal*,
                                                UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysDestroy(FlowSys*);

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

/* Test */
// UCFD_EXTERN void Exportarr(FlowSys, UCFDInt, UCFDReal*);
// UCFD_EXTERN void Exportdiag(FlowSys, UCFDInt, UCFDReal*);
