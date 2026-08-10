#pragma once

#include "ucfdtypes.h"

typedef struct _FlowSys *FlowSys;

UCFD_EXTERN ucfd_status_t UCFDFlowSysCreate(FlowSys*, UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt, UCFDInt, UCFDInt,
                                            UCFDInt*, UCFDInt*, UCFDInt*, UCFDInt8*,
                                            UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetRANS(FlowSys*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysSetElement(FlowSys*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*,
                                                UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDFlowSysDestroy(FlowSys*);

UCFD_EXTERN ucfd_status_t UCFDLUSGS_Pack(FlowSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_Update(FlowSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSPrepare(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSPrepare(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSLowerSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSUpperSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSLowerSweep(FlowSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSUpperSweep(FlowSys, UCFDReal);

UCFD_EXTERN void TestPack(FlowSys, UCFDInt);
UCFD_EXTERN void Exportranku(FlowSys, UCFDReal*);
UCFD_EXTERN void Exportrankdu(FlowSys, UCFDReal*);
UCFD_EXTERN void Exportrankdiag(FlowSys, UCFDReal*);
