#pragma once

#include <float.h>
#include "ucfdtypes.h"

typedef struct _UCFDSys *UCFDSys;

/* Matrix inverse parameter */
#if !defined(PIVOT_TOL)
    #define PIVOT_TOL DBL_EPSILON
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Common functions
 * - Create/set system and set elements for each method
 * - SysCreate
 *      - Basic 
 */
UCFD_EXTERN ucfd_status_t UCFDSysCreate(UCFDSys*, UCFDInt, UCFDInt, UCFDInt,
                                        UCFDInt, UCFDInt, UCFDInt, UCFDInt*,
                                        UCFDInt*, UCFDInt*, UCFDInt8*, UCFDReal*,
                                        UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDColoredSysCreate(UCFDSys*, UCFDInt, UCFDInt, UCFDInt,
                                               UCFDInt, UCFDInt, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDSysDestroy(UCFDSys*);

/* LU-SGS */
UCFD_EXTERN ucfd_status_t UCFDLUSGS_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_Update(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSPrepare(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSPrepare(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSLowerSweep(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSLowerSweep(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_NSUpperSweep(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDLUSGS_RANSUpperSweep(UCFDSys, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDSysSetLUSGS(UCFDSys*, UCFDReal*, UCFDReal*, UCFDReal*,
                                          UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDEleSetLUSGS(UCFDSys*, UCFDInt, UCFDInt,
                                          UCFDInt, UCFDInt*, UCFDReal*,
                                          UCFDReal*, UCFDReal*, UCFDReal*);

/* Block LU-SGS */
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Pack(UCFDSys, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_KWSST_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_SA_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Update(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_SubResidual(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSPrepare(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSPrepare(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSLowerSweep(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSLowerSweep(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_NSUpperSweep(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_RANSUpperSweep(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDBLUSGS_Reset(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDSysSetBLUSGS(UCFDSys*, UCFDReal*, UCFDReal*, UCFDReal*,
                                           UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDEleSetBLUSGS(UCFDSys*, UCFDInt, UCFDInt, UCFDInt, UCFDInt*,
                                           UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*,
                                           UCFDReal*, UCFDReal*);

/* Colored BLU-SGS */
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_Pack(UCFDSys, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_KWSST_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_SA_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_Update(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_SubResidual(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_NSPrepare(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_RANSPrepare(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_NSLowerSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_RANSLowerSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_NSUpperSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_RANSUpperSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCBLUSGS_Reset(UCFDSys);
UCFD_EXTERN ucfd_status_t UCFDColoredSysSetBLUSGS(UCFDSys*, UCFDReal*, UCFDReal*,
                                                  UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);
UCFD_EXTERN ucfd_status_t UCFDColoredEleSetBLUSGS(UCFDSys*, UCFDInt, UCFDInt, UCFDInt,
                                                  UCFDInt, UCFDInt, UCFDInt*, UCFDInt*,
                                                  UCFDInt*, UCFDInt*, UCFDInt*, UCFDReal*,
                                                  UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*,
                                                  UCFDReal*, UCFDReal*, UCFDReal*, UCFDReal*);

/* CUDA BLU-SGS */
#if defined(__CUDACC__)
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Pack(UCFDSys, UCFDInt, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_KWSST_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_SA_Pack(UCFDSys, UCFDInt, UCFDReal, UCFDReal);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Update(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_SubResidual(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSPrepare(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSPrepare(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSLowerSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_NSUpperSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSLowerSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_RANSUpperSweep(UCFDSys, UCFDInt);
UCFD_EXTERN ucfd_status_t UCFDCUBLUSGS_Reset(UCFDSys);
#endif

#if defined(__cplusplus)
}
#endif
