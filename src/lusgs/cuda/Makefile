# LU-SGS cuda library Makefile

# =================================
# Configuration Variables
# =================================
# Default values
NVCC ?= nvcc
INT64 ?= 0
FLOAT32 ?= 0

# Get configuration from `Makefile.inc`
include $(UCFD_PATH)/Makefile.inc

# =================================
# Derived Settings
# =================================
# Integer type (Default is int32)
ifeq ($(INT64), 1)
  CUFLAGS += -DUCFD_INT64
  $(info 64-bit integer type enabled)
endif

# Real type (Default is float64)
ifeq ($(FLOAT32), 1)
  CUFLAGS += -DUCFD_FLOAT32
  $(info 32-bit float type enabled)
endif

# Basic flow variables
ifeq ($(NVARS),)
  $(error `NVARS` must be defined)
else
  CUFLAGS += -DNVARS=$(NVARS)
endif

ifeq ($(NFVARS),)
  $(error `NFVARS` must be defined)
else
  CUFLAGS += -DNFVARS=$(NFVARS)
endif

ifeq ($(NTURBVARS),)
  $(info `NTURBVARS` is not defined, Auto setting...)
  TEMP = $(shell expr $(NVARS) - $(NFVARS))
  ifeq ($(TEMP), 0)
    $(info N-S equations mode; Do not use rans-prefixed functions)
    NTURBVARS = 1
  else
    NTURBVARS = $TEMP
  endif
  CUFLAGS += -DNTURBVARS=$(NTURBVARS)
else ifeq ($(filter $(NTURBVARS), 0 1 2),)
  $(error `NTURBVARS` must be defined as 0, 1 or 2)
else
  CUFLAGS += -DNTURBVARS=$(NTURBVARS)
endif

ifeq ($(NDIMS),)
  $(error `NDIMS` must be defined)
else
  CUFLAGS += -DNDIMS=$(NDIMS)
endif

# Constants
ifneq ($(GAMMA),)
  CUFLAGS += -DGAMMA=$(GAMMA)
endif

ifneq ($(PMIN),)
  CUFLAGS += -DPMIN=$(PMIN)
endif

ifneq ($(BETAST),)
  CUFLAGS += -DBETAST=$(BETAST)
endif

ifneq ($(TPB),)
  CUFLAGS += -DTPB=$(TPB)
else
  $(info Default threads per block=128)
  CUFLAGS += -DTPB=128
endif

# =================================
# Build Targets
# =================================
SLIB = liblusgs_cuda.a
DLIB = liblusgs_cuda.so

OBJDIR = ./obj
INCDIR = ./include
IFLAGS += -I$(INCDIR)

CLUSGS = $(OBJDIR)/coloredlusgs_cuda.o
FLUX = $(OBJDIR)/flux_cuda.o


# Build object files
$(FLUX) : flux_cuda.cu $(INCDIR)/flux_cuda.cuh $(CONFIG) $(TYPES)
	$(NVCC) $(CUFLAGS) -rdc=true -Xcompiler -fPIC $(IFLAGS) -c $< -o $@

$(CLUSGS) : coloredlusgs_cuda.cu $(INCDIR)/coloredlusgs_cuda.cuh $(INCDIR)/flux_cuda.cuh $(CONFIG) $(TYPES)
	$(NVCC) $(CUFLAGS) -rdc=true -Xcompiler -fPIC $(IFLAGS) -c $< -o $@

# Build library files
$(SLIB) : $(CLUSGS) $(FLUX)
	ar rcs $@ $^

$(DLIB) : $(CLUSGS) $(FLUX)
	$(NVCC) -rdc=true -Xcompiler -fPIC -shared $^ -o $@

all : $(SLIB) $(DLIB)
	@echo "All cuda libraries are created"

static : $(SLIB)
	@echo "Static cuda library is created"

dynamic : $(DLIB)
	@echo "Dynamic cuda library is created"

.PHONY : clean
clean :
	rm -rf $(OBJDIR)
	rm -f $(SLIB) $(DLIB)