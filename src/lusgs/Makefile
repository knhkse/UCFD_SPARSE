# LU-SGS library Makefile

# =================================
# Configuration Variables
# =================================
# Default values
CC ?= gcc
INT64 ?= 0
FLOAT32 ?= 0

# Get configuration from `Makefile.inc`
include $(UCFD_PATH)/Makefile.inc

# =================================
# Derived Settings
# =================================
# Integer type (Default is int32)
ifeq ($(INT64), 1)
  CFLAGS += -DUCFD_INT64
  $(info 64-bit integer type enabled)
endif

# Real type (Default is float64)
ifeq ($(FLOAT32), 1)
  CFLAGS += -DUCFD_FLOAT32
  $(info 32-bit float type enabled)
endif

# Basic flow variables
ifeq ($(NVARS),)
  $(error `NVARS` must be defined)
else
  CFLAGS += -DNVARS=$(NVARS)
endif

ifeq ($(NFVARS),)
  $(error `NFVARS` must be defined)
else
  CFLAGS += -DNFVARS=$(NFVARS)
endif

ifeq ($(NDIMS),)
  $(error `NDIMS` must be defined)
else
  CFLAGS += -DNDIMS=$(NDIMS)
endif

ifeq ($(NTURBVARS),)
  $(info `NTURBVARS` is not defined, default value is 1.)
  $(info Do not rans_* functions; it may expect incorrect results.)
  CFLAGS += -DNTURBVARS=1
else ifeq ($(filter $(NTURBVARS), 1 2),)
  $(error `NTURBVARS` must be defined as 1 or 2)
else
  CFLAGS += -DNTURBVARS=$(NTURBVARS)
endif

# Constants
ifneq ($(GAMMA),)
  CFLAGS += -DGAMMA=$(GAMMA)
endif

ifneq ($(PMIN),)
  CFLAGS += -DPMIN=$(PMIN)
endif

ifneq ($(BETAST),)
  CFLAGS += -DBETAST=$(BETAST)
endif


# =================================
# Build Targets
# =================================
SLIB = liblusgs.a
SCLIB = libclusgs.a
DLIB = liblusgs.so
DCLIB = libclusgs.so

OBJDIR = ./obj
INCDIR = ./include
IFLAGS += -I$(INCDIR)

LUSGS = $(OBJDIR)/lusgs.o
BLUSGS = $(OBJDIR)/blusgs.o
CLUSGS = $(OBJDIR)/coloredlusgs.o
CBLUSGS = $(OBJDIR)/coloredblusgs.o
FLUX = $(OBJDIR)/flux.o
INVERSE = $(OBJDIR)/inverse.o

# Build object(.o) files
$(FLUX) : flux.c $(INCDIR)/flux.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(INVERSE) : inverse.c $(INCDIR)/inverse.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(LUSGS) : lusgs.c $(INCDIR)/lusgs.h $(INCDIR)/flux.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(CLUSGS) : coloredlusgs.c $(INCDIR)/coloredlusgs.h $(INCDIR)/flux.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(BLUSGS) : blusgs.c $(INCDIR)/blusgs.h $(INCDIR)/flux.h $(INVERSEH) $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(CBLUSGS) : coloredblusgs.c $(INCDIR)/coloredblusgs.h $(INCDIR)/flux.h $(INVERSEH) $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@


# Build library files
$(SLIB) : $(LUSGS) $(BLUSGS) $(FLUX) $(INVERSE)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(SCLIB) : $(CLUSGS) $(CBLUSGS) $(FLUX) $(INVERSE)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(LUSGS) $(BLUSGS) $(FLUX) $(INVERSE)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DCLIB) : $(CLUSGS) $(CBLUSGS) $(FLUX) $(INVERSE)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib


all : $(SLIB) $(SCLIB) $(DLIB) $(DCLIB)
	@echo "All libraries are created"

static : $(SLIB) $(SCLIB)
	@echo "Static libraries are created"

dynamic : $(DLIB) $(DCLIB)
	@echo "Dynamic libraries are created"

.PHONY : clean
clean :
	rm -rf $(OBJDIR)
	rm -f $(SLIB) $(SCLIB) $(DLIB) $(DCLIB)
