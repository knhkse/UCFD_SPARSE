# LU-SGS type methods library
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

# Resolve paths from this Makefile, including when invoked with make -f.
THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCDIR := $(patsubst %/,%,$(dir $(THIS_MAKEFILE)))
INCDIR := $(UCFDPATH)/include

# Get configuration from `Makefile.inc`
include $(UCFDPATH)/Makefile.inc

# Default compiler - Recommend Intel compiler
CC ?= mpiicx
ifeq ($(origin CXX), default)
CXX := mpicxx
endif
CXX ?= mpicxx
NVCC ?= nvcc
OPTFLAGS ?= -O3 -ffast-math
NVCC_OPTFLAGS ?= -O3 --use_fast_math
BUILD_DIR ?= $(abspath $(SRCDIR)/../../build)
CFLAGS := -Wall $(OPTFLAGS) -fPIC -fvisibility=hidden
CUDA_HOME ?= $(shell nvcc_path=$$(command -v $(NVCC)); if [ -n "$$nvcc_path" ]; then dirname "$$(dirname "$$(readlink -f "$$nvcc_path")")"; fi)
CUDA_ARCH ?= sm_80
NVCCFLAGS ?= $(NVCC_OPTFLAGS) -std=c++14 -arch=$(CUDA_ARCH)
CUDA_LIBDIR ?= $(CUDA_HOME)/lib64
CUDA_LDLIBS ?= -L"$(CUDA_LIBDIR)" -lcudart

# Multi-threading flag
ifeq ($(USE_OMP), 1)
	ifneq ($(findstring icx,$(CC))$(findstring icc,$(CC)),)
		OMPFLAG := -qopenmp
	else
		OMPFLAG := -fopenmp
	endif
	CFLAGS += $(OMPFLAG)
	MACROS += -DUSE_OMP
endif

# Additional flags
ifeq ($(USE_INT64), 1)
	MACROS += -DUCFD_INT64
endif

ifeq ($(USE_FLOAT32), 1)
	MACROS += -DUCFD_FLOAT32
endif

BUILDDIR := $(BUILD_DIR)/lusgs
CPPFLAGS := $(MACROS) -I$(INCDIR)
C_SRCS := $(sort $(shell find $(SRCDIR) -type f -name '*.c'))
CU_SRCS := $(sort $(shell find $(SRCDIR) -type f -name '*.cu' ! -name 'test*.cu'))

# Replace the absolute source prefix instead of prepending to it.
C_OBJS := $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(C_SRCS))
CU_OBJS := $(patsubst $(SRCDIR)/%.cu,$(BUILDDIR)/%.cu.o,$(CU_SRCS))
OBJS := $(C_OBJS) $(CU_OBJS)
DEPS := $(OBJS:.o=.d)
SLIB := $(BUILDDIR)/libucfdlusgs.a
DLIB := $(BUILDDIR)/libucfdlusgs.so

# The source's directory supplies its individual include path.
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(THIS_MAKEFILE)
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) $(CPPFLAGS) -I"$(dir $<)" -MMD -MP -c "$<" -o "$@"

# Keep CUDA objects distinct from C objects with the same source basename.
$(BUILDDIR)/%.cu.o: $(SRCDIR)/%.cu $(THIS_MAKEFILE)
	@mkdir -p "$(@D)"
	$(NVCC) $(NVCCFLAGS) -allow-unsupported-compiler -rdc=true -dlto -ccbin $(CXX) \
		-Xcompiler=-Wall,-Wextra,-fPIC,-fvisibility=hidden $(if $(OMPFLAG),-Xcompiler=$(OMPFLAG)) \
		$(CPPFLAGS) -I"$(dir $<)" -MMD -MP -c "$<" -o "$@"

CUDA_DLINK_OBJ := $(BUILDDIR)/cuda_dlink.o

$(CUDA_DLINK_OBJ): $(CU_OBJS) $(THIS_MAKEFILE)
	@mkdir -p "$(@D)"
	$(NVCC) $(NVCCFLAGS) -allow-unsupported-compiler -ccbin $(CXX) \
		-dlink -dlto -Xcompiler=-fPIC $(CU_OBJS) -o "$@"

$(SLIB): $(OBJS) $(THIS_MAKEFILE)
	$(RM) "$@"
	$(AR) rcs "$@" $(OBJS)
	@mkdir -p "$(UCFDPATH)/lib"
	cp $(SLIB) $(UCFDPATH)/lib/$(notdir $(SLIB))

$(DLIB): $(OBJS) $(CUDA_DLINK_OBJ) $(THIS_MAKEFILE)
	$(CXX) $(LDFLAGS) $(OMPFLAG) -shared -o "$@" \
		$(OBJS) $(CUDA_DLINK_OBJ) $(CUDA_LDLIBS) $(LDLIBS)
	@mkdir -p "$(UCFDPATH)/lib"
	cp $(DLIB) $(UCFDPATH)/lib/$(notdir $(DLIB))


.PHONY: all static dynamic clean
all : $(SLIB) $(DLIB)

static : $(SLIB)

dynamic: $(DLIB)

clean:
	$(RM) $(OBJS) $(DEPS) $(SLIB) $(DLIB) $(CUDA_DLINK_OBJ)

# auto-generated header dependencies (from -MMD -MP)
-include $(DEPS)
