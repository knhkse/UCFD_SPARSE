# Krylov subspace methods library
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

# Resolve paths from this Makefile, including when invoked with make -f.
THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCDIR := $(patsubst %/,%,$(dir $(THIS_MAKEFILE)))
INCDIR := $(UCFDPATH)/include

# Get configuration from `Makefile.inc`
include $(UCFDPATH)/Makefile.inc

# Default compiler
CC ?= mpicc
# Override GNU make's built-in C++ compiler default for MPI headers/linking.
ifeq ($(origin CXX), default)
CXX := mpicxx
endif
CXX ?= mpicxx
NVCC ?= nvcc
C_OPTFLAGS ?= -O3 -ffast-math
NVCC_OPTFLAGS ?= -O3 --use_fast_math
BUILD_DIR ?= $(abspath $(SRCDIR)/../../build)
CFLAGS := -Wall $(C_OPTFLAGS) -fPIC -fvisibility=hidden
CUDA_HOME ?= $(shell nvcc_path=$$(command -v $(NVCC)); if [ -n "$$nvcc_path" ]; then dirname "$$(dirname "$$(readlink -f "$$nvcc_path")")"; fi)
CUDA_ARCH ?= sm_80
NVCCFLAGS ?= $(NVCC_OPTFLAGS) -std=c++14 -arch=$(CUDA_ARCH)
CUDA_LIBDIR ?= $(CUDA_HOME)/lib64
CUDA_LDLIBS ?= -L"$(CUDA_LIBDIR)" -lcudart -lcublas -lcusparse

# Multi-threading flag
ifeq ($(USE_OMP), 1)
	ifneq ($(findstring icx,$(CC))$(findstring iicc,$(CC)),)
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

# Default MKL flags with static linking
ifeq ($(USE_MKL), 1)
	MKL_CPPFLAGS := -DUSE_MKL -I"$(MKLROOT)/include"
	MKL_LIBDIR := $(MKLROOT)/lib
	MKLFLAGS ?= -Wl,--start-group "$(MKL_LIBDIR)/libmkl_intel_lp64.a" \
		"$(MKL_LIBDIR)/libmkl_intel_thread.a" "$(MKL_LIBDIR)/libmkl_core.a" \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
endif


BUILDDIR := $(BUILD_DIR)/krylov
COMMON_CPPFLAGS := $(MACROS) $(addprefix -I,$(addprefix $(SRCDIR)/,utils impl)) -I$(INCDIR)

# Discover library sources only; test programs are separate executables.
SOURCE_DIRS := $(addprefix $(SRCDIR)/,matrix precon solver utils)
C_SRCS := $(sort $(shell find $(SOURCE_DIRS) -type f -name '*.c'))
CU_SRCS := $(sort $(shell find $(SOURCE_DIRS) -type f -name '*.cu' ! -name 'test*.cu'))

# Replace the absolute source prefix instead of prepending to it.
C_OBJS := $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(C_SRCS))
CU_OBJS := $(patsubst $(SRCDIR)/%.cu,$(BUILDDIR)/%.cu.o,$(CU_SRCS))
OBJS := $(C_OBJS) $(CU_OBJS)
DEPS := $(OBJS:.o=.d)
SLIB := $(BUILDDIR)/libucfdkrylov.a
DLIB := $(BUILDDIR)/libucfdkrylov.so

# Only matrix files include Intel MKL
MKL_SRCS := $(wildcard $(SRCDIR)/matrix/*.c)
MKL_OBJS := $(patsubst $(SRCDIR)/matrix/%.c,$(BUILDDIR)/matrix/%.o,$(MKL_SRCS))

$(MKL_OBJS): PRIVATE_MKL_CPPFLAGS = $(MKL_CPPFLAGS)

# The source's directory supplies its individual include path.
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(THIS_MAKEFILE) | check-mkl
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) $(COMMON_CPPFLAGS) -I"$(dir $<)" $(PRIVATE_MKL_CPPFLAGS) -MMD -MP -c "$<" -o "$@"

# Keep CUDA objects distinct from C objects with the same source basename.
$(BUILDDIR)/%.cu.o: $(SRCDIR)/%.cu $(THIS_MAKEFILE) | check-mkl
	@mkdir -p "$(@D)"
	$(NVCC) $(NVCCFLAGS) -allow-unsupported-compiler -rdc=true -dlto -ccbin $(CXX) \
		-Xcompiler=-Wall,-Wextra,-fPIC,-fvisibility=hidden $(if $(OMPFLAG),-Xcompiler=$(OMPFLAG)) \
		$(COMMON_CPPFLAGS) -I"$(SRCDIR)/matrix" -I"$(dir $<)" \
		-MMD -MP -c "$<" -o "$@"

# CUDA dlink object
CUDA_DLINK_OBJ := $(BUILDDIR)/cuda_dlink.o

$(CUDA_DLINK_OBJ): $(CU_OBJS) $(THIS_MAKEFILE)
	$(NVCC) $(NVCCFLAGS) -allow-unsupported-compiler -ccbin $(CXX) \
		-dlink -dlto -Xcompiler=-fPIC \
		$(CU_OBJS) -o "$@"

$(SLIB): $(OBJS) $(THIS_MAKEFILE) | check-mkl
	$(RM) "$@"
	$(AR) rcs "$@" $(OBJS)
	@mkdir -p "$(UCFDPATH)/lib"
	cp $(SLIB) $(UCFDPATH)/lib/$(notdir $(SLIB))

$(DLIB): $(OBJS) $(CUDA_DLINK_OBJ) $(THIS_MAKEFILE) | check-mkl
	$(CXX) $(LDFLAGS) $(OMPFLAG) -shared -o "$@" \
	$(OBJS) $(CUDA_DLINK_OBJ) $(MKLFLAGS) $(CUDA_LDLIBS) $(LDLIBS)
	@mkdir -p "$(UCFDPATH)/lib"
	cp $(DLIB) $(UCFDPATH)/lib/$(notdir $(DLIB))


.PHONY: all static dynamic clean check-mkl
all : $(SLIB) $(DLIB)

static : $(SLIB)

dynamic: $(DLIB)

clean:
	$(RM) $(OBJS) $(DEPS) $(SLIB) $(DLIB) $(CUDA_DLINK_OBJ)

check-mkl:
ifeq ($(USE_MKL), 1)
	@test -n "$(strip $(MKLROOT))" || { echo "MKLROOT must be set when USE_MKL=1" >&2; exit 1; }
endif

# auto-generated header dependencies (from -MMD -MP)
-include $(DEPS)
