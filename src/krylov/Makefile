# Krylov subspace methods library
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

# Resolve paths from this Makefile, including when invoked with make -f.
THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
SRCDIR := $(patsubst %/,%,$(dir $(THIS_MAKEFILE)))
COMDIR := $(UCFDPATH)/src/common
INCDIR := $(UCFDPATH)/include

# Get configuration from `Makefile.inc`
include $(UCFDPATH)/Makefile.inc

# Default compiler - Recommend Intel compiler
CC ?= mpiicx
OPTFLAGS ?= -O3 -ffast-math
BUILD_DIR ?= $(abspath $(SRCDIR)/../../build)
CFLAGS := -Wall $(OPTFLAGS) -fPIC -fvisibility=hidden

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

# Default MKL flags with static linking
ifeq ($(USE_MKL), 1)
	MKL_CPPFLAGS := -DUSE_MKL -I"$(MKLROOT)/include"
	MKL_LIBDIR := $(MKLROOT)/lib
	MKLFLAGS ?= -Wl,--start-group "$(MKL_LIBDIR)/libmkl_intel_lp64.a" \
		"$(MKL_LIBDIR)/libmkl_intel_thread.a" "$(MKL_LIBDIR)/libmkl_core.a" \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl
endif


BUILDDIR := $(BUILD_DIR)/krylov
COMMON_CPPFLAGS := $(MACROS) $(addprefix -I,$(addprefix $(SRCDIR)/,utils impl)) -I$(COMDIR) -I$(INCDIR)

# Discover library sources only; root-level test programs are separate executables.
SOURCE_DIRS := $(addprefix $(SRCDIR)/,matrix precon solver utils)
SRCS := $(sort $(shell find $(SOURCE_DIRS) -type f -name '*.c'))

# Replace the absolute source prefix instead of prepending to it.
OBJS := $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)
SLIB := $(BUILDDIR)/libucfdkrylov.a
DLIB := $(BUILDDIR)/libucfdkrylov.so

# Only matrix files include Intel MKL
# MKL_OBJS := $(addprefix $(BUILDDIR)/matrix/,sparsemat.o mpicsrmat.o mpibsrmat.o)
MKL_SRCS := $(wildcard $(SRCDIR)/matrix/*.c)
MKL_OBJS := $(patsubst $(SRCDIR)/matrix/%.c,$(BUILDDIR)/matrix/%.o,$(MKL_SRCS))

$(MKL_OBJS): PRIVATE_MKL_CPPFLAGS = $(MKL_CPPFLAGS)

# The source's directory supplies its individual include path.
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(THIS_MAKEFILE) | check-mkl
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) $(COMMON_CPPFLAGS) -I"$(dir $<)" $(PRIVATE_MKL_CPPFLAGS) -MMD -MP -c "$<" -o "$@"

$(SLIB): $(OBJS) $(THIS_MAKEFILE) | check-mkl
	$(RM) "$@"
	$(AR) rcs "$@" $(OBJS)
	cp $(SLIB) $(UCFDPATH)/lib/$(notdir $(SLIB))

$(DLIB): $(OBJS) $(THIS_MAKEFILE) | check-mkl
	$(CC) $(LDFLAGS) -shared -o "$@" $(OBJS) $(MKLFLAGS) $(LDLIBS)
	cp $(DLIB) $(UCFDPATH)/lib/$(notdir $(DLIB))


.PHONY: all static dynamic clean check-mkl
all : $(SLIB) $(DLIB)

static : $(SLIB)

dynamic: $(DLIB)

clean:
	$(RM) $(OBJS) $(DEPS) $(SLIB) $(DLIB)

check-mkl:
ifeq ($(USE_MKL), 1)
	@test -n "$(strip $(MKLROOT))" || { echo "MKLROOT must be set when USE_MKL=1" >&2; exit 1; }
endif

# auto-generated header dependencies (from -MMD -MP)
-include $(DEPS)
