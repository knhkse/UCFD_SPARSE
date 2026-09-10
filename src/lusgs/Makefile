# LU-SGS type methods library
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

BUILDDIR := $(BUILD_DIR)/lusgs
CPPFLAGS := $(MACROS) -I$(COMDIR) -I$(INCDIR)
SRCS := $(sort $(shell find $(SRCDIR) -type f -name '*.c'))

# Replace the absolute source prefix instead of prepending to it.
OBJS := $(patsubst $(SRCDIR)/%.c,$(BUILDDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)
SLIB := $(BUILDDIR)/libucfdlusgs.a
DLIB := $(BUILDDIR)/libucfdlusgs.so

# The source's directory supplies its individual include path.
$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(THIS_MAKEFILE)
	@mkdir -p "$(@D)"
	$(CC) $(CFLAGS) $(CPPFLAGS) -I"$(dir $<)" -MMD -MP -c "$<" -o "$@"

$(SLIB): $(OBJS) $(THIS_MAKEFILE)
	$(RM) "$@"
	$(AR) rcs "$@" $(OBJS)
	cp $(SLIB) $(UCFDPATH)/lib/$(notdir $(SLIB))

$(DLIB): $(OBJS) $(THIS_MAKEFILE)
	$(CC) $(LDFLAGS) -shared -o "$@" $(OBJS)
	cp $(DLIB) $(UCFDPATH)/lib/$(notdir $(DLIB))


.PHONY: all static dynamic clean
all : $(SLIB) $(DLIB)

static : $(SLIB)

dynamic: $(DLIB)

clean:
	$(RM) $(OBJS) $(DEPS) $(SLIB) $(DLIB)

# auto-generated header dependencies (from -MMD -MP)
-include $(DEPS)
