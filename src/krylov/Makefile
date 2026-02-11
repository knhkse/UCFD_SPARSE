# Krylov subspace methods library Makefile

# =================================
# Configuration Variables
# =================================
CC ?= icx
include $(UCFD_PATH)/Makefile.inc

# =================================
# Derived Settings
# =================================
ifeq ($(MKLROOT),)
  $(error Intel MKL root directory must be defined)
else
  IFLAGS += -I$(MKLROOT)/include
endif

ifeq ($(MKLFLAGS),)
  $(error MKL flags must be defined)
endif

# Block size : Required
ifeq ($(BLOCK),)
  $(error `BLOCK` must be defined)
else
  CFLAGS += -DBLOCK=$(BLOCK)
endif

ifeq ($(EPS),)
  $(info Machine epsilon undefined::Default value is 2.22e-16)
else
  CFLAGS += -DEPS=$(EPS)
endif

# =================================
# Build Targets
# =================================
SLIB = libkrylov.a
DLIB = libkrylov.so

OBJDIR = ./obj
INCDIR = ./include
IFLAGS += -I$(INCDIR)

KRYLOV = $(OBJDIR)/krylov.o
PRECON = $(OBJDIR)/precon.o
INVERSE = $(OBJDIR)/inverse.o

# Build object(.o) files
$(INVERSE) : inverse.c $(INCDIR)/inverse.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(KRYLOV) : krylov.c $(INCDIR)/krylov.h $(INCDIR)/precon.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

$(PRECON) : precon.c $(INCDIR)/precon.h $(CONFIG) $(TYPES)
	$(CC) $(CFLAGS) -fPIC $(IFLAGS) -c $< -o $@

# Build library files
$(SLIB) : $(KRYLOV) $(PRECON) $(INVERSE)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(KRYLOV) $(PRECON) $(INVERSE)
	$(CC) $(CFLAGS) -fPIC -shared $(MKLFLAGS) -o $@ $^
	cp $@ $(UCFD_PATH)/lib


all : $(SLIB) $(DLIB)
	@echo All libraries are created

static : $(SLIB)
	@echo Static library is created

dynamic : $(DLIB)
	@echo Dynamic library is created

.PHONY : clean
clean :
	rm -rf $(OBJDIR)
	rm -f $(SLIB) $(DLIB)