include $(UCFD_PATH)/Makefile.inc

MKL_REQD := all dynamic
define ERR_MSG
MKLFlagsUndefinedError ***
`MKLFLAGS` must be defined in `Makefile.inc` when you want to make shared library.
Run `Prep.sh` bash script with appropriate flags.
For more information, refer to the document
endef
ifneq (,$(filter $(MKL_REQD),$(MAKECMDGOALS)))
  ifeq ($(origin MKLFLAGS),undefined)
    $(error $(ERR_MSG))
  endif
endif

SLIB = libkrylov.a
DLIB = libkrylov.so

OBJDIR = ./obj
INCDIR = ./include

KRYLOV = $(OBJDIR)/krylov.o
PRECON = $(OBJDIR)/precon.o
TYPES = types.h

$(KRYLOV) : krylov.c krylov.h precon.h $(INCDIR)/$(TYPES)
	$(CC) $(CFLAGS) -fPIC -I${MKLROOT}/include -I$(INCDIR) -c $< -o $@

$(PRECON) : precon.c precon.h $(INCDIR)/$(TYPES)
	$(CC) $(CFLAGS) -fPIC -I${MKLROOT}/include -I$(INCDIR) -c $< -o $@

$(SLIB) : $(KRYLOV) $(PRECON)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(KRYLOV) $(PRECON)
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