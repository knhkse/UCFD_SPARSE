include $(UCFD_PATH)/Makefile.inc

CC := mpiicx

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

DIR = ./utils
OBJDIR = ./obj

LIB = $(UCFD_PATH)/lib/liblusgs.a
CLIB = $(UCFD_PATH)/lib/libclusgs.a
MPIEXE = mpi3d.x
OMPEXE = omp3d.x

IFLAGS += -I$(DIR)
IFLAGS += -I$(UCFD_PATH)/src/lusgs/include

MPIOBJ = $(OBJDIR)/mpi3d.o
OMPOBJ = $(OBJDIR)/omp3d.o
ARRYAS = $(OBJDIR)/arrays.o
PBUTIL = $(OBJDIR)/pbutils.o
INPUT = $(OBJDIR)/readinput.o

OBJS = $(ARRYAS) $(PBUTIL) $(INPUT)

$(OBJDIR)/%.o : $(DIR)/%.c $(DIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@ -I$(DIR)

$(MPIOBJ) : mpi3d.c $(OBJS)
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

$(OMPOBJ) : omp3d.c $(OBJS)
	$(CC) $(CFLAGS) -c $< -o $@ $(IFLAGS)

%(MPIEXE) : %(MPIOBJ) %(OBJS)
%(OMPEXE) : %(OMPOBJ) %(OBJS)

all : $(OBJS) $(MPIOBJ) $(OMPOBJ)
	$(CC) $(CFLAGS) $(OBJS) $(MPIOBJ) -o $(MPIEXE) $(LIB)
	$(CC) $(CFLAGS) $(OBJS) $(OMPOBJ) -o $(OMPEXE) $(CLIB)
	mv $(MPIEXE) $(UCFD_PATH)/run
	mv $(OMPEXE) $(UCFD_PATH)/run

.PHONY : clean
clean :
	rm -rf $(OBJDIR) $(UCFD_PATH)/run/*.x