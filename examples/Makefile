include $(UCFD_PATH)/Makefile.inc

LIB = $(UCFD_PATH)/lib/liblusgs.a

DIR = ./utils
OBJDIR = ./obj

MPIEXE = mpi3d.x
OMPEXE = omp3d.x

MPIOBJ = $(OBJDIR)/mpi3d.o
OMPOBJ = $(OBJDIR)/omp3d.o
ARRYAS = $(OBJDIR)/arrays.o
PBUTIL = $(OBJDIR)/pbutils.o
QUEUES = $(OBJDIR)/queue.o
RINPUT = $(OBJDIR)/readinput.o

OBJS = $(ARRYAS) $(PBUTIL) $(QUEUES) $(RINPUT)

$(OBJDIR)/%.o : $(DIR)/%.c $(DIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@ -I$(DIR)

$(MPIOBJ) : mpi3d.c $(OBJS)
	$(CC) $(CFLAGS) -c $< -lm -o $@ -I$(DIR) -I$(UCFD_PATH)/src

$(OMPOBJ) : omp3d.c $(OBJS)
	$(CC) $(CFLAGS) -lm -c $< -o $@ -I$(DIR) -I$(UCFD_PATH)/src

%(MPIEXE) : %(MPIOBJ) %(OBJS)
%(OMPEXE) : %(OMPOBJ) %(OBJS)

all : $(OBJS) $(MPIOBJ) $(OMPOBJ)
	$(CC) $(CFLAGS) $(OBJS) $(MPIOBJ) -o $(MPIEXE) $(LIB)
	$(CC) $(CFLAGS) $(OBJS) $(OMPOBJ) -o $(OMPEXE) $(LIB)
	mv $(MPIEXE) $(UCFD_PATH)/run
	mv $(OMPEXE) $(UCFD_PATH)/run


.PHONY : clean
clean :
	rm -rf $(OBJDIR) $(UCFD_PATH)/run/*.x $(UCFD_PATH)/run/*.txt
	rm -rf *.txt