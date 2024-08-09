include $(UCFD_PATH)/Makefile.inc

SLIB = liblusgs.a
DLIB = liblusgs.so

OBJDIR = ./obj
INCDIR = ./include

LUSGS = $(OBJDIR)/lusgs.o
FLUX = $(OBJDIR)/flux.o
CLUSGS = $(OBJDIR)/coloredlusgs.o
CONST = const.h

$(FLUX) : flux.c flux.h $(INCDIR)/$(CONST)
	$(CC) $(CFLAGS) -fPIC -imacros $(INCDIR)/$(CONST) -c $< -o $@

$(LUSGS) : lusgs.c lusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(CLUSGS) : coloredlusgs.c coloredlusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(SLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	ar qc $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib

all : $(SLIB) $(DLIB)

static : $(SLIB)

dynamic : $(DLIB)

.PHONY : clean
clean :
	rm -rf $(UCFD_PATH)/lib $(OBJDIR)
	rm -f $(SLIB) $(DLIB)