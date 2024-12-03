include $(UCFD_PATH)/Makefile.inc

SLIB = liblusgs.a
SBLIB = libblusgs.a
DLIB = liblusgs.so
DBLIB = libblusgs.so

OBJDIR = ./obj
INCDIR = ./include

LUSGS = $(OBJDIR)/lusgs.o
BLUSGS = $(OBJDIR)/blusgs.o
FLUX = $(OBJDIR)/flux.o
CLUSGS = $(OBJDIR)/coloredlusgs.o
CBLUSGS = $(OBJDIR)/coloredblusgs.o
INVERSE = $(OBJDIR)/inverse.o
CONST = const.h

$(FLUX) : flux.c flux.h $(INCDIR)/$(CONST)
	$(CC) $(CFLAGS) -fPIC -imacros $(INCDIR)/$(CONST) -c $< -o $@

$(INVERSE) : inverse.c inverse.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(LUSGS) : lusgs.c lusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(CLUSGS) : coloredlusgs.c coloredlusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(BLUSGS) : blusgs.c blusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(CBLUSGS) : coloredblusgs.c coloredblusgs.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(SLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	ar qc $@ $^
	cp $@ $(UCFD_PATH)/lib

$(SBLIB) : $(FLUX) $(INVERSE) $(BLUSGS) $(CBLUSGS)
	ar qc $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DBLIB) : $(FLUX) $(INVERSE) $(BLUSGS) $(CBLUSGS)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib

all : $(SLIB) $(SBLIB) $(DLIB) $(DBLIB)

static : $(SLIB) $(SBLIB)

dynamic : $(DLIB) $(DBLIB)

.PHONY : clean
clean :
	rm -rf $(UCFD_PATH)/lib $(OBJDIR)
	rm -f $(SLIB) $(SBLIB) $(DLIB) $(DBLIB)