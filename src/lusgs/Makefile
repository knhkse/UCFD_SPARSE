include $(UCFD_PATH)/Makefile.inc

SLIB = liblusgs.a
SBLIB = libblusgs.a
DLIB = liblusgs.so
DBLIB = libblusgs.so

OBJDIR = ./obj
INCDIR = ./include

LUSGS = $(OBJDIR)/lusgs.o
BLUSGS = $(OBJDIR)/blusgs.o
CLUSGS = $(OBJDIR)/coloredlusgs.o
CBLUSGS = $(OBJDIR)/coloredblusgs.o
FLUX = $(OBJDIR)/flux.o
INVERSE = $(OBJDIR)/inverse.o
CONST = const.h

# Object files
$(FLUX) : flux.c flux.h $(INCDIR)/$(CONST)
	$(CC) $(CFLAGS) -fPIC -I$(INCDIR) -c $< -o $@

$(INVERSE) : inverse.c inverse.h
	$(CC) $(CFLAGS) -fPIC -c $< -o $@

$(LUSGS) : lusgs.c lusgs.h flux.h
	$(CC) $(CFLAGS) -fPIC -I$(INCDIR) -c $< -o $@

$(CLUSGS) : coloredlusgs.c coloredlusgs.h flux.h inverse.h
	$(CC) $(CFLAGS) -fPIC -I$(INCDIR) -c $< -o $@

$(BLUSGS) : blusgs.c blusgs.h
	$(CC) $(CFLAGS) -fPIC -I$(INCDIR) -c $< -o $@

$(CBLUSGS) : coloredblusgs.c coloredblusgs.h flux.h inverse.h
	$(CC) $(CFLAGS) -fPIC -I$(INCDIR) -c $< -o $@


# Library files
$(SLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(SBLIB) : $(FLUX) $(INVERSE) $(BLUSGS) $(CBLUSGS)
	ar rcs $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DLIB) : $(FLUX) $(LUSGS) $(CLUSGS)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib

$(DBLIB) : $(FLUX) $(INVERSE) $(BLUSGS) $(CBLUSGS)
	$(CC) -fPIC -shared -o $@ $^
	cp $@ $(UCFD_PATH)/lib


all : $(SLIB) $(SBLIB) $(DLIB) $(DBLIB)
	@echo All libraries are created

static : $(SLIB) $(SBLIB)
	@echo Static libraries are created

dynamic : $(DLIB) $(DBLIB)
	@echo Dynamic libraries are created

.PHONY : clean
clean :
	rm -rf $(OBJDIR)
	rm -f $(SLIB) $(SBLIB) $(DLIB) $(DBLIB)