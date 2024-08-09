include Makefile.inc

lib :
	mkdir -p lib
	cd src; mkdir -p obj; make all

example :
	cd examples; mkdir -p obj; make all

all :
	mkdir -p lib
	cd src; mkdir -p obj; make all
	cd examples; mkdir -p obj; make all

static :
	mkdir -p lib
	cd src; mkdir -p obj; make static

dynamic :
	mkdir -p lib
	cd src; mkdir -p obj; make dynamic
	

.PHONY : clean
clean :
	cd src; make clean
	cd examples; make clean