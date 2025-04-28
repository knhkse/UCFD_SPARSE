include Makefile.inc

lib :
	mkdir -p lib
	cd src/krylov; mkdir -p obj; make all
	cd src/lusgs; mkdir -p obj; make all

lusgs :
	mkdir -p lib
	cd src/lusgs; mkdir -p obj; make all

example :
	cd examples; mkdir -p obj; make all

all :
	mkdir -p lib
	cd src; make all
	cd examples; mkdir -p obj; make all

.PHONY : clean
clean :
	cd src/krylov; make clean
	cd src/lusgs; make clean
	rm -rf $(UCFD_PATH)/lib
#	cd examples; make clean