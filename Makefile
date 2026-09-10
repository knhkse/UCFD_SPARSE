include Makefile.inc

# TBU
# CUDA ?= 0

# ifeq ($(CUDA),1)
#   LUSGS_CUDA_BUILD = cd src/lusgs/cuda; mkdir -p obj; make all
# endif

lib :
	mkdir -p lib
	make -C src/lusgs all
	make -C src/krylov all

lusgs :
	mkdir -p lib
	make -C src/lusgs all

krylov :
	mkdir -p lib
	make -C src/krylov all

example :
	mkdir -p examples/obj
	make -C examples

all :
	mkdir -p lib
	make -C src/lusgs all
	make -C src/krylov all
	cd examples; mkdir -p obj; make all


.PHONY : clean
clean :
	make -C src/krylov clean
	make -C src/lusgs clean
	make -C examples clean
	rm -rf $(UCFDPATH)/lib
