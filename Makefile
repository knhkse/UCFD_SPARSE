include Makefile.inc

CUDA ?= 0

ifeq ($(CUDA),1)
  LUSGS_CUDA_BUILD = cd src/lusgs/cuda; mkdir -p obj; make all
endif

lib :
	mkdir -p lib
	cd src/krylov; mkdir -p obj; make all
	cd src/lusgs; mkdir -p obj; make all
	$(LUSGS_CUDA_BUILD)

lusgs :
	mkdir -p lib
	cd src/lusgs; mkdir -p obj; make all
	$(LUSGS_CUDA_BUILD)

krylov :
	mkdir -p lib
	cd src/krylov; mkdir -p obj; make all

example :
	cd examples; mkdir -p obj; make all

all :
	mkdir -p lib
	cd src/krylov; mkdir -p obj; make all
	cd src/lusgs; mkdir -p obj; make all
	$(LUSGS_CUDA_BUILD)
	cd examples; mkdir -p obj; make all

.PHONY : clean
clean :
	cd src/krylov; make clean
	cd src/lusgs; make clean
	rm -rf $(UCFD_PATH)/lib
	cd examples; make clean
