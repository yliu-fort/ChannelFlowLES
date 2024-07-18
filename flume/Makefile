CFLAGS += -O2 -I../src/
SOURCE := main.c

ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
INST = $(shell sudo ln -sf $(ROOT_DIR)/bin/* /usr/local/bin/)
UNINST = $(shell ls $(ROOT_DIR)/bin/* | xargs -n 1 find -L /usr/local/bin -samefile | xargs sudo rm -f )

all:
	make clean
	mkdir -p bin
	make mpi3d

2d: main.c
	rm -f _main.c
	make clean
	$(QCC) $(CFLAGS) -events -D_MPI=0 -o bin/$@ $^ -lm

3d: main.c
	rm -f _main.c
	make clean
	$(QCC) $(CFLAGS) -D_MPI=0 -grid=octree -o bin/$@ $^ -lm

mpi2d: main.c
	rm -f _main.c
	make clean
	$(QCC) -source $(CFLAGS) -D_MPI=1 $<
	mpicc -I$(BASILISK) $(CFLAGS) -o bin/$@ _main.c -lm

mpi3d: main.c
	rm -f _main.c
	make clean
	$(QCC) -source $(CFLAGS) -D_MPI=1 -grid=octree $<
	mpicc -I$(BASILISK) $(CFLAGS) -o bin/$@ _main.c -lm

winstall:
	echo $(ROOT_DIR)
	echo $(INST)

wuninstall:
	echo $(ROOT_DIR)
	echo $(UNINST)

allclean:
	rm -f *.exec _*.c 
	rm -rf bin .qcc*
	make clean
	
include $(BASILISK)/Makefile.defs

