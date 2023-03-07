LIBS = -fopenmp -ftree-vectorize -labpoa -lz -lm -lz -lhts

DEPS = str_utilities.h bloom_filter.h

CC=g++

OBJ = str_utilities.o bloom_filter.o strtoolkit.o 

HTS?=install

CONDA_LIB= -I/.mounts/labs/simpsonlab/sw/miniconda3/envs/strline_dup/lib/ -L/.mounts/labs/simpsonlab/sw/miniconda3/envs/strline_dup/lib/

ABPOA_INCLUDE= -I./abPOA/include -I./abPOA -I./abPOA/src -L./abPOA/lib 

#ifeq ($(HTS), install)
HTS_INCLUDE=-I./htslib -I./htslib/htslib -L./htslib/htslib -L./htslib
#else
#	HTS_INCLUDE=
#endif	

#all: librealign.so

#clean:
#	rm librealign.so ReAlign.o

# If htslib needs to be compiled
#htslib/libhts.so:
#	cd htslib && make CC='gcc' || exit 255

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -fPIC $(ABPOA_INCLUDE) $(HTS_INCLUDE) $(LIBS)

strtoolkit: $(OBJ)
	$(CC) -g -Wall -Wextra -std=c++17 -O2 -o $@ $^ -fPIC $(ABPOA_INCLUDE) $(HTS_INCLUDE) $(LIBS) $(CFLAGS)

