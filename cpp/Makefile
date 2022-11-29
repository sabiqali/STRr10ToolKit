LIBS = -fopenmp -ftree-vectorize -lz -lm -llzma -lbz2 -lz -lhts

DEPS = str_utilities.h

CC=gcc

OBJ = str_utilities.o strtoolkit.o 

HTS?=install

ifeq ($(HTS), install)
	HTS_INCLUDE=-I./htslib -I./htslib/htslib -L./htslib/htslib -L./htslib
else
	HTS_INCLUDE=
endif	

#all: librealign.so

#clean:
#	rm librealign.so ReAlign.o

# If htslib needs to be compiled
#htslib/libhts.so:
#	cd htslib && make CC='gcc' || exit 255

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(HTS_INCLUDE) $(LIBS)

strtoolkit: $(OBJ)
	$(CC) -g -Wall -Wextra -std=c++17 -O2 -o $@ $^ $(HTS_INCLUDE) $(LIBS) $(CFLAGS)

