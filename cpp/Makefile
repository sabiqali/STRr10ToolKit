LIBS = -lz -lm -llzma -lbz2 -lpthread -lcurl -lssl -lcrypto -lz -lhts -ftree-vectorize

DEPS = str_utilities.h

CC=g++

OBJ = str_utilities.o STRr10ToolKit.o 

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

STRr10ToolKit: $(OBJ)
	$(CC) -shared -g -Wall -Wextra -std=c++14 -fPIC -o $@ $^ $(HTS_INCLUDE) $(LIBS) $(CFLAGS)

