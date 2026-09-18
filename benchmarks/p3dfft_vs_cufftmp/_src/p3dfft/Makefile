FC = mpif90
P3DFFT_DIR = $(HOME)/work/.local/p3dfft.2
FFTW_DIR = /work/opt/local/aarch64/apps/nvidia/25.9/fftw/3.3.10
INCLUDES = -I$(P3DFFT_DIR)/include
CPPFLAGS = -DHAVE_CONFIG_H -I.
FFLAGS = -DMEASURE -DSTRIDE1 -DFFTW $(INCLUDES) -g
LIBS = $(P3DFFT_DIR)/lib/libp3dfft.a $(FFTW_DIR)/lib/libfftw3.a

all: test_rand_f.x

test_rand_f.x: driver_rand.o
	$(FC) $(FFLAGS) -o test_rand_f.x driver_rand.o $(LIBS)

driver_rand.o: driver_rand.F90
	$(FC) $(CPPFLAGS) $(INCLUDES) $(FFLAGS) -c driver_rand.F90

clean:
	rm -f driver_rand.o test_rand_f.x

.PHONY: all clean
