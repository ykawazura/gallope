#=====================================================================
# cuFFTMp subcommunicator smoke test build
#   Reuses /work/gr96/o07001/p3dfft_vs_cufftmp/cufftmp-slab/Makefile
#   (link settings proven on 256 nodes / 2048^3); only the exe name changed.
#   Requires: module load nvidia nv-hpcx
#=====================================================================
f90 := mpif90

# compile flags (Grace-Hopper unified memory)
CFLAGS  = -O3 -Mfree -fast -Mextend -Mpreprocess -Minform=warn
CFLAGS += -Minfo=accel -cuda -cudalib=cufftmp -gpu=mem:unified:nomanagedalloc

# link flags
LFLAGS  = -cuda -gpu=mem:unified:nomanagedalloc
LFLAGS += -L${NVHPC_ROOT}/compilers/lib -lnvhpcwrapcufft -lnvhpcwrapcufftmp
LFLAGS += -L${NVHPC_ROOT}/comm_libs/nvshmem/lib -lnvshmem
LFLAGS += -L${NVHPC_ROOT}/math_libs/lib64/ -lcufftMp
LFLAGS += -L${NVHPC_ROOT}/comm_libs/mpi/lib/ -lmpi

exe = cufftmp_subcomm

all : $(exe)

$(exe) : $(exe).o
	$(f90) -o $@ $< $(LFLAGS)

%.o : %.f90
	$(f90) -c $< -o $@ $(CFLAGS)

.PHONY: clean
clean:
	@rm -rf *.mod *.o $(exe)
