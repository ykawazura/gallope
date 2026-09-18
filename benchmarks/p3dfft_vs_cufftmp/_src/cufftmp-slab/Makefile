f90   := mpif90

# コンパイルフラグ
CFLAGS  = -O3 -Mfree -fast -Mextend -Mpreprocess -Minform=warn
CFLAGS += -Minfo=accel -cuda -cudalib=cufftmp -gpu=mem:unified:nomanagedalloc

# リンクフラグ（-cudaのみ、コンパイルオプションは除外）
LFLAGS  = -cuda -gpu=mem:unified:nomanagedalloc
LFLAGS += -L${NVHPC_ROOT}/compilers/lib -lnvhpcwrapcufft -lnvhpcwrapcufftmp
LFLAGS += -L${NVHPC_ROOT}/comm_libs/nvshmem/lib -lnvshmem
LFLAGS += -L${NVHPC_ROOT}/math_libs/lib64/ -lcufftMp
LFLAGS += -L${NVHPC_ROOT}/comm_libs/mpi/lib/ -lmpi

exe = cufftmp_r2c

all : $(exe)

$(exe) : $(exe).o
	$(f90) -o $@ $< $(LFLAGS)

%.o : %.f90
	$(f90) -c $< -o $@ $(CFLAGS)

.PHONY: clean
clean:
	@rm -rf *.mod *.o $(exe)

run: $(exe)
	LD_LIBRARY_PATH="${NVSHMEM_LIB}:${CUFFT_LIB}:${LD_LIBRARY_PATH}" mpirun -n 2 $(exe)
