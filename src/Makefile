#=======================================================================
# Makefile for gallope modified by Y. Kawazura @ Jun 24 2024
# Basically, don't change this file!
#=======================================================================
include ../Makefile.in
include ../arch/$(arch).in

MODEL_DIR = model/$(MODEL)

FLAGS += $(NETCDF_INC)
FLAGS += -Minfo=accel -cuda -cudalib=cufftmp 
LINKER := -L${NVHPC_ROOT}/compilers/lib -lnvhpcwrapcufft -L${NVHPC_ROOT}/math_libs/lib64/ -lcufftMp -lmpi $(NETCDF_LIB)

ifeq ($(DEBG),yes)
	FLAGS += -DDEBG
endif

### List of files for 'MHD_INCOMP'
SRC_MHD_INCOMP = \
			time_stamp.F90 file.F90 \
			$(MODEL_DIR)/params.F90 mp.F90 mpiio.F90 \
			terminate.F90 grid.F90 cuFFTmp.F90 time.F90 utils.F90 dealias.F90 \
			shearing_box.F90 \
			$(MODEL_DIR)/fields.F90 \
			$(MODEL_DIR)/force.F90 \
			$(MODEL_DIR)/$(MODEL)_specific.F90 \
			$(MODEL_DIR)/io.F90 $(MODEL_DIR)/diagnostics.F90 \
			$(MODEL_DIR)/advance.F90 \
			gallope.F90 

### List of files for 'RMHD'
SRC_RMHD = $(SRC_MHD_INCOMP)

### List of files for 'KRMHD'
SRC_KRMHD = $(SRC_MHD_INCOMP)

OBJ = $(SRC_$(MODEL):.F90=.o)

gallope : $(OBJ)
	$(F90) -o ../$@ $+ $(FLAGS) $(LINKER) -acc 

%.o : %.F90
	$(F90) -c $< -o $@ $(FLAGS) -D$(MODEL) -D_MODEL_=\"$(MODEL)\" -acc 


.PHONY: print clean

print:
	@echo OBJ = $(OBJ)
	@echo LIBS = $(LIBS)
	@echo F90FLAGS = $(F90FLAGS)
	@echo LDFLAGS = $(LDFLAGS)

clean:
	rm -f *.o *.mod $(MODEL_DIR)/*.o $(MODEL_DIR)/*.mod gallope
