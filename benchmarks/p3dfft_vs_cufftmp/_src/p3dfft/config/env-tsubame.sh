#=======================================================================
# Run-time module environment for the P3DFFT benchmark on TSUBAME4.0.
#
# Same module set as calliope/batch/env-tsubame.sh -- ~/work/.local/p3dfft was
# built with it and libp3dfft.a carries gfortran .mod files, so nothing else can
# link against it.
#
# Two ordering traps:
#   * openmpi must be named explicitly.  The site default is openmpi/5.0.10-gcc,
#     but the gcc builds of FFTW/HDF5/netCDF all sit under an "openmpi5.0.2"
#     prefix.
#   * cuda/12.x must come BEFORE openmpi.  openmpi/5.0.2-gcc is CUDA-aware and
#     its libmpi.so.40 needs libcudart.so.12, but its `prereq cuda` resolves to
#     the default cuda/13.1.1 (libcudart.so.13), so anything linked against MPI
#     dies with "libcudart.so.12: cannot open shared object file".
#=======================================================================
. /etc/profile.d/modules.sh

module purge
module load cuda/12.3.2
module load openmpi/5.0.2-gcc
module load fftw/3.3.10-gcc
