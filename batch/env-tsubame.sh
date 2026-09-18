#=======================================================================
# Module environment for gallope on TSUBAME4.0 (Institute of Science Tokyo)
#
# Required at RUN time (LD_LIBRARY_PATH); batch/job.sh-tsubame sources it.
# Not needed to build -- arch/tsubame.in pins the same toolchain as absolute
# paths, so 'make' works in a bare shell. Keep the two in sync.
#
#   source batch/env-tsubame.sh
#
# nvhpc must be loaded FIRST: the hdf5-parallel modulefile pulls in
# openmpi/5.0.2-nvhpc, whose 'prereq nvhpc' would otherwise be satisfied by
# whatever nvhpc happens to be around. netCDF/HDF5 exist only as nvhpc24.1
# builds, but they are C libraries and their .mod files are readable by the
# 26.1 nvfortran, so the mismatch in the module names is only cosmetic --
# see arch/tsubame.in for why the compiler must be 26.1.
#=======================================================================
. /etc/profile.d/modules.sh

module purge
module load nvhpc/26.1_cuda13.1
module load openmpi/5.0.2-nvhpc
module load hdf5-parallel/1.14.3/nvhpc24.1
module load netcdf-parallel/4.9.2/nvhpc24.1

#-----------------------------------------------------------------------
# GDRCopy, used by NVSHMEM (and hence cuFFTMp) for inter-node transfers.
#
# The gdrdrv kernel module is loaded on the compute nodes and libgdrapi.so.2
# is installed, but under /usr/local/lib, which is not in the ldconfig search
# path -- so NVSHMEM's dlopen of it fails and it silently disables GDRCopy.
# Harmless on a single node, where the remote transport is never used.
#-----------------------------------------------------------------------
export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
