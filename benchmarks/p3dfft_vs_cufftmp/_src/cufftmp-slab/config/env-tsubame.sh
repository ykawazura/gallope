#=======================================================================
# Run-time module environment for the cuFFTMp benchmarks on TSUBAME4.0.
#
# Not needed to build: config/Makefile.tsubame pins the same toolchain as
# absolute paths, so `make` works in a bare shell.  Keep the two in sync.
#
# nvhpc must be loaded FIRST, before openmpi/5.0.2-nvhpc, so that the latter's
# `prereq nvhpc` is satisfied by 26.1 and not by whatever happens to be around.
#
# GDRCopy: NVSHMEM (and hence cuFFTMp) uses libgdrapi.so.2 for its inter-node
# path.  The gdrdrv kernel module is loaded on the compute nodes and the library
# is installed, but under /usr/local/lib, which is not on the ldconfig search
# path -- so NVSHMEM's dlopen fails and it silently drops GDRCopy.  Harmless on
# one node, where no remote transport is opened at all; from 2 nodes on it is
# precisely the path being measured, so put it on LD_LIBRARY_PATH.
#=======================================================================
. /etc/profile.d/modules.sh

module purge
module load nvhpc/26.1_cuda13.1
module load openmpi/5.0.2-nvhpc

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
