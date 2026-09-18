#!/bin/sh
#=======================================================================
# Per-rank GPU binding for gallope on TSUBAME4.0 (4 x H100 per node).
#
# init_cuda (src/cuFFTmp.F90) already derives device_id from the intra-node
# MPI rank and calls cudaSetDevice, but that is too late for OpenACC: the
# NVHPC OpenACC runtime initialises itself -- and creates the device-side
# descriptors of every "!$acc declare create" module array -- before the
# Fortran main program starts. Those descriptors are then only ever updated
# on the device the runtime picked at start-up (device 0). A rank running on
# any other GPU therefore launches kernels whose dope vectors are still zero,
# and the first !$acc region dies with CUDA error 700 even though the host
# side reports the arrays as present with valid device pointers.
#
# ACC_DEVICE_NUM is read at that start-up, so setting it here -- before the
# executable is exec'd -- is what actually pins the OpenACC device. It must
# agree with the device_id computed in init_cuda, i.e. the intra-node rank.
#
# Only needed where a node holds several GPUs; on 1-GPU-per-node machines
# (Miyabi) the local rank is always 0 and this is a no-op.
#=======================================================================
export ACC_DEVICE_NUM=${OMPI_COMM_WORLD_LOCAL_RANK:-0}
exec "$@"
