#!/bin/bash
# Rebuild the fft micro-benchmark jobs, SPLIT into two PBS jobs per (node,trial):
#   fftgpu : cuFFTMp slab/pencil, launched EXACTLY like the working Gallope runs
#            -> select=N:mpiprocs=1  +  plain `mpirun ./exe`
#            (the old select=72 + `--map-by ppr:1:node` hangs NVSHMEM at node>=8)
#   fftcpu : P3DFFT (O3 + old-lib), flat-MPI 72 ranks/node  -> select=N:mpiprocs=72
# The per-sub-run dirs (slab-1024/, p3dfftO3-1024/, ...) already exist under
# runs/fft/nodeN-t and hold the exe symlinks; the new jobs just cd into them.
set -e
P=/work/gr96/o07001/bench-202607
SS=/work/gr96/o07001/p3dfft_vs_cufftmp
OLD=$P/results-p3dfft-oldlib

hdr_gpu () { # $1=nodes $2=walltime
cat <<EOF
#!/bin/sh -l
#PBS -q regular-g
#PBS -l select=$1:mpiprocs=1
#PBS -l walltime=$2
#PBS -W group_list=gr96
#PBS -j oe
module load nvidia/25.9
module load nv-hpcx
echo "===== ENV (R1-4) ====="; date; module list 2>&1 | tail -5
echo "unique nodes:" \$(sort -u \$PBS_NODEFILE | wc -l)
echo "======================"
EOF
}
hdr_cpu () { # $1=nodes $2=walltime
cat <<EOF
#!/bin/sh -l
#PBS -q regular-g
#PBS -l select=$1:mpiprocs=72
#PBS -l walltime=$2
#PBS -W group_list=gr96
#PBS -j oe
module load nvidia/25.9
module load nv-hpcx
module load hdf5
module load netcdf
module load netcdf-fortran
echo "===== ENV (R1-4) ====="; date; module list 2>&1 | tail -5
echo "unique nodes:" \$(sort -u \$PBS_NODEFILE | wc -l)
echo "======================"
EOF
}
# one sub-run block.  $1=subdir $2=launcher $3=dest
blk () {
cat <<EOF

echo '##### $1' \$(date)
cd $SRC/$1
$2 ./exe > stdout 2>&1
cp stdout '$3'
EOF
}

for d in $P/runs/fft/node*-*; do
  tag=$(basename "$d")                 # e.g. node16-1
  n=$(echo "$tag" | sed -E 's/node([0-9]+)-.*/\1/')
  SRC=$P/runs/fft/$tag                 # existing sub-run dirs live here
  gpu=$P/runs/fftgpu/$tag; cpu=$P/runs/fftcpu/$tag
  mkdir -p "$gpu" "$cpu"

  # ---- GPU job (cuFFTMp) : gallope-style launch ----
  { hdr_gpu "$n" 00:10:00
    blk slab-1024   "mpirun"  "$SS/cufftmp-slab/results-202607/miyabi/strong_scaling/1024^3/$tag"
    blk slab-2048   "mpirun"  "$SS/cufftmp-slab/results-202607/miyabi/strong_scaling/2048^3/$tag"
    blk pencil-1024 "mpirun"  "$SS/cufftmp-pencil/results-202607/miyabi/strong_scaling/1024^3/$tag"
    blk pencil-2048 "mpirun"  "$SS/cufftmp-pencil/results-202607/miyabi/strong_scaling/2048^3/$tag"
  } > "$gpu/job.pbs"

  # ---- CPU job (P3DFFT) : flat-MPI 72/node ----
  { hdr_cpu "$n" 00:20:00
    blk p3dfftO3-1024 "mpirun" "$SS/p3dfft/results-202607/miyabi/strong_scaling/1024^3/$tag"
    blk p3dfftg-1024  "mpirun" "$OLD/1024^3/$tag"
    blk p3dfftO3-2048 "mpirun" "$SS/p3dfft/results-202607/miyabi/strong_scaling/2048^3/$tag"
    blk p3dfftg-2048  "mpirun" "$OLD/2048^3/$tag"
  } > "$cpu/job.pbs"
done
echo "generated $(ls -d $P/runs/fftgpu/node*-* | wc -l) gpu + $(ls -d $P/runs/fftcpu/node*-* | wc -l) cpu jobs"
