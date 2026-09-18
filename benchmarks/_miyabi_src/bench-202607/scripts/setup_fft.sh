#!/bin/bash
set -euo pipefail
ROOT=/work/gr96/o07001
BENCH=$ROOT/bench-202607
BIN=$BENCH/bin
PV=$ROOT/p3dfft_vs_cufftmp
RUNS=$BENCH/runs/fft

[ -d "$RUNS" ] && rm -r "$RUNS"
mkdir -p "$RUNS" "$BENCH/results-p3dfft-oldlib"
for k in cufftmp-slab cufftmp-pencil p3dfft; do
  mkdir -p "$PV/$k/results-202607/miyabi/strong_scaling/1024^3" \
           "$PV/$k/results-202607/miyabi/strong_scaling/2048^3"
done
mkdir -p "$BENCH/results-p3dfft-oldlib/1024^3" "$BENCH/results-p3dfft-oldlib/2048^3"

# node-count -> benchmark list.  order: cheap/robust first, riskiest last
list_for () {
case $1 in
  1)   echo "slab:1024 pencil:1024 p3dfftO3:1024 p3dfftg:1024" ;;
  2)   echo "slab:1024 p3dfftO3:1024 p3dfftg:1024" ;;
  4)   echo "slab:1024 slab:2048 pencil:1024 p3dfftO3:1024 p3dfftg:1024" ;;
  8)   echo "slab:1024 slab:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
  16)  echo "slab:1024 slab:2048 pencil:1024 pencil:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
  32)  echo "slab:1024 slab:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
  64)  echo "slab:1024 slab:2048 pencil:1024 pencil:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
  128) echo "slab:1024 slab:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
  256) echo "slab:1024 slab:2048 pencil:1024 pencil:2048 p3dfftO3:1024 p3dfftg:1024 p3dfftO3:2048 p3dfftg:2048" ;;
esac
}

binof ()  { case $1 in slab) echo "$BIN/fft-slab-$2";; pencil) echo "$BIN/fft-pencil-$2";;
                       p3dfftO3) echo "$BIN/fft-p3dfft-$2-O3";; p3dfftg) echo "$BIN/fft-p3dfft-$2-g";; esac; }
destof () { case $1 in slab) echo "$PV/cufftmp-slab/results-202607/miyabi/strong_scaling/$2^3";;
                       pencil) echo "$PV/cufftmp-pencil/results-202607/miyabi/strong_scaling/$2^3";;
                       p3dfftO3) echo "$PV/p3dfft/results-202607/miyabi/strong_scaling/$2^3";;
                       p3dfftg) echo "$BENCH/results-p3dfft-oldlib/$2^3";; esac; }

for n in 1 2 4 8 16 32 64 128 256; do
for t in 1 2 3; do
  d="$RUNS/node$n-$t"; mkdir -p "$d"
  { cat <<H
#!/bin/sh -l
#PBS -q regular-g
#PBS -l select=$n:mpiprocs=72
#PBS -l walltime=00:30:00
#PBS -W group_list=gr96
#PBS -j oe

module load nvidia/25.9
module load nv-hpcx
module load hdf5
module load netcdf
module load netcdf-fortran

echo "===== ENVIRONMENT (for R1-4) ====="
date
module list 2>&1
nvfortran --version 2>&1 | head -3
mpif90 -show 2>&1 | head -1
echo "unique nodes: \$(sort -u \$PBS_NODEFILE | wc -l)"
echo "=================================="
H
    for item in $(list_for $n); do
      kind=${item%%:*}; sz=${item##*:}
      sub="$kind-$sz"
      mkdir -p "$d/$sub"
      ln -sfn "$(binof $kind $sz)" "$d/$sub/exe"
      case $kind in
        slab|pencil) mp="mpirun -n $n --map-by ppr:1:node ./exe" ;;
        *)           mp="mpirun ./exe" ;;
      esac
      echo ""
      echo "echo '##### $sub' \$(date)"
      echo "cd \${PBS_O_WORKDIR}/$sub"
      echo "$mp > stdout 2>&1"
      echo "cp stdout '$(destof $kind $sz)/node$n-$t'"
    done
  } > "$d/job.pbs"
done
done

echo "=== generated $(find "$RUNS" -name job.pbs | wc -l) fft jobs"
echo "=== sample node16-1:"; cat "$RUNS/node16-1/job.pbs"
