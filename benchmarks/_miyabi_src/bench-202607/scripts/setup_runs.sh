#!/bin/bash
set -euo pipefail
ROOT=/work/gr96/o07001
BENCH=$ROOT/bench-202607
BIN=$BENCH/bin
C=$ROOT/calliope_performance
G=$ROOT/gallope_performance
RUNS=$BENCH/runs

# freeze the inputs
cp -p "$C/calliope.in" "$BENCH/inputs/calliope.in"
cp -p "$G/gallope.in"  "$BENCH/inputs/gallope.in"

mkdir -p "$C/results-202607/1024^3" "$G/results-202607/1024^3" "$G/results-202607/1024^3-unified" "$BENCH/results-phaseA"
[ -d "$RUNS" ] && rm -r "$RUNS"
mkdir -p "$RUNS"

hdr () {  # $1 nodes  $2 mpiprocs  $3 walltime
cat <<H
#!/bin/sh -l
#PBS -q regular-g
#PBS -l select=$1:mpiprocs=$2
#PBS -l walltime=$3
#PBS -W group_list=gr96
#PBS -j oe

module load nvidia
module load nv-hpcx
module load hdf5
module load netcdf
module load netcdf-fortran

echo "===== ENVIRONMENT (for R1-4) ====="
date
echo "--- module list"; module list 2>&1
echo "--- nvfortran"; nvfortran --version 2>&1 | head -3
echo "--- mpirun";    mpirun --version 2>&1 | head -2
echo "--- nodes";     sort -u \$PBS_NODEFILE | wc -l
echo "=================================="
H
}

mkrun () {  # $1 rundir  $2 binary  $3 input-name
  mkdir -p "$1"/{out2d,out3d,restart}
  ln -sf "$2" "$1/$(basename "$3" .in)"
  cp -p "$BENCH/inputs/$3" "$1/$3"
}

############### Phase A : flag screening, node64, A/B/C/D in ONE allocation
for t in 1 2 3; do
  d="$RUNS/phaseA/t$t"
  for v in A B C D; do mkrun "$d/$v" "$BIN/calliope-$v" calliope.in; done
  { hdr 64 72 00:40:00
    echo "for v in A B C D; do"
    echo "  cd \${PBS_O_WORKDIR}/\$v"
    echo "  echo \"##### variant \$v  \$(date)\""
    echo "  mpirun ./calliope <calliope.in >calliope.out.std"
    echo "  cp calliope.out.std $BENCH/results-phaseA/t$t-\$v"
    echo "done"
  } > "$d/job.pbs"
done

############### Phase B : production
# calliope: binary decided after Phase A -> placeholder symlink 'calliope' set later
for n in 16 32 64 128; do
  case $n in 16) w=00:30:00;; 32) w=00:25:00;; *) w=00:20:00;; esac
  for t in 1 2 3; do
    d="$RUNS/calliope/node$n-$t"
    mkrun "$d" "$BIN/calliope-PICK" calliope.in
    { hdr "$n" 72 "$w"
      echo "cd \${PBS_O_WORKDIR}"
      echo "mpirun ./calliope <calliope.in >calliope.out.std"
      echo "cp calliope.out.std '$C/results-202607/1024^3/node$n-$t'"
    } > "$d/job.pbs"
  done
done

# gallope
for tag in sep uni; do
  case $tag in sep) out="$G/results-202607/1024^3";; uni) out="$G/results-202607/1024^3-unified";; esac
  for n in 16 32 64 128 256; do
    for t in 1 2 3; do
      d="$RUNS/gallope-$tag/node$n-$t"
      mkrun "$d" "$BIN/gallope-$tag" gallope.in
      { hdr "$n" 1 00:15:00
        echo "cd \${PBS_O_WORKDIR}"
        echo "mpirun ./gallope <gallope.in >gallope.out.std"
        echo "cp gallope.out.std '$out/node$n-$t'"
      } > "$d/job.pbs"
    done
  done
done

echo "=== generated:"
find "$RUNS" -name job.pbs | wc -l
echo "=== example (phaseA/t1):"; cat "$RUNS/phaseA/t1/job.pbs"
echo "=== example (gallope-uni/node256-1):"; cat "$RUNS/gallope-uni/node256-1/job.pbs"
