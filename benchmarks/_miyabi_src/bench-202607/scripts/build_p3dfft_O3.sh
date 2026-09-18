#!/bin/bash
set -euo pipefail
module load nvidia nv-hpcx hdf5 netcdf netcdf-fortran

ROOT=/work/gr96/o07001
BENCH=$ROOT/bench-202607
FFTW=/work/opt/local/aarch64/apps/nvidia/25.9/fftw/3.3.10
PREFIX=$ROOT/.local/p3dfft.2-O3

echo "=== nvfortran version"; nvfortran --version
echo "=== cleaning"
[ -d "$BENCH/p3dfft-src" ] && rm -r "$BENCH/p3dfft-src"
cp -a "$ROOT/p3dfft_vs_cufftmp/p3dfft/src" "$BENCH/p3dfft-src"
cd "$BENCH/p3dfft-src"
make distclean >/dev/null 2>&1 || true

echo "=== configure"
./configure \
  --enable-openmp --enable-stride1 --enable-fftw \
  --with-fftw="$FFTW" \
  FC=mpif90 CC=mpicc \
  FCFLAGS="-O3 -fast" CFLAGS="-O3 -fast" \
  --prefix="$PREFIX"

echo "=== effective flags"
grep -nE '^(FCFLAGS|CFLAGS|FDFLAGS) =' build/Makefile

echo "=== make"
make
echo "=== install"
[ -d "$PREFIX" ] && rm -r "$PREFIX"
make install

echo "=== SIMD check (new)"
printf 'vector-FP: '; objdump -d "$PREFIX/lib/libp3dfft.a" | grep -cE '^[[:space:]]+[0-9a-f]+:.*[[:space:]](fmla|fmul|fadd|fsub|fdiv)[[:space:]]+v[0-9]+\.' || true
printf 'total    : '; objdump -d "$PREFIX/lib/libp3dfft.a" | grep -cE '^[[:space:]]+[0-9a-f]+:' || true
echo "=== SIMD check (old, -g)"
printf 'vector-FP: '; objdump -d "$ROOT/.local/p3dfft.2/lib/libp3dfft.a" | grep -cE '^[[:space:]]+[0-9a-f]+:.*[[:space:]](fmla|fmul|fadd|fsub|fdiv)[[:space:]]+v[0-9]+\.' || true
printf 'total    : '; objdump -d "$ROOT/.local/p3dfft.2/lib/libp3dfft.a" | grep -cE '^[[:space:]]+[0-9a-f]+:' || true
echo "=== DONE"
