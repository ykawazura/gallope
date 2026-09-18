#!/bin/bash
set -euo pipefail
module load nvidia nv-hpcx hdf5 netcdf netcdf-fortran

ROOT=/work/gr96/o07001
BENCH=$ROOT/bench-202607
BIN=$BENCH/bin
PV=$ROOT/p3dfft_vs_cufftmp
FFTW=/work/opt/local/aarch64/apps/nvidia/25.9/fftw/3.3.10
WORK=$BENCH/fftbuild

[ -d "$WORK" ] && rm -r "$WORK"
mkdir -p "$WORK"

################ cuFFTMp slab / pencil ################
for kind in slab pencil; do
  case $kind in
    slab)   src=$PV/cufftmp-slab/cufftmp_r2c.f90;            base=cufftmp_r2c ;;
    pencil) src=$PV/cufftmp-pencil/cufftmp_r2c_c2r_pencils.f90; base=cufftmp_r2c_c2r_pencils ;;
  esac
  for N in 1024 2048; do
    d=$WORK/$kind-$N; mkdir -p "$d"
    sed -E "s/^([[:space:]]*n[xyz][[:space:]]*=[[:space:]]*)2048[[:space:]]*$/\1$N/" "$src" > "$d/$base.f90"
    echo "--- $kind-$N params:"; grep -nE '^[[:space:]]*(nx|ny|nz|nrepeat)[[:space:]]*=' "$d/$base.f90"
    cp -p "$PV/cufftmp-$kind/Makefile" "$d/Makefile"
    ( cd "$d" && make > "$BENCH/logs/build_fft-$kind-$N.log" 2>&1 )
    cp -p "$d/$base" "$BIN/fft-$kind-$N"
  done
done

################ P3DFFT driver (two library variants) ################
for LIBTAG in g O3; do
  case $LIBTAG in
    g)  P3D=$ROOT/.local/p3dfft.2    ; FOPT="-g" ;;
    O3) P3D=$ROOT/.local/p3dfft.2-O3 ; FOPT="-O3 -fast" ;;
  esac
  for N in 1024 2048; do
    d=$WORK/p3dfft-$N-$LIBTAG; mkdir -p "$d"
    sed -E "s/^([[:space:]]*n[xyz][[:space:]]+=[[:space:]]*)2048[[:space:]]*$/\1$N/" "$PV/p3dfft/driver_rand.F90" > "$d/driver_rand.F90"
    echo "--- p3dfft-$N-$LIBTAG params:"; grep -nE '^[[:space:]]*(nx|ny|nz|ndim|n)[[:space:]]+=' "$d/driver_rand.F90" | head -6
    ( cd "$d" && \
      mpif90 -DHAVE_CONFIG_H -I. -I"$P3D/include" -DMEASURE -DSTRIDE1 -DFFTW $FOPT -c driver_rand.F90 -o driver_rand.o && \
      mpif90 $FOPT -o test_rand_f.x driver_rand.o "$P3D/lib/libp3dfft.a" "$FFTW/lib/libfftw3.a" ) \
      > "$BENCH/logs/build_fft-p3dfft-$N-$LIBTAG.log" 2>&1
    cp -p "$d/test_rand_f.x" "$BIN/fft-p3dfft-$N-$LIBTAG"
  done
done

echo "=== BUILT"
ls -l "$BIN" | grep '^-.*fft-'
