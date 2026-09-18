#!/bin/bash
set -euo pipefail
module load nvidia nv-hpcx hdf5 netcdf netcdf-fortran

ROOT=/work/gr96/o07001
BENCH=$ROOT/bench-202607
BIN=$BENCH/bin
C=$ROOT/calliope_performance
G=$ROOT/gallope_performance

mkdir -p "$BIN"

simd () {  # $1 = binary
  local v t
  v=$(objdump -d "$1" | grep -cE '^[[:space:]]+[0-9a-f]+:.*[[:space:]](fmla|fmul|fadd|fsub|fdiv)[[:space:]]+v[0-9]+\.' || true)
  t=$(objdump -d "$1" | grep -cE '^[[:space:]]+[0-9a-f]+:' || true)
  echo "    SIMD-FP=$v  total=$t"
}

############################  CALLIOPE  ############################
cp -p "$C/arch/miyabi.in" "$BENCH/miyabi.in.calliope.orig"

CHEAD=$(sed -n '1,12p' "$BENCH/miyabi.in.calliope.orig" | grep -v '^P3DFFT_HOME')

build_calliope () {   # $1 tag  $2 p3dfft_home  $3 f90flags  $4 ldflags
  echo "=================== calliope-$1"
  { echo "FC = mpif90"
    echo ""
    echo "P3DFFT_HOME=$2"
    echo 'P3DFFT_INC=-I$(P3DFFT_HOME)/include/'
    echo 'P3DFFT_LIB=-L$(P3DFFT_HOME)/lib/ -lp3dfft'
    echo ""
    echo "FFTW3_HOME=/work/opt/local/aarch64/apps/nvidia/25.9/fftw/3.3.10"
    echo 'FFTW3_INC=-I$(FFTW3_HOME)/include/'
    echo 'FFTW3_LIB=-L$(FFTW3_HOME)/lib/ -lfftw3'
    echo ""
    echo 'NETCDF_INC=-I$(NETCDF_DIR)/include -I$(NETCDF_FORTRAN_DIR)/include'
    echo 'NETCDF_LIB=-L$(NETCDF_DIR)/lib -L$(NETCDF_FORTRAN_DIR)/lib -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz'
    echo ""
    echo "F90FLAGS = $3"
    echo "LDFLAGS = $4"
  } > "$C/arch/miyabi.in"
  cp -p "$C/arch/miyabi.in" "$BENCH/bin/arch-calliope-$1.in"
  ( cd "$C" && make clean >/dev/null 2>&1 && make > "$BENCH/logs/build_calliope-$1.log" 2>&1 )
  cp -p "$C/calliope" "$BIN/calliope-$1"
  simd "$BIN/calliope-$1"
  grep -icE 'error|warning' "$BENCH/logs/build_calliope-$1.log" | sed "s|^|    err+warn lines: |" || true
}

build_calliope A '$(HOME)/work/.local/p3dfft.2'    '-Mpreprocess -Mextend'                    ''
build_calliope B '$(HOME)/work/.local/p3dfft.2'    '-O3 -fast -Mpreprocess -Mextend'          '-O3 -fast'
build_calliope C '$(HOME)/work/.local/p3dfft.2-O3' '-O3 -fast -Mpreprocess -Mextend'          '-O3 -fast'
build_calliope D '$(HOME)/work/.local/p3dfft.2-O3' '-O3 -fast -Mnocache_align -Mpreprocess -Mextend' '-O3 -fast -Mnocache_align'

cp -p "$BENCH/miyabi.in.calliope.orig" "$C/arch/miyabi.in"
echo "calliope arch restored"

############################  GALLOPE  ############################
cp -p "$G/arch/miyabi.in" "$BENCH/miyabi.in.gallope.orig"

build_gallope () {   # $1 tag  $2 extra line (may be empty)
  echo "=================== gallope-$1"
  { echo "F90 = mpif90"
    echo ""
    echo 'NETCDF_INC=-I$(NETCDF_DIR)/include -I$(NETCDF_FORTRAN_DIR)/include'
    echo 'NETCDF_LIB=-L$(NETCDF_DIR)/lib -L$(NETCDF_FORTRAN_DIR)/lib -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz'
    echo ""
    echo "FLAGS  = -O3 -Mfree -fast -Mextend -Mpreprocess -Minform=warn -Minfo=accel -cuda -cudalib=cufftmp"
    [ -n "$2" ] && echo "$2"
    echo 'LINKER = -L${NVHPC_ROOT}/compilers/lib -lnvhpcwrapcufft -L${NVHPC_ROOT}/math_libs/lib64/ -lcufftMp -lmpi'
  } > "$G/arch/miyabi.in"
  cp -p "$G/arch/miyabi.in" "$BENCH/bin/arch-gallope-$1.in"
  ( cd "$G" && make clean >/dev/null 2>&1 && make > "$BENCH/logs/build_gallope-$1.log" 2>&1 )
  cp -p "$G/gallope" "$BIN/gallope-$1"
  ls -l "$BIN/gallope-$1"
  grep -c 'mem:unified' "$BENCH/bin/arch-gallope-$1.in" | sed "s|^|    unified-line: |" || true
}

build_gallope sep ''
build_gallope uni 'FLAGS  += -gpu=mem:unified'

cp -p "$BENCH/miyabi.in.gallope.orig" "$G/arch/miyabi.in"
echo "gallope arch restored"

echo "=================== BUILD SUMMARY"
ls -l "$BIN"
echo "=== ALL DONE"
