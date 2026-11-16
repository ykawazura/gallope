#!/bin/bash
#=======================================================================
# Gate G1: cuFFTMp with a size-2 dimension. Makes the run directories
# runs/isoth/g1/<nlx>x<nly>x<nlz>_np<P> and prints the CASES string for
# tests/run_cases-tsubame.sh. Run from the repository root after building
# bin/gallope-ISOTH_FFT (tests/MHD_COMP_ISOTH/Makefile.fft):
#
#   CASES=$(./tests/MHD_COMP_ISOTH/setup_g1.sh) RESOURCE=node_h NAME=g1 \
#     ./tests/run_cases-tsubame.sh
#   grep -H "PASS\|FAIL" runs/isoth/g1/*/out.std
#
# Only &box_parameters of input_example/MHD_COMP_ISOTH.in is changed.
#=======================================================================
set -e
SRC=input_example/MHD_COMP_ISOTH.in
cases=""
for n in "2 8 8" "8 2 8" "8 8 2" "6 6 6" "2 2 8" "2 8 2" "8 2 2" "16 32 64"; do
    read -r nx ny nz <<< "${n}"
    for p in 1 2; do
        d=runs/isoth/g1/${nx}x${ny}x${nz}_np${p}
        mkdir -p ${d}
        sed -e "s/^\( *nlx *=\).*/\1 ${nx}/" -e "s/^\( *nly *=\).*/\1 ${ny}/" \
            -e "s/^\( *nlz *=\).*/\1 ${nz}/" ${SRC} > ${d}/gallope.in
        cases="${cases:+${cases}+}${d}:${p}:ISOTH_FFT"
    done
done
echo "${cases}"
