#!/bin/sh
#=====================================================================
#  Run compare.py over every pair produced by job.pbs-gallope and
#  job.pbs-calliope, and collect the output in one file.
#
#  Usage:  ./compare_all.sh GALLOPE_DIR CALLIOPE_DIR [OUT]
#
#  GALLOPE_DIR   directory holding run_<series>_n<nstep>_g<ngpu>/
#  CALLIOPE_DIR  directory holding run_<series>_n<nstep>_c<np>/
#  OUT           output file (default compare_all.out)
#
#  Two comparisons are made for each (series, nstep):
#    gallope 4 GPU  vs  calliope        the code-to-code difference
#    gallope 1 GPU  vs  gallope 4 GPU   the summation-order control
#  The second must be small compared with the first, otherwise the
#  code-to-code difference is just a decomposition artefact.
#=====================================================================
set -u
G=${1:?GALLOPE_DIR}
C=${2:?CALLIOPE_DIR}
OUT=${3:-compare_all.out}
HERE=$(cd "$(dirname "$0")" && pwd)
NP=16
GRID="--nlx 128 --nly 128 --nlz 128 --lx 1.0 --ly 1.0 --lz 1.0"

: > "${OUT}"
for series in noshear shear; do
  # The solenoidal check is on kx + q*tsc*ky, so q must match the run.
  case ${series} in
    noshear) steps="1 10 100 1000"      ; Q="--q 0.0" ;;
    shear)   steps="1 10 100 1000 2000" ; Q="--q 1.5" ;;
  esac
  for n in ${steps}; do
    g4="${G}/run_${series}_n${n}_g4"
    g1="${G}/run_${series}_n${n}_g1"
    cc="${C}/run_${series}_n${n}_c${NP}"

    if [ -d "${cc}/restart" ] && [ -d "${g4}/restart" ]; then
      { echo "############ ${series}  nstep=${n}  gallope(4 GPU) vs calliope"
        python3 "${HERE}/compare.py" --a "${g4}" --a-code gallope \
                                     --b "${cc}" --b-code calliope ${GRID} ${Q}
        echo; } >> "${OUT}" 2>&1
    else
      echo "############ ${series}  nstep=${n}  SKIPPED (missing ${cc} or ${g4})" >> "${OUT}"
    fi

    if [ -d "${g1}/restart" ] && [ -d "${g4}/restart" ]; then
      { echo "############ ${series}  nstep=${n}  gallope 1 GPU vs 4 GPU [control]"
        python3 "${HERE}/compare.py" --a "${g1}" --a-code gallope \
                                     --b "${g4}" --b-code gallope ${GRID} ${Q}
        echo; } >> "${OUT}" 2>&1
    else
      echo "############ ${series}  nstep=${n}  control SKIPPED" >> "${OUT}"
    fi
  done
done

echo "wrote ${OUT}"
grep -c '^############' "${OUT}"
grep -n 'SKIPPED\|MISMATCH\|Traceback\|SystemExit' "${OUT}" || echo "no skips or gate failures"
