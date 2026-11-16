#!/bin/bash
#$ -cwd
#$ -j y
#
# Calliope's judgments of the MHD_COMP_ISOTH folders on the Gallope runs, and the comparison with the s13d
# outputs of Calliope, as one Grid Engine job on CPUs (the analyses read tens of GB of out3d/ and restart/;
# a login node stops a process after 10 minutes). The commands and what they judge are listed in
# validation_tsubame.md.
#
# Run this script in RUNROOT, do not qsub it: outside a job it submits itself (as job-tsubame.sh does):
#
#   cd $RUNROOT
#   FOLDERS="Euler1D MHD1D" CALLIOPE_ROOT=~/work/calliope_1d924a9 PYTHON=~/work/venv-isoth/bin/python \
#       REF=/gs/bs/tgh-26IDS/ua07691/s13d $GALLOPE/tests/MHD_COMP_ISOTH/judge-tsubame.sh
#
#   FOLDERS        folders judged, of Euler1D MHD1D Smooth_wave OT2D Turb3D            [all five]
#   CALLIOPE_ROOT  Calliope repository; its tests/MHD_COMP_ISOTH/runs must be a          required
#                  symlink to RUNROOT (the analysis scripts find some runs through it)
#   PYTHON         python3 with numpy, scipy and matplotlib                             [python3]
#   REF            the s13d outputs of Calliope, REF/<folder>/<case> (MHD1D, OT2D,      [none: not compared]
#                  Turb3D); compare_s13d_all.py compares the cases that have one
#   GALLOPE        Gallope repository                                                   [../.. of this script]
#   RESOURCE       CPU resource type of the job                                         [cpu_16]
#   H_RT           wall-time limit of the job, H:MM:SS                                  [3:00:00]
#   TMAX           seconds after which a single analysis command is stopped             [5400]
#   NAME           job name                                                             [isoth_judge]
#   TSUBAME_GROUP  group charged                                                        [tgh-26IDS]
#
# Output: RUNROOT/judge/<folder>.log with every command, its output and its exit status, the tables
# judge/euler1d.md (judge_stage1.py) and judge/<folder>_s13d.md (compare_s13d_all.py table --md), and the files
# that the analysis scripts write under RUNROOT (judge.json, summary.json, ...). Never runs
# analyze_stage15.py watch, which rewrites exit_status.

FOLDERS=${FOLDERS:-Euler1D MHD1D Smooth_wave OT2D Turb3D}
GALLOPE=${GALLOPE:-$(cd "$(dirname "$0")/../.." && pwd)}
RESOURCE=${RESOURCE:-cpu_16}
H_RT=${H_RT:-3:00:00}
TMAX=${TMAX:-5400}
TSUBAME_GROUP=${TSUBAME_GROUP:-tgh-26IDS}
REF=${REF:-none}

if [ -z "${JOB_ID}" ]; then     # not under the scheduler: check the setting and submit
  bad=0
  FOLDERS=$(echo ${FOLDERS} | tr ': ' '  ')
  for f in ${FOLDERS}; do
    case ${f} in Euler1D|MHD1D|Smooth_wave|OT2D|Turb3D) ;; *) echo "judge-tsubame.sh: unknown folder ${f}" >&2; bad=1 ;; esac
  done
  S=${CALLIOPE_ROOT}/tests/MHD_COMP_ISOTH/stages
  [ -n "${CALLIOPE_ROOT}" ] && [ -f ${S}/analyze_stage4.py ] || { echo "judge-tsubame.sh: set CALLIOPE_ROOT" >&2; exit 1; }
  CALLIOPE_ROOT=$(cd "${CALLIOPE_ROOT}" && pwd)
  [ "$(cd ${CALLIOPE_ROOT}/tests/MHD_COMP_ISOTH/runs 2>/dev/null && pwd -P)" = "$(pwd -P)" ] || {
    echo "judge-tsubame.sh: ${CALLIOPE_ROOT}/tests/MHD_COMP_ISOTH/runs must link to RUNROOT (${PWD})" >&2; bad=1; }
  PYTHON=$(command -v "${PYTHON:-python3}")
  "${PYTHON}" -c 'import numpy, scipy, matplotlib' 2>/dev/null || {
    echo "PYTHON (${PYTHON:-python3}) cannot import numpy, scipy and matplotlib" >&2; bad=1; }
  [ "${REF}" = none ] || REF=$(cd "${REF}" && pwd) || bad=1
  for v in "${CALLIOPE_ROOT}" "${PYTHON}" "${REF}" "${GALLOPE}"; do
    case "${v}" in *,*|*" "*) echo "judge-tsubame.sh: qsub -v cannot pass ${v} (comma or space)" >&2; bad=1 ;; esac
  done
  [ ${bad} -eq 0 ] || exit 1
  exec qsub -g ${TSUBAME_GROUP} -l ${RESOURCE}=1 -l h_rt=${H_RT} -N ${NAME:-isoth_judge} -cwd -j y \
       -v FOLDERS=$(echo ${FOLDERS} | tr ' ' ':'),CALLIOPE_ROOT=${CALLIOPE_ROOT},PYTHON=${PYTHON},REF=${REF},GALLOPE=${GALLOPE},TMAX=${TMAX} "$0"
fi

set -u
R=$(pwd)
S=${CALLIOPE_ROOT}/tests/MHD_COMP_ISOTH/stages
G=${GALLOPE}/tests/MHD_COMP_ISOTH
PY=${PYTHON}
J=${R}/judge
mkdir -p ${J}
export MPLBACKEND=Agg OMP_NUM_THREADS=$(nproc)
cd ${CALLIOPE_ROOT}

# run CMD ...: the command and its output, then its exit status
run() { echo "===== $(date '+%F %T') $*"; timeout ${TMAX} "$@" 2>&1; echo "===== rc=$?"; }

# s13d FOLDER LIST [FILTER]: compare the cases of LIST (with a file under REF/FOLDER) with their s13d outputs
s13d() {
  local f=$1 list=$2
  [ "${REF}" = none ] && { echo "===== no REF: ${f} not compared with s13d"; return; }
  (cd ${R} && for c in $(cat ${list}); do [ -d ${REF}/${f}/${c} ] && echo ${c}; done) > ${J}/${f,,}_s13d.txt
  echo "===== $(wc -l < ${J}/${f,,}_s13d.txt) of $(wc -l < ${list}) cases have an s13d output"
  run ${PY} ${G}/compare_s13d_all.py run ${J}/${f,,}_s13d.txt --runroot ${R} --ref ${REF}/${f}
  ${PY} ${G}/compare_s13d_all.py table ${J}/${f,,}_s13d.txt --runroot ${R} --md > ${J}/${f,,}_s13d.md 2>&1
  tail -3 ${J}/${f,,}_s13d.md
}

judge() {
  local f=$1
  ${PY} ${G}/make_cases.py list --runroot ${R} --folder ${f} --status all --long > ${J}/${f,,}_cases.txt
  awk '{print $1}' ${J}/${f,,}_cases.txt > ${J}/${f,,}_list.txt
  echo "===== ${f}: $(wc -l < ${J}/${f,,}_list.txt) cases;" \
       "exit_status: $(awk '{$1 = $2 = $3 = $4 = ""; print}' ${J}/${f,,}_cases.txt | sort | uniq -c | xargs)"
  case ${f} in
    Euler1D)          # REPORT_stage1: shock, smoke, scan, decomp, selftest
      run ${PY} ${S}/judge_stage1.py ${R}/stage1_es > ${J}/euler1d.md
      grep -c FAIL ${J}/euler1d.md | sed 's/^/===== lines with FAIL in judge\/euler1d.md: /' ;;
    MHD1D)            # REPORT_stage14 D1, D4; REPORT_stage16 CV1; D2 and the rest against s13d
      run ${PY} ${G}/check_mhd1d.py --calliope ${CALLIOPE_ROOT} --runroot ${R}
      s13d MHD1D ${J}/mhd1d_list.txt ;;
    Smooth_wave)      # REPORT_stage4: Alfven waves, cp diagnostic, planar waves
      for c in alfven cpdiag planar; do run ${PY} ${S}/analyze_stage4.py ${c} ${R}/stage4; done ;;
    OT2D)             # REPORT_stage14 O1, O2 (analyze_stage4.py ot); O3, O4 against s13d
      run ${PY} ${S}/analyze_stage4.py ot ${R}/stage14
      if [ "${REF}" != none ]; then
        # analyze_stage14.py ot compares <root>/ot with the earlier runs <root>/old/ot (both with the summary.json
        # of analyze_stage4.py ot); here <root>/ot are the Gallope runs and <root>/old/ot the s13d runs
        mkdir -p ${J}/ot_s13d/old
        ln -sfn ${R}/stage14/ot ${J}/ot_s13d/ot
        ln -sfn ${REF}/OT2D/stage14/ot ${J}/ot_s13d/old/ot
        run ${PY} ${S}/analyze_stage14.py ot ${J}/ot_s13d
      fi
      s13d OT2D ${J}/ot2d_list.txt ;;
    Turb3D)           # REPORT_stage4 (turb), 5, 6, 7, 8, 13, 15, 16
      run ${PY} ${S}/analyze_stage4.py turb ${R}/stage4
      local seams
      seams=$(cd ${R} && ls -d stage4/turb/*/chunk_1 stage5/*/*/chunk_1 stage6q/*/*/chunk_1 stage8/*/*/chunk_1 \
              stage13d/*/*/chunk_1 stage15/*/*/chunk_1 stage16/hf*/*/chunk_1 2>/dev/null | sed 's|/chunk_1||')
      [ -n "${seams}" ] && run ${PY} ${S}/check_seams.py $(for c in ${seams}; do echo ${R}/${c}; done)
      mkdir -p ${R}/stage5/ana
      for c in res ms; do run ${PY} ${S}/analyze_stage5.py ${c} ${R} --json ${R}/stage5/ana/${c}.json; done
      for c in track turb; do run ${PY} ${S}/analyze_stage6.py ${c} ${R}/stage6q; done
      for N in 64 128 256; do for c in tail hole; do run ${PY} ${S}/analyze_stage7.py ${c} ${R} ${N}; done; done
      run ${PY} ${S}/analyze_stage7.py judge ${R}
      for c in rhomin tail hole judge; do run ${PY} ${S}/analyze_stage8.py ${c} ${R}; done
      for c in rhomin fields instab judge; do run ${PY} ${S}/analyze_stage13.py ${c} ${R}; done
      run ${PY} ${S}/analyze_stage15.py pre ${R}/stage15
      for b in a b c; do run ${PY} ${S}/analyze_stage15.py runs ${b} ${R}/stage15; done
      for c in spec judge; do run ${PY} ${S}/analyze_stage15.py ${c} ${R}/stage15; done
      for b in hf64 hf128 hf256; do run ${PY} ${S}/analyze_stage16.py hf ${b}; done
      for c in spec_hf rejudge judge1; do run ${PY} ${S}/analyze_stage16.py ${c}; done
      s13d Turb3D ${J}/turb3d_list.txt ;;
  esac
}

echo "=== $(date '+%F %T')  job ${JOB_ID}: ${FOLDERS//:/ } in ${R}, Calliope ${CALLIOPE_ROOT}" \
     "($(git -C ${CALLIOPE_ROOT} rev-parse --short HEAD 2>/dev/null || echo 'no git')), Gallope" \
     "$(git -C ${GALLOPE} rev-parse --short HEAD 2>/dev/null || echo 'no git'), REF ${REF}"
for f in ${FOLDERS//:/ }; do
  echo "=== $(date '+%T') ${f} -> ${J}/${f}.log"
  judge ${f} > ${J}/${f}.log 2>&1
  grep -E '^===== rc=[^0]|^===== [0-9]+ of' ${J}/${f}.log | sed "s|^|    ${f}: |"
done
echo "=== $(date '+%F %T') done"
