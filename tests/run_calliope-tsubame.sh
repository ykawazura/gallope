#!/bin/bash
#=======================================================================
# Run Calliope reference cases for the Gallope <-> Calliope comparison (V8)
# in a single TSUBAME4.0 CPU job.
#
#   CASES="runs/c_a:16+runs/c_b:16" CALLIOPE_ROOT=$HOME/work/calliope_develop \
#     RESOURCE=cpu_16 ./tests/run_calliope-tsubame.sh
#
# Submit from the Gallope repository root. Each case is <dir>:<nranks> and
# runs inside <dir> as
#   mpirun -n <nranks> ${CALLIOPE_ROOT}/calliope calliope   (reads calliope.in)
# writing <dir>/out.std. The environment is Calliope's own
# ${CALLIOPE_ROOT}/batch/env-tsubame.sh (gcc + openmpi), not Gallope's.
#
#   RESOURCE       cpu_16 | cpu_40 | cpu_80 | cpu_160          [cpu_16]
#   H_RT           wall-time limit                            [0:10:00]
#   NAME           job name                                   [ctest]
#=======================================================================
#$ -cwd
#$ -j y

TSUBAME_GROUP=tgh-26IDS
RESOURCE=${RESOURCE:-cpu_16}
H_RT=${H_RT:-0:10:00}
CALLIOPE_ROOT=${CALLIOPE_ROOT:-$HOME/work/calliope_develop}

if [ -z "${JOB_ID}" ]; then     # not under AGE yet -> submit ourselves
    exec qsub -g ${TSUBAME_GROUP} -l ${RESOURCE}=1 -l h_rt=${H_RT} -N ${NAME:-ctest} \
              -v "CASES=${CASES},CALLIOPE_ROOT=${CALLIOPE_ROOT}" "$0"
fi

. ${CALLIOPE_ROOT}/batch/env-tsubame.sh
export OMP_NUM_THREADS=1

for c in ${CASES//+/ }; do
    IFS=: read -r dir np <<< "${c}"
    ( cd "${dir}" || exit 1
      mkdir -p restart out2d out3d
      t0=$(date +%s)
      mpirun -n ${np} --bind-to core --mca coll_hcoll_enable 0 \
          ${CALLIOPE_ROOT}/calliope calliope > out.std 2>&1
      rc=$?
      echo "[case] ${dir} np=${np} calliope rc=${rc} wall=$(( $(date +%s) - t0 ))s" )
done
echo "CASES_DONE"
