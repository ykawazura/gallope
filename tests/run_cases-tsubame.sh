#!/bin/bash
#=======================================================================
# Run several gallope cases one after another in a single TSUBAME4.0 job.
#
#   CASES="runs/a:1:RMHD+runs/b:4:RMHD" RESOURCE=node_f ./tests/run_cases-tsubame.sh
#
# Submit from the repository root. Each case is  <dir>:<nranks>:<MODEL>  and
# runs inside <dir> as
#   mpirun -n <nranks> bin/gallope-<MODEL>   (reads <dir>/gallope.in)
# writing <dir>/out.std. Cases are separated by '+', because qsub -v splits
# its argument at commas. Binaries are expected at bin/gallope-<MODEL>.
#
#   RESOURCE  node_f (4 GPUs) | node_h (2) | node_q (1)       [node_q]
#   NODES     number of nodes                                [1]
#   H_RT      wall-time limit                                [0:10:00]
#   NAME      job name                                       [gtest]
#=======================================================================
#$ -cwd
#$ -j y

TSUBAME_GROUP=tgh-26IDS
RESOURCE=${RESOURCE:-node_q}
NODES=${NODES:-1}
H_RT=${H_RT:-0:10:00}

if [ -z "${JOB_ID}" ]; then     # not under AGE yet -> submit ourselves
    exec qsub -g ${TSUBAME_GROUP} -l ${RESOURCE}=${NODES} -l h_rt=${H_RT} \
              -N ${NAME:-gtest} -v "CASES=${CASES},RESOURCE=${RESOURCE},NODES=${NODES}" "$0"
fi

case ${RESOURCE} in
    node_f) GPN=4 ;;
    node_h) GPN=2 ;;
    node_q) GPN=1 ;;
esac

. ./batch/env-tsubame.sh
ROOT=${PWD}

for c in ${CASES//+/ }; do
    IFS=: read -r dir np model <<< "${c}"
    ppn=$(( np < GPN ? np : GPN ))
    ( cd "${dir}" || exit 1
      mkdir -p restart out2d out3d
      t0=$(date +%s)
      mpirun -n ${np} -npernode ${ppn} -x LD_LIBRARY_PATH --mca coll_hcoll_enable 0 \
          ${ROOT}/batch/gpu-bind-tsubame.sh ${ROOT}/bin/gallope-${model} > out.std 2>&1
      rc=$?
      echo "[case] ${dir} np=${np} ${model} rc=${rc} wall=$(( $(date +%s) - t0 ))s" )
done
echo "CASES_DONE"
