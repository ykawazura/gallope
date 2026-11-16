#!/bin/bash
#$ -cwd
#$ -j y
#
# Runs the MHD_COMP_ISOTH cases of Calliope with Gallope on TSUBAME4.0: the port of
# tests/MHD_COMP_ISOTH/stages/job.pbs-ps-fes-dev of Calliope (PBS, CPUs) to Altair Grid Engine and GPUs.
#
# Runs every case listed in $CASEFILE (run directories relative to the directory the script is started from,
# e.g. the output of "make_cases.py list" started in the run root). Each directory must contain calliope.in
# and may contain a file "ngpu" holding the number of GPUs (= MPI ranks, default 1; make_cases.py setup writes
# it; nproc, the ranks of Calliope, is not used). The binary is $GALLOPE/bin/gallope-MHD_COMP_ISOTH.
#
# Run this script, do not qsub it: outside a job it submits itself (the group of "qsub -g" is accepted only on
# the command line, and qsub -v cannot pass a list with spaces or commas, hence CASEFILE):
#
#   cd $RUNROOT
#   make_cases.py list --runroot . --folder OT2D --status todo > ot2d.txt
#   CASEFILE=ot2d.txt RESOURCE=node_q H_RT=2:00:00 PYTHON=~/work/venv-isoth/bin/python \
#       $GALLOPE/tests/MHD_COMP_ISOTH/job-tsubame.sh
#
#   CASEFILE       file with the cases (whitespace separated)                          required
#   RESOURCE       node_q (1 GPU) | node_h (2 GPUs) | node_f (4 GPUs), one node          [node_q]
#   H_RT           wall-time limit of the job, H:MM:SS                                  [0:30:00]
#   BUDGET         minutes after the start of the job at which the runs stop            [H_RT - 3]
#   OVERSUB        cases that share a GPU                                               [1]
#                  (keep 1: without MPS, which users cannot start, cases that share a GPU are
#                  time-sliced; 4 small 1D cases on one GPU ran 8.4x slower in total, see README)
#   GALLOPE        Gallope repository (bin/, batch/, tests/MHD_COMP_ISOTH/)             [../.. of this script]
#   PYTHON         python3 with numpy, for gallope_layout.py                            [python3]
#   NAME           job name                                                             [isoth]
#   TSUBAME_GROUP  group charged                                                        [tgh-26IDS]
#
# The cases run concurrently on the GPUs of the node: each case starts, in the order listed (make_cases.py
# list lists them by decreasing cost), as soon as ngpu GPUs each run fewer than OVERSUB cases (the least used
# GPUs are taken). Its ranks see only these GPUs (CUDA_VISIBLE_DEVICES; rank r on the r-th of them, through
# batch/gpu-bind-tsubame.sh) and run unbound to cores. The GPUs of the job are those of CUDA_VISIBLE_DEVICES
# set by the scheduler, or 0 .. (GPUs of RESOURCE) - 1. A case with ngpu above the GPUs of the job is skipped.
#
# As in the PBS job: a run stops BUDGET minutes after the start of the job: the input calliope_run.in is
# calliope.in with max_wall_time set to the time left, the run saves its restart files when it is reached, and
# exit_status reads "incomplete" (also after a stop by the file "terminate", see &terminate).  Such a case
# continues from its restart files while at least 3 minutes are left, and when the same cases are submitted
# again: its earlier output (including out1d, out2d, out3d, *.out.nc, *.modes.out, layout.log) goes to
# chunk_<k>/ (the analysis scripts of Calliope concatenate chunk_1, chunk_2, ... and the case directory), and
# calliope_run.in has init_type = 'restart' and nstep = the steps left.  The steps left are computed from dt
# of calliope.in, so the chain is valid only for a fixed dt: a run whose dt was changed by reset_dt gets
# exit_status "dt_changed" (or "dt_stop" when reset_method = 'stop' stopped it) and is not continued.
# A case with a file "t_end" (one number) runs up to that time with the dt that reset_dt leaves (e.g.
# reset_method = 'decrement'): a continuation takes dt of the last step from the last row of cfl.dat (written
# at every step in full precision; its t + dt must equal t of restart/time.dat) and nstep =
# int((t_end - t)/dt + 0.5), a run that ends without error before t_end - dt/2 (its dt was reduced) is
# "incomplete", and a change of dt does not give "dt_changed".  When calliope.in also sets t_end
# (&time_parameters, e.g. with reset_method = 'track', which raises dt as well), the run stops at t_end
# itself and a continuation keeps nstep of calliope.in (an upper bound) in place of the steps left.
# A case with any other exit_status is skipped; remove exit_status (and chunk_*) to run it from the start.
# A case with a directory ic/ (initial condition in the restart format, from make_cases.py ic) starts
# from it: ic/* is copied to restart/ when the case runs from the start.
#
# Differences from the PBS job besides the scheduler and the GPUs:
#   calliope_run.in  is "translate calliope.in" (below) with the same edits of max_wall_time, init_type, nstep
#                    and dt. translate drops &mpi_settings and every "pruned = ..." (Gallope has neither: it
#                    dealiases by the 2/3 rule only, so a pruned = T case runs as pruned = F), and sets
#                    restart_dir = './restart_g/' (comments are kept).
#   restart files    Gallope reads and writes restart_g/ in its own layout. When a case with ic/ runs from the
#                    start, ic/ is converted into restart_g/ (gallope_layout.py to-gallope); after every run
#                    restart_g/ is converted into restart/ (to-calliope), which the chain, status_of and the
#                    analysis read, so restart/ always holds the Calliope layout.
#   out3d/           is converted in place to the Calliope layout after every run (gallope_layout.py
#                    out3d-to-calliope, after exit_status is written). A conversion that the end of the job
#                    interrupted is finished when the case is submitted again (make_cases.py list --status todo
#                    lists it). The conversions log to layout.log; a failed one gives exit_status
#                    "exit layout: ...".
#   exe, fftw_estimate, nproc of a case are not used (make_cases.py setup notes the first two in calliope_extra).
#
# DRYRUN=1 runs on any machine without submitting, loading modules or changing a file: for each case it prints
# the GPUs it would get, what it would do and the difference of calliope_run.in from calliope.in.
#   DRYRUN=1 CASEFILE=list.txt RESOURCE=node_f bash job-tsubame.sh
# The translation alone (without the edits of the job) is printed by
#   bash job-tsubame.sh translate calliope.in

# value of "name = value" in calliope.in (Fortran d exponent as e)
nml() { sed -n "s/^ *$1 *= *\([^ ]*\).*/\1/p" $2 | head -1 | sed 's/[dD]/e/'; }

# nlx nly nlz of calliope.in
grid() { echo $(nml nlx $1) $(nml nly $1) $(nml nlz $1); }

# the Gallope input of the Calliope input $1, on stdout (see the header; one group per line range "&group" .. "/")
translate() {
  awk '
  BEGIN { q = "\047"; rdv = q "./restart_g/" q
          rre = "restart_dir[ \t]*=[ \t]*(" q "[^" q "]*" q "|\"[^\"]*\"|[^ \t,]*)" }
  {
    line = $0; code = line; com = ""; inq = ""
    for (i = 1; i <= length(line); i++) {          # the code ends at the first ! outside quotes
      ch = substr(line, i, 1)
      if (inq != "") { if (ch == inq) inq = "" }
      else if (ch == q || ch == "\"") inq = ch
      else if (ch == "!") { code = substr(line, 1, i - 1); com = substr(line, i); break }
    }
    lc = tolower(code)
    if (grp == "" && lc ~ /^[ \t]*&[a-z_0-9]+/) {
      grp = lc; sub(/^[ \t]*&/, "", grp); sub(/[^a-z_0-9].*$/, "", grp)
      if (grp == "restart_parameters") { nrp++; krd = 0 }
      if (grp == "terminate") { nterm++; kmwt = 0 }
      if (grp != "mpi_settings") print line
      next
    }
    if (grp != "" && lc ~ /^[ \t]*\//) {
      if (grp == "restart_parameters" && !krd) print "  restart_dir = " rdv
      if (grp == "terminate" && !kmwt) print "  max_wall_time = 0.0d0"
      if (grp != "mpi_settings") print line
      grp = ""
      next
    }
    if (grp == "mpi_settings") next
    n = 0
    while (match(lc, /(^|[ \t,])pruned[ \t]*=[ \t]*[^ \t,\/]+[ \t]*,?/)) {
      code = substr(code, 1, RSTART - 1) " " substr(code, RSTART + RLENGTH)
      lc = substr(lc, 1, RSTART - 1) " " substr(lc, RSTART + RLENGTH)
      n++
    }
    if (n && code ~ /^[ \t]*$/ && com == "") next
    if (grp == "restart_parameters" && match(lc, rre)) {
      code = substr(code, 1, RSTART - 1) "restart_dir = " rdv substr(code, RSTART + RLENGTH)
      krd = 1; n++
    }
    if (grp == "terminate" && lc ~ /(^|[ \t,])max_wall_time[ \t]*=/) kmwt = 1
    print (n ? code com : line)
  }
  END {
    if (!nrp) { print "&restart_parameters"; print "  restart_dir = " rdv; print "/" }
    if (!nterm) { print "&terminate"; print "  max_wall_time = 0.0d0"; print "/" }
  }' "$1"
}

if [ "$1" = translate ]; then
  [ -f "$2" ] || { echo "job-tsubame.sh translate: no file $2" >&2; exit 1; }
  translate "$2"
  exit
fi

TSUBAME_GROUP=${TSUBAME_GROUP:-tgh-26IDS}
RESOURCE=${RESOURCE:-node_q}
H_RT=${H_RT:-0:30:00}
OVERSUB=${OVERSUB:-1}
GALLOPE=${GALLOPE:-$(cd "$(dirname "$0")/../.." && pwd)}
case ${RESOURCE} in
  node_f) GPN=4 ;;
  node_h) GPN=2 ;;
  node_q) GPN=1 ;;
  *) echo "job-tsubame.sh: RESOURCE must be node_q, node_h or node_f (one node)" >&2; exit 1 ;;
esac
# minutes of H_RT ([[H:]M:]S)
BUDGET=${BUDGET:-$(( $(echo ${H_RT} | awk -F: '{s = 0; for (i = 1; i <= NF; i++) s = 60*s + $i; printf "%d", s/60}') - 3 ))}

if [ -z "${JOB_ID}" ] && [ -z "${DRYRUN}" ]; then     # not under the scheduler: check the cases and submit
  [ -n "${CASEFILE}" ] && [ -f "${CASEFILE}" ] || { echo "job-tsubame.sh: set CASEFILE to a file with the cases" >&2; exit 1; }
  CASEFILE=$(cd "$(dirname "${CASEFILE}")" && pwd)/$(basename "${CASEFILE}")
  GALLOPE=$(cd "${GALLOPE}" && pwd)
  PYTHON=$(command -v "${PYTHON:-python3}")
  bad=0
  [ -x ${GALLOPE}/bin/gallope-MHD_COMP_ISOTH ] || { echo "no ${GALLOPE}/bin/gallope-MHD_COMP_ISOTH" >&2; bad=1; }
  "${PYTHON}" -c 'import numpy' 2>/dev/null || { echo "PYTHON (${PYTHON:-python3}) cannot import numpy" >&2; bad=1; }
  [ ${BUDGET} -ge 4 ] || { echo "BUDGET = ${BUDGET} min: no run would start (3 min must be left)" >&2; bad=1; }
  for c in $(cat ${CASEFILE}); do
    if [ ! -f ${c}/calliope.in ]; then echo "no ${c}/calliope.in (the cases are relative to ${PWD})" >&2; bad=1; continue; fi
    [ -s ${c}/exit_status ] && continue
    # the inputs that make_cases.py ic writes (profile1d.dat, ic/)
    for f in $(awk '{print $2}' ${c}/generated.md5 2>/dev/null); do
      [ -f ${c}/${f} ] || { echo "no ${c}/${f}: run make_cases.py ic" >&2; bad=1; break; }
    done
  done
  for v in "${CASEFILE}" "${GALLOPE}" "${PYTHON}"; do
    case "${v}" in *,*|*" "*) echo "job-tsubame.sh: qsub -v cannot pass ${v} (comma or space)" >&2; bad=1 ;; esac
  done
  [ ${bad} -eq 0 ] || exit 1
  exec qsub -g ${TSUBAME_GROUP} -l ${RESOURCE}=1 -l h_rt=${H_RT} -N ${NAME:-isoth} -cwd -j y \
       -v CASEFILE=${CASEFILE},GALLOPE=${GALLOPE},PYTHON=${PYTHON},RESOURCE=${RESOURCE},OVERSUB=${OVERSUB},BUDGET=${BUDGET} "$0"
fi

[ -z "${DRYRUN}" ] && . ${GALLOPE}/batch/env-tsubame.sh
export OMP_NUM_THREADS=1
EXE=${GALLOPE}/bin/gallope-MHD_COMP_ISOTH
BIND=${GALLOPE}/batch/gpu-bind-tsubame.sh
MPIRUN="mpirun --bind-to none -x LD_LIBRARY_PATH -x CUDA_VISIBLE_DEVICES --mca coll_hcoll_enable 0"
PY=${PYTHON:-python3}
LAYOUT=${GALLOPE}/tests/MHD_COMP_ISOTH/gallope_layout.py
T_END=$(( $(date +%s) + ${BUDGET:-27} * 60 ))

[ -n "${CASES}" ] || [ -f "${CASEFILE}" ] || { echo "job-tsubame.sh: set CASEFILE (or CASES with DRYRUN)" >&2; exit 1; }
CASES=${CASES:-$(cat ${CASEFILE})}

# GNU time if it is there (not BSD time, which has no -f), else the same line from date
if /usr/bin/time -f "wall %e s" true > /dev/null 2>&1; then
  timed() { /usr/bin/time -f "wall %e s" "$@"; }
else
  timed() { local s=$(date +%s) r; "$@"; r=$?; echo "wall $(( $(date +%s) - s )) s" >&2; return ${r}; }
fi

# exit_status of the run in the current directory ($1: exit code of mpirun, $2: t_end or empty)
status_of() {
  if grep -q 'This run is stopped' calliope.out.std; then echo dt_stop
  elif [ -z "$2" ] && grep -q -e 'dt is decreased' -e 'dt is increased' calliope.out.std; then echo dt_changed
  elif grep -q -e 'Wall time has exceeded' -e 'Termination file has been detected' calliope.out.std; then
    echo incomplete
  elif [ $1 -eq 0 ] && [ -n "$2" ] && awk -v t=$(sed -n 2p restart/time.dat | awk '{print $1}') \
         -v d=$(tail -1 cfl.dat | awk '{print $2}') -v e=$2 'BEGIN{exit !(t < e - 0.5*d)}'; then
    echo incomplete
  else echo "exit $1"; fi
}

# sed script "sub" of calliope_run.in of case $1 with $2 s left; with $3 (cfl.dat of the last run) the
# continuation from restart/ of an "incomplete" case (tend: t_end of the case or empty). Sets err on failure.
make_sub() {
  local c=$1 tt t1 d1 rest nleft
  err=
  sub="s/^\( *max_wall_time *=\).*/\1 $(awk -v s=$2 'BEGIN{printf "%.4f", (s - 60)/3600}')d0/"
  [ -z "$3" ] && return
  tt=$(sed -n 2p ${c}/restart/time.dat | awk '{print $1}')
  if [ -n "${tend}" ]; then
    # t and dt of the last step, which ended at tt
    read t1 d1 rest < <(tail -1 $3)
    if ! awk -v t=${tt} -v r=${t1} -v d=${d1} 'BEGIN{e = t - (r + d); exit !((e < 0 ? -e : e) <= 1e-12*t)}'; then
      err="exit t_end: t + dt of the last row of cfl.dat is not t of restart/time.dat"
      return
    fi
    if [ -n "$(nml t_end ${c}/calliope.in)" ]; then
      nleft=$(nml nstep ${c}/calliope.in)
    else
      nleft=$(awk -v t=${tt} -v d=${d1} -v e=${tend} 'BEGIN{printf "%d", (e - t)/d + 0.5}')
    fi
    sub="${sub};s/^\( *dt *=\).*/\1 ${d1}/"
  else
    nleft=$(awk -v t=${tt} -v d=$(nml dt ${c}/calliope.in) -v n=$(nml nstep ${c}/calliope.in) \
            'BEGIN{printf "%d", n - int(t/d + 0.5)}')
  fi
  sub="${sub};s/^\( *init_type *=\).*/\1 'restart'/;s/^\( *nstep *=\).*/\1 ${nleft}/"
}

# true if out3d/ of case $1 holds a file of the Gallope layout (not converted since the last run)
out3d_pending() {
  local f
  for f in rho mx my mz bx by bz; do
    [ -f $1/out3d/${f}.dat.gallope ] && return 0
    [ -f $1/out3d/${f}.dat ] && ! grep -qx ${f} $1/out3d/.calliope_layout 2>/dev/null && return 0
  done
  return 1
}

# run gallope_layout.py "$@" in the current directory, logged to layout.log
layout() {
  echo "=== $(date '+%F %T')  gallope_layout.py $*" >> layout.log
  ${PY} ${LAYOUT} "$@" >> layout.log 2>&1
}

# finish the conversion of out3d/ of case $1 that the end of an earlier job interrupted (its run had ended)
finish_layout() {
  local c=$1 st
  st=$(cat ${c}/exit_status 2>/dev/null)
  case "${st}" in ''|"exit layout"*) return ;; esac
  out3d_pending ${c} || return
  if [ -n "${DRYRUN}" ]; then
    echo "=== ${c}  (dry run) would finish: ${PY} ${LAYOUT} out3d-to-calliope out3d --n $(grid ${c}/calliope.in)"
    return
  fi
  ( cd ${c} && layout out3d-to-calliope out3d --n $(grid calliope.in) ) ||
    echo "exit layout: out3d-to-calliope out3d failed (run: ${st}; see layout.log)" > ${c}/exit_status
  echo "=== ${c}  out3d/ conversion finished: $(cat ${c}/exit_status)"
}

run_once() {
  c=$1
  np=$(cat ${c}/ngpu 2>/dev/null || echo 1)
  st=$(cat ${c}/exit_status 2>/dev/null)
  tend=$(cat ${c}/t_end 2>/dev/null)
  n3="$(grid ${c}/calliope.in)"
  left=$(( T_END - $(date +%s) ))
  if [ ${left} -lt 180 ]; then
    echo "=== ${c}  not started (${left} s left)"
    return
  fi
  if [ "${st}" = "incomplete" ]; then
    k=1; while [ -d ${c}/chunk_${k} ]; do k=$(( k + 1 )); done
    mkdir ${c}/chunk_${k}
    ( cd ${c} && for f in conservation.dat budget.dat cfl.dat calliope.out.std calliope_run.in exit_status \
                          out1d out2d out3d *.out.nc *.modes.out layout.log; do [ -e "${f}" ] && mv "${f}" chunk_${k}/; done )
    make_sub ${c} ${left} ${c}/chunk_${k}/cfl.dat
    if [ -n "${err}" ]; then
      echo "${err}" > ${c}/exit_status
      return
    fi
  elif [ -d ${c}/chunk_1 ]; then
    echo "exit chunk_1 exists" > ${c}/exit_status
    return
  else
    make_sub ${c} ${left}
    if [ -d ${c}/ic ]; then
      mkdir -p ${c}/restart && cp ${c}/ic/* ${c}/restart/
      if ! ( cd ${c} && layout to-gallope ic restart_g --n ${n3} ); then
        echo "exit layout: to-gallope ic restart_g failed (see layout.log)" > ${c}/exit_status
        echo "=== ${c}  $(cat ${c}/exit_status)"
        return
      fi
    fi
  fi
  translate ${c}/calliope.in | sed -e "${sub}" > ${c}/calliope_run.in
  t0="$(date '+%F %T')"
  # out3d/ of this run starts in the Gallope layout: no file of an earlier conversion may be taken for it
  ( cd ${c} && mkdir -p out1d out2d out3d restart restart_g && \
    rm -f out3d/.calliope_layout out3d/*.dat.gallope out3d/*.dat.tmp && \
    timed ${MPIRUN} -n ${np} ${BIND} ${EXE} calliope_run > calliope.out.std 2>&1 ; \
    rc=$?
    # restart/ (Calliope layout) from restart_g/, before status_of reads restart/time.dat
    if [ -f restart_g/time.dat ] && ! layout to-calliope restart_g restart --n ${n3}; then
      echo "exit layout: to-calliope restart_g restart failed (run: exit ${rc}; see layout.log)" > exit_status
      exit
    fi
    status_of ${rc} "${tend}" > exit_status
    if out3d_pending . && ! layout out3d-to-calliope out3d --n ${n3}; then
      echo "exit layout: out3d-to-calliope out3d failed (run: $(cat exit_status); see layout.log)" > exit_status
    fi )
  echo "=== ${c}  (ngpu = ${np}, GPUs ${CUDA_VISIBLE_DEVICES}, $(basename ${EXE}), $(cat ${c}/exit_status), chunk $(( $(ls -d ${c}/chunk_* 2>/dev/null | wc -l) + 1 )))  ${t0}
$(tail -3 ${c}/calliope.out.std)"
}

run_case() {
  finish_layout $1
  st=$(cat $1/exit_status 2>/dev/null)
  while [ -z "${st}" ] || [ "${st}" = "incomplete" ]; do
    [ $(( T_END - $(date +%s) )) -lt 180 ] && return
    run_once $1
    st=$(cat $1/exit_status 2>/dev/null)
  done
}

# what run_once would do for case $1 on the GPUs $2, without changing a file
plan_once() {
  c=$1
  np=$(cat ${c}/ngpu 2>/dev/null || echo 1)
  st=$(cat ${c}/exit_status 2>/dev/null)
  tend=$(cat ${c}/t_end 2>/dev/null)
  n3="$(grid ${c}/calliope.in)"
  left=$(( T_END - $(date +%s) ))
  before=
  if [ "${st}" = "incomplete" ]; then
    k=1; while [ -d ${c}/chunk_${k} ]; do k=$(( k + 1 )); done
    make_sub ${c} ${left} ${c}/cfl.dat
    what="continues from restart/ at t = $(sed -n 2p ${c}/restart/time.dat | awk '{print $1}'), the last output to chunk_${k}/"
  elif [ -d ${c}/chunk_1 ]; then
    err="exit chunk_1 exists"
  else
    make_sub ${c} ${left}
    what="runs from the start${tend:+ up to t_end = ${tend}}"
    [ -d ${c}/ic ] && before="cp ic/* restart/; ${PY} ${LAYOUT} to-gallope ic restart_g --n ${n3}"
  fi
  if [ -n "${err}" ]; then
    echo "=== ${c}  (dry run) would stop: ${err}"
    return
  fi
  echo "=== ${c}  (dry run: ngpu = ${np}, CUDA_VISIBLE_DEVICES=$2, ${what})"
  [ -n "${before}" ] && echo "    before: ${before}"
  echo "    run:    ${MPIRUN} -n ${np} ${BIND} ${EXE} calliope_run"
  echo "    after:  ${PY} ${LAYOUT} to-calliope restart_g restart --n ${n3}"
  echo "            ${PY} ${LAYOUT} out3d-to-calliope out3d --n ${n3}"
  echo "    calliope_run.in = translate calliope.in | sed -e \"${sub}\"; diff calliope.in calliope_run.in:"
  diff ${c}/calliope.in <(translate ${c}/calliope.in | sed -e "${sub}") | sed 's/^/      /'
}

# GPUs of the job: gid[0 .. NG-1]; used[i]: cases running on gid[i]; gpusof[pid]: positions of the GPUs of a case
NG=0
for g in $(echo ${CUDA_VISIBLE_DEVICES} | tr , ' '); do gid[NG]=${g}; NG=$(( NG + 1 )); done
if [ ${NG} -eq 0 ]; then
  while [ ${NG} -lt ${GPN} ]; do gid[NG]=${NG}; NG=$(( NG + 1 )); done
fi
i=0; while [ ${i} -lt ${NG} ]; do used[i]=0; i=$(( i + 1 )); done
echo "=== $(date '+%F %T')  ${RESOURCE}, GPUs ${gid[*]}, OVERSUB = ${OVERSUB}, BUDGET = ${BUDGET} min, ${EXE}${DRYRUN:+ (dry run)}"

# positions of $1 GPUs that run fewer than OVERSUB cases, the least used first; empty if they do not fit
pick() {
  local u=0 i n=0 out=
  while [ ${u} -lt ${OVERSUB} ] && [ ${n} -lt $1 ]; do
    i=0
    while [ ${i} -lt ${NG} ] && [ ${n} -lt $1 ]; do
      [ ${used[i]} -eq ${u} ] && { out="${out} ${i}"; n=$(( n + 1 )); }
      i=$(( i + 1 ))
    done
    u=$(( u + 1 ))
  done
  [ ${n} -eq $1 ] && echo $(echo ${out} | tr ' ' '\n' | sort -n)
}

for c in ${CASES}; do
  np=$(cat ${c}/ngpu 2>/dev/null || echo 1)
  st=$(cat ${c}/exit_status 2>/dev/null)
  if [ -n "${st}" ] && [ "${st}" != "incomplete" ]; then
    echo "=== ${c}  skipped: ${st}"
    # a conversion of out3d/ that an earlier job left unfinished runs on the CPUs, beside the cases on the GPUs
    if [ -n "${DRYRUN}" ]; then finish_layout ${c}; else finish_layout ${c} & fi
    continue
  fi
  if [ ${np} -gt ${NG} ]; then
    echo "=== ${c}  skipped: ngpu = ${np}, the job has ${NG} GPUs"
    continue
  fi
  while :; do
    if [ -z "${DRYRUN}" ]; then
      i=0; while [ ${i} -lt ${NG} ]; do used[i]=0; i=$(( i + 1 )); done
      for p in $(jobs -rp); do
        for i in ${gpusof[$p]}; do used[i]=$(( used[i] + 1 )); done
      done
    fi
    sel=$(pick ${np})
    [ -n "${sel}" ] && break
    if [ -n "${DRYRUN}" ]; then
      echo "    (dry run) the next case waits for running cases to end; all GPUs are taken as free again"
      i=0; while [ ${i} -lt ${NG} ]; do used[i]=0; i=$(( i + 1 )); done
      continue
    fi
    # not wait -n, which misses a case that ended before it was called but after jobs -rp listed it
    sleep 2
  done
  ids=
  for i in ${sel}; do ids=${ids:+${ids},}${gid[i]}; done
  if [ -n "${DRYRUN}" ]; then
    plan_once ${c} ${ids}
    for i in ${sel}; do used[i]=$(( used[i] + 1 )); done
  else
    ( export CUDA_VISIBLE_DEVICES=${ids}; run_case ${c} ) &
    gpusof[$!]="${sel}"
  fi
done
wait
echo "=== all done $(date '+%F %T')"
