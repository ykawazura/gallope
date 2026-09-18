#!/bin/bash
# Kill calliope/gallope jobs that hang BEFORE mpirun ever starts.
#
# Observed failure (2026-07-25, calliope node128 = 9216 ranks): the job sits
# RUNNING but "<rundir>/*.out.std" is never created, i.e. the shell never
# reached the `mpirun ... >*.out.std` line.  Left alone it holds 128 nodes and a
# RUN slot until the 20-min walltime kill, starving the queued node256 jobs.
#
# Liveness signal that does NOT depend on Fortran stdout buffering:
#   a healthy attempt reaches mpirun within seconds and the `>` redirect
#   truncates/creates *.out.std, so its mtime >= the job's start time.
#   If, THRESH seconds after start, *.out.std is absent OR older than start
#   (a stale file from a previous attempt), mpirun never launched -> qdel.
# Only absence/staleness triggers a kill; an existing empty file (output merely
# buffered) is left alone.  fft jobs have many sub-runs -> left to walltime.
BENCH=/work/gr96/o07001/bench-202607
CLOG=$BENCH/logs/campaign.log
WLOG=$BENCH/logs/watchdog.log
THRESH=360        # 6 min; healthy init prints within ~1.3 min even at node128
say(){ echo "$(date '+%m-%d %H:%M:%S') $*" >> "$WLOG"; }
say "watchdog start THRESH=${THRESH}s"
idle=0
while :; do
  running=$(qstat 2>/dev/null | awk '$1 ~ /^[0-9]+$/ {print $1}')
  if [ -z "$running" ]; then
    idle=$((idle+1))
    if grep -qE 'campaign complete|need manual attention' "$CLOG" 2>/dev/null && [ $idle -ge 3 ]; then
      say "queue empty and campaign finished -- exit"; break
    fi
    sleep 120; continue
  fi
  idle=0
  now=$(date +%s)
  for j in $running; do
    line=$(grep -E "qsub ${j}\.[^ ]+  <- " "$CLOG" 2>/dev/null | tail -1)
    [ -n "$line" ] || continue
    rel=${line##*<- }
    case "$rel" in
      calliope/*)  f=$BENCH/runs/$rel/calliope.out.std ;;
      gallope-*/*) f=$BENCH/runs/$rel/gallope.out.std ;;
      *) continue ;;
    esac
    st=$(qstat -f "$j" 2>/dev/null | awk -F'= ' '/ stime = /{print $2}')
    [ -n "$st" ] || continue
    sts=$(date -d "$st" +%s 2>/dev/null) || continue
    age=$((now - sts))
    [ "$age" -gt "$THRESH" ] || continue
    fmt=$(stat -c %Y "$f" 2>/dev/null || echo 0)
    if [ "${fmt:-0}" -lt "$sts" ]; then
      qdel "$j" 2>/dev/null && say "QDEL $j age=${age}s stdout-not-written <- $rel"
    fi
  done
  sleep 60
done
