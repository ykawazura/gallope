#!/bin/bash
# usage: submit_driver.sh <list-file>   (one run-directory per line)
BENCH=/work/gr96/o07001/bench-202607
LOG=$BENCH/logs/submit.log
MAXQ=28
while read -r d; do
  [ -z "$d" ] && continue
  [ -f "$d/job.pbs" ] || { echo "MISSING $d/job.pbs" >> "$LOG"; continue; }
  while :; do
    nq=$(qstat 2>/dev/null | awk '$1 ~ /^[0-9]+$/' | wc -l)
    [ "$nq" -lt "$MAXQ" ] && break
    sleep 60
  done
  out=$( cd "$d" && qsub job.pbs 2>&1 )
  echo "$(date +%m-%d\ %H:%M:%S)  $out  $d" >> "$LOG"
  sleep 3
done < "$1"
echo "$(date) === ALL SUBMITTED from $1" >> "$LOG"
