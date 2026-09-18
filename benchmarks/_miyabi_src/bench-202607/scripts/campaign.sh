#!/bin/bash
# Submit every run directory in $1, wait for the queue to drain, then resubmit
# any run whose destination files are missing / empty / lack a success marker.
# Repeat up to $MAXROUND times.  Node-hours are cheap compared with losing a
# data point: at 64 nodes roughly 1 in 4 runs dies at startup with
# "Bus error / invalid address alignment", so retrying is mandatory.

LIST=${1:?usage: campaign.sh <listfile>}
BENCH=/work/gr96/o07001/bench-202607
LOG=$BENCH/logs/campaign.log
MAXQ=28          # queue ACCEPT limit is 32
MAXROUND=6
mkdir -p "$BENCH/logs"

say () { echo "$(date '+%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }

nq () { qstat 2>/dev/null | awk '$1 ~ /^[0-9]+$/' | wc -l; }

# destinations this job is supposed to produce
dests () { grep -oE "^cp [^ ]+ '[^']+'" "$1/job.pbs" | sed "s/.*'\(.*\)'\$/\1/"; }

ok () {                       # ok <rundir> -> 0 if every destination is good
  # Completion markers, printed only at the very end of a full run:
  #   Calliope / Gallope stdout : "# of steps advanced   100"
  #   FFT micro-benchmark stdout: "cpu time per loop ="
  # A startup crash (bus error ~16 s in) never reaches either line.
  local d f n=0
  while read -r f; do
    n=$((n+1))
    [ -s "$f" ] || return 1
    grep -qE '# of steps advanced|cpu time per loop' "$f" || return 1
  done < <(dests "$1")
  [ "$n" -gt 0 ]
}

for round in $(seq 1 $MAXROUND); do
  todo=()
  while read -r d; do
    [ -n "$d" ] || continue
    ok "$d" || todo+=("$d")
  done < "$LIST"

  if [ ${#todo[@]} -eq 0 ]; then
    say "ROUND $round: nothing left to do -- campaign complete"
    exit 0
  fi
  say "ROUND $round: submitting ${#todo[@]} job(s)"

  for d in "${todo[@]}"; do
    while [ "$(nq)" -ge $MAXQ ]; do sleep 60; done
    out=$( cd "$d" && qsub job.pbs 2>&1 )
    say "  qsub $out  <- ${d#$BENCH/runs/}"
    sleep 3
  done

  say "ROUND $round: all submitted, waiting for the queue to drain"
  while [ "$(nq)" -gt 0 ]; do sleep 120; done
  say "ROUND $round: queue empty"
done

say "reached MAXROUND=$MAXROUND; remaining failures need manual attention"
