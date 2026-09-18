#!/bin/bash
# Wait for the jobs the (killed) buggy driver already put in PBS to finish
# writing their result files, THEN run the fixed campaign driver.  Starting the
# fixed driver while those jobs are still in flight would double-submit them,
# because ok() keys on the destination file, which an in-flight job hasn't
# written yet.
BENCH=/work/gr96/o07001/bench-202607
LOG=$BENCH/logs/campaign.log
echo "$(date '+%m-%d %H:%M:%S') RECOVER: waiting for pre-existing queue to drain" >> "$LOG"
while :; do
  nq=$(qstat 2>/dev/null | awk '$1 ~ /^[0-9]+$/' | wc -l)
  [ "$nq" -eq 0 ] && break
  sleep 120
done
echo "$(date '+%m-%d %H:%M:%S') RECOVER: queue empty, launching fixed driver" >> "$LOG"
exec bash "$BENCH/scripts/campaign.sh" "$BENCH/lists/all.txt"
