#!/usr/bin/env bash
## launch.sh -- run ON klebsiella. Fans replicates across independent chunks.
##
##   bash launch.sh STUDY [NCHUNK] [REP_PER_CHUNK] [EXTRA_ENV...]
##
##   bash launch.sh coverage 12 10
##     -> 120 replicates of study_coverage.R across 12 processes
##
##   bash launch.sh selection 16 6 NSUBJ=400 SHAPE=4
##     -> 96 replicates of study_selection.R, extra settings passed through
##
## Each chunk gets its own TAG, so its log and its results cannot collide with
## another's, and seeds are disjoint by construction: chunk i covers
## SEED0 = BASE + (i-1)*REP, so no seed is ever used twice. Every replicate
## writes its own rds, so a chunk that dies costs one replicate and a rerun
## skips whatever is already on disk.
##
## saemix-copula is single-threaded pure R, so one chunk is one core and there
## is no compile cache to thrash -- no stagger is needed, unlike the harness
## this came from. Keep NCHUNK below the free core count; the box has 64 and
## other work usually holds some.
set -eu
STUDY=${1:?usage: launch.sh STUDY [NCHUNK] [REP] [KEY=VALUE...]}
NCHUNK=${2:-12}; REP=${3:-10}
shift 3 2>/dev/null || shift $#
BASE=${SEED_BASE:-1000}
RDIR=~/saemix_copula_study
SCRIPT="study_${STUDY}.R"

cd "$RDIR"
[ -f "$SCRIPT" ] || { echo "no such study: $SCRIPT"; exit 1; }
mkdir -p out "out_${STUDY}"

free=$(( $(nproc) - $(cut -d' ' -f1 /proc/loadavg | cut -d. -f1) ))
echo "cores $(nproc), roughly $free free, launching $NCHUNK chunks"
[ "$NCHUNK" -gt "$free" ] && echo "  WARNING: asking for more chunks than free cores"

echo "study $STUDY: $NCHUNK chunks x $REP reps = $((NCHUNK*REP)) replicates"
for i in $(seq 1 "$NCHUNK"); do
  TAG=$(printf "c%02d" "$i")
  SEED0=$(( BASE + (i-1)*REP ))
  env PKG="$RDIR/lib" REPS="$REP" SEED0="$SEED0" TAG="$TAG" \
      OUTDIR="$RDIR/out_${STUDY}" "$@" \
      nohup Rscript "$SCRIPT" > "out/${STUDY}_${TAG}.log" 2>&1 &
  echo "  $TAG  seeds $((SEED0+1))..$((SEED0+REP))  pid $!"
done
echo "started $(date +%Y-%m-%d\ %H:%M)"
echo "watch with:  bash watch.sh $STUDY"
