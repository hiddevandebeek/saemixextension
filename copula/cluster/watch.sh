#!/usr/bin/env bash
## watch.sh -- poll klebsiella for study progress. Needs NM_PASS.
##
##   NM_PASS=... bash watch.sh [STUDY]
##
## Two traps inherited from the harness this came from, both of which produce a
## confident and wrong "everything died":
##
##   R runs as `--file=study_coverage.R`, so `pgrep -f "Rscript study"` matches
##   nothing. Match on `file=` instead.
##
##   There are THREE states, not two. A chunk that has not finished its first
##   replicate has written no rds at all, which is not the same as having
##   failed. Report started, finished and failed separately and never infer
##   death from an absence.
set -u
STUDY=${1:-coverage}
: "${NM_PASS:?set NM_PASS}"
HOST='klebsiella.lacdr.leidenuniv.nl'; USER_='beekh'; RDIR='saemix_copula_study'
AP=$(mktemp); printf '#!/bin/sh\necho "$NM_PASS"\n' > "$AP"; chmod +x "$AP"
export SSH_ASKPASS="$AP" SSH_ASKPASS_REQUIRE=force DISPLAY=:0
trap 'rm -f "$AP"' EXIT
SSHO="-o StrictHostKeyChecking=accept-new -o ConnectTimeout=25 \
-o PreferredAuthentications=password -o PubkeyAuthentication=no \
-o NumberOfPasswordPrompts=1 -o BatchMode=no"

ssh $SSHO "$USER_@$HOST" "cd ~/$RDIR 2>/dev/null || exit 0
  live=\$(pgrep -fc 'file=study_${STUDY}.R' || echo 0)
  done_=\$(ls out_${STUDY}/*.rds 2>/dev/null | wc -l)
  part=\$(ls out_${STUDY}/*.part 2>/dev/null | wc -l)
  echo \"processes alive: \$live   replicates written: \$done_   in flight: \$part\"
  echo '--- per chunk, last line ---'
  for f in out/${STUDY}_c*.log; do
    [ -e \"\$f\" ] || continue
    printf '%-28s %s\n' \"\$(basename \$f)\" \"\$(tail -1 \"\$f\")\"
  done
  echo '--- failed replicates ---'
  grep -h 'FAILED' out/${STUDY}_c*.log 2>/dev/null | head -10 || true
  n=\$(grep -hc 'FAILED' out/${STUDY}_c*.log 2>/dev/null | paste -sd+ | bc 2>/dev/null || echo 0)
  echo \"(total failures: \${n:-0})\"
  echo '--- R errors, if any ---'
  grep -l -iE 'error|halted' out/${STUDY}_c*.log 2>/dev/null || echo '(none)'
" </dev/null 2>&1 | grep -vE "Permanently added"
