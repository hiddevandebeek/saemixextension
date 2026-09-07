#!/usr/bin/env bash
## deploy.sh -- put saemix-copula and the study scripts on klebsiella.
##
##   NM_PASS=... bash deploy.sh
##
## Tars the WORKING TREE rather than `git archive HEAD`, because the changes
## that matter are usually uncommitted -- on this package the whole point of a
## deploy is to run something that has just been edited.
##
## Unlike the admixr2 harness this came from, saemix-copula has no src/: it is
## pure R. So there is no C++ to rebuild, nothing to exclude to stop a Linux
## link picking up Windows objects, and no compile cache to thrash. Install is
## just untar and R CMD INSTALL.
##
## Excluded: .git, and the whole of copula/. That directory is 1.8 GB of
## experiment results, examples and design notes, and R CMD INSTALL does not
## look at it -- it is not a package directory. The study scripts that live
## under copula/cluster are copied separately below. Including it made the
## first deploy a 1.8 GB upload of things the cluster has no use for.
set -eu
: "${NM_PASS:?set NM_PASS first}"
HOST='klebsiella.lacdr.leidenuniv.nl'; USER_='beekh'; RDIR='saemix_copula_study'
SRC="${SAEMIX_SRC:-C:/package/saemix-copula}"
HERE=$(cd "$(dirname "$0")" && pwd)

AP=$(mktemp); printf '#!/bin/sh\necho "$NM_PASS"\n' > "$AP"; chmod +x "$AP"
export SSH_ASKPASS="$AP" SSH_ASKPASS_REQUIRE=force DISPLAY=:0
trap 'rm -f "$AP"' EXIT
SSHO="-o StrictHostKeyChecking=accept-new -o ConnectTimeout=25 \
-o PreferredAuthentications=password -o PubkeyAuthentication=no \
-o NumberOfPasswordPrompts=1 -o BatchMode=no"

echo ">> packing working tree"
TAR=$(mktemp -u).tar.gz
tar czf "$TAR" -C "$SRC" \
  --exclude='.git' --exclude='.claude' --exclude='.Rproj.user' \
  --exclude='copula' --exclude='testsaemix4' \
  .
echo "   $(du -h "$TAR" | cut -f1)"

echo ">> uploading"
ssh $SSHO "$USER_@$HOST" "mkdir -p ~/$RDIR/pkg ~/$RDIR/out" </dev/null 2>&1 \
  | grep -vE "Permanently added" || true
scp $SSHO "$TAR" "$USER_@$HOST:~/$RDIR/saemix-copula.tar.gz" </dev/null 2>&1 \
  | grep -vE "Permanently added" || true
## The study scripts, plus the experiment's own functions.R -- the margin
## selection study is that file's `combined_run_replicate`, and copula/ is
## excluded from the tarball above, so it has to travel separately.
scp $SSHO "$HERE"/study_*.R "$HERE"/combine_*.R "$HERE"/launch.sh \
  "$SRC/copula/experiments/combined-natural-frem-study/functions.R" \
  "$SRC/copula/experiments/combined-natural-frem-study/functions_testfirst.R" \
  "$SRC/copula/experiments/combined-natural-frem-study/summarize.R" \
  "$USER_@$HOST:~/$RDIR/" </dev/null 2>&1 | grep -vE "Permanently added" || true
rm -f "$TAR"

echo ">> installing into a private library"
## A private library, not the system one: this is a fork of a package that may
## already be installed for other work, and silently replacing it would be
## rude and hard to notice.
ssh $SSHO "$USER_@$HOST" "set -eu
  cd ~/$RDIR
  mkdir -p lib
  rm -rf pkg && mkdir -p pkg && tar xzf saemix-copula.tar.gz -C pkg
  R CMD INSTALL --library=\$HOME/$RDIR/lib pkg 2>&1 | tail -15
" </dev/null 2>&1 | grep -vE "Permanently added"

echo ">> verifying the study entry points exist"
ssh $SSHO "$USER_@$HOST" "Rscript -e '
  suppressMessages(library(saemix, lib.loc = \"~/$RDIR/lib\"))
  need <- c(\"copulaPopulation\",\"copulaGaussianRvineFromCor\",
            \"copulaNaturalMarginLognormal\",\"copulaNaturalMarginGeneralizedGamma\",
            \"copulaStandardErrors\",\"copulaGet\",\"copulaPosteriorEtaDraws\",
            \"copulaNaturalPosteriorData\",\"copulaNaturalMarginStart\",
            \"copulaNaturalMarginFamilies\")
  miss <- need[!vapply(need, exists, TRUE, where = asNamespace(\"saemix\"))]
  cat(if (length(miss)) paste(\"MISSING:\", paste(miss, collapse = \", \")) else
      \"all study entry points present\", \"\n\")
  cat(\"version:\", as.character(packageVersion(\"saemix\")), \"\n\")'" \
  </dev/null 2>&1 | grep -vE "Permanently added"

echo ">> deploy done"
