#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

IMAGE=ghcr.io/precimed/gsa-mixer@sha256:5bf54ddd6f7f81b93eeb5450b1a6dc809d1fcf926b3a815d6d514f36ce04f51c
MIXER_PY=/tools/mixer/precimed/mixer.py
FIGURES_PY=/tools/mixer/precimed/mixer_figures.py
REF=ref/mixer/reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC
THREADS=16

usage() {
  echo "Usage:"
  echo "  bash scripts/38_run_mixer_task.sh univariate TRAIT REP"
  echo "  bash scripts/38_run_mixer_task.sh combine-univariate TRAIT"
  echo "  bash scripts/38_run_mixer_task.sh bivariate SLEEP_TRAIT NON_SLEEP_TRAIT REP"
  echo "  bash scripts/38_run_mixer_task.sh combine-bivariate SLEEP_TRAIT NON_SLEEP_TRAIT"
}

[ "$#" -ge 2 ] || { usage >&2; exit 2; }
MODE=$1
shift

python3 scripts/35_mixer_preflight.py >/dev/null
if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  echo "ERROR: pinned MiXeR image is absent; inspect then explicitly run scripts/37_pull_mixer_image.sh --pull" >&2
  exit 1
fi

run_container() {
  docker run --rm --platform linux/amd64 \
    -v "$ROOT:/home" -w /home "$IMAGE" python "$@"
}

validate_trait() {
  trait=$1
  domain=${2:-any}
  awk -F '\t' -v target="$trait" -v wanted="$domain" '
    NR == 1 { for (i=1; i<=NF; i++) { if ($i=="trait_id") t=i; if ($i=="domain") d=i } ; next }
    $t == target && (wanted=="any" || (wanted=="sleep" && $d=="sleep") || (wanted=="non_sleep" && $d!="sleep")) { found=1 }
    END { exit(found ? 0 : 1) }
  ' config/analysis_panel.tsv || { echo "ERROR: invalid $domain trait: $trait" >&2; exit 1; }
}

common=(--ld-file "${REF}.@.run4.ld" --bim-file "${REF}.@.bim" --threads "$THREADS" --exclude-ranges MHC)

case "$MODE" in
  univariate)
    [ "$#" -eq 2 ] || { usage >&2; exit 2; }
    trait=$1
    rep=$2
    validate_trait "$trait"
    case "$rep" in ''|*[!0-9]*) echo "ERROR: REP must be 1-20" >&2; exit 2 ;; esac
    [ "$rep" -ge 1 ] && [ "$rep" -le 20 ] || { echo "ERROR: REP must be 1-20" >&2; exit 2; }
    mkdir -p results/mixer/univariate
    prefix="results/mixer/univariate/${trait}"
    extract="${REF}.prune_maf0p05_rand2M_r2p8.rep${rep}.snps"
    seed=$((1000 + rep))
    run_container "$MIXER_PY" fit1 "${common[@]}" --extract "$extract" --seed "$seed" \
      --trait1-file "data/mixer/${trait}.sumstats.gz" --out "${prefix}.fit.rep${rep}"
    run_container "$MIXER_PY" test1 "${common[@]}" --seed "$seed" \
      --trait1-file "data/mixer/${trait}.sumstats.gz" \
      --load-params "${prefix}.fit.rep${rep}.json" --out "${prefix}.test.rep${rep}"
    ;;
  combine-univariate)
    [ "$#" -eq 1 ] || { usage >&2; exit 2; }
    trait=$1
    validate_trait "$trait"
    prefix="results/mixer/univariate/${trait}"
    for rep in $(seq 1 20); do
      [ -s "${prefix}.fit.rep${rep}.json" ] && [ -s "${prefix}.test.rep${rep}.json" ] || {
        echo "ERROR: incomplete 20-replicate univariate family for $trait" >&2
        exit 1
      }
    done
    run_container "$FIGURES_PY" combine --json "${prefix}.fit.rep@.json" --out "${prefix}.fit"
    run_container "$FIGURES_PY" combine --json "${prefix}.test.rep@.json" --out "${prefix}.test"
    run_container "$FIGURES_PY" one --json "${prefix}.fit.json" --out "${prefix}.fit.summary" \
      --trait1 "$trait" --statistic mean std --ext svg
    ;;
  bivariate)
    [ "$#" -eq 3 ] || { usage >&2; exit 2; }
    sleep=$1
    disease=$2
    rep=$3
    validate_trait "$sleep" sleep
    validate_trait "$disease" non_sleep
    case "$rep" in ''|*[!0-9]*) echo "ERROR: REP must be 1-20" >&2; exit 2 ;; esac
    [ "$rep" -ge 1 ] && [ "$rep" -le 20 ] || { echo "ERROR: REP must be 1-20" >&2; exit 2; }
    python3 scripts/39_collate_mixer.py --univariate-only
    python3 - "$sleep" "$disease" <<'PY'
import csv, sys
rows = list(csv.DictReader(open("results/tables/mixer_univariate.tsv"), delimiter="\t"))
status = {row["trait_id"]: row["bivariate_eligibility"] for row in rows}
missing = [trait for trait in sys.argv[1:] if status.get(trait) != "ELIGIBLE"]
if missing:
    raise SystemExit("ERROR: bivariate MiXeR requires eligible univariate models: " + ", ".join(missing))
PY
    mkdir -p results/mixer/bivariate
    prefix="results/mixer/bivariate/${sleep}_vs_${disease}"
    first="results/mixer/univariate/${sleep}"
    second="results/mixer/univariate/${disease}"
    extract="${REF}.prune_maf0p05_rand2M_r2p8.rep${rep}.snps"
    seed=$((1000 + rep))
    run_container "$MIXER_PY" fit2 "${common[@]}" --extract "$extract" --seed "$seed" \
      --trait1-file "data/mixer/${sleep}.sumstats.gz" --trait2-file "data/mixer/${disease}.sumstats.gz" \
      --trait1-params "${first}.fit.rep${rep}.json" --trait2-params "${second}.fit.rep${rep}.json" \
      --out "${prefix}.fit.rep${rep}"
    run_container "$MIXER_PY" test2 "${common[@]}" --seed "$seed" \
      --trait1-file "data/mixer/${sleep}.sumstats.gz" --trait2-file "data/mixer/${disease}.sumstats.gz" \
      --load-params "${prefix}.fit.rep${rep}.json" --out "${prefix}.test.rep${rep}"
    ;;
  combine-bivariate)
    [ "$#" -eq 2 ] || { usage >&2; exit 2; }
    sleep=$1
    disease=$2
    validate_trait "$sleep" sleep
    validate_trait "$disease" non_sleep
    prefix="results/mixer/bivariate/${sleep}_vs_${disease}"
    for rep in $(seq 1 20); do
      [ -s "${prefix}.fit.rep${rep}.json" ] && [ -s "${prefix}.test.rep${rep}.json" ] || {
        echo "ERROR: incomplete 20-replicate bivariate family for ${sleep}/${disease}" >&2
        exit 1
      }
    done
    run_container "$FIGURES_PY" combine --json "${prefix}.fit.rep@.json" --out "${prefix}.fit"
    run_container "$FIGURES_PY" combine --json "${prefix}.test.rep@.json" --out "${prefix}.test"
    run_container "$FIGURES_PY" two --json-fit "${prefix}.fit.json" --json-test "${prefix}.test.json" \
      --out "$prefix" --trait1 "$sleep" --trait2 "$disease" --statistic mean std --ext svg
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
