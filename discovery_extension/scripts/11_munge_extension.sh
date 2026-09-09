#!/usr/bin/env bash
# Munging is extension-only; it never rewrites core harmonized or munged files.
set -euo pipefail
cd "$(dirname "$0")/../.."
source scripts/_common.sh

[ "$#" -gt 0 ] || die "usage: bash discovery_extension/scripts/11_munge_extension.sh --all|TRAIT [TRAIT ...]"
python3 discovery_extension/scripts/08_validate_source_contracts.py
require_ldsc

if [ "$1" = "--all" ]; then
  [ "$#" -eq 1 ] || die "--all cannot be combined with trait IDs"
  traits=$(awk -F'\t' 'NR > 1 {print $1}' discovery_extension/config/candidate_traits.tsv)
else
  traits="$*"
fi

mkdir -p discovery_extension/data/munged
for trait in $traits; do
  awk -F'\t' -v wanted="$trait" 'NR > 1 && $1 == wanted {found=1} END {exit !found}' \
    discovery_extension/config/candidate_traits.tsv || die "$trait is not in the locked extension panel"
  harmonized="discovery_extension/data/harmonized/$trait.txt.gz"
  require_file "$harmonized"
  echo "==> extension munge $trait"
  "$LDSC_PYTHON" "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "$harmonized" \
    --merge-alleles ref/w_hm3.snplist \
    --snp SNP --a1 A1 --a2 A2 --frq FRQ --p P --N-col N \
    --signed-sumstats BETA,0 --chunksize 500000 \
    --out "discovery_extension/data/munged/$trait"
  require_file "discovery_extension/data/munged/$trait.sumstats.gz"
  require_file "discovery_extension/data/munged/$trait.log"
done

echo "Extension munging complete; inspect every extension log before h2."
