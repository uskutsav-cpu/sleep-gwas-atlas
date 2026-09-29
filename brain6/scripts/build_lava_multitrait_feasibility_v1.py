#!/usr/bin/env python3
"""Pre-compute immutable lower bounds for any power-optimized LAVA rescue family."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CANONICAL = ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv"
SCREEN = ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/summary.json"
FEASIBILITY = ROOT / "brain6/results/power_optimized_sensitivity_v1/sensitivity_family_feasibility_v1.json"
OUT = ROOT / "brain6/results/lava_multitrait_feasibility_v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    f = json.loads(FEASIBILITY.read_text())
    s = json.loads(SCREEN.read_text())
    counts = f["canonical_not_run_by_trait"]
    ceiling = f["family"]["maximum_not_run"]
    total = f["family"]["canonical_not_run"]
    if sum(counts.values()) != total or (total,ceiling) != (3720,873):
        raise ValueError("frozen LAVA failure accounting changed")
    if sha(CANONICAL) != f["source_sha256"]["brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv"]:
        raise ValueError("canonical LAVA trait receipt changed")
    if sha(SCREEN) != f["source_sha256"]["brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/summary.json"]:
        raise ValueError("continuous-duration screen receipt changed")
    if s["rows"] != 2495 or s["not_run"] != 876:
        raise ValueError("observed alternative screen count changed")
    traits = tuple(sorted(counts))
    success = []
    for n in range(1,len(traits)+1):
        for subset in itertools.combinations(traits,n):
            lower = total - sum(counts[t] for t in subset)
            if lower <= ceiling:
                success.append({"traits_replaced":list(subset),"replacement_count":n,
                                "best_possible_not_run_if_all_replacements_perfect":lower})
    minimum = min(x["replacement_count"] for x in success)
    no_insomnia = [x for x in success if "insomnia" not in x["traits_replaced"]]
    no_insomnia_minimum = min(x["replacement_count"] for x in no_insomnia)
    result = {
        "analysis_id":"brain6-lava-multitrait-feasibility-v1",
        "status":"PROVEN_INFEASIBLE_FOR_CURRENT_SINGLE_TRAIT_REPLACEMENT",
        "canonical_family_cells":f["family"]["cells"],"canonical_not_run":total,
        "frozen_maximum_not_run":ceiling,"canonical_not_run_by_trait":counts,
        "minimum_number_of_perfect_trait_replacements":minimum,
        "perfect_replacement_sets_at_minimum":[x for x in success if x["replacement_count"]==minimum],
        "minimum_number_of_perfect_replacements_without_insomnia":no_insomnia_minimum,
        "perfect_replacement_sets_at_minimum_without_insomnia":[x for x in no_insomnia if x["replacement_count"]==no_insomnia_minimum],
        "observed_continuous_duration_not_run":s["not_run"],
        "observed_continuous_duration_alone_exceeds_family_ceiling_by":s["not_run"]-ceiling,
        "observed_continuous_duration_best_possible_whole_family_not_run_even_if_every_other_trait_perfect":s["not_run"],
        "source_sha256":{str(p.relative_to(ROOT)):sha(p) for p in (CANONICAL,SCREEN,FEASIBILITY)},
        "interpretation":"These are arithmetic lower bounds, not predicted experimental results. The observed continuous-duration substitute cannot pass the unchanged seven-trait 5% ceiling even if all other six traits had zero NOT_RUN cells. Any viable separate rescue must use a stronger long-sleep representation or an explicitly different prospective family/QC protocol; the canonical v3 failure remains unchanged."
    }
    OUT.mkdir(parents=True,exist_ok=False)
    (OUT / "lower_bounds.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:result[k] for k in ("status","minimum_number_of_perfect_trait_replacements","minimum_number_of_perfect_replacements_without_insomnia","observed_continuous_duration_alone_exceeds_family_ceiling_by")}))


if __name__ == "__main__":
    main()
