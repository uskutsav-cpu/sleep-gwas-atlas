#!/usr/bin/env python3
"""Exercise full 100-trait h2 and 12x98 rg collation in isolation."""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MARKER = "SYNTHETIC EXTENSION LDSC OUTPUT - NOT REAL RESULTS"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    panel = read_tsv(ROOT / "discovery_extension/config/candidate_traits.tsv")
    sleeps = [
        row["trait_id"] for row in read_tsv(ROOT / "config/analysis_panel.tsv")
        if row["domain"] == "sleep"
    ]
    ids = [row["extension_trait_id"] for row in panel]
    with tempfile.TemporaryDirectory(prefix="synthetic_extension_ldsc_") as temporary:
        work = Path(temporary)
        h2_logs, rg_logs = work / "h2", work / "rg"
        h2_logs.mkdir()
        rg_logs.mkdir()
        h2_out = work / "h2.tsv"
        rg_out = work / "rg.tsv"
        universe_out = work / "universe.tsv"

        for index, trait_id in enumerate(ids):
            h2, se, intercept = 0.2, 0.02, 1.0
            if index == 0:
                se = 0.1
            if index == 1:
                intercept = 1.3
            (h2_logs / f"h2_{trait_id}.log").write_text(
                f"{MARKER}\nTotal Observed scale h2: {h2} ({se})\n"
                f"Intercept: {intercept} (0.01)\nMean Chi^2: 1.2\nLambda GC: 1.1\nRatio: 0.1\n"
            )
        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
                "--mode", "h2", "--logdir", str(h2_logs), "--out", str(h2_out),
                "--provenance-out", str(work / "h2.json"),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--core", str(ROOT / "config/analysis_panel.tsv"),
            ], check=True,
        )
        h2_rows = read_tsv(h2_out)
        if len(h2_rows) != 100:
            raise SystemExit(f"ERROR: synthetic h2 rows={len(h2_rows)} expected 100")
        if sum(row["primary_rg_eligibility"] == "PRIMARY_PASS" for row in h2_rows) != 98:
            raise SystemExit("ERROR: synthetic h2 gate did not yield exactly 98 primary traits")

        passing = ids[2:]
        header = "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se"
        for sleep_index, sleep in enumerate(sleeps):
            p_values = [1e-8 if sleep_index == 0 and index == 0 else 0.5 for index in range(len(passing))]
            lines = [MARKER]
            lines.extend(f"P: {p_value}" for p_value in p_values)
            lines.extend(["Summary of Genetic Correlation Results", header])
            for index, (trait_id, p_value) in enumerate(zip(passing, p_values)):
                rg = 0.2 if sleep_index == 0 and index == 0 else 0.01
                lines.append(
                    f"data/munged/{sleep}.sumstats.gz "
                    f"discovery_extension/data/munged/{trait_id}.sumstats.gz "
                    f"{rg} 0.03 {rg/0.03:.5g} {p_value} 0.2 0.02 1.0 0.01 0.01 0.001"
                )
            lines.append("")
            (rg_logs / f"rg_{sleep}.log").write_text("\n".join(lines))

        subprocess.run(
            [
                "python3", str(ROOT / "discovery_extension/scripts/13_collate_extension_ldsc.py"),
                "--mode", "rg", "--logdir", str(rg_logs), "--h2", str(h2_out),
                "--out", str(rg_out), "--pair-universe-out", str(universe_out),
                "--provenance-out", str(work / "rg.json"),
                "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
                "--core", str(ROOT / "config/analysis_panel.tsv"),
            ], check=True,
        )
        rg_rows, universe = read_tsv(rg_out), read_tsv(universe_out)
        if len(rg_rows) != 1176 or len(universe) != 1200:
            raise SystemExit(f"ERROR: rg rows={len(rg_rows)} universe={len(universe)}")
        if sum(row["pair_status"] == "H2_FAILED_PRIMARY_EXCLUSION" for row in universe) != 24:
            raise SystemExit("ERROR: synthetic pair universe did not preserve 24 h2 exclusions")
        top = next(row for row in rg_rows if row["sleep_trait"] == sleeps[0] and row["extension_trait_id"] == passing[0])
        if float(top["extension_fdr"]) >= 0.05 or top["initial_screen_status"] != "PAIR_NOVELTY_AUDIT_REQUIRED":
            raise SystemExit(f"ERROR: synthetic isolated FDR/priority failed: {top}")
        figure = work / "screen.png"
        subprocess.run(
            [
                str(ROOT / ".venv/bin/python"),
                str(ROOT / "discovery_extension/scripts/15_plot_extension.py"),
                "--rg", str(rg_out), "--out", str(figure),
            ], check=True,
        )
        if not figure.is_file() or figure.stat().st_size < 10000 or not figure.with_suffix(".pdf").is_file():
            raise SystemExit("ERROR: synthetic discovery figure was not created")
    print("EXTENSION_LDSC_COLLATION_SYNTHETIC_OK h2=100 primary_pairs=1176 universe=1200 figure=true isolated=true")


if __name__ == "__main__":
    main()
