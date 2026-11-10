#!/usr/bin/env python3
"""Synthetic integration checks for research scorers, not score validation."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

root = Path(__file__).resolve().parent.parent / "outputs"
cases = [("score_15_research.py", "unvalidated_15_variant_score_prototype.csv"),
         ("score_22_research.py", "unvalidated_22_variant_score_prototype.csv"),
         ("shared_score_22_research.py", "unvalidated_22_shared_direction_score.csv")]
fields = ["sample_id", "variant", "dosage_allele", "other_allele", "dosage"]
for program, model in cases:
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        with (root / model).open(newline="") as handle:
            variants = list(csv.DictReader(handle))
        for reverse in (False, True):
            with (folder / ("reversed.csv" if reverse else "direct.csv")).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for item in variants:
                    a, b = item["counted_allele"], item["other_allele"]
                    writer.writerow({"sample_id": "synthetic", "variant": item["variant"],
                                     "dosage_allele": b if reverse else a,
                                     "other_allele": a if reverse else b,
                                     "dosage": "0.75" if reverse else "1.25"})
        base = [sys.executable, str(root / program), "--dosages"]
        direct = subprocess.run(base + [str(folder / "direct.csv"), "--output", str(folder / "direct.out")], capture_output=True, text=True)
        rejected = subprocess.run(base + [str(folder / "reversed.csv"), "--output", str(folder / "rejected.out")], capture_output=True, text=True)
        allowed = subprocess.run(base + [str(folder / "reversed.csv"), "--output", str(folder / "allowed.out"), "--allow-verified-reversal"], capture_output=True, text=True)
        assert direct.returncode == 0, direct.stderr
        assert rejected.returncode != 0 and "Reversed dosage" in rejected.stderr
        assert not (folder / "rejected.out").exists()
        assert allowed.returncode == 0, allowed.stderr
        assert (folder / "direct.out").read_text() == (folder / "allowed.out").read_text()
        print(f"{program}: direct accepted; reversal rejected by default; verified reversal matches direct")
