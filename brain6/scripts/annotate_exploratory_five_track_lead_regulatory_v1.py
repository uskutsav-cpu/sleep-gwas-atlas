#!/usr/bin/env python3
"""Cache GRCh37 Ensembl regulatory overlap at every five-track lead coordinate."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SIGNED = ROOT / "brain6/results/exploratory_five_track_v1"
B_ANNOTATION = SIGNED / "b_annotation"
MEMBERS = ROOT / "brain6/results/loci/five_track_cross_pair_region_reconciliation_v1/pair_candidate_members.tsv"
OUT = ROOT / "brain6/results/exploratory_five_track_v6/lead_regulatory_context"
API = "https://grch37.rest.ensembl.org/overlap/region/human"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path: Path, records: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    if (OUT / "provenance.json").exists():
        raise FileExistsError("Versioned regulatory lookup already complete")
    source_receipt = json.loads((B_ANNOTATION / "provenance.json").read_text())
    for name, expected in source_receipt["output_sha256"].items():
        if sha(B_ANNOTATION / name) != expected:
            raise ValueError(f"Track B annotation hash mismatch: {name}")
    b_raw = json.loads((B_ANNOTATION / "raw_official_api_responses.json").read_text())
    b_by_snp = {r["snp"]: r for r in b_raw}
    if len(b_by_snp) != 6:
        raise ValueError("Expected six cached protected Track B leads")
    signed_paths = [SIGNED / "b_signed_ld/candidate_signed_ld.tsv", SIGNED / "fourpair_signed_ld/candidate_signed_ld.tsv"]
    leads = [r for path in signed_paths for r in read(path) if r["candidate_status"] == "LEAD"]
    if len(leads) != 27 or len({(r["CHR"], r["BP"], r["SNP"]) for r in leads}) != 26:
        raise ValueError("Expected 27 pair-specific lead rows, 26 unique variants")
    region = {r["candidate_locus_id"]: r["region_group"] for r in read(MEMBERS)}
    if {r["locus_id"] for r in leads} != set(region):
        raise ValueError("Protected 25 candidate loci changed")
    OUT.mkdir(parents=True, exist_ok=True)
    raw_dir = OUT / "raw_new_queries"
    raw_dir.mkdir(exist_ok=True)
    raw_hashes = {}
    responses = {}
    for chrom, bp, snp in sorted({(r["CHR"], r["BP"], r["SNP"]) for r in leads},
                                 key=lambda x: (int(x[0]), int(x[1]), x[2])):
        key = (chrom, bp, snp)
        url = f"{API}/{chrom}:{bp}-{bp}?feature=regulatory"
        if snp in b_by_snp:
            previous = b_by_snp[snp]
            if str(previous["chromosome"]) != chrom or str(previous["bp_grch37"]) != bp or previous["regulatory_url"] != url:
                raise ValueError(f"Protected Track B response identity changed: {snp}")
            records = previous["regulatory_response"]
            origin = "REUSED_PROTECTED_B_SNAPSHOT"
        else:
            path = raw_dir / f"{snp}_{chrom}_{bp}.json"
            if path.exists():
                previous = json.loads(path.read_text())
                if previous["url"] != url or previous["snp"] != snp:
                    raise ValueError(f"Cached query identity changed: {path}")
                records = previous["response"]
            else:
                result = subprocess.run(
                    ["curl", "--fail", "--silent", "--show-error", "--location", "--retry", "3",
                     "--retry-delay", "2", "--max-time", "30", "--header", "Accept: application/json", url],
                    check=True, capture_output=True, text=True,
                )
                records = json.loads(result.stdout)
                if not isinstance(records, list):
                    raise ValueError(f"Official Ensembl response not a list: {snp}")
                with path.open("x") as stream:
                    json.dump({"url": url, "snp": snp, "chromosome": chrom, "bp_grch37": int(bp),
                               "response": records}, stream, indent=2, sort_keys=True)
                    stream.write("\n")
                time.sleep(0.3)
            raw_hashes[str(path.relative_to(OUT))] = sha(path)
            origin = "NEW_OFFICIAL_ENSEMBL_QUERY"
        if not isinstance(records, list):
            raise ValueError(f"Invalid regulatory records for {snp}")
        for record in records:
            if int(record["start"]) > int(bp) or int(record["end"]) < int(bp):
                raise ValueError(f"Regulatory feature does not overlap query coordinate: {snp}")
        responses[key] = (records, origin, url)
        print(f"{snp}: {len(records)} regulatory features ({origin})", flush=True)

    summary = []
    features = []
    for lead in leads:
        records, origin, url = responses[(lead["CHR"], lead["BP"], lead["SNP"])]
        base = {"region_group": region[lead["locus_id"]], "pair_id": lead["pair_id"],
                "locus_id": lead["locus_id"], "lead_snp": lead["SNP"],
                "chr_grch37": lead["CHR"], "bp_grch37": lead["BP"],
                "source_url": url, "source_origin": origin}
        summary.append(base | {"n_regulatory_features_overlapping_exact_lead": str(len(records)),
                               "interpretation": "EXPLORATORY_COORDINATE_OVERLAP_ONLY",
                               "analysis_label": "EXPLORATORY", "lava_local_rg": "BLOCKED_LAVA"})
        for feature in records:
            features.append(base | {"feature_id": str(feature.get("id", "")),
                                    "feature_type": str(feature.get("feature_type", "")),
                                    "description": str(feature.get("description", "")),
                                    "feature_start": str(feature["start"]), "feature_end": str(feature["end"]),
                                    "interpretation": "EXPLORATORY_COORDINATE_OVERLAP_ONLY"})
    if len(summary) != 27:
        raise ValueError("Lead annotation coverage incomplete")
    write(OUT / "lead_regulatory_summary.tsv", summary, list(summary[0]))
    write(OUT / "lead_regulatory_features.tsv", features,
          list(features[0]) if features else list(summary[0])[:8] + ["feature_id", "feature_type", "description", "feature_start", "feature_end", "interpretation"])
    report = ["# Five-track lead regulatory overlap — exploratory", "",
              "All 27 pair-specific lead rows (26 unique variants) were queried at their exact",
              "GRCh37 coordinates using the official Ensembl GRCh37 regulatory overlap endpoint.",
              "Six protected Track B responses were reused byte-for-byte; the other 20 unique",
              "queries are individually cached and hashed. A returned feature means coordinate",
              "overlap only, not regulatory function, variant causality, or shared-trait support.",
              "An empty response cannot exclude regulation by other variants in the locus.", "",
              f"Exact-lead regulatory feature rows: {len(features)}.",
              f"Lead rows with at least one overlap: {sum(int(x['n_regulatory_features_overlapping_exact_lead']) > 0 for x in summary)} / 27.",
              "All LAVA local-rg and final tiers remain BLOCKED_LAVA.", ""]
    (OUT / "report.md").write_text("\n".join(report))
    receipt = {"analysis_label": "EXPLORATORY", "assembly": "GRCh37", "lava_local_rg": "BLOCKED_LAVA",
               "official_api_documentation": "https://grch37.rest.ensembl.org/documentation/info/overlap_region",
               "protected_b_receipt_sha256": sha(B_ANNOTATION / "provenance.json"),
               "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in signed_paths + [MEMBERS]},
               "raw_query_sha256": raw_hashes,
               "output_sha256": {p.name: sha(p) for p in (
                   OUT / "lead_regulatory_summary.tsv", OUT / "lead_regulatory_features.tsv", OUT / "report.md")}}
    (OUT / "provenance.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
