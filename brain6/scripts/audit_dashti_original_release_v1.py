#!/usr/bin/env python3
"""Stream the original Dashti long-sleep ZIP and inventory its native fields."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = Path("/Volumes/Extreme SSD/brain6-work/atlas-raw/.archives/longsumstats.txt.zip")
README = ROOT / "brain6/results/lava_confirmatory_source_triage_v1/source_docs/Saxena_fullUKBB_Longsleep_summary_stats_README"
CONTRACT = ROOT / "brain6/results/lava_longsleep_source_rescue_v1/phenotype_contract_v1.md"
FREEZE = CONTRACT.with_suffix(".freeze.json")
OUTPUT = ROOT / "brain6/results/lava_longsleep_source_rescue_v1/original_release_inventory_v1.json"
EXPECTED = ("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "INFO",
            "BETA_LONGSLEEP", "SE_LONGSLEEP", "P_LONGSLEEP")
MISSING = {b"", b"NA", b"NAN", b"NULL", b"NONE", b"."}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    frozen = json.loads(FREEZE.read_text())
    if digest(CONTRACT) != frozen["contract_sha256"]:
        raise ValueError("Phenotype contract changed after freeze")
    if digest(ARCHIVE) != frozen["source_archive_sha256"] or digest(README) != frozen["source_readme_sha256"]:
        raise ValueError("Original source or README changed")
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)

    rows = 0
    malformed = 0
    missing = Counter()
    chromosome = Counter()
    id_form = Counter()
    member_hash = hashlib.sha256()
    with ZipFile(ARCHIVE) as z:
        infos = z.infolist()
        if len(infos) != 1 or infos[0].filename != "longsumstats.txt":
            raise ValueError("Unexpected archive member list")
        info = infos[0]
        with z.open(info) as stream:
            header_line = stream.readline()
            member_hash.update(header_line)
            header = tuple(field.decode("utf-8") for field in header_line.rstrip(b"\r\n").split(b"\t"))
            if header != EXPECTED:
                raise ValueError(f"Unexpected source header: {header}")
            for line in stream:
                member_hash.update(line)
                rows += 1
                fields = line.rstrip(b"\r\n").split(b"\t")
                if len(fields) != len(header):
                    malformed += 1
                    continue
                for field, value in zip(header, fields):
                    if value.strip().upper() in MISSING:
                        missing[field] += 1
                chromosome[fields[1].decode("ascii", "replace")] += 1
                marker = fields[0]
                id_form["rsid" if marker.startswith(b"rs") else "other"] += 1
                if rows % 3_000_000 == 0:
                    print(f"source_rows_scanned={rows}", flush=True)

    cases, controls = 34_184, 305_742
    total = cases + controls
    neff = 4 * cases * controls / total
    result = {
        "analysis_id": "dashti_original_release_inventory_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_portal": "https://www.kp4cd.org/node/235",
        "source_zip_url": "https://personal.broadinstitute.org/ryank/longsumstats.txt.zip",
        "portal_release_date": "UNKNOWN_NOT_STATED_ON_INSPECTED_PORTAL",
        "archive_path": str(ARCHIVE),
        "archive_bytes": ARCHIVE.stat().st_size,
        "archive_sha256": frozen["source_archive_sha256"],
        "member_name": info.filename,
        "member_bytes": info.file_size,
        "member_crc32": f"{info.CRC:08x}",
        "member_sha256": member_hash.hexdigest(),
        "readme_path": str(README),
        "readme_sha256": frozen["source_readme_sha256"],
        "phenotype_contract_sha256": frozen["contract_sha256"],
        "header_exact": header_line.decode("utf-8").rstrip("\r\n"),
        "fields": list(header),
        "data_rows": rows,
        "malformed_width_rows": malformed,
        "missingness_by_field": {field: missing[field] for field in header},
        "chromosome_rows": dict(sorted(chromosome.items(), key=lambda x: x[0])),
        "id_forms": dict(id_form),
        "native_n_or_case_control_fields_present": False,
        "study_level_cases": cases,
        "study_level_controls": controls,
        "study_level_total_n": total,
        "balanced_equivalent_neff_formula": "4*cases*controls/(cases+controls)",
        "balanced_equivalent_neff": neff,
    }
    with OUTPUT.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps({"rows": rows, "malformed": malformed,
                      "member_sha256": result["member_sha256"], "neff": neff}, sort_keys=True))


if __name__ == "__main__":
    main()
