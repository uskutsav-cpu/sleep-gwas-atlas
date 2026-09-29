"""Build a descriptive comparison to Jia et al.'s published 7 x 3 rg table."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE_MAP = ROOT / "brain6/results/global/brain6_72_locked.tsv"
OUTPUT = ROOT / "brain6/results/novelty/jia_2025_sleep_psychiatric_rg_context.tsv"
PROVENANCE = OUTPUT.with_suffix(".provenance.json")
SOURCE_DOI = "10.1093/sleep/zsae209"
SOURCE_URL = "https://academic.oup.com/sleep/article/48/1/zsae209/7750695"

# Table 1 values transcribed from the publisher's article page. The source page
# identifies these as LDSC, HDL, and GPA results and flags significant LDSC
# sample-overlap intercepts with an asterisk.
SOURCE_TSV = """source_sleep_trait\tbrain6_sleep_trait\tsource_disorder\tbrain6_disorder\tldsc_rg\tldsc_se\tldsc_p\tldsc_intercept\tldsc_intercept_se\tintercept_overlap_flag\thdl_rg\thdl_se\thdl_p\tpm11\tpar\tgpa_p
Insomnia\tinsomnia\tBP\tbipolar\t0.1100\t0.0261\t2.49e-05\t-0.0100\t0.0061\tfalse\t0.0833\t0.0196\t2.11e-05\t0.167\t0.433\t<1e-300
Insomnia\tinsomnia\tMDD\tmdd\t0.4435\t0.0255\t7.88e-68\t0.0856\t0.0064\ttrue\t0.4495\t0.0226\t6.06e-88\t0.217\t0.764\t<1e-300
Insomnia\tinsomnia\tSCZ\tscz\t0.0365\t0.0232\t1.15e-01\t0.0024\t0.0065\tfalse\t0.0255\t0.0163\t1.17e-01\t0.166\t0.407\t<1e-300
Chronotype\tchronotype\tBP\tbipolar\t-0.0534\t0.0215\t1.30e-02\t0.0021\t0.0065\tfalse\t-0.0545\t0.0113\t1.47e-06\t0.164\t0.414\t<1e-300
Chronotype\tchronotype\tMDD\tmdd\t-0.0811\t0.0219\t2.00e-04\t-0.0364\t0.0069\ttrue\t-0.0892\t0.0172\t2.07e-07\t0.156\t0.429\t<1e-300
Chronotype\tchronotype\tSCZ\tscz\t-0.1298\t0.0189\t6.51e-12\t-0.0012\t0.0070\tfalse\t-0.1224\t0.0138\t7.42e-19\t0.173\t0.424\t<1e-300
Nap\tnapping\tBP\tbipolar\t0.1353\t0.0232\t5.29e-09\t0.0087\t0.0061\tfalse\t0.1371\t0.0147\t9.43e-21\t0.154\t0.425\t<1e-300
Nap\tnapping\tMDD\tmdd\t0.1950\t0.0213\t5.90e-20\t0.0832\t0.0065\ttrue\t0.2207\t0.0211\t1.16e-25\t0.176\t0.516\t<1e-300
Nap\tnapping\tSCZ\tscz\t0.1700\t0.0197\t5.76e-18\t-0.0045\t0.0073\tfalse\t0.1389\t0.0123\t1.70e-29\t0.161\t0.419\t<1e-300
EDS\tsleepiness\tBP\tbipolar\t0.1236\t0.0249\t6.83e-07\t0.0031\t0.0059\tfalse\t0.1193\t0.0171\t2.68e-12\t0.151\t0.428\t<1e-300
EDS\tsleepiness\tMDD\tmdd\t0.1668\t0.0253\t4.17e-11\t0.0642\t0.0061\ttrue\t0.1758\t0.0220\t1.24e-15\t0.156\t0.476\t<1e-300
EDS\tsleepiness\tSCZ\tscz\t0.0919\t0.0222\t3.58e-05\t-0.0018\t0.0069\tfalse\t0.0702\t0.0154\t5.38e-06\t0.155\t0.419\t<1e-300
Duration\tsleepdur\tBP\tbipolar\t0.1122\t0.0230\t1.07e-06\t-0.0034\t0.0061\tfalse\t0.1093\t0.0142\t1.16e-14\t0.162\t0.441\t<1e-300
Duration\tsleepdur\tMDD\tmdd\t-0.1095\t0.0238\t4.16e-06\t-0.0127\t0.0065\tfalse\t-0.1236\t0.0224\t3.24e-08\t0.158\t0.462\t<1e-300
Duration\tsleepdur\tSCZ\tscz\t0.1498\t0.0201\t8.34e-14\t-0.0057\t0.0070\tfalse\t0.1329\t0.0128\t3.30e-25\t0.166\t0.432\t<1e-300
Long duration\tlongsleep\tBP\tbipolar\t0.2315\t0.0344\t1.62e-11\t-0.0037\t0.0052\tfalse\t0.2103\t0.0234\t2.16e-19\t0.103\t0.325\t<1e-300
Long duration\tlongsleep\tMDD\tmdd\t0.2966\t0.0331\t3.27e-19\t0.0343\t0.0055\ttrue\t0.2752\t0.0344\t1.17e-15\t0.158\t0.525\t<1e-300
Long duration\tlongsleep\tSCZ\tscz\t0.2853\t0.0299\t1.42e-21\t-0.0032\t0.0062\tfalse\t0.2503\t0.0180\t7.02e-44\t0.110\t0.333\t<1e-300
Short duration\tshortsleep\tBP\tbipolar\t-0.0053\t0.0252\t8.32e-01\t-0.0005\t0.0058\tfalse\t-0.0151\t0.0136\t2.64e-01\t0.163\t0.428\t<1e-300
Short duration\tshortsleep\tMDD\tmdd\t0.3084\t0.0256\t1.59e-33\t0.0435\t0.0065\ttrue\t0.3047\t0.0214\t4.81e-46\t0.185\t0.601\t<1e-300
Short duration\tshortsleep\tSCZ\tscz\t-0.0286\t0.0234\t2.20e-01\t0.0054\t0.0064\tfalse\t-0.0259\t0.0145\t7.54e-02\t0.173\t0.438\t<1e-300
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_rows() -> list[dict[str, str]]:
    source_rows = list(csv.DictReader(io.StringIO(SOURCE_TSV), delimiter="\t"))
    with SOURCE_MAP.open(encoding="utf-8", newline="") as stream:
        global_rows = list(csv.DictReader(stream, delimiter="\t"))
    global_by_key = {(r["sleep_trait"], r["brain_disorder"]): r for r in global_rows}
    out = []
    for source in source_rows:
        key = (source["brain6_sleep_trait"], source["brain6_disorder"])
        current = global_by_key.get(key)
        if current is None:
            raise ValueError(f"Published pair is missing from locked Brain6 map: {key}")
        row = dict(source)
        row.update({
            "brain6_rg": current["rg"],
            "brain6_se": current["se"],
            "brain6_p": current["p"],
            "brain6_original_396_family_q": current["original_BH_FDR_q"],
            "brain6_original_family_significant": current["significance_under_original_396_family"],
            "rg_difference_published_minus_brain6": format(float(source["ldsc_rg"]) - float(current["rg"]), ".12g"),
            "direction_concordant": str(float(source["ldsc_rg"]) * float(current["rg"]) > 0).lower(),
            "comparison_independence": "NOT_ESTABLISHED; European GWAS; source IDs/participant overlap not fully verified",
            "interpretation": "Published descriptive rg comparison; not independent replication and not a new Brain6 test",
        })
        out.append(row)
    if len(out) != 21 or len({(r["brain6_sleep_trait"], r["brain6_disorder"]) for r in out}) != 21:
        raise ValueError("Expected 21 unique published pair rows")
    return out


def expected_bytes(rows: list[dict[str, str]]) -> bytes:
    fields = list(rows[0])
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def run(validate: bool) -> dict[str, object]:
    rows = build_rows()
    data = expected_bytes(rows)
    if validate:
        if not OUTPUT.is_file() or OUTPUT.read_bytes() != data:
            raise ValueError("Published-rg context does not rebuild byte-identically")
        provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
        expected = {
            "output_sha256": hashlib.sha256(data).hexdigest(),
            "global_map_sha256": sha256(SOURCE_MAP),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "row_count": 21,
            "source_doi": SOURCE_DOI,
            "status": "PASS_DESCRIPTIVE_PUBLISHED_RG_CONTEXT_NOT_INDEPENDENT_REPLICATION",
            "source_overlap_note": "The publisher marks six MDD-pair LDSC intercepts as significant sample overlap. Exact source-study identity and participant-level overlap with the locked Brain6 inputs are not fully verified.",
        }
        for key, value in expected.items():
            if provenance.get(key) != value:
                raise ValueError(f"Published-rg provenance mismatch for {key}")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        if OUTPUT.exists() or PROVENANCE.exists():
            raise FileExistsError("Refusing to overwrite published-rg context outputs")
        OUTPUT.write_bytes(data)
        provenance = {
            "schema_version": 1,
            "status": "PASS_DESCRIPTIVE_PUBLISHED_RG_CONTEXT_NOT_INDEPENDENT_REPLICATION",
            "source_title": "From single nucleotide variations to genes: identifying the genetic links between sleep and psychiatric disorders",
            "source_doi": SOURCE_DOI,
            "source_url": SOURCE_URL,
            "source_publication": "Sleep. 2025;48(1):zsae209. Published 2024-09-07.",
            "source_accessed_date": "2026-09-26",
            "source_table": "Table 1; LDSC, HDL, GPA columns",
            "source_scope": "7 sleep/circadian traits x 3 European-ancestry psychiatric disorders; 21 pairs",
            "source_overlap_note": "The publisher marks six MDD-pair LDSC intercepts as significant sample overlap. Exact source-study identity and participant-level overlap with the locked Brain6 inputs are not fully verified.",
            "interpretation": "Descriptive published comparison only; not an independent replication, new discovery family, or mechanism claim.",
            "row_count": len(rows),
            "global_map_path": str(SOURCE_MAP.relative_to(ROOT)),
            "global_map_sha256": sha256(SOURCE_MAP),
            "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "output_path": str(OUTPUT.relative_to(ROOT)),
            "output_sha256": sha256(OUTPUT),
        }
        PROVENANCE.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"status": provenance["status"], "rows": len(rows), "output_sha256": sha256(OUTPUT)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.validate), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
