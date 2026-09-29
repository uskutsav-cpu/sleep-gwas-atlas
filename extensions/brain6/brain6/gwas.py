"""Bounded-memory GWAS normalization and allele-aware joining.

Never changes genome build by renaming it. Build conversions must be supplied
as reviewed upstream artifacts and given new checksums. No effect sizes or
sample sizes are inferred from p-values or paper headlines.
"""
from __future__ import annotations
import collections
import csv
import math
import sqlite3
from pathlib import Path
from .artifacts import transaction, verify_artifact
from .io import (ContractError, check_hash, ensure_free, open_text, read_json,
                 require, safe_id, write_json, write_tsv)
from .stats import effective_n

FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N", "EAF", "INFO"]
PAIR_FIELDS = ["SNP", "CHR", "BP", "A1", "A2", "BETA1", "SE1", "P1", "N1", "EAF1",
               "BETA2", "SE2", "P2", "N2", "EAF2", "Z1", "Z2", "alignment"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def align(a1: str, a2: str, b1: str, b2: str) -> tuple[int, str] | None:
    """Sign needed to align B effects to A. Palindromes always fail closed."""
    a1, a2, b1, b2 = [x.upper() for x in (a1, a2, b1, b2)]
    if any(x not in {"A", "C", "G", "T"} for x in (a1, a2, b1, b2)):
        return None
    if a1 == a2 or b1 == b2 or {a1, a2} in ({"A", "T"}, {"C", "G"}):
        return None
    options = [(b1, b2, 1, "same"), (b2, b1, -1, "swapped"),
               (b1.translate(COMPLEMENT), b2.translate(COMPLEMENT), 1, "complement"),
               (b2.translate(COMPLEMENT), b1.translate(COMPLEMENT), -1, "complement_swapped")]
    for x, y, sign, label in options:
        if (x, y) == (a1, a2):
            return sign, label
    return None


def optional_float(value: str | None) -> float | None:
    if value is None or value.strip().lower() in {"", "na", "nan", ".", "null", "none"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("nonfinite")
    return result


def validate_source(source: dict, *, synthetic=False) -> None:
    for key in ["trait_id", "path", "sha256", "study_id", "phenotype_definition", "ancestry",
                "genome_build", "effect_scale", "n_semantics", "column_map"]:
        require(key in source and source[key] not in (None, "", "UNRESOLVED"), f"Missing source field: {key}")
    safe_id(source["trait_id"])
    require(source["genome_build"] in {"GRCh37", "GRCh38"}, "Unsupported/missing build")
    require(source["ancestry"] == "EUR", "This extension contract requires EUR; use a separate ancestry protocol")
    require(source["effect_scale"] in {"beta", "log_odds", "odds_ratio"}, "Unknown effect scale")
    require(source["n_semantics"] in {"total", "effective"}, "Sample-size semantics must be explicit")
    require(source.get("synthetic", False) == synthetic, "Synthetic/empirical source mix forbidden")
    if not synthetic:
        require(source.get("access_permitted") is True, "Source access/licensing not acknowledged")
        require(source.get("source_uri") not in (None,"","UNRESOLVED") and source.get("citation") not in (None,"","UNRESOLVED"), "Missing source provenance")
    cols = source["column_map"]
    require(set(["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]) <= set(cols),
            "Explicit SNP/coordinate/allele/effect/SE/P mapping required")
    require("N" in cols or source.get("sample_size") is not None or
            {"N_CASES", "N_CONTROLS"} <= set(cols), "No sample-size route")
    require(len(set(cols.values())) == len(cols), "Ambiguous source-column mapping")
    check_hash(source["path"], source["sha256"])


def normalize(source_path: str | Path, out_root: str | Path, relative: str,
              *, synthetic=False) -> Path:
    source = read_json(source_path)
    validate_source(source, synthetic=synthetic)
    policy = {"min_maf": 0.01, "min_info": 0.9, "require_eaf": True,
              "require_info": False, "min_n": 100, "batch_rows": 20000,
              "exclude_regions": [], **source.get("qc", {})}
    require(0 <= policy["min_maf"] < 0.5 and 0 <= policy["min_info"] <= 1,
            "Invalid MAF/INFO thresholds")
    ensure_free(out_root, 0 if synthetic else int(source.get("minimum_free_bytes", 2 * 2**30)))
    target = Path(out_root) / relative
    with transaction(out_root, relative, stage="normalize", inputs=[source_path, source["path"]],
                     parameters=policy, synthetic=synthetic) as (work, meta):
        db = sqlite3.connect(work / "variants.sqlite")
        db.execute("PRAGMA cache_size=-65536")
        db.execute("PRAGMA temp_store=FILE")
        db.execute("""CREATE TABLE variants (
          snp TEXT PRIMARY KEY, chr INTEGER NOT NULL, bp INTEGER NOT NULL,
          a1 TEXT NOT NULL, a2 TEXT NOT NULL, beta REAL NOT NULL, se REAL NOT NULL,
          p REAL NOT NULL, n REAL NOT NULL, eaf REAL, info REAL, seen INTEGER NOT NULL DEFAULT 1)
          WITHOUT ROWID""")
        qc = collections.Counter()
        with open_text(source["path"]) as f:
            reader = csv.DictReader(f, delimiter=source.get("delimiter", "\t"))
            require(reader.fieldnames is not None, "Empty GWAS input")
            require(len(set(reader.fieldnames)) == len(reader.fieldnames), "Duplicate source headers")
            cols = source["column_map"]
            require(set(cols.values()) <= set(reader.fieldnames), "Column map does not match source header")
            rows = []
            for line, raw in enumerate(reader, 2):
                qc["input_rows"] += 1
                require(None not in raw and all(v is not None for v in raw.values()),
                        f"Malformed source row {line}; fix source-specific materialization instead of silently skipping")
                row = {k: raw[v] for k, v in cols.items()}
                try:
                    snp = row["SNP"].strip()
                    if not snp or any(c.isspace() for c in snp):
                        qc["invalid_id"] += 1; continue
                    chrom = int(row["CHR"].lower().removeprefix("chr"))
                    bp = int(row["BP"])
                    if chrom not in range(1, 23) or bp < 1:
                        qc["non_autosomal_or_invalid_coordinate"] += 1; continue
                    a1, a2 = row["A1"].upper(), row["A2"].upper()
                    if align(a1, a2, a1, a2) is None:
                        qc["palindromic_or_non_biallelic_snp"] += 1; continue
                    beta, se, p = [optional_float(row[x]) for x in ("BETA", "SE", "P")]
                    if None in (beta, se, p) or se <= 0 or not 0 <= p <= 1:
                        qc["invalid_statistics"] += 1; continue
                    if source["effect_scale"] == "odds_ratio":
                        if beta <= 0:
                            qc["invalid_odds_ratio"] += 1; continue
                        require(source.get("se_scale") == "log_odds", "OR inputs require log-odds SEs")
                        beta = math.log(beta)
                    if "N" in row:
                        n = optional_float(row["N"])
                        if n is not None:
                            n *= float(source.get("n_multiplier", 1))
                    elif {"N_CASES", "N_CONTROLS"} <= set(row):
                        nc, nn = float(row["N_CASES"]), float(row["N_CONTROLS"])
                        if not (math.isfinite(nc) and math.isfinite(nn) and nc>0 and nn>0):
                            qc["invalid_case_control_counts"] += 1; continue
                        n = effective_n(nc, nn) if source["n_semantics"] == "effective" else nc + nn
                    else:
                        n = float(source["sample_size"])
                    if n is None or not math.isfinite(n) or n < policy["min_n"]:
                        qc["invalid_or_low_n"] += 1; continue
                    eaf, info = optional_float(row.get("EAF")), optional_float(row.get("INFO"))
                    if eaf is None:
                        qc["missing_eaf_seen"] += 1
                        if policy["require_eaf"]:
                            qc["missing_eaf_excluded"] += 1; continue
                    elif not 0 <= eaf <= 1 or min(eaf, 1 - eaf) < policy["min_maf"]:
                        qc["frequency_filter"] += 1; continue
                    if info is None:
                        qc["missing_info_seen"] += 1
                        if policy["require_info"]:
                            qc["missing_info_excluded"] += 1; continue
                    elif not 0 <= info <= 1 or info < policy["min_info"]:
                        qc["info_filter"] += 1; continue
                    if any(chrom == int(reg[0]) and int(reg[1]) <= bp <= int(reg[2])
                           for reg in policy["exclude_regions"]):
                        qc["excluded_region"] += 1; continue
                    # Duplicate rsIDs are all quarantined, not first-row-wins.
                    rows.append((snp, chrom, bp, a1, a2, beta, se, p, n, eaf, info))
                    if len(rows) >= policy["batch_rows"]:
                        _insert(db, rows); rows.clear()
                except (ValueError, OverflowError) as exc:
                    if isinstance(exc, ContractError):
                        raise
                    qc["parse_error"] += 1
            _insert(db, rows)
        qc["duplicate_ids"] = db.execute("SELECT count(*) FROM variants WHERE seen>1").fetchone()[0]
        db.execute("DELETE FROM variants WHERE seen>1")
        # Distinct IDs at the same coordinate+alleles are also ambiguous.
        db.execute("CREATE INDEX variant_position ON variants(chr,bp)")
        # Keep ambiguous-position filtering disk-backed; fetchall can exhaust RAM
        # on pathological multi-allelic sources.
        db.execute("CREATE TEMP TABLE ambiguous_positions AS SELECT chr,bp FROM variants GROUP BY chr,bp HAVING count(*)>1")
        qc["ambiguous_positions"] = db.execute("SELECT count(*) FROM ambiguous_positions").fetchone()[0]
        qc["ambiguous_position_rows"] = db.execute("SELECT count(*) FROM variants JOIN ambiguous_positions USING(chr,bp)").fetchone()[0]
        db.execute("DELETE FROM variants WHERE (chr,bp) IN (SELECT chr,bp FROM ambiguous_positions)")
        db.commit()
        count = db.execute("SELECT count(*) FROM variants").fetchone()[0]
        require(count > 0, "No variants passed QC")
        def export():
            for row in db.execute("SELECT snp,chr,bp,a1,a2,beta,se,p,n,eaf,info FROM variants ORDER BY chr,bp,snp"):
                yield dict(zip(FIELDS, ["NA" if x is None else x for x in row]))
        write_tsv(work / "sumstats.tsv.gz", FIELDS, export())
        qc["retained_rows"] = count
        meta["source"] = {k: source[k] for k in ["trait_id", "study_id", "phenotype_definition",
                                               "genome_build", "ancestry", "effect_scale", "n_semantics"]}
        meta["qc"] = dict(qc)
        meta["scientific_status"] = "PASS"
        meta["interpretation"] = "Source normalization QC only; not a hypothesis-test result"
        write_json(work / "qc.json", dict(qc))
        db.close()
    return target


def _insert(db, rows):
    db.executemany("""INSERT INTO variants(snp,chr,bp,a1,a2,beta,se,p,n,eaf,info)
      VALUES (?,?,?,?,?,?,?,?,?,?,?) ON CONFLICT(snp) DO UPDATE SET seen=seen+1""", rows)
    db.commit()


def join_pair(left: str | Path, right: str | Path, out_root: str | Path,
              relative: str, *, min_variants=1000000, min_overlap=0.8,
              max_eaf_difference=0.15, synthetic=False) -> Path:
    """Disk-backed SNP/coordinate join. B effects/frequencies are aligned to A."""
    left, right = Path(left), Path(right)
    lm, rm = verify_artifact(left), verify_artifact(right)
    require(lm["synthetic"] == rm["synthetic"] == synthetic, "Cannot mix synthetic/empirical inputs")
    require(lm["source"]["genome_build"] == rm["source"]["genome_build"], "Genome build mismatch; no implicit liftOver")
    require(lm["source"]["ancestry"] == rm["source"]["ancestry"], "Ancestry mismatch")
    require(0 <= min_overlap <= 1 and min_variants >= 1, "Invalid overlap threshold")
    require(0 <= max_eaf_difference <= 1, "Invalid allele-frequency difference threshold")
    params = dict(min_variants=min_variants, min_overlap=min_overlap,
                  max_eaf_difference=max_eaf_difference)
    with transaction(out_root, relative, stage="join_pair",
                     inputs=[left / "variants.sqlite", right / "variants.sqlite", left / "receipt.json", right / "receipt.json"],
                     parameters=params, synthetic=synthetic) as (work, meta):
        db = sqlite3.connect(f"{(left / 'variants.sqlite').resolve().as_uri()}?mode=ro", uri=True)
        db.execute("ATTACH DATABASE ? AS b", (f"{(right / 'variants.sqlite').resolve().as_uri()}?mode=ro",))
        query = """SELECT a.snp,a.chr,a.bp,a.a1,a.a2,a.beta,a.se,a.p,a.n,a.eaf,
                   b.chr,b.bp,b.a1,b.a2,b.beta,b.se,b.p,b.n,b.eaf
                   FROM main.variants a JOIN b.variants b ON a.snp=b.snp ORDER BY a.chr,a.bp,a.snp"""
        qc = collections.Counter()
        chromosomes = set()
        def rows():
            for r in db.execute(query):
                qc["common_ids"] += 1
                if r[1:3] != r[10:12]:
                    qc["coordinate_mismatch"] += 1; continue
                aligned = align(r[3], r[4], r[12], r[13])
                if aligned is None:
                    qc["allele_mismatch"] += 1; continue
                sign, label = aligned
                freq = r[18] if sign == 1 or r[18] is None else 1 - r[18]
                if r[9] is not None and freq is not None and abs(r[9] - freq) > max_eaf_difference:
                    qc["frequency_mismatch"] += 1; continue
                beta2 = sign * r[14]
                values = [*r[:10], beta2, r[15], r[16], r[17], freq,
                          r[5] / r[6], beta2 / r[15], label]
                qc[label] += 1; qc["retained_rows"] += 1
                chromosomes.add(r[1])
                yield dict(zip(PAIR_FIELDS, ["NA" if v is None else v for v in values]))
        write_tsv(work / "pair.tsv.gz", PAIR_FIELDS, rows())
        smaller = min(lm["qc"]["retained_rows"], rm["qc"]["retained_rows"])
        fraction = qc["retained_rows"] / smaller
        require(qc["retained_rows"] >= min_variants, "BLOCKED_BY_VARIANT_COUNT")
        require(fraction >= min_overlap, f"BLOCKED_BY_OVERLAP: {fraction:.3f}")
        if not synthetic:
            require(chromosomes == set(range(1, 23)), "Missing autosome coverage")
        meta.update(qc=dict(qc), overlap_fraction=fraction, chromosomes=sorted(chromosomes),
                    left=lm["source"], right=rm["source"], scientific_status="PASS")
        write_json(work / "qc.json", {**dict(qc), "overlap_fraction": fraction})
        db.close()
    return Path(out_root) / relative
