#!/usr/bin/env python3
"""
Phase 0: harmonize one GWAS summary statistics file into a standard schema.

Every QC decision below is a CHOICE, and each one is tagged [CDG3] where it
matches the Methods of Grotzinger/Werme et al. Nature 649:406-415 (2026), or
[LDSC] where it is required by munge_sumstats.py. You will be asked in your
meeting why you did each of these. Read them, don't just run them.

Output schema (tab-separated, gzipped):
    SNP  CHR  BP  A1  A2  FRQ  BETA  SE  P  N
where A1 is the EFFECT allele and BETA is on the log-odds scale for binary
traits.

Usage:
    python3 01_harmonize.py --trait insomnia --config config/traits.tsv \
        --infile data/raw/insomnia.txt.gz --outdir data/harmonized

    # With effective-N mode (CDG3 ascertainment correction):
    python3 01_harmonize.py --trait insomnia --config config/traits.tsv \
        --infile data/raw/insomnia.txt.gz --outdir data/harmonized \
        --neff-mode effective

    # With build sanity-checking:
    python3 01_harmonize.py --trait insomnia --config config/traits.tsv \
        --infile data/raw/insomnia.txt.gz --outdir data/harmonized \
        --build-check
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

# --- QC thresholds (edit here, they are logged to the report) ------------------
INFO_MIN = 0.6      # [CDG3] restrict to imputation INFO > 0.6 when available
MAF_MIN = 0.01      # [CDG3] minor allele frequency > 1%
MHC = ("6", 25_000_000, 34_000_000)   # [CDG3] MHC excluded from all sumstats
AMBIG = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
VALID = {"A", "C", "G", "T"}
NEFF_SNP_FRAC = 0.50  # [CDG3] drop SNPs where SNP-specific Neff < 50% of total

# hg38 sentinel positions — if >10% of chr6 SNPs have BP > 170Mb, suspect hg38
# (hg19 chr6 is ~171Mb, hg38 chr6 is ~171Mb but many datasets have different
# coordinate distributions)
HG38_SENTINEL_CHRS = {
    "1": 248_956_422,  # hg38 length
    "6": 170_805_979,
}

# Column aliases seen in the wild. Extend as you meet new files.
# Covers: standard, UKB, PGC, GIANT, GLGC, DIAMANTE, FinnGen, METAL output
ALIASES = {
    "SNP":  ["snp", "rsid", "rs_id", "rsids", "markername", "variant_id", "id",
             "marker", "snpid", "rsid_dbsnp", "varid", "snp_id", "#snpid"],
    "CHR":  ["chr", "chrom", "chromosome", "#chrom", "hg19chr", "chr_id",
             "chr_name", "#chr", "seqnames"],
    "BP":   ["bp", "pos", "position", "base_pair_location", "bp_hg19",
             "pos_hg19", "bpos", "chrompos", "pos_b37", "bp_grch37",
             "bp_grch38", "genpos", "start"],
    "A1":   ["a1", "effect_allele", "ea", "allele1", "tested_allele", "alt",
             "a1_effect", "risk_allele", "inc_allele", "effectallele",
             "reference_allele", "coded_allele"],
    "A2":   ["a2", "other_allele", "nea", "allele2", "non_effect_allele", "ref",
             "a2_other", "noneffect_allele", "noncoded_allele", "otherallele",
             "baseline_allele", "non_coded_allele"],
    "FRQ":  ["frq", "freq", "eaf", "effect_allele_frequency", "maf",
             "a1freq", "freq1", "freq_a1", "af_alt", "af_coded",
             "eaf_hapmap", "allelefreq", "af", "coded_af", "alt_af",
             "effect_allele_freq", "frq_a1", "a1_freq"],
    "BETA": ["beta", "effect", "b", "log_odds", "logor", "effect_size",
             "est", "all_inv_var_meta_beta", "frequentist_add_beta_1",
             "meta_beta", "gwas_beta", "effect_weight"],
    "OR":   ["or", "odds_ratio", "oddsratio", "or_random", "or_fixed"],
    "SE":   ["se", "standard_error", "stderr", "sebeta", "log_odds_se",
             "se_beta", "all_inv_var_meta_sebeta", "frequentist_add_se_1",
             "meta_se", "gwas_se", "se_effect"],
    "P":    ["p", "pval", "pvalue", "p_value", "p_bolt_lmm", "p-value",
             "p.value", "pval_nominal", "frequentist_add_pvalue",
             "all_inv_var_meta_p", "meta_p", "gwas_p", "p_dgc",
             "p_random", "p_fixed", "pval_meta"],
    "N":    ["n", "n_total", "samplesize", "n_complete_samples", "neff",
             "n_eff", "n_samples", "ntotal", "total_n", "nmiss",
             "num_samples", "weight", "n_analyzed"],
    "INFO": ["info", "imputation_info", "rsq", "r2", "imp_quality",
             "info_score", "impinfo", "info_type0", "metric"],
}


def log(msg):
    print(f"  {msg}", flush=True)


def sniff_sep(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        head = fh.readline()
    for sep, name in ((",", "comma"), ("\t", "tab"), (" ", "space")):
        if head.count(sep) >= 4:
            return (sep, name)
    return ("\t", "tab")


def map_columns(cols):
    """Return {standard_name: original_name} using the alias table."""
    lower = {c.lower().strip(): c for c in cols}
    found = {}
    for std, opts in ALIASES.items():
        for o in opts:
            if o in lower:
                found[std] = lower[o]
                break
    return found


def effective_n(ncase, ncontrol):
    """Sum of effective sample size, 4/(1/ncase + 1/ncontrol). [CDG3] uses
    effective N for the liability-scale ascertainment correction."""
    return 4.0 / (1.0 / float(ncase) + 1.0 / float(ncontrol))


def check_build_heuristic(df, expected_build):
    """Heuristic build detection: check whether BP coordinates are consistent
    with the expected genome build. Issues a WARNING, never blocks. [LDSC rule 4]"""
    if "CHR" not in df or "BP" not in df:
        return
    if expected_build not in ("hg19", "GRCh37"):
        return  # only check for hg19 expected

    # Check if any chr has positions beyond hg19 chromosome lengths
    chr6 = df[df["CHR"].astype(str).str.replace("chr", "") == "6"]
    if len(chr6) > 100:
        bp_vals = pd.to_numeric(chr6["BP"], errors="coerce").dropna()
        if len(bp_vals) > 0:
            max_bp = bp_vals.max()
            if max_bp > 175_000_000:
                log(f"WARNING: max BP on chr6 = {max_bp:,.0f} — "
                    f"expected build is {expected_build} but this looks like hg38. "
                    f"[Rule 4: Never silently lift over genome builds]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outdir", default="data/harmonized")
    ap.add_argument("--keep-ambiguous", action="store_true",
                    help="skip strand-ambiguous SNP removal (NOT recommended)")
    ap.add_argument("--neff-mode", choices=["total", "effective"],
                    default="effective",
                    help="'effective' uses 4/(1/Ncase+1/Nctrl) for binary traits "
                         "[CDG3]; 'total' uses raw N_total. Default: effective")
    ap.add_argument("--build-check", action="store_true",
                    help="run heuristic build detection against config build")
    args = ap.parse_args()

    cfg = pd.read_csv(args.config, sep="\t", dtype=str).set_index("trait_id")
    if args.trait not in cfg.index:
        sys.exit(f"ERROR: trait '{args.trait}' not in {args.config}")
    meta = cfg.loc[args.trait]

    sep, sepname = sniff_sep(args.infile)
    log(f"delimiter detected: {sepname}")
    df = pd.read_csv(args.infile, sep=sep, dtype=str, low_memory=False,
                     comment=None, engine="c")
    n0 = len(df)
    log(f"read {n0:,} rows, {len(df.columns)} columns")

    cmap = map_columns(df.columns)
    need = ["SNP", "A1", "A2", "P"]
    missing = [c for c in need if c not in cmap]
    if missing:
        sys.exit(f"ERROR: could not map required columns {missing}.\n"
                 f"       file has: {list(df.columns)}\n"
                 f"       -> add the right alias to ALIASES in this script.")
    if "BETA" not in cmap and "OR" not in cmap:
        sys.exit("ERROR: neither BETA nor OR present; cannot get effect size.")
    log(f"column map: { {k: v for k, v in cmap.items()} }")

    out = pd.DataFrame(index=df.index)
    for std, orig in cmap.items():
        out[std] = df[orig]

    # --- Build check (heuristic, non-blocking) --------------------------------
    if args.build_check:
        check_build_heuristic(out, meta.get("build", "hg19"))

    # --- effect size to log-odds / beta --------------------------------------
    if "BETA" in out:
        out["BETA"] = pd.to_numeric(out["BETA"], errors="coerce")
    else:
        orv = pd.to_numeric(out["OR"], errors="coerce")
        out["BETA"] = np.log(orv.where(orv > 0))
        log("converted OR -> log(OR)")
    out = out.drop(columns=[c for c in ["OR"] if c in out])

    for c in ["BP", "FRQ", "SE", "P", "N", "INFO"]:
        if c in out:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    # --- alleles --------------------------------------------------------------
    out["A1"] = out["A1"].str.upper().str.strip()
    out["A2"] = out["A2"].str.upper().str.strip()

    steps = []

    def drop(mask, reason):
        nonlocal out
        k = int(mask.sum())
        if k:
            out = out[~mask]
        steps.append((reason, k, len(out)))
        return k

    drop(~out["A1"].isin(VALID) | ~out["A2"].isin(VALID),
         "non-SNP / indel / multi-char allele")
    drop(out["A1"] == out["A2"], "A1 == A2")

    if not args.keep_ambiguous:
        amb = [(a, b) in AMBIG for a, b in zip(out["A1"], out["A2"])]
        drop(pd.Series(amb, index=out.index), "strand-ambiguous (A/T, C/G) [CDG3]")
    else:
        steps.append(("strand-ambiguous KEPT (--keep-ambiguous)", 0, len(out)))

    drop(out["P"].isna() | (out["P"] <= 0) | (out["P"] > 1), "P missing or out of (0,1]")
    drop(out["BETA"].isna() | ~np.isfinite(out["BETA"]), "BETA missing / non-finite")
    if "SE" in out:
        drop(out["SE"].isna() | (out["SE"] <= 0), "SE missing or <= 0")

    if "INFO" in out:
        drop(out["INFO"].notna() & (out["INFO"] <= INFO_MIN),
             f"INFO <= {INFO_MIN} [CDG3]")
    else:
        steps.append(("INFO column absent — filter skipped", 0, len(out)))

    if "FRQ" in out:
        maf = out["FRQ"].where(out["FRQ"] <= 0.5, 1 - out["FRQ"])
        drop(out["FRQ"].notna() & (maf <= MAF_MIN), f"MAF <= {MAF_MIN} [CDG3]")
    else:
        steps.append(("FRQ column absent — MAF filter skipped", 0, len(out)))

    # --- MHC ------------------------------------------------------------------
    if "CHR" in out and "BP" in out:
        chrs = out["CHR"].astype(str).str.replace("chr", "", case=False, regex=False)
        out["CHR"] = chrs
        inmhc = (chrs == MHC[0]) & out["BP"].between(MHC[1], MHC[2])
        drop(inmhc, f"in MHC chr{MHC[0]}:{MHC[1]}-{MHC[2]} [CDG3]")
    else:
        steps.append(("CHR/BP absent — MHC not excluded (FIX THIS)", 0, len(out)))

    drop(out["SNP"].duplicated(keep="first"), "duplicate SNP ID")

    # --- SNP-specific N_eff filtering [CDG3] ----------------------------------
    if "N" in out and out["N"].notna().any():
        n_median = out["N"].median()
        snp_neff_mask = out["N"].notna() & (out["N"] < NEFF_SNP_FRAC * n_median)
        drop(snp_neff_mask,
             f"SNP-specific N < {NEFF_SNP_FRAC*100:.0f}% of median N [CDG3]")

    # --- sample size & effective N -------------------------------------------
    calc_neff = None
    if meta["type"] == "binary" and str(meta.get("ncase", "")) not in ("NA", "", "nan", "None"):
        calc_neff = effective_n(meta["ncase"], meta["ncontrol"])

    if "N" not in out or out["N"].isna().all():
        if args.neff_mode == "effective" and calc_neff is not None:
            out["N"] = calc_neff
            log(f"N absent -> using effective N = {calc_neff:,.0f} [CDG3]")
        elif calc_neff is not None and args.neff_mode == "total":
            out["N"] = float(meta["n_total"])
            log(f"N absent -> using n_total = {float(meta['n_total']):,.0f} (--neff-mode total)")
        else:
            out["N"] = float(meta["n_total"])
            log(f"N absent -> using n_total = {float(meta['n_total']):,.0f}")

    for c in ["CHR", "BP", "FRQ", "INFO"]:
        if c not in out:
            out[c] = np.nan

    out = out[["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]]

    os.makedirs(args.outdir, exist_ok=True)
    dest = os.path.join(args.outdir, f"{args.trait}.harmonized.tsv.gz")
    out.to_csv(dest, sep="\t", index=False, na_rep="NA", compression="gzip")

    # --- QC report with summary statistics ------------------------------------
    rpt = os.path.join(args.outdir, f"{args.trait}.qc.txt")
    with open(rpt, "w") as fh:
        fh.write(f"trait\t{args.trait}\nsource\t{meta.get('source_note', '')}\n")
        fh.write(f"build\t{meta.get('build', 'hg19')}\n")
        fh.write(f"n_total\t{meta.get('n_total', 'NA')}\n")
        fh.write(f"neff_mode\t{args.neff_mode}\n")
        if calc_neff is not None:
            fh.write(f"n_eff\t{calc_neff:.2f}\n")
        fh.write(f"infile\t{args.infile}\nrows_in\t{n0}\n")
        fh.write("\nstep\tdropped\tremaining\n")
        for r, k, rem in steps:
            fh.write(f"{r}\t{k}\t{rem}\n")
        fh.write(f"\nrows_out\t{len(out)}\n")
        fh.write(f"pct_retained\t{100*len(out)/max(n0,1):.2f}\n")

        # Summary statistics
        fh.write("\n--- Summary Statistics ---\n")
        if "INFO" in out and out["INFO"].notna().any():
            fh.write(f"median_INFO\t{out['INFO'].median():.4f}\n")
        if "FRQ" in out and out["FRQ"].notna().any():
            maf_out = out["FRQ"].where(out["FRQ"] <= 0.5, 1 - out["FRQ"])
            fh.write(f"median_MAF\t{maf_out.median():.4f}\n")
        if "CHR" in out and out["CHR"].notna().any():
            fh.write("\nchr\tsnp_count\n")
            for ch, cnt in out["CHR"].value_counts().sort_index().items():
                fh.write(f"{ch}\t{cnt}\n")

    log(f"kept {len(out):,} / {n0:,} ({100*len(out)/max(n0,1):.1f}%)")
    log(f"wrote {dest}")
    log(f"wrote {rpt}")


if __name__ == "__main__":
    main()
