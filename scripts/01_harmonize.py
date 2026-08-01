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
    # Explicit rsID names FIRST: a file may carry a coordinate-style "SNP"
    # column AND a real rsID column (Watanabe 2022 insomnia does). LDSC
    # merges against HapMap3 by rsID, so the rsID column must win.
    "SNP":  ["rsid", "rs_id", "rsid_ukb", "rsid_dbsnp", "rsids", "snp",
             "markername", "variant_id", "id", "marker", "snpid",
             "varid", "snp_id", "#snpid", "rs"],
    "CHR":  ["chr", "chrom", "chromosome", "#chrom", "hg19chr", "chr_id",
             "chr_name", "#chr", "seqnames"],
    "BP":   ["bp", "pos", "position", "base_pair_location", "bp_hg19",
             "pos_hg19", "bpos", "chrompos", "pos_b37", "bp_grch37",
             "bp_grch38", "genpos", "start"],
    "A1":   ["a1", "effect_allele", "effect_all", "ea", "allele1", "tested_allele", "alt",
             "a1_effect", "risk_allele", "inc_allele", "effectallele",
             "reference_allele", "coded_allele"],
    "A2":   ["a2", "a0", "other_allele", "other_all", "nea", "allele2", "allele0", "non_effect_allele", "ref",
             "a2_other", "noneffect_allele", "noncoded_allele", "otherallele",
             "baseline_allele", "non_coded_allele"],
    "FRQ":  ["frq", "freq", "eaf", "effect_allele_frequency", "maf",
             "a1freq", "freq1", "effect_allele_freq", "eaf_ukb",
             "freq_tested_allele", "freq_tested_allele_in_hrs",
             "pooled_alt_af", "af_allele2", "all_meta_af"],
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
    "INFO": ["info", "info_ukb", "imputation_info", "rsq", "r2", "imp_quality",
             "info_score", "impinfo", "info_type0", "metric"],
}


def log(msg):
    print(f"  {msg}", flush=True)


def _open_text(path):
    """Open .gz, .zip or plain text uniformly.

    Several core sleep releases (Dashti 2019, Wang 2019) ship a single .txt
    inside a .zip. Reading them required no other change, so support the
    container rather than making the caller pre-extract 380 MB.
    """
    if path.endswith(".zip"):
        import zipfile, io
        z = zipfile.ZipFile(path)
        members = [n for n in z.namelist()
                   if not n.endswith("/") and not n.startswith("__MACOSX")]
        members.sort(key=lambda n: z.getinfo(n).file_size, reverse=True)
        return io.TextIOWrapper(z.open(members[0]), errors="replace")
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def zip_member(path):
    """Name of the member _open_text would pick, for the provenance ledger."""
    if not path.endswith(".zip"):
        return ""
    import zipfile
    z = zipfile.ZipFile(path)
    m = [n for n in z.namelist()
         if not n.endswith("/") and not n.startswith("__MACOSX")]
    m.sort(key=lambda n: z.getinfo(n).file_size, reverse=True)
    return m[0] if m else ""


def sniff_sep(path):
    with _open_text(path) as fh:
        head = fh.readline()
    for sep, name in ((",", "comma"), ("\t", "tab"), (" ", "space")):
        if head.count(sep) >= 4:
            return (sep, name)
    return ("\t", "tab")


# Effect/SE/P columns are often suffixed with the phenotype name
# (BETA_SLEEPDURATION, SE_SHORTSLEEP, P_LONGSLEEP). Exact aliases cannot
# enumerate those, so these standard names also accept a prefix match --
# but only these, because a prefix rule on e.g. "N" would match anything.
PREFIX_OK = {"BETA": ("beta_", "b_"), "SE": ("se_", "stderr_"),
             "P": ("p_", "pval_", "pvalue_"), "OR": ("or_",)}


def map_columns(cols):
    """Return {standard_name: original_name} using the alias table."""
    lower = {c.lower().strip(): c for c in cols}
    found = {}
    for std, opts in ALIASES.items():
        for o in opts:
            if o in lower:
                found[std] = lower[o]
                break
    for std, prefixes in PREFIX_OK.items():
        if std in found:
            continue
        for lc, orig in lower.items():
            if lc.startswith(prefixes):
                found[std] = orig
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



def _apply_filters(out, args, counts):
    """Apply every QC filter to one frame, accumulating counts by reason.

    Identical logic to the single-pass path; factored out so a chunked run
    produces the same ledger.
    """
    def drop(mask, reason):
        k = int(mask.sum())
        counts[reason] = counts.get(reason, 0) + k
        return out[~mask] if k else out

    out = drop(~out["A1"].isin(VALID) | ~out["A2"].isin(VALID),
               "non-SNP / indel / multi-char allele")
    out = drop(out["A1"] == out["A2"], "A1 == A2")
    if not args.keep_ambiguous:
        amb = pd.Series([(a, b) in AMBIG for a, b in zip(out["A1"], out["A2"])],
                        index=out.index)
        out = drop(amb, "strand-ambiguous (A/T, C/G) [CDG3]")
    out = drop(out["P"].isna() | (out["P"] <= 0) | (out["P"] > 1),
               "P missing or out of (0,1]")
    out = drop(out["BETA"].isna() | ~np.isfinite(out["BETA"]),
               "BETA missing / non-finite")
    if "SE" in out:
        out = drop(out["SE"].isna() | (out["SE"] <= 0), "SE missing or <= 0")
    if "INFO" in out:
        out = drop(out["INFO"].notna() & (out["INFO"] <= INFO_MIN),
                   f"INFO <= {INFO_MIN} [CDG3]")
    if "FRQ" in out:
        maf = out["FRQ"].where(out["FRQ"] <= 0.5, 1 - out["FRQ"])
        out = drop(out["FRQ"].notna() & (maf <= MAF_MIN), f"MAF <= {MAF_MIN} [CDG3]")
    if "CHR" in out and "BP" in out:
        out["CHR"] = out["CHR"].astype(str).str.replace("chr", "", case=False,
                                                        regex=False)
        out = drop((out["CHR"] == MHC[0]) & out["BP"].between(MHC[1], MHC[2]),
                   f"in MHC chr{MHC[0]}:{MHC[1]}-{MHC[2]} [CDG3]")
    out = drop(out["SNP"].duplicated(keep="first"), "duplicate SNP ID (within chunk)")
    return out


def run_chunked(args, meta, sep, cmap, strcols):
    """Stream the file in chunks, filtering each and appending to the output.

    Two documented differences from the single-pass path, both forced by not
    holding the file in memory:
      * duplicate rsIDs are removed WITHIN each chunk, not across the file.
        LDSC's munge_sumstats removes duplicate rs numbers itself, so this is
        covered downstream rather than skipped.
      * the SNP-specific effective-N filter uses the median N of the FIRST
        chunk as the reference, not the whole-file median.
    Both are recorded in the QC ledger so the run stays auditable.
    """
    import gzip as _gzip
    counts, n_in, n_out = {}, 0, 0
    dest = os.path.join(args.outdir, f"{args.trait}.harmonized.tsv.gz")
    os.makedirs(args.outdir, exist_ok=True)
    cols_final, n_median, wrote_header = None, None, False
    any_frq = False

    reader = pd.read_csv(_open_text(args.infile), sep=sep, usecols=list(cmap.values()),
                         dtype={c: str for c in strcols}, chunksize=args.chunksize,
                         low_memory=False, engine="c")
    fh = _gzip.open(dest, "wt")
    for chunk in reader:
        n_in += len(chunk)
        out = pd.DataFrame(index=chunk.index)
        for std, orig in cmap.items():
            out[std] = chunk[orig]
        if "BETA" in out:
            out["BETA"] = pd.to_numeric(out["BETA"], errors="coerce")
        else:
            orv = pd.to_numeric(out["OR"], errors="coerce")
            out["BETA"] = np.log(orv.where(orv > 0))
        out = out.drop(columns=[c for c in ["OR"] if c in out])
        for c in ["BP", "FRQ", "SE", "P", "N", "INFO"]:
            if c in out:
                out[c] = pd.to_numeric(out[c], errors="coerce")
        out["A1"] = out["A1"].str.upper().str.strip()
        out["A2"] = out["A2"].str.upper().str.strip()

        out = _apply_filters(out, args, counts)

        if "N" in out and out["N"].notna().any():
            if n_median is None:
                n_median = out["N"].median()
            m = out["N"].notna() & (out["N"] < NEFF_SNP_FRAC * n_median)
            k = int(m.sum())
            counts[f"SNP-specific N < {NEFF_SNP_FRAC*100:.0f}% of first-chunk median N [CDG3]"] = \
                counts.get(f"SNP-specific N < {NEFF_SNP_FRAC*100:.0f}% of first-chunk median N [CDG3]", 0) + k
            if k:
                out = out[~m]
        if "N" not in out or out["N"].isna().all():
            out["N"] = float(meta["n_total"])
        for c in ["CHR", "BP", "FRQ", "INFO"]:
            if c not in out:
                out[c] = np.nan
        if cols_final is None:
            cols_final = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
            if "FRQ" not in out or out["FRQ"].notna().sum() == 0:
                cols_final.remove("FRQ")
                counts["FRQ empty in first chunk — column omitted so munge does "
                       "not MAF-filter every SNP away"] = 0
        if "FRQ" in cols_final:
            any_frq = any_frq or bool(out["FRQ"].notna().any())
        out = out[cols_final]
        out.to_csv(fh, sep="\t", index=False, na_rep="NA", header=not wrote_header)
        wrote_header = True
        n_out += len(out)
        log(f"  chunk done: {n_in:,} read, {n_out:,} kept")
    fh.close()

    rpt = os.path.join(args.outdir, f"{args.trait}.qc.txt")
    with open(rpt, "w") as r:
        r.write(f"trait\t{args.trait}\nsource\t{meta['source_note']}\n")
        r.write(f"build\t{meta.get('build','hg19')}\nn_total\t{meta['n_total']}\n")
        r.write(f"mode\tCHUNKED (chunksize={args.chunksize})\n")
        r.write(f"infile\t{args.infile}\nrows_in\t{n_in}\n")
        r.write("\nstep\tdropped\n")
        for k, v in counts.items():
            r.write(f"{k}\t{v}\n")
        r.write(f"\nrows_out\t{n_out}\n")
        r.write(f"pct_retained\t{100*n_out/max(n_in,1):.2f}\n")
        r.write("\nNOTE\tduplicate rsIDs removed within chunks only; LDSC munge "
                "removes cross-file duplicates\n")
    log(f"kept {n_out:,} / {n_in:,} ({100*n_out/max(n_in,1):.1f}%)")
    log(f"wrote {dest}")
    log(f"wrote {rpt}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True)
    ap.add_argument("--config", default="config/traits.tsv")
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outdir", default="data/harmonized")
    ap.add_argument("--chunksize", type=int, default=0,
                    help="process the file in chunks of this many rows. "
                         "Needed for large releases: the in-memory path holds "
                         "the whole frame plus filtering copies, which on an "
                         "8 GB machine is killed outright for a 13M-row file. "
                         "0 (default) keeps the original single-pass path.")
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
    # Read the header alone first, map the columns we actually need, then read
    # only those. The previous approach pulled every column in as Python
    # strings, which for a 13M-row file is roughly 10 GB of RAM -- large
    # releases (mdd, t2d) simply died. Strings are kept only where the value is
    # categorical (IDs, alleles, chromosome); the rest are parsed as numbers.
    hdr = pd.read_csv(_open_text(args.infile), sep=sep, nrows=0, engine="c")
    pre = map_columns(hdr.columns)
    want = {v: k for k, v in pre.items()}
    strcols = {v for k, v in pre.items() if k in ("SNP", "A1", "A2", "CHR")}
    if args.chunksize:
        # Dispatch BEFORE the full read -- reading first and then chunking
        # would already have blown the memory we are trying to avoid.
        need0 = [c for c in ("SNP", "A1", "A2", "P") if c not in pre]
        if need0:
            sys.exit(f"ERROR: could not map required columns {need0} from header")
        if "BETA" not in pre and "OR" not in pre:
            sys.exit("ERROR: neither BETA nor OR present; cannot get effect size.")
        log(f"column map: { {k: v for k, v in pre.items()} }")
        log(f"chunked mode, {args.chunksize:,} rows per chunk")
        run_chunked(args, meta, sep, pre, strcols)
        return

    df = pd.read_csv(_open_text(args.infile), sep=sep, usecols=list(want),
                     dtype={c: str for c in strcols},
                     low_memory=False, comment=None, engine="c")
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

    # An all-missing FRQ column is worse than no column at all: LDSC's
    # munge_sumstats applies --maf-min to it, every value is NA, and the merge
    # silently drops EVERY SNP ("No objects to concatenate"). Emit the column
    # only when it carries data, and record the decision in the ledger.
    cols = ["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
    if "FRQ" not in out or out["FRQ"].notna().sum() == 0:
        cols.remove("FRQ")
        steps.append(("FRQ entirely missing — column omitted so that munge "
                      "does not MAF-filter every SNP away", 0, len(out)))
    out = out[cols]

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
