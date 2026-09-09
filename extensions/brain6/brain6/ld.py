"""LD-block mapping and signed-LD alignment; no SNP-thinning shortcuts."""
from __future__ import annotations
import bisect
import math
import sqlite3
from pathlib import Path
import numpy as np
from .artifacts import transaction, verify_artifact
from .gwas import align, PAIR_FIELDS
from .io import (check_hash, open_text, read_json, read_tsv, require, write_json, write_tsv)


def index_pair(pair_dir, root, relative):
    pair_dir = Path(pair_dir)
    receipt = verify_artifact(pair_dir)
    with transaction(root, relative, stage="index_pair", inputs=[pair_dir / "pair.tsv.gz"],
                     parameters={}, synthetic=receipt["synthetic"]) as (work, meta):
        con = sqlite3.connect(work / "pair.sqlite")
        cols = ",".join(f'"{c}" ' + ("TEXT" if c in {"SNP", "A1", "A2", "alignment"} else "REAL") for c in PAIR_FIELDS)
        con.execute(f"CREATE TABLE pair ({cols})")
        batch = []
        for r in read_tsv(pair_dir / "pair.tsv.gz", PAIR_FIELDS):
            batch.append([None if r[c] == "NA" else r[c] for c in PAIR_FIELDS])
            if len(batch) >= 20000:
                con.executemany(f"INSERT INTO pair VALUES ({','.join('?' for _ in PAIR_FIELDS)})", batch)
                con.commit(); batch.clear()
        con.executemany(f"INSERT INTO pair VALUES ({','.join('?' for _ in PAIR_FIELDS)})", batch)
        con.execute("CREATE UNIQUE INDEX pair_snp ON pair(SNP)")
        con.execute("CREATE INDEX pair_coordinate ON pair(CHR,BP)")
        con.commit(); con.close()
        meta["pair_metadata"] = receipt
    return Path(root) / relative


def read_blocks(path):
    blocks = {}
    names = set()
    for r in read_tsv(path, ["LOC", "CHR", "START", "STOP"]):
        ch, start, stop = int(r["CHR"]), int(r["START"]), int(r["STOP"])
        require(ch in range(1, 23) and 1 <= start <= stop, "Invalid block bounds")
        require(r["LOC"] not in names, "Duplicate locus IDs")
        names.add(r["LOC"])
        blocks.setdefault(ch, []).append((start, stop, r["LOC"]))
    for ch in blocks:
        blocks[ch].sort()
        require(all(blocks[ch][i-1][1] < blocks[ch][i][0] for i in range(1, len(blocks[ch]))),
                "LD blocks overlap; independent-locus counting would be ambiguous")
    return blocks


def map_loci(leads, blocks_path, out):
    """Map PLINK-clumped leads to intact, predeclared blocks; merge repeated blocks."""
    blocks = read_blocks(blocks_path)
    starts = {ch: [b[0] for b in rs] for ch, rs in blocks.items()}
    grouped = {}
    for r in read_tsv(leads, ["SNP", "CHR", "BP", "P"]):
        ch, pos = int(r["CHR"]), int(r["BP"])
        i = bisect.bisect_right(starts.get(ch, []), pos) - 1
        require(i >= 0 and pos <= blocks[ch][i][1], f"Lead {r['SNP']} not covered by locked LD blocks")
        lo, hi, name = blocks[ch][i]
        p = float(r["P"])
        require(math.isfinite(p) and 0 <= p <= 1, "Invalid lead P")
        grouped.setdefault((name, ch, lo, hi), []).append((p, r["SNP"]))
    rows = []
    for (name, ch, lo, hi), leadlist in sorted(grouped.items(), key=lambda x: (x[0][1], x[0][2])):
        leadlist.sort()
        rows.append(dict(locus_id=name, CHR=ch, START=lo, STOP=hi,
                         lead_snp=leadlist[0][1], lead_p=leadlist[0][0],
                         clump_leads=";".join(s for _, s in leadlist), n_clump_leads=len(leadlist)))
    if out is not None:
        write_tsv(out, ["locus_id", "CHR", "START", "STOP", "lead_snp", "lead_p", "clump_leads", "n_clump_leads"], rows)
    return rows


def validate_ld(matrix, *, tolerance=1e-6):
    matrix = np.asarray(matrix, dtype=float)
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1] and matrix.shape[0] >= 2,
            "LD matrix must be square with >=2 SNPs")
    require(np.isfinite(matrix).all(), "LD contains missing or nonfinite correlations")
    require(np.max(np.abs(matrix)) <= 1 + tolerance, "LD correlations outside [-1,1]")
    require(np.allclose(matrix, matrix.T, atol=tolerance, rtol=0), "LD not symmetric")
    require(np.allclose(np.diag(matrix), 1, atol=tolerance, rtol=0), "LD diagonal must equal one")
    eigen_min = float(np.linalg.eigvalsh(matrix)[0])
    require(eigen_min >= -tolerance, "LD not positive semidefinite; no automatic eigenvalue repair")
    return eigen_min


def prepare_locus(pair_index, ld_manifest, locus, root, relative, *, max_working_bytes=4*2**30,
                  min_coverage=.8, synthetic=False):
    """Create aligned locus table + signed LD matrix for both trait-specific fits.

    Reference manifest binds matrix bytes, BIM order, counted allele, build and
    ancestry. BOTH GWAS effects flip to the allele counted in LD; R stays fixed.
    """
    index = Path(pair_index)
    im = verify_artifact(index)
    require(im["synthetic"] == synthetic, "Synthetic/real mismatch")
    lm = read_json(ld_manifest)
    require(lm.get("kind") == "signed_r" and lm.get("counted_allele") == "BIM_A1", "Signed allele-count LD required, not r-squared")
    require(lm.get("synthetic", False) == synthetic, "Synthetic reference mismatch")
    pair_meta = im["pair_metadata"]
    for trait in ["left", "right"]:
        require(lm["genome_build"] == pair_meta[trait]["genome_build"] and
                lm["ancestry"] == pair_meta[trait]["ancestry"], "Reference build/ancestry mismatch")
    check_hash(lm["matrix_path"], lm["matrix_sha256"])
    check_hash(lm["bim_path"], lm["bim_sha256"])
    require(0 < min_coverage <= 1, "Invalid locus coverage minimum")
    bim = []
    with open_text(lm["bim_path"]) as f:
        for line in f:
            r = line.split()
            require(len(r) == 6, "Expected six-column PLINK BIM")
            bim.append((r[1], int(r[0]), int(r[3]), r[4], r[5]))
    require(len({r[0] for r in bim}) == len(bim), "Duplicate BIM SNP IDs")
    n = len(bim)
    require(8 * n * n * 8 <= max_working_bytes, "BLOCKED_BY_COMPUTE: full locus matrix exceeds working-memory budget; do not thin")
    R = np.loadtxt(lm["matrix_path"])
    require(R.shape == (n, n), "LD/BIM dimension mismatch")
    min_eig = validate_ld(R)
    ch, lo, hi = int(locus["CHR"]), int(locus["START"]), int(locus["STOP"])
    require(all(r[1] == ch and lo <= r[2] <= hi for r in bim), "Reference matrix is not restricted to the intact locus")
    con = sqlite3.connect(f"{(index / 'pair.sqlite').resolve().as_uri()}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    rows = {r["SNP"]: dict(r) for r in con.execute("SELECT * FROM pair WHERE CHR=? AND BP BETWEEN ? AND ?", (ch, lo, hi))}
    con.close()
    require(len(rows) >= 2, "Insufficient GWAS variants in locus")
    keep, aligned = [], []
    dropped = []
    for i, (snp, c, bp, a1, a2) in enumerate(bim):
        if snp not in rows:
            continue
        row = rows[snp].copy()
        match = align(a1, a2, row["A1"], row["A2"])
        if match is None or (int(row["CHR"]), int(row["BP"])) != (c, bp):
            dropped.append(snp); continue
        sign, _ = match
        row.update(A1=a1, A2=a2, CHR=c, BP=bp)
        for suffix in ["1", "2"]:
            row["BETA" + suffix] *= sign
            row["Z" + suffix] *= sign
            if sign == -1 and row["EAF" + suffix] is not None:
                row["EAF" + suffix] = 1 - row["EAF" + suffix]
        keep.append(i); aligned.append(row)
    coverage = len(aligned) / len(rows)
    require(coverage >= min_coverage and len(aligned) >= 2, f"BLOCKED_BY_LD_COVERAGE: {coverage:.3f}")
    if locus.get("lead_snp"):
        require(locus["lead_snp"] in {r["SNP"] for r in aligned}, "Lead SNP absent from usable LD")
    with transaction(root, relative, stage="prepare_locus",
                     inputs=[index / "pair.sqlite", ld_manifest, lm["matrix_path"], lm["bim_path"]],
                     parameters={"locus": locus, "min_coverage": min_coverage}, synthetic=synthetic) as (work, meta):
        np.savetxt(work / "ld.tsv", R[np.ix_(keep, keep)], delimiter="\t", fmt="%.12g")
        write_tsv(work / "locus.tsv", PAIR_FIELDS,
                  ({k: "NA" if r[k] is None else r[k] for k in PAIR_FIELDS} for r in aligned))
        write_json(work / "qc.json", {"coverage": coverage, "n_gwas": len(rows), "n_kept": len(aligned),
                                      "minimum_eigenvalue": min_eig, "allele_or_coordinate_failures": dropped,
                                      "ld_unavailable_snps": sorted(set(rows) - {r['SNP'] for r in aligned})})
        meta["coverage"] = coverage
    return Path(root) / relative
