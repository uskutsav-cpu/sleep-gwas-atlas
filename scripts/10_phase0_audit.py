#!/usr/bin/env python3
"""Strict Phase 0 audit.

The workflow brief requires "the Phase 0 strict audit" as a gate on promoting a
trait to CURATED and as a pre-commit check. No such audit existed in the
repository, so this is it.

The audit is deliberately unforgiving in one direction only: it never promotes
anything, it never edits anything, and it exits non-zero if any ISSUE is found.
A trait that cannot be proven ready stays TODO. That asymmetry is the point --
this project has already shipped one artefact whose numbers were invented and
whose inputs did not exist.

    python3 scripts/10_phase0_audit.py                  # audit everything
    python3 scripts/10_phase0_audit.py --block neuro    # only owned traits
    python3 scripts/10_phase0_audit.py --check-build data/raw/ibd.txt.gz

Exit status: 0 = clean, 1 = at least one ISSUE.
"""
import argparse
import csv
import gzip
import os
import re
import sys
import zipfile

TRAITS = "config/traits.tsv"
SOURCES = "config/public_gwas_sources.tsv"
ANCHORS = "config/build_anchors.tsv"

# The 21 trait IDs owned by this workstream.
OWNED = ["alz", "parkinson", "mdd", "scz", "bipolar", "adhd",
         "ibd", "crohn", "uc", "ra", "ms", "asthma",
         "bmi", "t2d", "ldl", "hdl", "triglycerides", "cad", "stroke",
         "atrial_fibrillation", "sbp"]

VALID_ACCESS = {"PUBLIC", "REGISTRATION", "CONTROLLED", "UNAVAILABLE"}
MIN_DELTA = 1000        # anchors closer than this do not discriminate builds
CITATION_RE = re.compile(r"PMID[:\s]?\d{6,8}|doi[:\s]|10\.\d{4,}/", re.I)

ISSUES, NOTES = [], []


def issue(tag, msg):
    ISSUES.append(f"[{tag}] {msg}")


def note(tag, msg):
    NOTES.append(f"[{tag}] {msg}")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------- build check
def load_anchors():
    out = []
    with open(ANCHORS) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("rsid"):
                continue
            parts = line.rstrip("\n").split("\t")
            rs, c, p19, p38 = parts[0], parts[1], parts[2], parts[3]
            al = parts[4] if len(parts) > 4 else ""
            if abs(int(p38) - int(p19)) >= MIN_DELTA:
                out.append((rs, c, int(p19), int(p38), al))
    return out


def open_any(path):
    # .tar.gz must be tested before .gz: gzip.open on a tarball yields the tar
    # container, whose first "line" is a 512-byte header, not data.
    if path.endswith((".tar.gz", ".tgz")):
        import tarfile
        t = tarfile.open(path)
        members = [m for m in t.getmembers()
                   if m.isfile() and not m.name.upper().endswith("README")]
        members.sort(key=lambda m: m.size, reverse=True)
        f = t.extractfile(members[0])
        return f, members[0].name
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    if path.endswith(".zip"):
        z = zipfile.ZipFile(path)
        names = [n for n in z.namelist() if not n.endswith("/")]
        return z.open(names[0]), names[0]
    return open(path, errors="replace")


def check_build(path, max_lines=4_000_000):
    """Decide GRCh37 vs GRCh38 from coordinates in the file itself.

    Not from the filename. Scans for anchor rsIDs and compares the recorded
    base-pair position against both builds. Returns (verdict, hits19, hits38).
    """
    anchors = load_anchors()
    by_rs = {a[0]: a for a in anchors}
    fh = open_any(path)
    member = None
    if isinstance(fh, tuple):
        fh, member = fh
        fh = (l.decode("utf8", "replace") for l in fh)
    header = next(fh).rstrip("\n")
    sep = "\t" if header.count("\t") >= 2 else (
        "," if header.count(",") >= 2 else None)
    cols = [c.strip().lower() for c in (header.split(sep) if sep else header.split())]
    def find(names):
        for i, c in enumerate(cols):
            if c in names:
                return i
        return None
    i_rs = find({"snp", "rsid", "rs_id", "variant_id", "markername", "snpid",
                 "marker", "id", "rs"})
    i_bp = find({"bp", "pos", "position", "base_pair_location", "bp_hg19",
                 "pos_hg19", "bpos", "chrompos", "genpos"})
    i_chr = find({"chr", "chrom", "chromosome", "#chrom", "chr_id", "#chr"})

    # Some releases carry no rsID and no separate position column: the marker
    # itself is the coordinate, e.g. IIBDGC de Lange 2017 uses
    # "1:100000012_G_T". That is still a perfectly good build signal, so parse
    # it rather than giving up. Matching is then on chr:pos, not rsID.
    coord_mode = False
    if i_bp is None and i_rs is not None:
        coord_mode = True
    if i_rs is None and i_bp is None:
        return ("INDETERMINATE: no usable marker or position column "
                f"(header: {cols[:10]})", 0, 0), member

    if coord_mode:
        # Position alone is NOT discriminating here: a densely imputed file
        # contains a variant at both builds' coordinates for any anchor, so
        # naive position matching returns a conflict every time. Require the
        # alleles to agree as well.
        pos19 = {f"{c}:{p19}": set(al.split("/")) for _, c, p19, _, al in anchors}
        pos38 = {f"{c}:{p38}": set(al.split("/")) for _, c, _, p38, al in anchors}
        h19 = h38 = seen = 0
        rxc = re.compile(r"^(?:chr)?(\d{1,2}|X|Y):(\d+)[_:]([ACGTacgt]+)[_:]([ACGTacgt]+)")
        for n, line in enumerate(fh):
            if n > max_lines:
                break
            f = line.split(sep) if sep else line.split()
            if len(f) <= i_rs:
                continue
            m = rxc.match(f[i_rs].strip())
            if not m:
                continue
            key = f"{m.group(1)}:{m.group(2)}"
            obs = {m.group(3).upper(), m.group(4).upper()}
            if key in pos19 and obs <= pos19[key]:
                h19 += 1; seen += 1
            elif key in pos38 and obs <= pos38[key]:
                h38 += 1; seen += 1
        if seen == 0:
            return ("INDETERMINATE: coordinate-keyed markers parsed, but no "
                    "anchor position matched either build", 0, 0), member
        if h19 and not h38:
            return (f"GRCh37 ({h19}/{seen} anchor coordinates matched hg19, "
                    f"parsed from chr:pos markers)", h19, h38), member
        if h38 and not h19:
            return (f"GRCh38 ({h38}/{seen} anchor coordinates matched hg38, "
                    f"parsed from chr:pos markers)", h19, h38), member
        return (f"CONFLICT: {h19} hg19 and {h38} hg38 coordinate matches",
                h19, h38), member
    # Preferred when the file carries explicit CHR and BP columns: match the
    # anchor on chromosome+position regardless of how the marker is named.
    # Needed for releases whose "MarkerName" is itself a coordinate
    # (e.g. "6:28571110:T"), where rsID matching can never fire.
    if i_chr is not None and i_bp is not None:
        # Alleles are required here for the same reason as in coordinate mode:
        # a dense file has a variant at both builds' positions, so position
        # alone produces a false CONFLICT.
        i_a1 = find({"a1", "allele1", "effect_allele", "ea", "alt",
                     "tested_allele", "effectallele"})
        i_a2 = find({"a2", "allele2", "other_allele", "nea", "ref",
                     "non_effect_allele", "otherallele"})
        p19 = {(c, p): set(al.split("/")) for _, c, p, _, al in anchors}
        p38 = {(c, p): set(al.split("/")) for _, c, _, p, al in anchors}
        h19 = h38 = seen = 0
        for n, line in enumerate(fh):
            if n > max_lines:
                break
            f = line.split(sep) if sep else line.split()
            if len(f) <= max(i_chr, i_bp):
                continue
            try:
                key = (str(f[i_chr]).strip().replace("chr", ""),
                       int(float(f[i_bp])))
            except ValueError:
                continue
            if key not in p19 and key not in p38:
                continue
            ok = True
            if i_a1 is not None and i_a2 is not None and len(f) > max(i_a1, i_a2):
                obs = {f[i_a1].strip().upper(), f[i_a2].strip().upper()}
                exp = p19.get(key) or p38.get(key)
                ok = obs <= exp
            if not ok:
                continue
            if key in p19:
                h19 += 1; seen += 1
            elif key in p38:
                h38 += 1; seen += 1
        if seen:
            if h19 and not h38:
                return (f"GRCh37 ({h19}/{seen} anchors matched hg19 on chr+pos)",
                        h19, h38), member
            if h38 and not h19:
                return (f"GRCh38 ({h38}/{seen} anchors matched hg38 on chr+pos)",
                        h19, h38), member
            return (f"CONFLICT: {h19} hg19 vs {h38} hg38 on chr+pos",
                    h19, h38), member

    h19 = h38 = seen = 0
    for n, line in enumerate(fh):
        if n > max_lines:
            break
        f = line.split(sep) if sep else line.split()
        if len(f) <= max(i_rs, i_bp):
            continue
        rs = f[i_rs].strip()
        if rs not in by_rs:
            continue
        try:
            bp = int(float(f[i_bp]))
        except ValueError:
            continue
        seen += 1
        _, _, p19, p38, _al = by_rs[rs]
        if bp == p19:
            h19 += 1
        elif bp == p38:
            h38 += 1
    if seen == 0:
        return ("INDETERMINATE: no anchor SNPs found in scanned region", 0, 0), member
    if h19 and not h38:
        return (f"GRCh37 ({h19}/{seen} anchors matched hg19)", h19, h38), member
    if h38 and not h19:
        return (f"GRCh38 ({h38}/{seen} anchors matched hg38)", h19, h38), member
    if h19 and h38:
        return (f"CONFLICT: {h19} hg19 and {h38} hg38 matches", h19, h38), member
    return (f"INDETERMINATE: {seen} anchors found, none matched either build",
            h19, h38), member


# --------------------------------------------------------------- audits
def audit_registry(traits, block_only):
    ids = [t["trait_id"] for t in traits]
    dupes = {i for i in ids if ids.count(i) > 1}
    if dupes:
        issue("REG", f"duplicate trait_id(s): {sorted(dupes)}")
    for t in traits:
        if len(t) != 12:
            issue("REG", f"{t['trait_id']}: expected 12 columns, got {len(t)}")

    # identical numeric fingerprints = same phenotype under two IDs
    seen = {}
    for t in traits:
        fp = (t["ncase"], t["ncontrol"], t["n_total"], t["type"])
        if fp[0] in ("NA", "") and fp[2] in ("NA", ""):
            continue
        seen.setdefault(fp, []).append(t["trait_id"])
    for fp, names in seen.items():
        if len(names) <= 1:
            continue
        owned = [n for n in names if n in OWNED]
        # Two different things share a fingerprint:
        #   (a) the SAME phenotype under two IDs -> a genuine duplicate row;
        #   (b) DIFFERENT phenotypes measured in one cohort (LDL/HDL/TG in the
        #       same 1.3M people) -> sample overlap, which is legitimate and is
        #       handled by the LDSC bivariate intercept, not by deletion.
        # Only (a) is an error. Detect it by a shared name stem.
        stems = {n.split("_")[0] for n in names}
        dup = len(stems) < len(names)
        if dup:
            (issue if owned else note)(
                "DUP", f"identical fingerprint {fp[:3]} shared by {names} -- "
                       f"same phenotype under two trait_ids")
        else:
            note("OVERLAP",
                 f"{names} share sample-size fingerprint {fp[:3]}: distinct "
                 f"phenotypes from one cohort. Legitimate, but they are NOT "
                 f"independent -- rely on the LDSC gcov intercept and do not "
                 f"present them as independent tests")


def audit_traits(traits, sources, block_only):
    smap = {s["trait_id"]: s for s in sources}
    tmap = {t["trait_id"]: t for t in traits}
    targets = OWNED if block_only else [t["trait_id"] for t in traits]

    for tid in targets:
        t = tmap.get(tid)
        if t is None:
            issue("MISSING", f"{tid}: owned trait absent from {TRAITS}")
            continue
        st = t["status"]

        if tid in OWNED:
            s = smap.get(tid)
            if s is None:
                issue("SRC", f"{tid}: owned trait has no row in {SOURCES}")
                continue
            if s["access"] not in VALID_ACCESS:
                issue("SRC", f"{tid}: access={s['access']!r} not one of "
                             f"{sorted(VALID_ACCESS)}")
            if s["pmid"] in ("", "NA", "UNRESOLVED") and st == "CURATED":
                issue("CITE", f"{tid}: CURATED without a resolved PMID")

            if st == "CURATED":
                # every promotion condition, enforced
                if s["access"] != "PUBLIC":
                    issue("GATE", f"{tid}: CURATED but access={s['access']}")
                if s["build_verified"] != "GRCh37":
                    issue("GATE", f"{tid}: CURATED but build_verified="
                                  f"{s['build_verified']!r} (must be GRCh37, "
                                  f"proven from the file)")
                if s["sha256"] in ("", "NA", "PENDING", "NOT_DOWNLOADED"):
                    issue("GATE", f"{tid}: CURATED but no SHA-256 recorded")
                local = os.path.join("data/raw", s["local_file"])
                if s["local_file"] in ("NA", "") or not os.path.exists(local):
                    issue("GATE", f"{tid}: CURATED but local file missing "
                                  f"({s['local_file']})")

        # A binary pop_prev needs a citation for THE PREVALENCE, which is a
        # different claim from the GWAS citation. A PMID in source_note cites
        # the study, not the epidemiology, so it must not satisfy this check.
        # The dedicated field lives in public_gwas_sources.tsv.
        if t["type"] == "binary" and t["pop_prev"] not in ("NA", "", "UNKNOWN"):
            s = smap.get(tid, {})
            cite = (s.get("pop_prev_citation") or "").strip()
            if not cite or cite in ("NA", "PENDING", "UNKNOWN"):
                (issue if tid in OWNED else note)(
                    "PREV", f"{tid}: pop_prev={t['pop_prev']} with no prevalence "
                            f"citation -- rule 1 requires a real source or UNKNOWN")
            elif not CITATION_RE.search(cite):
                (issue if tid in OWNED else note)(
                    "PREV", f"{tid}: pop_prev_citation={cite!r} carries no "
                            f"PMID/DOI -- not a verifiable source")

        if st not in ("TODO", "CURATED", "PASS", "DROP"):
            issue("STATUS", f"{tid}: unrecognised status {st!r}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", action="store_true",
                    help="restrict issue-level checks to the owned 21 traits")
    ap.add_argument("--check-build", metavar="FILE",
                    help="run the discriminating build check on one file and exit")
    a = ap.parse_args()

    if a.check_build:
        (verdict, h19, h38), member = check_build(a.check_build)
        print(f"{a.check_build}")
        if member:
            print(f"  archive member: {member}")
        print(f"  build verdict : {verdict}")
        sys.exit(0 if verdict.startswith("GRCh") else 1)

    traits = read_tsv(TRAITS)
    sources = read_tsv(SOURCES) if os.path.exists(SOURCES) else []
    if not sources:
        issue("SRC", f"{SOURCES} missing or empty")

    audit_registry(traits, a.block)
    audit_traits(traits, sources, a.block)

    print("=" * 72)
    print(f"PHASE 0 STRICT AUDIT   traits={len(traits)}  sources={len(sources)}"
          f"  scope={'owned block' if a.block else 'whole registry'}")
    print("=" * 72)
    counts = {}
    for t in traits:
        counts[t["status"]] = counts.get(t["status"], 0) + 1
    print("status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if sources:
        acc = {}
        for s in sources:
            acc[s["access"]] = acc.get(s["access"], 0) + 1
        print("access: " + ", ".join(f"{k}={v}" for k, v in sorted(acc.items())))
    print()
    if NOTES:
        print(f"NOTES ({len(NOTES)}) -- outside owned block or advisory:")
        for n in NOTES:
            print(f"  {n}")
        print()
    if ISSUES:
        print(f"ISSUES ({len(ISSUES)}):")
        for i in ISSUES:
            print(f"  {i}")
        print("\nAUDIT FAILED")
        sys.exit(1)
    print("no issues\nAUDIT PASSED")


if __name__ == "__main__":
    main()
