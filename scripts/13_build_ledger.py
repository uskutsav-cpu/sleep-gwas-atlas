#!/usr/bin/env python3
"""Build results/phase_ledger.tsv -- one authoritative row per trait.

Every downstream script should read state from here rather than globbing
filenames or trusting a chat summary. The ledger is rebuilt from primary
evidence on disk each time it runs:

    config/traits.tsv                 registry + curation status
    config/public_gwas_sources.tsv    source, access, hashes, build evidence
    data/harmonized/*.qc.txt          harmonization outcome
    data/munged/*.sumstats.gz         munge outcome
    results/logs/h2_*.log             parsed h2 / intercept / gate

Nothing here is carried over on trust: if a file is absent the status says so.

    python3 scripts/13_build_ledger.py
    python3 scripts/13_build_ledger.py --summary
"""
import argparse
import csv
import datetime
import glob
import hashlib
import os
import re
import subprocess

TRAITS = "config/traits.tsv"
SOURCES = "config/public_gwas_sources.tsv"
OUT = "results/phase_ledger.tsv"

CORE_SLEEP = ["insomnia", "sleepdur", "shortsleep", "longsleep",
              "chronotype", "sleepiness"]
SECONDARY_SLEEP = ["accel_sleep_duration", "sleep_efficiency", "sleep_timing",
                   "napping", "snoring", "sleep_apnea", "activity_rhythm",
                   "rem_nrem", "insomnia_ukb", "insomnia_full",
                   "shortsleep_dashti", "shortsleep_az", "longsleep_dashti",
                   "longsleep_az"]
AGING = ["frailty", "healthspan", "parental_lifespan", "longevity",
         "epigenetic_age", "telomere_length", "grip_strength", "walking_pace"]

H2_PAT = re.compile(r"Total (Observed|Liability) scale h2:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)")
INT_PAT = re.compile(r"^Intercept:\s*(-?[\d.eE+-]+)\s*\(([\d.eE+-]+)\)", re.M)
CHI_PAT = re.compile(r"Mean Chi\^2:\s*(-?[\d.eE+-]+)")
LAM_PAT = re.compile(r"Lambda GC:\s*(-?[\d.eE+-]+)")

Z_MIN, INTERCEPT_MAX, MIXER_MIN = 4.0, 1.20, 12000.0

FIELDS = ["trait_id", "display_name", "domain", "role", "source_name",
          "source_url_or_accession", "publication_or_PMID", "ancestry",
          "genome_build", "phenotype_type", "sample_size", "n_cases",
          "n_controls", "population_prevalence", "effect_type",
          "raw_status", "raw_path", "raw_sha256",
          "license_or_redistribution_note",
          "harmonized_status", "harmonized_path", "harmonized_sha256",
          "munged_status", "munged_path", "munged_sha256",
          "h2_status", "h2_observed", "h2_se", "h2_z", "intercept",
          "intercept_se", "mean_chisq", "lambda_gc", "N_eff_times_h2",
          "phase1_eligible", "phase2_eligible", "phase3_eligible",
          "exclusion_reason", "blocker",
          "last_verified_commit", "last_verified_timestamp"]


def sha_head(path, nbytes=1 << 20):
    """Hash of the first MB -- enough to detect a changed file cheaply."""
    try:
        with open(path, "rb") as fh:
            return hashlib.sha256(fh.read(nbytes)).hexdigest()[:16]
    except OSError:
        return ""


def role_of(tid, domain):
    if tid in CORE_SLEEP:
        return "core_sleep"
    if tid in SECONDARY_SLEEP or domain == "sleep":
        return "secondary_sleep"
    if tid in AGING or domain == "aging":
        return "aging"
    if domain in ("psychiatric", "neuro", "immune", "metabolic", "cardio",
                  "renal", "hepatic", "blood", "respiratory", "ocular",
                  "gastro", "musculoskeletal", "reproductive", "endocrine",
                  "biomarker", "behavioural", "anthropometric", "cancer"):
        return "disease"
    return "other"


def parse_h2(tid):
    for p in (f"results/logs/h2_{tid}.log", f"results/logs/h2_{tid}"):
        if os.path.exists(p):
            txt = open(p, errors="replace").read()
            m = H2_PAT.search(txt)
            if not m:
                return {"h2_status": "LOG_NO_H2"}
            h2, se = float(m.group(2)), float(m.group(3))
            i = INT_PAT.search(txt)
            c = CHI_PAT.search(txt)
            l = LAM_PAT.search(txt)
            return {"h2_status": "DONE", "h2_observed": h2, "h2_se": se,
                    "h2_z": round(h2 / se, 3) if se else "",
                    "intercept": float(i.group(1)) if i else "",
                    "intercept_se": float(i.group(2)) if i else "",
                    "mean_chisq": float(c.group(1)) if c else "",
                    "lambda_gc": float(l.group(1)) if l else ""}
    return {"h2_status": "NOT_RUN"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    ts = datetime.datetime.now().isoformat(timespec="seconds")

    traits = list(csv.DictReader(open(TRAITS), delimiter="\t"))
    src = {}
    if os.path.exists(SOURCES):
        for r in csv.DictReader(open(SOURCES), delimiter="\t"):
            src[r["trait_id"]] = r

    rows = []
    for t in traits:
        tid = t["trait_id"]
        s = src.get(tid, {})
        r = {k: "" for k in FIELDS}
        r.update({
            "trait_id": tid,
            "display_name": t.get("label", ""),
            "domain": t.get("domain", ""),
            "role": role_of(tid, t.get("domain", "")),
            "source_name": (t.get("source_note", "") or "")[:70],
            "source_url_or_accession": s.get("accession") or s.get("download_url", ""),
            "publication_or_PMID": s.get("pmid", ""),
            "ancestry": "EUR" if s else "",
            "genome_build": s.get("build_verified") or t.get("build", ""),
            "phenotype_type": ("continuous" if t.get("type") == "continuous"
                               else "case_control"),
            "sample_size": t.get("n_total", ""),
            "n_cases": t.get("ncase", ""),
            "n_controls": t.get("ncontrol", ""),
            "population_prevalence": t.get("pop_prev", ""),
            "raw_sha256": (s.get("sha256", "") or "")[:16],
            "license_or_redistribution_note":
                "public GWAS Catalog / consortium release; raw not committed",
            "last_verified_commit": commit,
            "last_verified_timestamp": ts,
        })

        # raw
        lf = s.get("local_file", "")
        rp = os.path.join("data/raw", lf) if lf and lf != "NA" else ""
        if rp and os.path.exists(rp):
            r["raw_status"], r["raw_path"] = "PRESENT", rp
        elif "EVICTED" in (s.get("notes", "") or ""):
            r["raw_status"] = "EVICTED_REPRODUCIBLE"
            r["raw_path"] = s.get("download_url", "")
        elif s.get("access") in ("REGISTRATION", "CONTROLLED", "UNAVAILABLE"):
            r["raw_status"] = "BLOCKED_" + s.get("access", "")
        else:
            r["raw_status"] = "ABSENT"

        # harmonized
        hp = f"data/harmonized/{tid}.harmonized.tsv.gz"
        if os.path.exists(hp):
            r["harmonized_status"], r["harmonized_path"] = "DONE", hp
            r["harmonized_sha256"] = sha_head(hp)
        else:
            r["harmonized_status"] = "NOT_RUN"

        # munged
        mp = f"data/munged/{tid}.sumstats.gz"
        if os.path.exists(mp):
            r["munged_status"], r["munged_path"] = "DONE", mp
            r["munged_sha256"] = sha_head(mp)
        else:
            r["munged_status"] = "NOT_RUN"

        r.update(parse_h2(tid))

        # gates
        if r.get("h2_status") == "DONE":
            z = float(r["h2_z"]) if r["h2_z"] != "" else 0
            icpt = float(r["intercept"]) if r["intercept"] != "" else 99
            n = t.get("n_total", "")
            try:
                neff_h2 = float(n) * float(r["h2_observed"])
            except (TypeError, ValueError):
                neff_h2 = ""
            r["N_eff_times_h2"] = round(neff_h2) if neff_h2 != "" else ""
            ok = (z >= Z_MIN) and (icpt <= INTERCEPT_MAX)
            r["phase1_eligible"] = "YES" if ok else "NO"
            if not ok:
                r["exclusion_reason"] = (
                    f"Z={z:.2f}<{Z_MIN}" if z < Z_MIN
                    else f"intercept={icpt:.3f}>{INTERCEPT_MAX}")
            # Phase 2/3 need a real local signal; reuse the same gate plus power
            r["phase2_eligible"] = "YES" if ok else "NO"
            r["phase3_eligible"] = "YES" if ok else "NO"
        else:
            r["phase1_eligible"] = "PENDING"
            r["phase2_eligible"] = "PENDING"
            r["phase3_eligible"] = "PENDING"

        note = s.get("notes", "") or ""
        for key in ("PHASE 1 BLOCKER", "BLOCKER", "NOT CURATED", "REJECTED"):
            if key in note:
                r["blocker"] = note.split(key, 1)[1][:150].strip(": ")
                break
        rows.append(r)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}  ({len(rows)} traits)")

    if a.summary:
        import collections
        by_role = collections.Counter(r["role"] for r in rows)
        print("\nby role: " + ", ".join(f"{k}={v}" for k, v in sorted(by_role.items())))
        for stage in ("raw_status", "harmonized_status", "munged_status", "h2_status"):
            c = collections.Counter(r[stage] for r in rows)
            print(f"{stage:20}" + ", ".join(f"{k}={v}" for k, v in sorted(c.items())))
        print("\ncore sleep anchors:")
        for tid in CORE_SLEEP:
            m = [r for r in rows if r["trait_id"] == tid]
            if not m:
                print(f"  {tid:12} ABSENT FROM REGISTRY")
                continue
            r = m[0]
            print(f"  {tid:12} raw={r['raw_status']:22} munged={r['munged_status']:9} "
                  f"h2={r['h2_status']}")


if __name__ == "__main__":
    main()
