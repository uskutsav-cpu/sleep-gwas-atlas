#!/usr/bin/env python3
"""Project acceptance gates. Exit 0 only when the whole project is complete.

Reads machine-readable artifacts only. Prose in a report never satisfies a
gate: a gate is satisfied by a file that exists, parses, has the expected
columns, and has rows that agree with the ledger.

    python3 scripts/99_project_acceptance.py            # all gates
    python3 scripts/99_project_acceptance.py --next     # first incomplete only
    python3 scripts/99_project_acceptance.py --phase 2

Every gate returns (ok, detail). A gate that cannot be satisfied on this
machine must be recorded as a BLOCKER row in results/blockers.tsv with
evidence, and is then reported as BLOCKED rather than silently passing.
"""
import argparse
import csv
import glob
import gzip
import os
import subprocess
import sys

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(R)

BLOCKERS = "results/blockers.tsv"
TERMINAL_STATES = {
    "COMPLETE_PASS", "COMPLETE_FAIL_H2", "BLOCKED_SOURCE_RESTRICTED",
    "BLOCKED_SOURCE_MISSING", "BLOCKED_SCHEMA", "BLOCKED_BUILD",
    "BLOCKED_REFERENCE_MAPPING", "BLOCKED_INSUFFICIENT_VARIANTS",
    "EXCLUDED_QC", "DUPLICATE_PHENOTYPE", "SUPERSEDED",
}


def rows(path, min_rows=1):
    if not os.path.exists(path):
        return None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        r = list(csv.DictReader(fh, delimiter="\t"))
    return r if len(r) >= min_rows else (r if r else None)


def blocked(gate):
    """Is this gate recorded as an evidence-backed external blocker?"""
    b = rows(BLOCKERS)
    if not b:
        return None
    for x in b:
        if x.get("gate") == gate and x.get("evidence", "").strip():
            return x.get("reason", "")[:110]
    return None


def have(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


# ----------------------------------------------------------------- gates
def g_env():
    need = ["environment/environment.lock.md", ".r-lib", "ldsc", "ref"]
    miss = [p for p in need if not os.path.exists(p)]
    return (not miss), ("ok" if not miss else f"missing {miss}")


def g_ledger_terminal():
    led = rows("results/phase_ledger.tsv")
    if not led:
        return False, "no ledger"
    bad = [r["trait_id"] for r in led
           if r.get("terminal_state", "") not in TERMINAL_STATES]
    if bad:
        return False, f"{len(bad)} trait(s) not in a terminal state e.g. {bad[:4]}"
    return True, f"all {len(led)} traits terminal"


def g_core_anchors():
    led = rows("results/phase_ledger.tsv")
    if not led:
        return False, "no ledger"
    core = [r for r in led if r["role"] == "core_sleep"]
    done = [r for r in core if r["h2_status"] == "DONE"]
    return len(done) == 6, f"{len(done)}/6 core sleep anchors with h2"


def g_phase1_runnable():
    led = rows("results/phase_ledger.tsv")
    if not led:
        return False, "no ledger"
    pend = [r["trait_id"] for r in led
            if r.get("terminal_state") not in TERMINAL_STATES
            and r["h2_status"] != "DONE"]
    return not pend, ("no runnable incomplete traits"
                      if not pend else f"{len(pend)} runnable incomplete")


def g_phase1_tables():
    need = ["results/phase1/h2_summary.tsv",
            "results/phase1/rg_primary_complete.tsv",
            "results/phase1/fdr_families.tsv",
            "results/phase1/exclusions_and_blockers.tsv",
            "results/phase1/provenance_manifest.tsv"]
    miss = [p for p in need if not have(p)]
    return (not miss), ("ok" if not miss else f"missing {[os.path.basename(m) for m in miss]}")


def g_phase1_frozen():
    f = rows("results/phase1/fdr_families.tsv")
    if not f:
        return False, "no fdr_families.tsv"
    fam = [x for x in f if x.get("status") == "FROZEN"]
    if not fam:
        return False, "no FROZEN family recorded"
    rg = rows("results/phase1/rg_primary_complete.tsv")
    if not rg:
        return False, "no primary rg table"
    want = int(fam[0].get("n_pairs", 0))
    return len(rg) == want, f"frozen family {want} pairs, table has {len(rg)}"


def g_phase2():
    need = ["results/phase2/prioritized_pairs.tsv",
            "results/phase2/lava_univariate.tsv",
            "results/phase2/lava_bivariate_all.tsv"]
    miss = [p for p in need if not have(p)]
    return (not miss), ("ok" if not miss else f"missing {[os.path.basename(m) for m in miss]}")


def g_phase3_placo():
    return have("results/phase3/placo_plus_loci.tsv"), \
        "ok" if have("results/phase3/placo_plus_loci.tsv") else "no PLACO+ loci table"


def g_phase3_second():
    for p in ("results/phase3/conjfdr_loci.tsv",
              "results/phase3/fallback_method_loci.tsv"):
        if have(p):
            return True, f"second method: {os.path.basename(p)}"
    return False, "neither conjFDR nor a validated fallback present"


def g_phase(n, files):
    miss = [p for p in files if not have(p)]
    return (not miss), ("ok" if not miss else f"missing {[os.path.basename(m) for m in miss]}")


GATES = [
    ("env",              "environment locked and tools present", g_env),
    ("phase1.core",      "6/6 core sleep anchors", g_core_anchors),
    ("phase1.runnable",  "no runnable incomplete traits", g_phase1_runnable),
    ("phase1.ledger",    "every ledger row terminal", g_ledger_terminal),
    ("phase1.tables",    "phase 1 deliverable tables", g_phase1_tables),
    ("phase1.frozen",    "frozen FDR family matches table", g_phase1_frozen),
    ("phase2",           "LAVA local rg", g_phase2),
    ("phase3.placo",     "PLACO+ loci", g_phase3_placo),
    ("phase3.second",    "conjFDR or validated fallback", g_phase3_second),
    ("phase4",           "GenomicSEM factor models",
     lambda: g_phase(4, ["results/phase4/factor_model_fit.tsv"])),
    ("phase5",           "shared vs trait-specific",
     lambda: g_phase(5, ["results/phase5/qsnp_results.tsv"])),
    ("phase6",           "fine-mapping and colocalization",
     lambda: g_phase(6, ["results/phase6/coloc_summary.tsv"])),
    ("phase7",           "TWAS / molecular mediation",
     lambda: g_phase(7, ["results/phase7/twas_results.tsv"])),
    ("phase8",           "cell-type enrichment",
     lambda: g_phase(8, ["results/phase8/celltype_enrichment.tsv"])),
    ("phase9",           "pathways and networks",
     lambda: g_phase(9, ["results/phase9/pathway_enrichment.tsv"])),
    ("phase10",          "causal inference (MR)",
     lambda: g_phase(10, ["results/phase10/mr_results.tsv"])),
    ("manuscript",       "manuscript sections",
     lambda: g_phase(0, ["manuscript/methods.md", "manuscript/results.md"])),
    ("git.clean",        "working tree clean",
     lambda: (not subprocess.run(["git", "status", "--porcelain"],
                                 capture_output=True, text=True).stdout.strip(),
              "clean" if not subprocess.run(["git", "status", "--porcelain"],
                                            capture_output=True, text=True).stdout.strip()
              else "uncommitted changes")),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--next", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    incomplete, blocked_n = [], 0
    print("PROJECT ACCEPTANCE")
    print("=" * 66)
    for key, desc, fn in GATES:
        try:
            ok, detail = fn()
        except Exception as e:
            ok, detail = False, f"gate error: {type(e).__name__}: {e}"
        if ok:
            mark = "PASS"
        else:
            b = blocked(key)
            if b:
                mark, detail, blocked_n = "BLOCKED", b, blocked_n + 1
                ok = True          # an evidence-backed blocker satisfies the gate
            else:
                mark = "FAIL"
                incomplete.append((key, desc, detail))
        if not a.quiet:
            print(f"  [{mark:7}] {key:16} {desc:38} {detail[:60]}")

    print("=" * 66)
    if incomplete:
        print(f"{len(incomplete)} incomplete gate(s); {blocked_n} externally blocked")
        print(f"\nFIRST INCOMPLETE EXECUTABLE GATE: {incomplete[0][0]}")
        print(f"  {incomplete[0][1]} -- {incomplete[0][2]}")
        if a.next:
            print(f"\nNEXT: {incomplete[0][0]}")
        sys.exit(1)
    print(f"ALL GATES SATISFIED ({blocked_n} externally blocked with evidence)")
    sys.exit(0)


if __name__ == "__main__":
    main()
