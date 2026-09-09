"""Generate auditable figures and status reports without fabricated results."""
from __future__ import annotations
import html
import json
from pathlib import Path
from .artifacts import verify_artifact
from .io import read_json, read_tsv, write_json, atomic_text


def report(root, out):
    root=Path(root)
    records=[]
    for p in sorted(root.rglob("receipt.json")):
        if any(x.startswith(".") for x in p.relative_to(root).parts):continue
        try:
            raw = read_json(p)
            if raw.get("status") == "FAILED":
                records.append({"path":str(p.parent.relative_to(root)),"stage":raw.get("stage","UNKNOWN"),
                    "status":"FAILED_DIAGNOSTIC_ONLY", "synthetic":raw.get("synthetic","UNKNOWN"),
                    "fingerprint":raw.get("fingerprint","")})
                continue
            r=verify_artifact(p.parent)
            records.append({"path":str(p.parent.relative_to(root)),"stage":r["stage"],
                "status":r.get("scientific_status","COMPUTED"),"synthetic":r["synthetic"],
                "fingerprint":r["fingerprint"]})
        except Exception as e:
            records.append({"path":str(p.parent.relative_to(root)),"stage":"UNKNOWN",
                            "status":"CORRUPTED_OR_INVALID","synthetic":"UNKNOWN","fingerprint":str(e)})
    with atomic_text(out) as f:
        f.write("# Brain6 execution report\n\n")
        f.write("Computed artifacts are not synonymous with accepted science. Synthetic tests are not empirical findings.\n\n")
        f.write("| Artifact | Stage | Status | Synthetic |\n|---|---|---|---|\n")
        for r in records:
            f.write(f"| {r['path']} | {r['stage']} | {r['status']} | {r['synthetic']} |\n")
        if not records:f.write("\nNo executed artifact receipts found. Production is not complete.\n")
        f.write("\n## Interpretation boundary\n\nGenetic correlation, pleiotropic association, shared statistical signals, gene links, and causal mechanisms are different evidence levels.\n")
    return records


def forest(ranking_file, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    rows=list(read_tsv(ranking_file,["pair_id","rg","se","eligible"]))
    rows=[r for r in rows if r["eligible"] in {"True","true","1"}]
    if not rows:raise ValueError("No eligible measured associations to plot")
    rows=rows[:24]
    x=np.array([float(r["rg"]) for r in rows]); se=np.array([float(r["se"]) for r in rows])
    fig,ax=plt.subplots(figsize=(10,max(4,len(rows)*.32)))
    ax.errorbar(x,np.arange(len(rows)),xerr=1.96*se,fmt="o",capsize=3)
    ax.set_yticks(range(len(rows)),[r["pair_id"].replace("__"," × ") for r in rows])
    ax.invert_yaxis();ax.axvline(0,linestyle="--")
    ax.set_xlabel("Genetic correlation (nominal 95% confidence interval)")
    ax.set_title(("SYNTHETIC — " if "SYNTHETIC" in str(ranking_file).upper() else "")+"Brain6: atlas associations — not independent replication")
    fig.tight_layout(); fig.savefig(out,dpi=180);plt.close(fig)
