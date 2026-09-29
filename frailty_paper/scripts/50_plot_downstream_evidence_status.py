#!/usr/bin/env python3
"""Render transparent eligibility/status panels for unsupported Figures 4 and 5."""
from __future__ import annotations
import argparse, csv, hashlib, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "frailty_paper/analysis"
BASE_INPUTS = [
    "frailty_paper/config/lava_frailty_sensitivity_v1.yaml",
    "frailty_paper/analysis/local_sharing_readiness_audit_2026-09-24.md",
    "frailty_paper/manifests/targeted_resource_audit.tsv",
    "frailty_paper/analysis/pqtl_access_recheck_2026-09-24.md",
    "frailty_paper/analysis/hlma_processed_data_portal_audit_2026-09-23.md",
    "frailty_paper/manifests/muscle_single_cell_resource_index.tsv",
]

def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()
def wrap(text, width=42):
    words=text.split(); lines=[]; line=""
    for w in words:
        if line and len(line)+len(w)+1>width: lines.append(line); line=w
        else: line=(line+" "+w).strip()
    if line: lines.append(line)
    return lines

def draw_figure(number, title, subtitle, cards, footer, stem):
    fig, ax=plt.subplots(figsize=(14,8)); fig.patch.set_facecolor("#F5F7FA"); ax.set_facecolor("#F5F7FA")
    ax.set_xlim(0,14); ax.set_ylim(0,8); ax.axis("off")
    ax.text(.48,7.48,f"Figure {number}",fontsize=12,weight="bold",color="#145C72")
    ax.text(.48,6.96,title,fontsize=21,weight="bold",color="#182230")
    ax.text(.48,6.56,subtitle,fontsize=10.5,color="#475467")
    gap=.24; n=len(cards); cols=2 if n>2 else n; rows=(n+cols-1)//cols
    x0=.48; ytop=6.05; total_w=13.04; cw=(total_w-gap*(cols-1))/cols; ch=2.06 if rows>1 else 2.28; vgap=.18
    for k,(heading,value,detail,evidence) in enumerate(cards):
        r,c=divmod(k,cols); x=x0+c*(cw+gap); y=ytop-(r+1)*ch-r*vgap
        ax.add_patch(FancyBboxPatch((x,y),cw,ch,boxstyle="round,pad=0.02,rounding_size=0.12",facecolor="white",edgecolor="#D9E0E8",linewidth=1.2))
        ax.text(x+.22,y+ch-.34,heading.upper(),fontsize=9,weight="bold",color="#667085",va="top")
        ax.text(x+.22,y+ch-.72,value,fontsize=15,weight="bold",color="#9A5A12" if any(t in value.lower() for t in ["failed","unavailable","none","not justified"]) else "#145C72",va="top")
        lines=wrap(detail, 66 if cw>6 else 38)
        ax.text(x+.22,y+ch-1.12,"\n".join(lines),fontsize=9.2,color="#344054",va="top",linespacing=1.35)
        evidence_lines=wrap("Evidence: "+evidence, 78 if cw>6 else 38)
        ax.text(x+.22,y+.13,"\n".join(evidence_lines),fontsize=7.0,color="#667085",va="bottom",linespacing=1.12)
    ax.add_patch(FancyBboxPatch((.48,.37),13.04,.72,boxstyle="round,pad=0.02,rounding_size=0.1",facecolor="#E9EEF3",edgecolor="#CBD5DF",linewidth=1))
    ax.text(.72,.73,footer,fontsize=9.2,color="#344054",va="center",wrap=True)
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    paths=[OUT/f"{stem}.png",OUT/f"{stem}.pdf",OUT/f"{stem}.svg"]
    fig.savefig(paths[0],dpi=300,bbox_inches="tight",facecolor=fig.get_facecolor())
    fig.savefig(paths[1],bbox_inches="tight",facecolor=fig.get_facecolor(),metadata={"Title":title,"Subject":"Evidence status only; no local or molecular inference"})
    fig.savefig(paths[2],bbox_inches="tight",facecolor=fig.get_facecolor())
    svg_text=paths[2].read_text(encoding="utf-8")
    paths[2].write_text("\n".join(line.rstrip() for line in svg_text.splitlines())+"\n",encoding="utf-8")
    plt.close(fig)
    return paths

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lava-diagnostic",type=Path,default=Path("frailty_paper/analysis/lava_interim_receipt_diagnostic_2026-09-24_1148.json"),help="Read-only interim audit JSON to show in Figure 4")
    parser.add_argument("--lava-full-audit",type=Path,help="Read-only full-family receipt audit JSON; when supplied, this replaces the older interim prefix in Figure 4")
    args=parser.parse_args()
    full_audit=None
    if args.lava_full_audit:
        audit_path=(ROOT/args.lava_full_audit if not args.lava_full_audit.is_absolute() else args.lava_full_audit).resolve()
        try: audit_rel=audit_path.relative_to(ROOT).as_posix()
        except ValueError: raise SystemExit("LAVA full audit must be inside the repository so its provenance path is portable")
        full_audit=json.loads(audit_path.read_text())
        lock_rel="frailty_paper/config/lava_frailty_sensitivity_v1.yaml"
        if full_audit.get("lock_sha256") != sha(ROOT/lock_rel) or full_audit.get("runner_state_lock_sha256") != full_audit.get("lock_sha256"):
            raise SystemExit("Full audit does not match the frozen LAVA analysis lock")
        if full_audit.get("mutates_run_or_receipts") is not False:
            raise SystemExit("Full audit is not certified read-only")
        if full_audit.get("expected_receipts") != 29940 or len(full_audit.get("process_failures_by_trait", {})) != 12:
            raise SystemExit("Full audit does not cover the locked 29,940-slot, 12-trait family")
        inputs=[lock_rel,audit_rel,*BASE_INPUTS]
    else:
        diagnostic=(ROOT/args.lava_diagnostic if not args.lava_diagnostic.is_absolute() else args.lava_diagnostic).resolve()
        try: diagnostic_rel=diagnostic.relative_to(ROOT).as_posix()
        except ValueError: raise SystemExit("LAVA diagnostic must be inside the repository so its provenance path is portable")
        inputs=["frailty_paper/config/lava_frailty_sensitivity_v1.yaml",diagnostic_rel,*BASE_INPUTS]
    for rel in inputs:
        if not (ROOT/rel).is_file(): raise SystemExit(f"Missing provenance input: {rel}")
    if full_audit:
        total=int(full_audit["expected_receipts"]); audited=int(full_audit["receipt_files_validated"])
        failed_traits=sum(not row["passes_1pct_gate"] for row in full_audit["process_failures_by_trait"].values())
        stamp=full_audit["audit_time_end_utc"].replace("T"," ")[:19]+" UTC"
        audit_ref=audit_rel.removeprefix("frailty_paper/")
        inventory_note=("receipt inventory advanced during scanning." if full_audit.get("receipt_file_list_changed_during_scan") else "receipt inventory and runner state were stable during scanning.")
        fig4=[
            ("Locked analysis", "RUNNING · 12 × 2,495", "The prespecified Frailty Index × sleep/circadian sensitivity family is still executing.",f"config/lava_frailty_sensitivity_v1.yaml; {audit_ref}"),
            ("Latest full-family audit", f"{audited:,}/{total:,} receipts", f"Read-only audit completed {stamp}; {inventory_note}",audit_ref),
            ("Locked 1% trait gates", f"FAILED · {failed_traits}/12", "All 12 trait-level failure lower bounds already exceed the frozen 1% allowance; no thresholds or exclusions were changed.",audit_ref),
            ("Interpretation", "NO LOCAL INFERENCE", "The family is incomplete and local estimates remain inadmissible. Failed gates do not establish absence of local sharing.",audit_ref),
        ]
        trait="all 12 traits"; flagged=failed_traits; denom=12; max_fail=0
        audit_generated=full_audit["audit_time_end_utc"]
    else:
        diag=json.loads(diagnostic.read_text())
        trait=str(diag["trait"]); audited=int(diag["receipts_audited"]); flagged=int(diag["loci_flagged_by_frozen_collator_rule"]); denom=int(diag["locus_family_denominator"]); max_fail=int(diag["locked_maximum_failed_loci"])
        if (denom,max_fail)!=(2495,24): raise SystemExit("Frozen LAVA denominator or failure limit changed; review before drawing")
        if diag["analysis_lock_sha256"] != sha(ROOT/inputs[0]): raise SystemExit("Diagnostic does not match the frozen LAVA analysis lock")
        if not diag["gate_already_exceeded"]: raise SystemExit("Interim gate is not failed; review panel claims before drawing")
        stamp=diag["generated_utc"].replace("T"," ")[:19]+" UTC"
        common_refs=f"config/lava_frailty_sensitivity_v1.yaml; {diagnostic_rel.removeprefix('frailty_paper/')}"
        fig4=[
            ("Locked analysis", "RUNNING · 12 × 2,495", "The prespecified Frailty Index × sleep/circadian local-sharing sensitivity family is still executing.",common_refs),
            (f"Latest audited {trait} prefix", f"{audited:,} receipts", f"As of {stamp}; {flagged:,}/{denom:,} loci ({flagged/denom:.2%}) are flagged under the frozen collator rule.",diagnostic_rel.removeprefix("frailty_paper/")),
            ("Completeness gate", f"FAILED · limit {max_fail}/{denom:,}", "The interim flagged count exceeds the unchanged 1% family limit. All process failures remain counted; no thresholds or exclusions were relaxed.",diagnostic_rel.removeprefix("frailty_paper/")),
            ("Interpretation", "NO LOCAL INFERENCE", "No family-level local result or example is interpretable at this checkpoint. Final 29,940-slot collation remains pending.","analysis/local_sharing_readiness_audit_2026-09-24.md; config/lava_frailty_sensitivity_v1.yaml"),
        ]
    local_ref=audit_ref if full_audit else diagnostic_rel.removeprefix("frailty_paper/")
    fig5=[
        ("Frailty-specific loci", "NONE ELIGIBLE", "No eligible shared-locus set is available to drive fine-mapping, trait–trait colocalization or downstream molecular queries.",f"analysis/local_sharing_readiness_audit_2026-09-24.md; {local_ref}"),
        ("Molecular-QTL resources", "INDEXED · NOT LOCUS-MATCHED", "eQTL Catalogue and BrainSCOPE are resource indexes; pQTL objects remain unacquired. No locus-specific QTL integration is reported.","manifests/targeted_resource_audit.tsv; analysis/pqtl_access_recheck_2026-09-24.md"),
        ("Muscle and cell resources", "METADATA PARTIAL", "HLMA donor/library supplements and muscle-resource indexes exist, but archive identity, donor mapping, build and reuse details remain incomplete.","analysis/hlma_processed_data_portal_audit_2026-09-23.md; manifests/muscle_single_cell_resource_index.tsv"),
        ("Biological interpretation", "NO GENE / TISSUE CLAIM", "There is no eligible locus-to-gene, protein, tissue or cell-type evidence in this package. This status panel reports evidence readiness, not biological absence.","analysis/local_sharing_readiness_audit_2026-09-24.md; analysis/pqtl_access_recheck_2026-09-24.md"),
    ]
    OUT.mkdir(parents=True,exist_ok=True)
    figures={
      "figure4_local_sharing_evidence_status":draw_figure(4,"Local sharing: interim evidence status","Locked sensitivity analysis · checkpoint only · no locus-level inference",fig4,"The run continues. Failed trait gates prohibit local-sharing inference; they do not show that local sharing is absent.","figure4_local_sharing_evidence_status"),
    }
    if full_audit:
        # Preserve Figure 5's prior rendering and hash: this update changes only
        # the LAVA evidence panel, while its molecular eligibility evidence is unchanged.
        fig5_paths=[OUT/f"figure5_molecular_evidence_status.{suffix}" for suffix in ("png","pdf","svg")]
        if not all(path.is_file() for path in fig5_paths): raise SystemExit("Existing Figure 5 outputs are required when refreshing Figure 4 only")
        figures["figure5_molecular_evidence_status"]=fig5_paths
    else:
        figures["figure5_molecular_evidence_status"]=draw_figure(5,"Molecular follow-up: upstream eligibility status","No eligible frailty-specific locus set · resource availability is not locus evidence",fig5,"No genes, proteins, tissues or cell types are inferred. Reassess only after eligible upstream loci and locus-matched molecular data become available.","figure5_molecular_evidence_status")
    source_path=OUT/"figures4_5_downstream_evidence_status_source_data.tsv"
    with source_path.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(["figure","panel","status","detail","evidence_files"])
        for stem,cards in [("Figure 4",fig4),("Figure 5",fig5)]:
            for i,(h,v,d,e) in enumerate(cards,1): w.writerow([stem,str(i),v,d,e])
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),"command":"python frailty_paper/scripts/50_plot_downstream_evidence_status.py "+" ".join(sys.argv[1:]),"script":"frailty_paper/scripts/50_plot_downstream_evidence_status.py","script_sha256":sha(Path(__file__).resolve()),"python":platform.python_version(),"matplotlib":matplotlib.__version__,"inputs_sha256":{rel:sha(ROOT/rel) for rel in inputs},"latest_lava_trait":trait,"latest_lava_receipts_audited":audited,"loci_flagged":flagged,"locus_denominator":denom,"locked_maximum_failures":max_fail,"audit_generated_utc":audit_generated if full_audit else diag["generated_utc"],"audit_type":"full_family" if full_audit else "interim_trait_prefix","figure5_rendered_this_run":not bool(full_audit),"figures":{stem:{p.relative_to(ROOT).as_posix():sha(p) for p in paths} for stem,paths in figures.items()},"source_data":{"frailty_paper/analysis/figures4_5_downstream_evidence_status_source_data.tsv":sha(source_path)},"interpretation_limits":["The FI × sleep LAVA family remains incomplete; failure gates preclude local inference.","A failed completeness gate does not establish absence of local sharing.","Figure 5 reports eligibility and resource status only; it contains no biological finding.","No gene, protein, tissue, or cell-type claim is made."]}
    out=OUT/"figures4_5_downstream_evidence_status.provenance.json";out.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"audit_receipts_validated":audited,"gates_failed":flagged,"source_data":str(source_path),"provenance":str(out),"figures":{k:[p.name for p in v] for k,v in figures.items()}}))
if __name__=="__main__":main()
