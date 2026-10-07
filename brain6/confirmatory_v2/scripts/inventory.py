#!/usr/bin/env python3
"""Inventory current artifacts and retained local data without following SSD links."""
import csv, hashlib, json, os, subprocess
from pathlib import Path
from audit_evidence import ROOT, AREA, sha, write_tsv

VALID={'COMPLETE_VALIDATED','COMPLETE_EXPLORATORY','FAILED_QC','NOT_RUN','BLOCKED_DATA','BLOCKED_ACCESS','BLOCKED_SOFTWARE','NEEDS_REVIEW','NOT_APPLICABLE'}

def main():
    paths=subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines()
    files=[]
    for name in paths:
        if not (name.startswith(('brain6/','extensions/brain6/','results/track_b/','.github/workflows/','scripts/','tests/','config/','docs/','environment/')) or name in {'README.md','CLAUDE.md','Snakefile'}): continue
        p=ROOT/name
        if not p.is_file(): status='BLOCKED_DATA'; digest=''; size=0; kind='MISSING'
        else:
            size=p.stat().st_size; digest=sha(p)
            with p.open('rb') as f: beginning=f.read(128)
            if beginning.startswith(b'version https://git-lfs.github.com/spec/v1'): status='BLOCKED_DATA';kind='LFS_POINTER'
            elif name.endswith(('.py','.R','.sh','.yaml','.json')) and ('scripts/' in name or '/tests/' in name or name.startswith(('scripts/','tests/','.github/'))):status='NEEDS_REVIEW';kind='IMPLEMENTATION_OR_TEST_NOT_EMPIRICAL'
            elif 'exploratory' in name or '/annotation/' in name or '/novelty/' in name:status='COMPLETE_EXPLORATORY';kind='EXPLORATORY_OR_SOURCE_CONTEXT'
            elif '/final_package_v1/' in name:status='COMPLETE_EXPLORATORY';kind='BOUNDED_PACKAGE_SEPARATE_HASH_AUDIT'
            else:status='NEEDS_REVIEW';kind='HISTORICAL_ARTIFACT_SEMANTICS_NOT_INFERRED_FROM_PRESENCE'
        files.append({'path':name,'bytes':size,'sha256':digest,'classification':status,'evidence_kind':kind})
    write_tsv(AREA/'qc/file_inventory.tsv',files)
    component_specs=[
        ('global_map','COMPLETE_VALIDATED','72 rows; 35 inherited q<0.05; original LDSC fits not rerun','brain6/results/global/brain6_72_locked.tsv'),
        ('placo_candidates','COMPLETE_EXPLORATORY','25 candidates/20 regions; no validated shared causal locus','brain6/paper/final_package_v1/BRAIN6_FINAL_CANDIDATES.tsv'),
        ('canonical_lava','FAILED_QC','13745 TESTED;3720 NOT_RUN;ceiling873; immutable nonpromotion','brain6/paper/final_package_v1/BRAIN6_FINAL_LAVA.tsv'),
        ('longsleep_source_rescue','BLOCKED_DATA','Exact model/per-SNP N mapping unavailable; perfect other repair still leaves1291','brain6/results/confirmatory_source_rescue_20260927/TERMINAL_BLOCKER_AND_RESTART.md'),
        ('secondary_supergnova','COMPLETE_EXPLORATORY','3304 estimated/8465 slots;3 FWER blocks;0 candidate overlaps;same source','brain6/results/brain6_alternative_local_validation_v1/genome_wide_blocks.tsv'),
        ('longsleep_supergnova','NOT_APPLICABLE','Three pairs/5079 slots prospectively method-inapplicable','brain6/results/brain6_alternative_local_validation_v1/source_eligibility_lock_v1.json'),
        ('legacy_ld_susie','FAILED_QC','32 LD matrices failed;no accepted multi-signal fine-mapping','brain6/results/brain6_exploratory_finemap_coloc_v1/ld_numeric_qc.tsv'),
        ('single_signal_abf','COMPLETE_EXPLORATORY','6 trait coloc;3 H4>=0.8 default;1 lower prior;not SuSiE','brain6/results/brain6_exploratory_finemap_coloc_v2/trait_coloc_25.tsv'),
        ('molecular_coloc','COMPLETE_EXPLORATORY','Four ADHD-component QTL tests;maxH4 0.0824;no pair-level molecular support','brain6/results/brain6_exploratory_functional_v1/coloc_results.tsv'),
        ('full20_enrichment','FAILED_QC','2 regions failed matched-control pool gate;no full-family p values','brain6/results/brain6_exploratory_enrichment_v3/matching_status_20.tsv'),
        ('subset18_enrichment','COMPLETE_EXPLORATORY','54 tissues/1680 pathways;none q<0.05;subset positional sensitivity','brain6/results/brain6_exploratory_enrichment_v3/run_summary.json'),
        ('independent_two_trait_replication','NOT_RUN','None established;reused sleep and documented cohort overlap','brain6/paper/final_package_v1/BRAIN6_FINAL_REPLICATION.tsv'),
        ('cell_resolved_enrichment','BLOCKED_DATA','No reviewed eligible set and tested background','brain6/paper/final_package_v1/BRAIN6_FINAL_TISSUE_CELLTYPE.tsv'),
        ('native_finemapping_runtime','BLOCKED_SOFTWARE','Installed R lacks susieR,coloc,LAVA;SSD runtime unavailable','brain6/confirmatory_v2/NEW_ANALYSIS_PROTOCOL.md'),
        ('manuscript_trackA','NOT_RUN','Validated novel key finding not established','brain6/paper/final_package_v1/BRAIN6_FINAL_READINESS.md'),
        ('human_scientific_review','NEEDS_REVIEW','Automated checks do not replace statistical-genetics collaborator','brain6/paper/REVIEWER_2_AUDIT.md')]
    components=[{'component':c,'classification':s,'quantitative_evidence_or_blocker':e,'source':p,'source_sha256':sha(ROOT/p)} for c,s,e,p in component_specs]
    assert all(x['classification'] in VALID for x in components)
    write_tsv(AREA/'COMPONENT_INVENTORY.tsv',components)
    local=[]
    parent=Path('/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03')
    for sub in ['data','ref','brain6','work/lava-canonical-v3-production','work/preparation-v1']:
        for base,ds,fs in os.walk(parent/sub,followlinks=False):
            ds[:]=[d for d in ds if d not in {'.venv','__pycache__','conda-envs'}]
            for name in fs:
                p=Path(base)/name
                if not p.exists(): status='BROKEN_EXTERNAL_SYMLINK';size=0
                else: status='PRESENT_NOT_REVALIDATED';size=p.stat().st_size
                # No giant raw-file hashing until a concrete admitted use requires it.
                local.append({'path':str(p),'bytes':size,'availability':status,'symlink_target':os.readlink(p) if p.is_symlink() else ''})
    write_tsv(AREA/'qc/retained_local_inventory.tsv',local)
    history=subprocess.check_output(['git','log','--all','--date=iso-strict','--format=%H\t%ad\t%s','--','brain6','extensions/brain6','results/track_b','.github/workflows'],cwd=ROOT,text=True)
    (AREA/'qc/relevant_commit_history.tsv').write_text('commit\tdate\tsubject\n'+history)
    for args,name in [(['branch','-avv'],'branches.txt'),(['worktree','list','--porcelain'],'worktrees.txt')]:
        (AREA/'qc'/name).write_text(subprocess.check_output(['git',*args],cwd=ROOT,text=True))
    print(json.dumps({'tracked_audited_files':len(files),'components':len(components),'retained_local_files':len(local),'external_ssd_mounted':Path('/Volumes/Extreme SSD').exists()}))

if __name__=='__main__':main()
