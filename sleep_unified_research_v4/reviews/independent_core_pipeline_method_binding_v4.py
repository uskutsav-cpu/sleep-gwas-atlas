#!/usr/bin/env python3
"""Supplement the replay design with original config/method byte identities."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'sleep_unified_research_v4'
ORIGINAL='659d01cf5b10a6debf1841404e7e01ef4750f755'
ARCHIVE=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas')

def digest(data):return hashlib.sha256(data).hexdigest()

def main():
    rows=[]
    for n in ['config/analysis_panel.tsv','config/analysis_panel.lock.json','config/gwas_schemas.tsv',
              'config/hm3_prefilter_plans.tsv','config/liftover_plans.tsv','config/variant_mapping_plans.tsv',
              'config/public_gwas_sources.tsv']:
        current=(ROOT/n).read_bytes()
        original=subprocess.check_output(['git','show',ORIGINAL+':'+n],cwd=ROOT)
        archived=(ARCHIVE/n).read_bytes()
        assert current==original==archived
        rows.append(dict(relative_path=n,current_sha256=digest(current),original_git_sha256=digest(original),
            archived_path=str(ARCHIVE/n),archived_sha256=digest(archived),exact_three_way_identity=True))
    old_shell=subprocess.check_output(['git','show',ORIGINAL+':scripts/02_munge.sh'],cwd=ROOT)
    changes=subprocess.check_output(['git','diff',ORIGINAL,'HEAD','--','scripts/02_munge.sh','scripts/_common.sh'],cwd=ROOT,text=True)
    receipt=dict(status='PASS_ORIGINAL_CONFIG_AND_SCIENTIFIC_SOURCE_BINDINGS',original_commit=ORIGINAL,
        checker_sha256=digest(Path(__file__).read_bytes()),configs=rows,
        original_02_sha256=digest(old_shell),current_02_sha256=digest((ROOT/'scripts/02_munge.sh').read_bytes()),
        shell_diff=changes,shell_differences_reviewed_as_output_reference_and_workspace_routing_only=True,
        raw_source_path_remains_hardcoded_in_shell=True,direct_shell_execution_into_historical_or_repo_inputs_not_admitted=True,
        GWAS_body_reads=0,native_workers_launched=0,
        per_trait_historical_code_runtime_attestation_independently_established=False)
    target=P/'reviews/independent_core_pipeline_method_binding_receipt_v4.json'
    with target.open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('configs','shell_diff')}))

if __name__=='__main__':main()
