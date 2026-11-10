"""Freeze source availability and component roles without inspecting test performance.

Only original finite-coefficient availability is an outcome eligibility rule.
No significance, activity, measured effect size, or precision threshold is used.
Reserved numeric coefficients are never written into development artifacts.
"""
import sys, math, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *

SOURCE=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/Siraj_unused_libraries'
DEST=ROOT/'outputs/sleep_frailty_ml_next'
CELL_MAP={'HepG2':('hepatocyte_tumor','HepG2'),'K562':('hematopoietic','K562'),'SKNSH':('neuronal_tumor','SK-N-SH')}

def number(value):
    try:return float(value)
    except (ValueError,TypeError):return float('nan')

def partition(component):
    # Whole-component, outcome-independent assignment, shared across cell lines.
    value=int(hashlib.sha256(('OL27-coverage-20261108|'+str(component)).encode()).hexdigest()[:16],16)%1000000
    return 0 if value<600000 else 1 if value<800000 else 2

def write_gzip_csv(path,rows):
    with gzip.open(path,'wt',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    started=time.monotonic();assert_previous_sources()
    target=DEST/'library_measurement_role_lock.json';assert not target.exists()
    guard=json.loads((SOURCE/'final_pretraining_guard_audit.json').read_text());assert guard['passed']
    retrieval=json.loads((SOURCE/'targeted_retrieval_receipt.json').read_text())
    roles={(r['library'],r['original_variant_GRCh37']):r for r in csv_rows(SOURCE/'pretraining_guarded_source_roles.csv')}
    designs={(r['library'],r['original_variant_GRCh37']):r for r in csv_rows(SOURCE/'full_reference_and_shared_core_roles.csv')}
    assert roles.keys()==designs.keys();available={};stats={};source_hashes={}
    manifest=[];development=[];all_components=collections.defaultdict(set)
    for r in roles.values():all_components[r['component_root']].add(r['role'])
    for component,assigned in all_components.items():
        assert not ('development_candidate' in assigned and any(r.startswith('reserved') for r in assigned))
    for item in retrieval['files']:
        library,cell=item['library'],item['cell'];path=ROOT/item['path'];assert sha(path)==item['sha256'];source_hashes[str(path.relative_to(ROOT))]=item['sha256']
        centered={};total=0
        with gzip.open(path,'rt') as handle:raw_rows=list(csv.DictReader(handle,delimiter='\t'))
        for r in raw_rows:
            total+=1
            if r['window']!='center' or r['strand']!='fwd' or r['haplotype']!='ref':continue
            identity='chr'+r['SNP'].removeprefix('chr');assert identity not in centered
            assert r['comb']==r['SNP']+'_center_fwd_ref'
            parts=r['SNP'].removeprefix('chr').split(':');assert len(parts)==4
            assert parts==[r['chr'].removeprefix('chr'),r['pos'],r['ref_allele'],r['alt_allele']]
            centered[identity]=r
        counts=collections.Counter(total_original_rows=total,unique_center_forward_rows=len(centered));availability=collections.Counter()
        selected=sorted(k for k in roles if k[0]==library)
        for key in selected:
            role=roles[key];design=designs[key];assert design['full_physical_reference_verified']=='True'
            raw=centered.get(key[1]);effect=number(raw.get('Log2Skew')) if raw else float('nan');se=number(raw.get('Skew_SE')) if raw else float('nan')
            finite=math.isfinite(effect);se_ok=math.isfinite(se) and se>0
            count_ok=bool(raw) and all(math.isfinite(number(raw.get(k))) and number(raw.get(k))>=0 for k in ['A_Ctrl_Mean','A_Exp_Mean','B_Ctrl_Mean','B_Exp_Mean'])
            p=partition(role['component_root']) if role['role']=='development_candidate' else -1
            availability[(role['role'],'designs')]+=1
            availability[(role['role'],'row_present')]+=int(raw is not None)
            availability[(role['role'],'finite_coefficient')]+=int(finite)
            availability[(role['role'],'positive_finite_SE')]+=int(se_ok)
            availability[(role['role'],'finite_nonnegative_source_count_means')]+=int(count_ok)
            entry={'library':library,'cell_file':cell,'cell':CELL_MAP[cell][0],'line':CELL_MAP[cell][1],
                   'original_variant_GRCh37':key[1],'variant_GRCh38':role['variant_GRCh38'],
                   'component_root':role['component_root'],'role':role['role'],'development_split':p,
                   'original_center_forward_row_present':raw is not None,'finite_original_magnitude_available':finite,
                   'positive_finite_original_SE':se_ok,'source_count_means_available_nonnegative':count_ok,
                   'magnitude_eligible':finite,'original_source_row_ID':raw['ID'] if raw else ''}
            manifest.append(entry)
            if role['role']=='development_candidate' and finite:
                development.append({**entry,'effect':abs(effect),'measurement_standard_error':se if se_ok else float('nan'),
                                    'signed_effect':float('nan'),'activity_reference':float('nan'),'activity_alternate':float('nan'),
                                    'reference':design['reference'],'alternate':design['alternate'],'changed_position':design['changed_position'],
                                    'chromosome':design['chromosome'],'position_hg38':design['position_hg38']})
        stats[library+'|'+cell]={**counts,'by_role':{role:{metric:value for (r,metric),value in availability.items() if r==role} for role in sorted(set(r for r,m in availability))}}
        print('MEASUREMENT AVAILABILITY CHECKED',library,cell,len(centered),flush=True)
    assert len(manifest)==len(roles)*3
    dev_components={r['component_root'] for r in development};reserved_components={r['component_root'] for r in manifest if r['role'].startswith('reserved')}
    assert not dev_components&reserved_components
    assert all(r['library']=='OL27' and r['role']=='development_candidate' for r in development)
    for r in development:
        assert len(r['reference'])==len(r['alternate'])==200
        assert [i for i,(a,b) in enumerate(zip(r['reference'],r['alternate'])) if a!=b]==[int(r['changed_position'])]
    write_gzip_csv(DEST/'library_measurement_availability_manifest.csv.gz',manifest)
    write_gzip_csv(DEST/'OL27_development_measurements.csv.gz',development)
    by_split=collections.Counter(str(r['development_split']) for r in development)
    split_components={i:{r['component_root'] for r in development if r['development_split']==i} for i in range(3)}
    assert not (split_components[0]&split_components[1] or split_components[0]&split_components[2] or split_components[1]&split_components[2])
    record={'passed':True,'source_files_sha256':source_hashes,'guard_sha256':sha(SOURCE/'final_pretraining_guard_audit.json'),
            'role_file_sha256':sha(SOURCE/'pretraining_guarded_source_roles.csv'),'physical_design_file_sha256':sha(SOURCE/'full_reference_and_shared_core_roles.csv'),
            'availability_manifest_sha256':sha(DEST/'library_measurement_availability_manifest.csv.gz'),
            'development_measurements_sha256':sha(DEST/'OL27_development_measurements.csv.gz'),
            'denominators':stats,'development_measurements':len(development),'development_physical_designs':len({r['original_variant_GRCh37'] for r in development}),
            'development_rows_by_split':dict(by_split),'development_components_by_split':{str(i):len(v) for i,v in split_components.items()},
            'no_reserved_components_in_development':True,'selection_policy':'Finite original Log2Skew only; abs(Log2Skew) is primary magnitude. No p-value, detected/activity, effect-size or SE threshold. Missing SE stays missing; signed and allele-activity targets stay unknown.',
            'partition_policy':'SHA256 OL27-coverage-20261108|component modulo 1,000,000; <600,000 training, <800,000 validation, otherwise calibration. All cells and connected physical designs stay together.',
            'assay_scope':'Original published paired GLM estimates; full source DNA verified. Count means are source aggregates, not independent raw barcode/repeat QC. Original GLM was not independently refit. Sign convention remains unproven.',
            'reserved_numeric_targets_exported_to_training':False,'reserved_model_predictions_or_performance_calculated':False,
            'not_independent_publication_or_donor_replication':True,'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-started}
    write_json(target,record);assert_previous_sources();print(json.dumps({k:record[k] for k in ['passed','development_measurements','development_physical_designs','development_rows_by_split','development_components_by_split']}),flush=True)

if __name__=='__main__':main()
