#!/usr/bin/env python3
"""Read-only small-artifact audit. No GWAS body reads, plots, or estimators.

Inventory scope deliberately excludes binary workbooks and large source bodies.
Field nonmissing counts do not authenticate upstream semantics or human approvals.
"""
import collections
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'sleep_unified_research_v4'
NA = {'', 'na', 'nan', 'none', 'null', 'n/a', 'not_documented', 'unsupplied', 'unknown', 'unresolved', 'pending', 'not_available'}

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(65536), b''):
            h.update(b)
    return h.hexdigest()

def rows(rel):
    with (ROOT / rel).open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def numeric(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None

def bh(ps):
    order = sorted(range(len(ps)), key=ps.__getitem__)
    q = [0.] * len(ps)
    previous = 1.
    for j in range(len(ps)-1, -1, -1):
        i = order[j]
        previous = min(previous, ps[i] * len(ps) / (j+1))
        q[i] = previous
    return q

def png_info(path):
    with path.open('rb') as f:
        assert f.read(8) == b'\x89PNG\r\n\x1a\n'
        result = {}
        while True:
            header = f.read(8)
            if not header:
                break
            n, kind = struct.unpack('>I4s', header)
            data = f.read(n)
            f.read(4)
            if kind == b'IHDR':
                result['width_px'], result['height_px'] = struct.unpack('>II', data[:8])
            elif kind == b'pHYs':
                x,y,unit = struct.unpack('>IIB', data)
                result.update({'physical_unit_meter': unit == 1, 'dpi_x': x * .0254 if unit == 1 else None, 'dpi_y': y * .0254 if unit == 1 else None})
            elif kind == b'IEND':
                break
        return result

def main():
    table_paths = set()
    for base in ['sleep_unified_research_v1', 'sleep_unified_research_v2', 'sleep_unified_research_v3', 'discovery_extension/sleep_submission_evidence_v1/tables']:
        table_paths.update((ROOT/base).rglob('*.tsv'))
    for base in ['config', 'discovery_extension/config', 'discovery_extension/sleep_submission_evidence_v1/sources']:
        table_paths.update((ROOT/base).glob('*.tsv'))
    inventory = []
    for p in sorted(table_paths):
        if p.stat().st_size > 20*1024*1024:
            continue
        with p.open(newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            count = 0
            present = collections.Counter()
            numeric_present = collections.Counter()
            values = collections.defaultdict(collections.Counter)
            for row in reader:
                count += 1
                for k,v in row.items():
                    if k is not None:
                        if v is not None and str(v).strip().lower() not in NA:
                            present[k] += 1
                        if numeric(v) is not None:
                            numeric_present[k] += 1
                        if len(values[k]) < 30 or v in values[k]:
                            values[k][v] += 1
            inventory.append({'path': str(p.relative_to(ROOT)), 'sha256': digest(p), 'bytes': p.stat().st_size, 'rows': count,
                              'columns': reader.fieldnames, 'nonmissing_by_column': {k:present[k] for k in reader.fieldnames or []},
                              'finite_numeric_by_column': {k:numeric_present[k] for k in reader.fieldnames or []},
                              'low_cardinality_values': {k:dict(v) for k,v in values.items() if len(v)<30}})
    (OUT/'statistical_validation/independent_table_inventory_v4.json').write_text(json.dumps({'missing_tokens': sorted(NA), 'scope': 'small TSV files; column presence is not scientific verification', 'tables': inventory}, indent=2)+'\n')

    v1 = ROOT/'sleep_unified_research_v1'
    manifest = json.loads((v1/'manifests/figure_manifest_v1.json').read_text())
    source_hash_checks = [{'path':str((v1/p).relative_to(ROOT)), 'expected':h, 'observed':digest(v1/p), 'pass':digest(v1/p)==h} for p,h in manifest['input_sha256'].items()]
    figures = []
    entries = [(v1, f['figure'], f['outputs'], manifest['output_sha256'], f['source_tables']) for f in manifest['figures']]
    v3 = ROOT/'sleep_unified_research_v3'
    qm = json.loads((v3/'manifests/mvp_QC_figure_receipt_v3.json').read_text())
    entries.append((v3,'MVP_source_preprocessing_QC_v3',list(qm['outputs']),qm['outputs'],list(qm['tables'])))
    for p,h in qm['tables'].items():
        source_hash_checks.append({'path':str((v3/p).relative_to(ROOT)), 'expected':h, 'observed':digest(v3/p), 'pass':digest(v3/p)==h})
    for p,h in [(v3/'logs/mvp_preprocessing_worker_receipt_v3.json',qm['source_receipt_sha256']), (v3/'scripts/24_make_mvp_qc_figure.py',qm['script_sha256'])]:
        source_hash_checks.append({'path':str(p.relative_to(ROOT)), 'expected':h,'observed':digest(p),'pass':digest(p)==h})
    for package, stem, paths, expected, sources in entries:
        f={'id':stem,'sources':[str((package/p).relative_to(ROOT)) for p in sources],'exports':[]}
        for rel in paths:
            p=package/rel
            e={'path':str(p.relative_to(ROOT)), 'sha256':digest(p),'expected_sha256':expected[rel], 'hash_pass':digest(p)==expected[rel],'bytes':p.stat().st_size,'format':p.suffix[1:]}
            if p.suffix == '.png':
                e.update(png_info(p))
            elif p.suffix == '.svg':
                tree=ET.parse(p)
                root=tree.getroot()
                e.update({'width':root.attrib.get('width'),'height':root.attrib.get('height'),'viewBox':root.attrib.get('viewBox'),
                          'text_elements':len(list(root.iter('{http://www.w3.org/2000/svg}text'))),
                          'image_elements':len(list(root.iter('{http://www.w3.org/2000/svg}image'))),
                          'title_elements':len(list(root.iter('{http://www.w3.org/2000/svg}title'))),
                          'desc_elements':len(list(root.iter('{http://www.w3.org/2000/svg}desc')))})
            else:
                info=subprocess.check_output(['/opt/homebrew/bin/pdfinfo',str(p)],text=True)
                img=subprocess.check_output(['/opt/homebrew/bin/pdfimages','-list',str(p)],text=True)
                e['pdfinfo']={line.split(':',1)[0]:line.split(':',1)[1].strip() for line in info.splitlines() if ':' in line}
                e['pdfimages_list']=img.splitlines()
                image_rows=[line for line in img.splitlines()[2:] if line.strip()]
                e['embedded_image_count']=len(image_rows)
                e['vector_status']='mixed vector text/paths and embedded raster' if image_rows else 'vector text/paths; no embedded images listed by Poppler'
            f['exports'].append(e)
        figures.append(f)
    (OUT/'qc/independent_figure_export_probe_v4.json').write_text(json.dumps({'figures':figures, 'source_hash_checks':source_hash_checks,'generation_script_hashes':{str(p.relative_to(ROOT)):digest(p) for p in [v1/'scripts/08_make_evidence_figures.py',v3/'scripts/24_make_mvp_qc_figure.py']}},indent=2)+'\n')

    core=rows('sleep_unified_research_v1/tables/original_core_396.tsv')
    ext=rows('sleep_unified_research_v1/tables/original_extension_1200.tsv')
    rep=rows('sleep_unified_research_v1/REPLICATION_RESULTS.tsv')
    check={}
    for name,data,field,keys in [('core',core,'fdr',['sleep_trait','disease_trait']),('extension',ext,'extension_fdr',['sleep_trait','extension_trait_id'])]:
        q=bh([float(r['p']) for r in data])
        check[name]={'rows':len(data),'unique_pairs':len({tuple(r[k] for k in keys) for r in data}),
                     'bh_positive_count':sum(x<.05 for x in q),'reported_positive_count':sum(float(r[field])<.05 for r in data),
                     'bh_max_abs_discrepancy':max(abs(a-float(r[field])) for a,r in zip(q,data)),
                     'rg_range':[min(float(r['rg']) for r in data),max(float(r['rg']) for r in data)],
                     'finite_rg_se_p':sum(all(numeric(r[k]) is not None for k in ['rg','se','p']) for r in data)}
    check['core']['primary_positive_count']=sum(float(r['fdr'])<.05 and r['analysis_tier']=='PRIMARY_PHASE1' for r in core)
    check['core']['sensitivity_positive_count']=sum(float(r['fdr'])<.05 and r['analysis_tier']!='PRIMARY_PHASE1' for r in core)
    check['core']['tier_counts']=dict(collections.Counter(r['analysis_tier'] for r in core))
    estimated=[r for r in rep if numeric(r['replication_rg']) is not None]
    hz=[float(r['heterogeneity_z']) for r in estimated]
    check['validation']={'locked_rows':len(rep),'unique_pairs':len({r['pair_id'] for r in rep}),'class_counts':dict(collections.Counter(r['replication_class'] for r in rep)),
                         'estimated_rows':len(estimated),'both_trait_independence_values':dict(collections.Counter(r['both_trait_independent_replication'] for r in rep)),
                         'sleep_independence_values':dict(collections.Counter(r['sleep_input_independent'] for r in estimated)),
                         'heterogeneity_z_range':[min(hz),max(hz)],'histogram_outside_minus5_plus5':sum(x < -5 or x > 5 for x in hz),
                         'nominal_heterogeneity_41':sum(float(r['heterogeneity_p'])<.05 for r in estimated),
                         'nominal_heterogeneity_23':sum(float(r['heterogeneity_p'])<.05 and r['replication_class']=='REPLICATED' for r in estimated),
                         'alpha_expected':.05/217,'alpha_max_abs_error':max(abs(float(r['replication_alpha'])-.05/217) for r in estimated),
                         'cov0_z_max_abs_discrepancy':max(abs((float(r['replication_rg'])-float(r['discovery_rg']))/math.hypot(float(r['discovery_se']),float(r['replication_se']))-float(r['heterogeneity_z'])) for r in estimated)}
    outcomes={r['extension_trait_id']:r['phenotype_name'] for r in ext}
    reverse=collections.defaultdict(list)
    for k,v in outcomes.items(): reverse[v].append(k)
    check['extension_duplicate_display_labels']={k:v for k,v in reverse.items() if len(v)>1}
    lipid=[r for r in core if r['disease_trait'] in {'ldl','hdl','triglycerides'}]
    check['lipid_pairwise_intercepts']={'rows':len(lipid),'above_1p2':sum(float(r['h2_int'])>1.2 for r in lipid), 'per_trait':{t:{'min':min(float(r['h2_int']) for r in lipid if r['disease_trait']==t),'max':max(float(r['h2_int']) for r in lipid if r['disease_trait']==t)} for t in ['ldl','hdl','triglycerides']}}
    jobs=json.loads((v1/'manifests/native_reproduction_jobs_v1.json').read_text())['jobs']
    executions=[]
    for base in ['sleep_unified_research_v1','sleep_unified_research_v2','sleep_unified_research_v3']:
        executions.extend((ROOT/base).rglob('*.execution_receipt.json'))
    check['native_campaign']={'planned_jobs':len(jobs),'by_stage':dict(collections.Counter(j['stage'] for j in jobs)),
                              'by_kind':dict(collections.Counter(j['kind'] for j in jobs)),'campaign_execution_receipt_count':len(executions),
                              'campaign_execution_receipt_paths':[str(p.relative_to(ROOT)) for p in executions],
                              'processed_input_rg_pilot_full_precision_capture_count':len(list((v1/'native').rglob('*.full_precision.json'))),
                              'full_raw_to_estimator_core_pairs':0,'full_raw_to_estimator_extension_pairs':0,'full_raw_to_estimator_validation_pairs':0,
                              'important_scope':'A single processed-input pilot is not a completed 33-pair core campaign job.'}
    stages=[]
    for rel in ['sleep_unified_research_v2/logs/'+t+'_native_prefilter_receipt_v2.json' for t in ['hdl','ldl','triglycerides']]+['sleep_unified_research_v3/logs/'+t+'_materialization_receipt_v3.json' for t in ['breast_cancer','ovarian_cancer']]+['sleep_unified_research_v3/logs/mvp_preprocessing_execution_receipt_v3.json']:
        p=ROOT/rel; j=json.loads(p.read_text())
        stages.append({'path':rel,'sha256':digest(p),'status':j['status'],'exit_code':j['exit_code'],'scientific_scope':'preprocessing only','h2_or_rg_estimated':j.get('h2_or_rg_estimated',False),'input_hashes_unchanged':j.get('input_hashes_unchanged')})
    check['successful_preprocessing_stage_units']=stages
    test=json.loads((v3/'manifests/portable_test_receipt_v3.json').read_text())
    check['portable_tests']={'test_count':test['test_count'],'scope':test['scope'],'log_bindings':[{'package':t['package'],'pass':t['pass'],'expected_tests':t['expected_tests'],'log_hash_pass':digest(v3/'logs'/('portable_'+t['package']+'_tests_v3.log'))==t['log_sha256']} for t in test['tests']]}
    check['mvp_counts']=json.loads((v3/'logs/mvp_preprocessing_worker_receipt_v3.json').read_text())['counts']
    check['input_hash_checks_all_pass']=all(c['pass'] for c in source_hash_checks)
    check['export_hash_checks_all_pass']=all(e['hash_pass'] for f in figures for e in f['exports'])
    check['reviewer_code_sha256']=digest(Path(__file__))
    check['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (OUT/'statistical_validation/independent_numeric_table_checks_v4.json').write_text(json.dumps(check,indent=2)+'\n')
    print(json.dumps({'tables':len(inventory),'figures':len(figures),'checks':check},indent=2))

if __name__ == '__main__':
    main()
