#!/usr/bin/env python3
"""Regenerate conservative pair-level prior-art crosswalks from retained evidence.
Does not search automatically, estimate priority, or modify historical results.
Requires openpyxl only to inspect source XLSX tables. Run from repository root.
"""
from pathlib import Path
import csv, hashlib, json, re, collections, sys, subprocess
B = Path('discovery_extension')
O = B / 'sleep_submission_evidence_v1'
D = '2026-10-08'

if '--verify-canonical' in sys.argv:
    def read(p):
        with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
    a=read(O/'tables/NOVELTY_MASTER_1200.tsv'); b=read(O/'tables/NOVELTY_REPLICATED_23.tsv'); c=read(O/'tables/novelty_crosswalk.tsv')
    assert len(a)==1200 and len({r['pair_id'] for r in a})==1200 and len(b)==23
    assert {r['pair_id'] for r in b}=={r['pair_id'] for r in a if r['historical_replication_class']=='REPLICATED'}
    assert all(r['adequately_supported_novel_result']=='NO' for r in a)
    counts=collections.Counter(r['pair_id'] for r in c)
    assert all(int(r['direct_comparator_count'])==counts[r['pair_id']] for r in a)
    print('Canonical source-free novelty verification passed:1200 pairs,23 positives,zero certified novel results.')
    sys.exit(0)

if '--acquire' in sys.argv:
    # Official URLs only. Existing bytes never silently replaced; hash drift fails.
    with (O/'sources/literature_acquisition_manifest.tsv').open(newline='') as f:
        manifest=list(csv.DictReader(f,delimiter='\t'))
    for row in manifest:
        if row['audit_required']=='MANUAL' and not Path(row['path']).exists():
            raise RuntimeError('Retrieve official historical supplement manually and verify hash: '+row['path']+' from '+row['url'])
        if row['audit_required']=='MANUAL':
            if hashlib.sha256(Path(row['path']).read_bytes()).hexdigest()!=row['sha256']:
                raise RuntimeError('Historical supplement hash differs: '+row['path'])
            continue
        if row['audit_required']!='YES': continue
        p=Path(row['path'])
        if not p.exists():
            p.parent.mkdir(parents=True,exist_ok=True)
            temp=p.with_suffix(p.suffix+'.download')
            subprocess.run(['curl','--fail','--location','--max-time','60','--output',str(temp),row['url']],check=True)
            if hashlib.sha256(temp.read_bytes()).hexdigest()!=row['sha256']:
                raise RuntimeError('Source hash differs; inspect local download before using: '+str(temp))
            temp.replace(p)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            raise RuntimeError('Cached source hash differs: '+str(p))
    from pypdf import PdfReader
    text_path=O/'sources/literature_sun2022_gerd.txt'
    extracted='\n'.join(page.extract_text() for page in PdfReader(O/'sources/literature_sun2022_gerd.pdf').pages)
    expected='64a9c40818f9042be518df15deb00b200a8a2a678fa65147d13a7c24000176ec'
    if hashlib.sha256(extracted.encode()).hexdigest()!=expected:
        raise RuntimeError('PDF extraction drift; inspect pypdf version/text before using derived evidence.')
    text_path.write_text(extracted)
    print('Required official-source bytes verified. Historical extracted inventory TSVs remain repository prerequisites.')
    sys.exit(0)

from openpyxl import load_workbook

def tsv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', str(s).lower()).strip()

def xrows(path, sheet):
    w = load_workbook(path, read_only=True, data_only=True)
    rows = list(w[sheet].values); w.close()
    return rows

G = tsv(B / 'results/ldsc/extension_rg_matrix.tsv')
H = {x['pair_id']: x for x in tsv(B / 'results/novelty/extension_novelty_audit.tsv')}
R = {x['pair_id']: x for x in tsv(B / 'results/replication/replication_results.tsv')}
P = {x['extension_trait_id']: x for x in tsv(B / 'results/panukbb_novelty_prescreen.tsv')}
I = tsv(B / 'results/prior_sleep_screen_inventory.tsv')
assert len(G) == 1200 and len({(r['sleep_trait'],r['extension_trait_id']) for r in G}) == 1200
assert len(H) == 603 and len(R) == 217

# Audit source worksheets directly rather than relying on historical coverage labels.
M = {}
workbook = B / 'provenance/prior_screens/morrison_2024_supplementary_tables_1_15.xlsx'
w = load_workbook(workbook, read_only=True, data_only=True)
extract_audit=[]
for f in sorted((B/'provenance/prior_screens/extracted').glob('morrison*LDSC.tsv')):
    extracted=tsv(f); number=re.search(r'_S(\d+)_', f.name).group(1)
    candidate=[s for s in w if re.match(rf'(S|Table\s*S?){number}(\D|$)', s.title)]
    if not candidate:
        candidate=[s for s in w if number in s.title.split(' ')[:2]]
    # Explicit names verified by inspection, not inferred row-count matches.
    if len(candidate)!=1:
        raise AssertionError((number,w.sheetnames))
    sheet=candidate[0]; rows=list(sheet.values)
    header=next(i for i,r in enumerate(rows) if 'phenotype' in r and 'rg' in r)
    cols=list(rows[header]); raw=[dict(zip(cols,r)) for r in rows[header+1:] if r[cols.index('phenotype')] is not None]
    assert len(raw)==len(extracted)==1403
    for a,b in zip(raw,extracted):
        assert str(a['phenotype'])==b['phenotype']
        for k in ['rg','se','p']:
            try: assert abs(float(a[k])-float(b[k])) <= max(1e-300,abs(float(a[k]))*1e-13)
            except (TypeError,ValueError): assert str(a[k])==str(b[k]) or str(a[k])=='None'
    extract_audit.append(dict(file=str(f),sheet=sheet.title,rows=len(raw),status='SOURCE_WORKBOOK_REEXTRACTION_VERIFIED',sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    if number in ['13','14']:
        key='sleepdur' if number=='13' else 'insomnia'
        M[key]=[(x, str(f),i+2) for i,x in enumerate(extracted)]
w.close()
# Independently inspect other historical workbooks; SD17 is TWAS, not rg evidence.
for filename,sheet,header_index,key in [
 ('goodman_2025_supplementary_data_1_27.xlsx','SD24',4,1),
 ('dashti_2019_supplementary_data_18.xlsx','Supplementary Data 18',3,1)]:
    rows=xrows(B/'provenance/prior_screens'/filename,sheet)
    data=[r for r in rows[header_index+1:] if len(r)>key and r[key] is not None]
    extract_audit.append(dict(file=str(B/'provenance/prior_screens'/filename),sheet=sheet,rows=len(data),status='SOURCE_WORKBOOK_CONTENT_INSPECTED',sha256=hashlib.sha256((B/'provenance/prior_screens'/filename).read_bytes()).hexdigest()))
    assert len(data)==(380 if 'goodman' in filename else 45)

# Full Nature 2026 SF5; blocks b/c are selected reported rows, NOT full 527-test families.
S=[]
for i,row in enumerate(xrows(O/'sources/literature_sleepchart_supplement.xlsx','SF5')[3:],4):
    for offset,block in [(0,'a_original_sleepchart'),(8,'b_Dashti_2019'),(16,'c_Austin_Zimmerman')]:
        if row[offset]:
            S.append(dict(block=block,excel_row=i,sleep=str(row[offset]),outcome=str(row[offset+1]),dataset=str(row[offset+2]),rg=row[offset+3],se=row[offset+4],z=row[offset+5],p=row[offset+6]))
write(O/'sources/literature_sleepchart_SF5_extracted.tsv',S)
U=[]
for i,row in enumerate(xrows(O/'sources/literature_multiorgan_MOESM12.xlsx','Sheet1')[2:],3):
    if row[0] and row[1]:
        U.append(dict(excel_row=i,sleep=row[0],outcome=row[1],rg=row[2],se=row[3],p=row[4],q=row[5],category=row[6]))
write(O/'sources/literature_multiorgan_SD11_extracted.tsv',U)
V={}
for sid,sleep,col,doi in [('wang2019','sleepiness',2,'10.1038/s41467-019-11456-7'),('dashti2021_napping','napping',1,'10.1038/s41467-020-20585-3')]:
    path=O/('sources/literature_'+sid+'_rg.xlsx')
    data=[(i,r) for i,r in enumerate(xrows(path,'Sheet1')[4:],5) if r[col]]
    assert len(data)==(232 if sid=='wang2019' else 257)
    V[sleep]=[(sid,doi,i,r,col) for i,r in data]
    extract_audit.append(dict(file=str(path),sheet='Sheet1',rows=len(data),status='SOURCE_WORKBOOK_CONTENT_INSPECTED',sha256=hashlib.sha256(path.read_bytes()).hexdigest()))

# These explicit aliases encode sufficiently related disease definitions, never byte-identical GWAS.
aliases={
 'Noninfectious gastroenteritis':['Noninfectious colitis NAS','Diagnoses - main ICD10: K52 Other non-infective gastro-enteritis and colitis'],
 'Abdominal pain':['Pain type(s) experienced in last month: Stomach or abdominal pain','Stomach/abdominal pain for 3+ months'],
 'Diaphragmatic hernia':['Diagnoses - main ICD10: K44 Diaphragmatic hernia'],
}
FG={'K21 Gastro-oesophageal reflux disease':'K11_REFLUX','J44 Other chronic obstructive pulmonary disease':'J10_COPD','M51 Other intervertebral disk disorders':'M13_INTERVERTEB','Diaphragmatic hernia':'K11_DIAHER','Noninfectious gastroenteritis':'K11_OTHENTERCOL','K80 Cholelithiasis':'K11_CHOLELITH'}
U_sleep={'insomnia':'Insomnia - Jansen et al., 2019','shortsleep':'Short sleep duration - Dashti et al., 2019','longsleep':'Long sleep duration - Dashti et al., 2019','sleepdur':'Sleep duration - Dashti et al., 2019','snoring':'Snoring - Campos et al., 2020','napping':'Daytime nap - Dashti et al., 2021','chronotype':'Chronotype - Jones et al., 2019a','accel_sleep_duration':'Accelerometer based sleep duration - Jones et al., 2019'}
# Sun 2022 Table 2 is extracted directly from retained PDF text, both Neale/FinnGen columns.
T=[]
for line in (O/'sources/literature_sun2022_gerd.txt').read_text().splitlines():
    m=re.match(r'^(Insomnia|Short sleep|Long sleep|Daytime sleepiness|Daytime napping)\s+(.+)$',line)
    if m:
        vals=m.group(2).split()
        if len(vals)==8:
            try:
                nums=[float(v.replace('−','-')) for v in vals]
                for off,cohort in [(0,'Neale'),(4,'FinnGen_R6')]:
                    T.append(dict(sleep=m.group(1),external='GORD_'+cohort,rg=nums[off],se=nums[off+1],z=nums[off+2],p=nums[off+3]))
            except ValueError: pass
assert len(T)==10
write(O/'sources/literature_sun2022_table2_extracted.tsv',T)
T_sleep={'insomnia':'Insomnia','shortsleep':'Short sleep','longsleep':'Long sleep','sleepiness':'Daytime sleepiness','napping':'Daytime napping'}
E=[]; C=[]
for g in G:
    pid=g['sleep_trait']+'__'+g['extension_trait_id']; h=H.get(pid,{})
    name=g['phenotype_name']; sleep=g['sleep_trait']; targets={norm(name)}|{norm(x) for x in aliases.get(name,[])}
    ev=[]
    for sid,doi,i,row,col in V.get(sleep,[]):
        if norm(row[col])==norm(name):
            off=6 if sid=='wang2019' else 4
            ev.append(dict(pair_id=pid,source_id=sid,location='sources/literature_'+sid+'_rg.xlsx#Sheet1:'+str(i),prior_sleep=sleep,prior_external=row[col],prior_rg=row[off],prior_se=row[off+1],prior_p=row[off+2],prior_rg_reoriented=row[off],comparability='SAME_SLEEP_CONSTRUCT_OUTCOME_PRIOR_GWAS;SOURCE_SUBSETS_UNRESOLVED',orientation='SAME_DIRECTION',doi=doi))
    for row,f,line in M.get(sleep,[]):
        stripped=re.sub(r'^Diagnoses - main ICD10:\s*','',row['phenotype'])
        # Categorical labels without coding are ambiguous; only full explicitly aliased labels match.
        if norm(row['phenotype']) in targets or norm(stripped) in targets:
            ev.append(dict(pair_id=pid,source_id='morrison_2024',location=f+':'+str(line),prior_sleep=row['sleep'],prior_external=row['phenotype'],prior_rg=row['rg'],prior_se=row['se'],prior_p=row['p'],prior_rg_reoriented=str(-float(row['rg'])) if sleep=='insomnia' else row['rg'],comparability='SINGLE_INDICATOR_SLEEP_SIMILAR_OUTCOME_DIFFERENT_GWAS_RELEASE',orientation='REVERSE_CODED_HEALTHY_NON_INSOMNIA' if sleep=='insomnia' else 'DURATION_SAME_DIRECTION',doi='10.1093/sleep/zsad320'))
    for row in S:
        ss='short' if sleep=='shortsleep' else 'long' if sleep=='longsleep' else None
        if ss and row['sleep']==ss and row['outcome']==FG.get(name):
            ev.append(dict(pair_id=pid,source_id='sleepchart_2026',location='sources/literature_sleepchart_supplement.xlsx#SF5:'+str(row['excel_row'])+':'+row['block'],prior_sleep=row['sleep'],prior_external=row['outcome'],prior_rg=row['rg'],prior_se=row['se'],prior_p=row['p'],prior_rg_reoriented=row['rg'],comparability='SAME_DASHTI_SLEEP_GWAS_OUTCOME_FINNGEN_EARLIER_RELEASE' if row['block']=='b_Dashti_2019' else 'RELATED_AUSTIN_ZIMMERMAN_SLEEP_GWAS_OUTCOME_FINNGEN',orientation='SAME_DIRECTION',doi='10.1038/s41586-026-10524-5'))
    if 'chronic obstructive pulmonary' in name.lower() and sleep in U_sleep:
        for row in U:
            if row['sleep']==U_sleep[sleep] and str(row['outcome']).startswith('Chronic obstructive pulmonary disease'):
                ev.append(dict(pair_id=pid,source_id='multiorgan_2026',location='sources/literature_multiorgan_MOESM12.xlsx#Sheet1:'+str(row['excel_row']),prior_sleep=row['sleep'],prior_external=row['outcome'],prior_rg=row['rg'],prior_se=row['se'],prior_p=row['p'],prior_rg_reoriented=row['rg'],comparability='PUBLISHED_SLEEP_CONSTRUCT_COPD_DIFFERENT_SOURCE_CHECK_23ANDME_SUBSET',orientation='SAME_DIRECTION',doi='10.1038/s43856-026-01656-w'))
    if name=='K21 Gastro-oesophageal reflux disease' and sleep in T_sleep:
        for row in T:
            if row['sleep']==T_sleep[sleep]:
                ev.append(dict(pair_id=pid,source_id='sun_2022_gerd',location='sources/literature_sun2022_gerd.pdf#Table2',prior_sleep=row['sleep'],prior_external=row['external'],prior_rg=row['rg'],prior_se=row['se'],prior_p=row['p'],prior_rg_reoriented=row['rg'],comparability='SAME_SLEEP_CONSTRUCT_DIFFERENT_GWAS_SOURCE;FINNGEN_R6_PRIOR_OUTCOME_VALIDATION',orientation='SAME_DIRECTION',doi='10.3389/fnut.2022.1009122'))
    for prior in ev:
        prior['prior_nominal_p_lt_005']=str(float(prior['prior_p'])<.05)
        prior['prior_family_corrected_support']='NOT_CERTIFIED_BY_THIS_CROSSWALK;CONSULT_ORIGINAL_FAMILY'
    E.extend(ev)
    coverage=P[g['extension_trait_id']]['prior_sleep_screen_coverage']
    if ev: cls='PREVIOUSLY_PUBLISHED_SUBSTANTIALLY_SIMILAR'; status='DIRECT_RESULT_TABLE_COMPARATOR_FOUND'
    elif coverage and coverage not in ['NO_EXACT_LOCAL_INVENTORY_MATCH','NA','NONE','']:
        cls='PUBLISHED_ASSOCIATION_DIFFERENT_SLEEP_CONSTRUCT_OR_UNRESOLVED';status='BROAD_SCREEN_EXTERNAL_COVERAGE_NOT_EXACT_PAIR'
    else: cls='LITERATURE_COMPARISON_UNRESOLVED';status='NO_COMPARATOR_IN_AUDITED_LOCAL_TABLES'
    contextual=[]
    if name=='Tobacco use disorder' and sleep in ['insomnia','shortsleep','longsleep','sleepdur','chronotype']: contextual=['10.1093/ntr/nty230','10.1016/j.drugalcdep.2020.108151','10.1016/j.drugalcdep.2020.108313']
    elif name=='Abdominal pain' and sleep in ['insomnia','shortsleep','longsleep','sleepdur','napping','chronotype']:contextual=['https://pmc.ncbi.nlm.nih.gov/articles/PMC8271146/','https://pmc.ncbi.nlm.nih.gov/articles/PMC11347297/']
    if name=='K21 Gastro-oesophageal reflux disease' and sleep=='snoring':contextual=['https://pmc.ncbi.nlm.nih.gov/articles/PMC10716286/']
    if contextual and not ev:cls='PUBLISHED_ASSOCIATION_SUBSTANTIALLY_SIMILAR_GENETIC_OR_MR_ONLY'
    if sleep=='insomnia' and 'Gastro-oesophageal reflux' in name:contextual+=['10.1016/j.jpsychires.2024.02.030']
    C.append(dict(pair_id=pid,sleep_trait=sleep,extension_trait_id=g['extension_trait_id'],phenotype_name=name,phenotype_domain=g['phenotype_domain'],rg=g['rg'],se=g['se'],p=g['p'],extension_fdr=g['extension_fdr'],historical_novelty_class=h.get('novelty_class','NOT_HISTORICALLY_PAIR_SEARCHED'),historical_search_date=h.get('search_date','NOT_SEARCHED'),historical_audit_status=h.get('audit_status','NOT_SEARCHED'),current_novelty_class=cls,current_pair_evidence_status=status,local_full_supplement_review_date=D,current_pair_specific_database_search='NOT_RERUN_ALL_1200;FOCUSED_DOMAIN_SEARCHES_ONLY',literature_coverage='BOUNDED_PRIMARY_TABLE_AND_DOMAIN_AUDIT;NOT_SYSTEMATIC_EXHAUSTIVE_REVIEW',prior_broad_screen_external_coverage=coverage,direct_comparator_count=len(ev),direct_comparator_sources=';'.join(sorted({x['source_id'] for x in ev})) or 'NONE_IN_AUDITED_TABLES',prior_dois=';'.join(sorted({x['doi'] for x in ev})),contextual_sources=';'.join(contextual),new_cohort_replication_of_prior_relationship='OUTCOME_SIDE_VALIDATION_OR_PHENOTYPE_SENSITIVITY' if pid in R and R[pid]['replication_class']=='REPLICATED' and ev else 'NOT_ESTABLISHED_BY_THIS_AUDIT',potential_underreported_relationship='UNRESOLVED;NO_PRIORITY_CLAIM',new_quantitative_characterization='VERSIONED_ESTIMATES_AND_SOURCE_COMPARISON',adequately_supported_novel_result='NO',historical_replication_class=R.get(pid,{}).get('replication_class','NOT_IN_217_FAMILY'),remaining_review='Confirm exact source subsets, historical phenotype coding, full 2020/2021 screen supplements and current indexing; author/specialist review required.'))
rep=[x for x in C if x['historical_replication_class']=='REPLICATED'];assert len(rep)==23
write(O/'tables/NOVELTY_MASTER_1200.tsv',C);write(O/'tables/NOVELTY_REPLICATED_23.tsv',rep);write(O/'tables/novelty_crosswalk.tsv',E)
receipt=dict(access_date=D,universe=len(C),historical_searched_positive_pairs=len(H),historical_unsearched_pairs=len(C)-len(H),replicated_pairs=len(rep),full_morrison_ldsc_rows_verified=sum(x['rows'] for x in extract_audit if x['status']=='SOURCE_WORKBOOK_REEXTRACTION_VERIFIED'),sleepchart_sf5_block_counts=dict(collections.Counter(x['block'] for x in S)),multiorgan_sd11_rows=len(U),pair_class_counts=dict(collections.Counter(x['current_novelty_class'] for x in C)),replicated_class_counts=dict(collections.Counter(x['current_novelty_class'] for x in rep)),supported_novel_results=0,crosswalk_rows=len(E),qualified_limits='Current complete-pair database searches not executed; exact source reuse unresolved except explicitly identified supplement labels.',source_workbook_checks=extract_audit)
(O/'logs').mkdir(exist_ok=True)
(O/'logs/literature_audit_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
files=list((B/'provenance/prior_screens').glob('*'))+list((B/'provenance/prior_screens/extracted').glob('*'))+list((O/'sources').glob('literature*'))+list((O/'sources').glob('journal*'))
write(O/'sources/literature_source_hashes.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),access_date=D,verification='LOCAL_RETAINED_BYTES;SEE_FETCH_LEDGER_FOR_HTTP_STATUS') for p in sorted(files) if p.is_file() and p.name!='literature_source_hashes.tsv'])
print(json.dumps({k:v for k,v in receipt.items() if k!='source_workbook_checks'},indent=2))
