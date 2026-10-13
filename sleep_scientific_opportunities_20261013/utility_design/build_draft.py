#!/usr/bin/env python3
"""Build six NEW DRAFT multirow items. No solvers, downloads, fits or controller.

Public output is authored tasks/methods/hashes only. Source numeric packets/gold
are private on SSD. Existing benchmark and human instrument are read-only.
"""
import argparse,csv,gzip,hashlib,json,math,random
from collections import Counter
from pathlib import Path

SEED=20261013
RG_FIELDS=['analysis_family','family_size','sleep_trait','sleep_source_id','outcome_trait','discovery_outcome_trait_id','outcome_source_id','rg','se','p','ci_lower_95','ci_upper_95','adjusted_p_or_q','adjustment_method']
H2_FIELDS=['analysis_family','trait_id','source_id','h2_observed','se_observed','h2_reported','se_reported','scale','intercept','intercept_se','sample_prevalence_argument','population_prevalence_argument','liability_factor']
BANNED={'result_class','threshold_pass','pairwise_h2_z_ge_4_both','pairwise_h2_intercept_le_1_2_both','fully_independent_two_trait_replication','interpretation','directional_interpretation_status','record_use_scope','causal_interpretation','validation_scope','sleep_source_reused','difference_test','h2_z','h2_z_ge_4','intercept_le_1_2','adjudication_note','match_class'}
TABLES={'core':'core_rg.tsv','extension':'extension_rg.tsv','validation':'external_validation_rg.tsv'}
POLICY={
 'primary_screen':'For this instrument only: a candidate primary marginal finding requires all of frozen-family adjusted P/q <0.05, both exact-source h2/SE >=4, both exact-source intercept <=1.2, and no unresolved exact-source interpretation contract. This is admissibility under preserved atlas rules, not causal evidence or raw-source certification.',
 'uncertainty':'Marginal intervals are separate estimate uncertainties. A difference test or difference interval requires calibrated joint sampling covariance, which is not supplied in these source snapshots. Sort/order or subtract estimates only descriptively.',
 'source_equivalence':'Publication PMID, public sample total, ancestry and similar labels are bibliographic checks. They do not authenticate source release/export, ascertainment/coding/effect direction or estimator/reference compatibility.',
 'source_validity':'These are frozen fitted estimates under historical source/N conventions. Authenticated artifact identity is not independently certified raw-GWAS validity. Do not diagnose confounding solely from an intercept or describe a fitted h2 total as a certified variance fraction.',
 'overlap':'An identical sleep GWAS source is reused evidence, not an independently remeasured sleep GWAS. Exact participant/relative intersections between different sources are not supplied.',
}
SYSTEM='Use only this complete packet. Resolve the natural scientific request against ALL candidate rows. Query wording is a requested estimand, not evidence that a record matches it. Return the closed JSON schema, selecting packet-local evidence IDs and source facts. Use NOT_IDENTIFIED_FROM_PACKET for unsupported positive claims; never infer a biological null from absence or a hold. Do not browse, inspect gold, use other-condition answers or generate answers by code. The packet contains no pre-resolved task decisions. Source-contract facts and general admissibility rules are evidence, not semantic certification.'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def key(r):return tuple(r[x] for x in ('analysis_family','sleep_trait','outcome_trait'))
def num(x):return float(x) if x not in ('',None,'NA','UNKNOWN','UNRESOLVED') else None
def select_fields(r,fields):return {f:r[f] for f in fields}
def walk_names(x):
 if isinstance(x,dict):return set(x)|set().union(*(walk_names(v) for v in x.values()))
 if isinstance(x,list):return set().union(*(walk_names(v) for v in x)) if x else set()
 return set()
def shuffled(xs,tag):
 xs=list(xs);random.Random(int(hashlib.sha256(f'{SEED}|{tag}'.encode()).hexdigest(),16)).shuffle(xs);return xs

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--release',type=Path,required=True);ap.add_argument('--hga',type=Path,required=True);ap.add_argument('--public',type=Path,required=True);ap.add_argument('--private',type=Path,required=True);a=ap.parse_args()
 old=a.repo/'sleep_interpretation_validation_20261012/benchmark';oldtasks=json.loads((old/'final_tasks.json').read_text())+json.loads((old/'development_tasks.json').read_text());exclude={tuple(t['record_key']) for t in oldtasks};pilot=a.repo/'sleep_unified_research_v8/utility_study/tasks.json';pilot_tasks=json.loads(pilot.read_text());examples=a.repo/'sleep_unified_research_v8/resource_tools/scientific_query_examples.json';example=json.loads(examples.read_text());exclude.update(tuple(x['record_key']) for x in example.values() if 'record_key' in x);exclude.update({('validation','insomnia','finngen_r13_K11_CHOLELITH'),('extension','insomnia','panukbb_phecode__785__both_sexes__na__na')})
 # Private HUMAN scoring key is never discovered or read.
 sources={f:read(a.release/f) for f in list(TABLES.values())+['sleep_constructs.tsv','outcome_phenotypes.tsv','heritability.tsv']};rg={key(r):r for f in TABLES.values() for r in sources[f]};sleep={r['trait_id']:r for r in sources['sleep_constructs.tsv']};outcome={(r['analysis_family'],r['trait_id']):r for r in sources['outcome_phenotypes.tsv']};h2={(r['analysis_family'],r['trait_id']):r for r in sources['heritability.tsv']};registry=read(a.repo/'config/public_gwas_sources.tsv');registry_by_id={r['source_id']:r for r in registry}
 dh=json.loads((a.release/'source_direction_hold.json').read_text());uh=json.loads((a.release/'source_uncertainty_holds.json').read_text());hold_keys={tuple(x['record_key']) for x in dh['affected_core_records']}|{tuple(x['record_key']) for x in uh['affected_records']};held_source_ids=sorted({rg[k]['outcome_source_id'] for k in hold_keys});contracts=[{'source_id':sid,'unresolved_contract':('effect allele/disease direction' if sid!='gwas_catalog_GCST90479148' else 'deposited CI confidence level/construction'),'observed_basis':dh['basis'] if sid!='gwas_catalog_GCST90479148' else uh['basis'],'qualification_date':dh['date'] if sid!='gwas_catalog_GCST90479148' else uh['date'],'analysis_family_scope':'core' if sid!='gwas_catalog_GCST90479148' else 'validation'} for sid in held_source_ids]
 # Fixed semantic scenario recipes are chosen for source/construct and prespecified
 # gate contrasts, never effect direction, magnitude, P/q or prior model errors.
 recipes={
  'D01':[('core',s,'breast_cancer') for s in ['sleepiness','napping','sleepdur','accel_sleep_duration']],
  'D02':[('core','snoring',o) for o in ['bmi','t2d','melanoma','uc','ra']],
  'D03':[('core',s,'lung_cancer') for s in ['sleepdur','accel_sleep_duration','chronotype']],
  'D04':[('core','insomnia','mdd')],
  'D05':[('validation',s,o) for s,o in [('snoring','finngen_r13_K11_REFLUX'),('napping','finngen_r13_K11_REFLUX'),('snoring','finngen_r13_K11_DIAHER')]],
  'D06':[('core','chronotype',o) for o in ['t2d','melanoma','ra']],
 }
 for k in recipes['D05'][:]:
  row=rg[k];recipes['D05'].append(('extension',k[1],row['discovery_outcome_trait_id']))
 selected=[k for ks in recipes.values() for k in ks];assert all(k in rg for k in selected);assert len(selected)==len(set(selected));assert not set(selected)&exclude, 'Exact prior record overlap'
 with gzip.open(a.hga/'gwasATLAS_v20191115.txt.gz','rt') as f:hga_rows=list(csv.DictReader(f,delimiter='\t'))
 hga_meta={r['id']:r for r in hga_rows}
 annotation=next(r for r in read(a.repo/'sleep_unified_research_v7/tables/cross_resource_match_annotations.tsv') if r['atlas_record_id']=='core|insomnia|mdd' and r['comparator']=='human_gwas_atlas');candidate_pairs=[tuple(x.rsplit(':',1)[1].split('|')) for x in annotation['candidate_result_refs'].split(';') if x]
 insomnia_pmid=registry_by_id[sleep['insomnia']['source_id']]['pmid'];mdd_pmid=outcome['core','mdd']['PMID'];groups={2:[],1:[],0:[]}
 for p in candidate_pairs:
  count=int(hga_meta[p[0]]['PMID']==insomnia_pmid)+int(hga_meta[p[1]]['PMID']==mdd_pmid);groups[count].append(p)
 hga_pairs=[]
 for count,n in [(2,1),(1,3),(0,2)]:
  hga_pairs += sorted(groups[count],key=lambda p:hashlib.sha256(f'{SEED}|bibliographic|{p}'.encode()).hexdigest())[:n]
 assert len(hga_pairs)==6
 hga_gc={};wanted={tuple(sorted(p)) for p in hga_pairs}
 with gzip.open(a.hga/'gwasATLAS_v20191115_GC.txt.gz','rt') as f:
  for line,row in enumerate(csv.DictReader(f,delimiter='\t'),2):
   k=tuple(sorted((row['id1'],row['id2'])))
   if k in wanted:hga_gc[k]=(row,line)
 assert wanted<=hga_gc.keys()
 tasks=[];golds=[];packets={c:[] for c in ['B','C']};provenance={};all_packet_rg_keys=set()
 def ptr(file,row,keys):
  line=next(i for i,r in enumerate(sources[file],2) if r==row);return {'file':file,'line':line,'row_key':{k:row[k] for k in keys},'sha256':sha(a.release/file)}
 for tid,ks in recipes.items():
  ks=shuffled(ks,tid+'-candidate-order');rows=[rg[k] for k in ks];id_by_key={k:f'{tid}-R{i:02d}' for i,k in enumerate(ks,1)};unique_sleep=sorted({k[1] for k in ks});unique_outcome=sorted({(k[0],k[2]) for k in ks});h2keys=sorted({(k[0],k[1]) for k in ks}|{(k[0],k[2]) for k in ks});valid_h2=[k for k in h2keys if k in h2]
  sleep_rows=[sleep[t] for t in unique_sleep];outcome_rows=[outcome[k] for k in unique_outcome];h2_rows=[h2[k] for k in valid_h2]
  factual_tables={'estimates':[{'evidence_id':id_by_key[key(r)],**select_fields(r,RG_FIELDS)} for r in rows],'sleep_constructs':sleep_rows,'outcome_phenotypes':outcome_rows,'heritability':[select_fields(r,H2_FIELDS) for r in h2_rows]}
  common={'admissibility_policy':POLICY,'source_contract_evidence':contracts,'resource_scope':'Frozen core396, extension1200, historical external217. No new fits. Draft interpretation experiment; candidate/source clearance pending.'}
  if tid=='D04':
   # Target PMID is an existing registry fact, not a guessed publication identity.
   factual_tables['sleep_publication_registry']=[{'source_id':sleep['insomnia']['source_id'],'PMID':insomnia_pmid}]
   factual_tables['hga_metadata']=[hga_meta[x] for x in sorted({i for p in hga_pairs for i in p},key=int)]
   factual_tables['hga_gc']=[{'evidence_id':f'{tid}-H{i:02d}',**hga_gc[tuple(sorted(p))][0]} for i,p in enumerate(shuffled(hga_pairs,'hga_order'),1)]
   common['hga_ordinary_documentation']={'reference':'https://atlas.ctglab.nl/documentation','section':'3. Database features','verified_on':'2026-10-12','paraphrase':'id identifies a GWAS. Trait retains original naming; uniqTrait harmonizes labels without proving identical phenotype definitions. PMID, population and public N describe the source. This contemporary documentation is not asserted to be frozen2019 documentation.'}
  # C performs entity joins only; it adds no task verdict, filter, sorted effect or
  # pre-computed numerical contrast. Every candidate survives unchanged.
  cards=[]
  for r in factual_tables['estimates']:
   k=key(r);sp=sleep[k[1]];op=outcome[k[0],k[2]];assert sp['source_id']==r['sleep_source_id'] and op['source_id']==r['outcome_source_id'];card={'estimate':r,'sleep_construct':sp,'outcome_phenotype':op}
   for role,t in [('sleep',k[1]),('outcome',k[2])]:
    hk=(k[0],t)
    if hk in h2:
     assert h2[hk]['source_id']==(r['sleep_source_id'] if role=='sleep' else r['outcome_source_id']);card[role+'_heritability']=select_fields(h2[hk],H2_FIELDS)
   cards.append(card);all_packet_rg_keys.add(k)
  c_evidence={'source_cards':cards,'common':common}
  extras={k:v for k,v in factual_tables.items() if k not in ['estimates','sleep_constructs','outcome_phenotypes','heritability']}
  if 'hga_gc' in extras:
   hmeta={x['id']:x for x in extras['hga_metadata']}
   c_evidence['conventional_source_cards']=[{'gc':x,'id1_metadata':hmeta[x['id1']],'id2_metadata':hmeta[x['id2']]} for x in extras['hga_gc']]
   c_extras={k:v for k,v in extras.items() if k not in ['hga_gc','hga_metadata']}
  else:c_extras=extras
  if c_extras:c_evidence['other_factual_tables']=c_extras
  b_evidence={'separate_factual_tables':factual_tables,'common':common}
  assert not walk_names(b_evidence)&BANNED and not walk_names(c_evidence)&BANNED
  # Exact reconstruction, including unique-table ordering, certifies parity.
  reconstructed={'estimates':[x['estimate'] for x in cards],'sleep_constructs':[next(x['sleep_construct'] for x in cards if x['sleep_construct']['trait_id']==s) for s in unique_sleep],'outcome_phenotypes':[next(x['outcome_phenotype'] for x in cards if (x['outcome_phenotype']['analysis_family'],x['outcome_phenotype']['trait_id'])==ok) for ok in unique_outcome],'heritability':[]}
  for hk in valid_h2:
   reconstructed['heritability'].append(next(x[role+'_heritability'] for x in cards for role in ['sleep','outcome'] if role+'_heritability' in x and (x[role+'_heritability']['analysis_family'],x[role+'_heritability']['trait_id'])==hk))
  reconstructed.update(c_evidence.get('other_factual_tables',{}))
  if 'conventional_source_cards' in c_evidence:
   hc=c_evidence['conventional_source_cards'];hm={x[role]['id']:x[role] for x in hc for role in ['id1_metadata','id2_metadata']}
   reconstructed['hga_gc']=[x['gc'] for x in hc];reconstructed['hga_metadata']=[hm[i] for i in sorted(hm,key=int)]
  assert reconstructed==factual_tables
  answer={};schema={};reason=[]
  if tid=='D01':
   prompt='A researcher asks for the marginal genetic association between unintended daytime dozing and breast-cancer case status. Which candidate estimates that requested sleep construct? Which supplied candidate is the reported daytime-nap-frequency alternative rather than an interchangeable measurement? Select both evidence IDs and state whether the packet establishes construct equivalence.'
   answer={'requested_measure_evidence_id':id_by_key['core','sleepiness','breast_cancer'],'alternative_napping_evidence_id':id_by_key['core','napping','breast_cancer'],'construct_equivalence':'NOT_IDENTIFIED_FROM_PACKET'};schema={'requested_measure_evidence_id':'evidence_id','alternative_napping_evidence_id':'evidence_id','construct_equivalence':['ESTABLISHED','NOT_IDENTIFIED_FROM_PACKET']};reason=['Resolve measurement roles against wording/coding, not a similar sleep label or a shared cohort. Different definitions do not by themselves quantify measurement equivalence.']
  elif tid=='D02':
   prompt='A figure would highlight primary marginal associations for partner/relative/friend-reported snoring across these five outcomes. Apply the preserved family, exact-source h2/intercept and unresolved-contract requirements together. For each candidate give all applicable blocking reasons; list every admissible evidence ID. Treat empty signed fields as unavailable, not zero.'
   decisions={};admissible=[]
   for r in rows:
    k=key(r);block=[];sh=h2[k[0],k[1]];oh=h2[k[0],k[2]]
    if num(r['adjusted_p_or_q'])>=.05:block.append('FROZEN_ADJUSTMENT_NOT_PASS')
    if min(float(z['h2_observed'])/float(z['se_observed']) for z in [sh,oh])<4:block.append('H2_Z_GATE_FAIL')
    if max(float(z['intercept']) for z in [sh,oh])>1.2:block.append('INTERCEPT_GATE_FAIL')
    if k in hold_keys:block.append('UNRESOLVED_SOURCE_CONTRACT')
    decisions[id_by_key[k]]=sorted(block)
    if not block:admissible.append(id_by_key[k])
   answer={'admissible_evidence_ids':sorted(admissible),'blocking_reasons_by_evidence':decisions};schema={'admissible_evidence_ids':'evidence_id_set','blocking_reasons_by_evidence':{eid:['FROZEN_ADJUSTMENT_NOT_PASS','H2_Z_GATE_FAIL','INTERCEPT_GATE_FAIL','UNRESOLVED_SOURCE_CONTRACT'] for eid in decisions}};reason=['Gold combines independent frozen-family, source-identity, numerical-QC and source-contract facts. q or a nominal interval alone is insufficient.']
  elif tid=='D03':
   prompt='For lung cancer, compare the 24-hour self-reported duration measure including naps with the device measure of detected sleep episodes inside the main sleep-period window. Identify both estimate rows, compute device minus reported rg descriptively, order those two rg values from smaller to larger, and decide whether the supplied evidence establishes a calibrated between-measurement difference.'
   reported=('core','sleepdur','lung_cancer');device=('core','accel_sleep_duration','lung_cancer');answer={'reported_duration_evidence_id':id_by_key[reported],'device_duration_evidence_id':id_by_key[device],'descriptive_device_minus_reported_rg':float(rg[device]['rg'])-float(rg[reported]['rg']),'ascending_rg_evidence_ids':[id_by_key[k] for k in sorted([reported,device],key=lambda k:float(rg[k]['rg']))],'calibrated_difference':'NOT_IDENTIFIED_FROM_PACKET'};schema={'reported_duration_evidence_id':'evidence_id','device_duration_evidence_id':'evidence_id','descriptive_device_minus_reported_rg':'number','ascending_rg_evidence_ids':'ordered_evidence_ids','calibrated_difference':['ESTABLISHED','NOT_IDENTIFIED_FROM_PACKET']};reason=['Numeric order/subtraction are recomposed from matched rows; no supplied statement gives their answer. Marginal uncertainty lacks joint contrast covariance.']
  elif tid=='D04':
   prompt='A reviewer asks whether the Jansen public UKB-only insomnia versus Howard public depression result is simply the same analysis already in Human GWAS ATLAS. Among the six conventional candidate pairs, identify every pair whose two publications and exact public sample totals match the atlas sources. Then state whether those matched facts authenticate identical source release, phenotype/coding and estimator/reference contracts.'
   fit=[];tests={};hm={x['id']:x for x in factual_tables['hga_metadata']}
   for pair in factual_tables['hga_gc']:
    s,o=hm[pair['id1']],hm[pair['id2']]
    # Resolve roles from source PMID rather than file id1/id2 orientation.
    if s['PMID']==mdd_pmid and o['PMID']==insomnia_pmid:s,o=o,s
    pmids=any(ss['PMID']==insomnia_pmid and oo['PMID']==mdd_pmid for ss,oo in [(s,o),(o,s)]);counts=any(num(ss['N'])==num(sleep['insomnia']['N']) and num(oo['N'])==num(outcome['core','mdd']['sample_size']) for ss,oo in [(s,o),(o,s)]);tests[pair['evidence_id']]={'both_pmids_match':pmids,'both_public_N_match':counts}
    if pmids and counts:fit.append(pair['evidence_id'])
   answer={'bibliographic_and_N_match_evidence_ids':sorted(fit),'identity_checks_by_evidence':tests,'exact_estimand_equivalence':'NOT_IDENTIFIED_FROM_PACKET'};schema={'bibliographic_and_N_match_evidence_ids':'evidence_id_set','identity_checks_by_evidence':{eid:{'both_pmids_match':'boolean','both_public_N_match':'boolean'} for eid in tests},'exact_estimand_equivalence':['ESTABLISHED','NOT_IDENTIFIED_FROM_PACKET']};reason=['Positive bibliographic matching is checked per pair, not assigned blanket abstention. Metadata agreement does not authenticate analysis/estimand equivalence.']
  elif tid=='D05':
   prompt='For snoring reported through a partner/relative/friend complaint and gastroesophageal reflux, find the discovery estimate and its Finnish-register outcome-side check among these records. Identify which additional validation row reuses that reflux outcome source for a different sleep construct. From the actual source identities, classify whether the target check reuses the sleep GWAS or documents independent remeasurement of both traits.'
   target=('validation','snoring','finngen_r13_K11_REFLUX');discovery=('extension','snoring',rg[target]['discovery_outcome_trait_id']);other=('validation','napping','finngen_r13_K11_REFLUX');assert rg[target]['sleep_source_id']==rg[discovery]['sleep_source_id'];assert rg[target]['outcome_source_id']!=rg[discovery]['outcome_source_id'];answer={'discovery_evidence_id':id_by_key[discovery],'outcome_check_evidence_id':id_by_key[target],'different_sleep_same_reflux_source_evidence_id':id_by_key[other],'two_trait_confirmation_state':'REUSED_SLEEP_SOURCE'};schema={'discovery_evidence_id':'evidence_id','outcome_check_evidence_id':'evidence_id','different_sleep_same_reflux_source_evidence_id':'evidence_id','two_trait_confirmation_state':['REUSED_SLEEP_SOURCE','INDEPENDENT_REMEASUREMENT_DOCUMENTED','NOT_IDENTIFIED_FROM_PACKET']};reason=['The response binds roles across two families and distinguishes outcome replacement from sleep-source replacement without a supplied independence flag.']
  else:
   prompt='A colleague proposes that using the reported h2 scale instead of the observed scale could rescue every chronotype–outcome candidate for the numerical primary QC gate. Match each outcome to its exact h2 record, calculate h2/SE on both supplied scales, and classify the pairwise gate state using BOTH traits. Do not call a high intercept a proven causal diagnosis of confounding.'
   decisions={}
   for r in rows:
    k=key(r);oh=h2[k[0],k[2]];sh=h2[k[0],k[1]];zo=float(oh['h2_observed'])/float(oh['se_observed']);zr=float(oh['h2_reported'])/float(oh['se_reported']);zpass=min(zo,float(sh['h2_observed'])/float(sh['se_observed']))>=4;ipass=max(float(oh['intercept']),float(sh['intercept']))<=1.2;state='BOTH_PASS' if zpass and ipass else 'H2_Z_GATE_FAIL' if not zpass and ipass else 'INTERCEPT_GATE_FAIL' if zpass else 'BOTH_FAIL';assert math.isclose(zo,zr,rel_tol=1e-10);decisions[id_by_key[k]]={'outcome_Z_observed':zo,'outcome_Z_reported':zr,'pairwise_gate_state':state}
   answer={'numerical_checks_by_evidence':decisions,'scale_change_rescues_all':'NOT_SUPPORTED_BY_NUMERICAL_CHECKS'};schema={'numerical_checks_by_evidence':{eid:{'outcome_Z_observed':'number','outcome_Z_reported':'number','pairwise_gate_state':['BOTH_PASS','H2_Z_GATE_FAIL','INTERCEPT_GATE_FAIL','BOTH_FAIL']} for eid in decisions},'scale_change_rescues_all':['SUPPORTED_BY_NUMERICAL_CHECKS','NOT_SUPPORTED_BY_NUMERICAL_CHECKS']};reason=['The two h2/SE ratios are calculated from original scale-specific values; no h2-Z or gate flag is supplied. Scale transformation does not change this preserved numerical Z eligibility.']
  evidence_ids=[x['evidence_id'] for x in factual_tables['estimates']]+[x['evidence_id'] for x in factual_tables.get('hga_gc',[])]
  def jsonschema(contract,path=()):
   if isinstance(contract,dict):return {'type':'object','additionalProperties':False,'required':list(contract),'properties':{f:jsonschema(v,path+(f,)) for f,v in contract.items()}}
   if isinstance(contract,list):
    return {'type':'array','uniqueItems':True,'items':{'type':'string','enum':contract}} if 'blocking_reasons_by_evidence' in path else {'type':'string','enum':contract}
   if contract in ['evidence_id_set','ordered_evidence_ids']:return {'type':'array','uniqueItems':True,'items':{'type':'string','enum':evidence_ids}}
   if contract=='evidence_id':return {'type':'string','enum':evidence_ids}
   return {'type':contract}
  closed_schema=jsonschema(schema)
  task={'response_schema':closed_schema,'task_id':tid,'status':'DRAFT_NOT_FOR_EXECUTION','prompt':prompt,'response_contract':schema,'operation':{'D01':'construct_role_resolution','D02':'joint_admissibility_screen','D03':'descriptive_recomposition_and_contrast_scope','D04':'cross_resource_identity_binding','D05':'discovery_validation_role_and_source_dependence','D06':'scale_specific_numeric_QC_reasoning'}[tid],'candidate_estimate_count':len(rows),'candidate_order_seed':SEED};tasks.append(task)
  for c,e in [('B',b_evidence),('C',c_evidence)]:
   packet={'task_id':tid,'condition':c,'status':'DRAFT_NOT_FOR_EXECUTION','system_prompt':SYSTEM,'user_prompt':prompt,'response_contract':schema,'response_schema':closed_schema,'evidence':e};packets[c].append(packet);dump(a.private/'packets'/c/f'{tid}.json',packet)
  row_prov=[]
  for r in rows:row_prov.append(ptr(TABLES[r['analysis_family']],r,['analysis_family','sleep_trait','outcome_trait']))
  for file,rr,keys in [('sleep_constructs.tsv',sleep_rows,['trait_id','source_id']),('outcome_phenotypes.tsv',outcome_rows,['analysis_family','trait_id','source_id']),('heritability.tsv',h2_rows,['analysis_family','trait_id','source_id'])]:
   row_prov.extend(ptr(file,r,keys) for r in rr)
  if tid=='D04':
   row_prov.append({'file':'config/public_gwas_sources.tsv','line':next(i for i,r in enumerate(registry,2) if r['source_id']==sleep['insomnia']['source_id']),'field':'pmid','sha256':sha(a.repo/'config/public_gwas_sources.tsv')})
   row_prov.extend({'file':'gwasATLAS_v20191115.txt.gz','line':next(i for i,r in enumerate(hga_rows,2) if r['id']==mid),'metadata_id':mid,'sha256':sha(a.hga/'gwasATLAS_v20191115.txt.gz')} for mid in sorted({i for p in hga_pairs for i in p},key=int))
   row_prov.extend({'file':'gwasATLAS_v20191115_GC.txt.gz','line':hga_gc[tuple(sorted(p))][1],'pair_ids':list(p),'sha256':sha(a.hga/'gwasATLAS_v20191115_GC.txt.gz')} for p in hga_pairs)
  field_prov={'estimates':RG_FIELDS,'sleep_constructs':['trait_id','source_id','construct','wording_device_and_coding','reliability_ascertainment','N','sample_overlap'],'outcome_phenotypes':['analysis_family','trait_id','source_id','phenotype_definition','source_release_label','PMID','sample_size'],'heritability':H2_FIELDS}
  contract_prov=[{'file':f,'json_fields':['basis','affected_core_records' if f=='source_direction_hold.json' else 'affected_records'],'sha256':sha(a.release/f)} for f in ['source_direction_hold.json','source_uncertainty_holds.json']]
  gold={'source_fields_used':field_prov,'source_contract_provenance':contract_prov,'task_id':tid,'status':'DRAFT_PENDING_INDEPENDENT_SCIENTIFIC_ADJUDICATION','condition_answers':{'B':answer,'C':answer},'reasoning_basis':reason,'row_field_provenance':row_prov,'numeric_tolerance':{'absolute':1e-6,'relative':1e-5},'answer_annotation_fields_removed':sorted(BANNED),'facts_reused_across_source_clusters':True};golds.append(gold)
 assert len(tasks)==6 and len(golds)==6 and len(packets['B'])==6 and len(packets['C'])==6
 assert not all_packet_rg_keys&exclude
 dump(a.public/'DRAFT_TASKS.json',tasks);dump(a.private/'DRAFT_GOLD_NOT_FOR_SOLVERS.json',{'status':'DRAFT_NOT_FROZEN_NEVER_GIVE_TO_SOLVERS','items':golds})
 for c in ['B','C']:dump(a.private/f'DRAFT_condition_{c}_inputs.json',{'status':'DRAFT_NOT_FOR_EXECUTION','items':packets[c]})
 manifest={'status':'DRAFT_REVIEW_REQUIRED_NO_SOLVERS','date':'2026-10-13','seed':SEED,'draft_items':6,'conditions':['B','C'],'authentic_new_data':False,'existing_data_reused':True,'exact_previous_benchmark_pair_overlap':0,'excluded_prior_benchmark_record_count':len(oldtasks),'pilot_public_task_file_sha256':sha(pilot),'pilot_private_key_read':False,'exact_pilot_prompt_overlap':sum(t['prompt']==p['prompt'] for t in tasks for p in pilot_tasks),'previous_frozen_instrument_modified':False,'human_instrument_modified':False,'all_new_estimate_record_keys':[list(k) for k in sorted(all_packet_rg_keys)],'B_C_full_factual_object_parity':'6/6 exact reconstruction','pre_resolved_verdicts_suppressed':sorted(BANNED),'source_fact_overlap_limitation':'The12sleep-source dictionary and qualification policies were already exposed; record/question disjointness does not establish independent new-source generalization.','selection_basis':'Fixed construct/source/gate scenarios before new responses; no sign/magnitude/P/q or earlier model error used to rank/select rows. D04 candidates stratified by bibliographic matches, not rg.','source_release_hashes':{f:sha(a.release/f) for f in sources}|{f:sha(a.release/f) for f in ['source_direction_hold.json','source_uncertainty_holds.json']},'hga_hashes':{f:sha(a.hga/f) for f in ['gwasATLAS_v20191115.txt.gz','gwasATLAS_v20191115_GC.txt.gz']},'public_tasks_sha256':sha(a.public/'DRAFT_TASKS.json'),'private_gold_sha256':sha(a.private/'DRAFT_GOLD_NOT_FOR_SOLVERS.json'),'private_packets_sha256':{f'{c}/{t["task_id"]}':sha(a.private/'packets'/c/f'{t["task_id"]}.json') for c in ['B','C'] for t in tasks},'max_packet_bytes':max((a.private/'packets'/c/f'{t["task_id"]}.json').stat().st_size for c in ['B','C'] for t in tasks),'model_responses_generated':0,'new_fits_or_downloads':0,'rights':'Unresolved comparator/source redistribution; numerical evidence/gold private on SSD; public authored prompts/code/metadata only.'}
 dump(a.public/'DRAFT_CONSTRUCTION_MANIFEST.json',manifest);dump(a.private/'DRAFT_CONSTRUCTION_MANIFEST.json',manifest);print(json.dumps({'draft_items':6,'new_estimate_records':len(all_packet_rg_keys),'parity':'6/6','max_packet_bytes':manifest['max_packet_bytes'],'exact_prior_record_overlap':0,'model_responses':0}))
if __name__=='__main__':main()
