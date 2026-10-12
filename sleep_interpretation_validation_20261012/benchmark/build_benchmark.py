#!/usr/bin/env python3
"""Construct tasks and private solver packets; no fitting, model calls, or source downloads."""
from __future__ import annotations
import argparse, csv, gzip, hashlib, importlib.util, json, random, shutil
from collections import Counter
from pathlib import Path

SEED = 20261012
STRATA = ['uncertainty', 'phenotype', 'overlap', 'provenance', 'testing_family', 'qc', 'direction', 'validation', 'causal', 'missing_and_contrast']
SCHEMA_VERSION = 'sleep-interpretation-60-v1'
TABLE_MAP = {'core':'core_rg.tsv', 'extension':'extension_rg.tsv', 'validation':'external_validation_rg.tsv'}
# Exact records used in the pre-existing human instrument and six V8 examples.
EXCLUDED = {('core','insomnia','frailty'), ('validation','insomnia','gwas_catalog_GCST90479148'), ('extension','insomnia','panukbb_phecode__785__both_sexes__na__na'), ('validation','insomnia','finngen_r13_K11_CHOLELITH'), ('core','accel_sleep_duration','hdl'), ('core','longsleep','melanoma'), ('core','insomnia','alz'), ('core','insomnia','ibd'), ('core','insomnia','does_not_exist')}
# Metadata-only pilot tasks are also excluded from exact question reuse. New phenotype
# tasks use other source constructs and ask different, pair-specific structured facts.
PHENOTYPE_ALLOWED = {'chronotype','sleepiness','napping','snoring','sleep_apnea'}
DOCS = {
 'resource':'Human GWAS ATLAS', 'snapshot':'v20191115',
 'snapshot_listing':'https://atlas.ctglab.nl/',
 'documentation_url':'https://atlas.ctglab.nl/documentation',
 'documentation_checked':'2026-10-12',
 'documentation_temporal_scope':'Contemporary ordinary documentation, not asserted to be a frozen 2019 document',
 'ordinary_documentation_paraphrase':'The id field identifies a database GWAS. Trait is its original label; uniqTrait harmonizes labels without certifying identical phenotype definitions. PMID, Population and N describe the source. Genetic correlations require eligible heritability estimates, predominantly European samples, explicit effect/non-effect alleles and heritability Z greater than 2. Correlation heatmap multiplicity uses the number of possible selected-GWAS pairs. The conventional snapshot has its own sources and coverage; its labels do not establish exact identity with another resource.',
 'documentation_sections':['3. Database features','5.3. Genetic correlation','7.2. Genetic correlation heatmap'],
 'cached_landing_page':'home.html from 2026-10-10 acquisition; original protected source files retained in place',
}
SYSTEM = ('Answer only using supplied evidence. Return one JSON object conforming to response_schema. '
          'Use null for factual fields unavailable for the exact requested source record. '
          'Do not substitute a similar trait, source, release, or estimator. Distinguish historical '
          'reproduction from cleared inference. A missing record does not establish biological absence. '
          'Do not use outside browsing, tools, previous task answers, private gold, or memory. '
          'status is record state: AVAILABLE means ESTIMATE_PRESENT without a hold; WITHHELD_DIRECTION means HISTORICAL_RECORD_WITH_INTERPRETATION_HOLD; WITHHELD_UNCERTAINTY means HISTORICAL_RECORD_WITH_UNCERTAINTY_HOLD; held rows still exist historically. NO_RELEASE_RECORD means NO_ESTIMATE_IN_THIS_RELEASE_TABLE; INSUFFICIENT_TARGET_SOURCE_EVIDENCE means this packet cannot authenticate the exact target row. '
          'Identifiers in the query are unverified selectors, not observed source evidence; source IDs must be copied from an authenticated evidence row. '
          'For unavailable metadata sentinels (blank, NA, UNKNOWN, UNRESOLVED), use null unless instructed to copy a sentinel exactly. '
          'Return no explanations outside the JSON object.')
STATUS_ENUM = ['AVAILABLE','WITHHELD_DIRECTION','WITHHELD_UNCERTAINTY','NO_RELEASE_RECORD','INSUFFICIENT_TARGET_SOURCE_EVIDENCE']

def dump(path,obj):
 path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+'\n')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def read(path):
 with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def key(r): return tuple(r[x] for x in ('analysis_family','sleep_trait','outcome_trait'))
def rank(k,tag): return hashlib.sha256((str(SEED)+'|'+tag+'|'+'|'.join(k)).encode()).hexdigest()
def num(v): return float(v) if v not in ('','NA',None) else None
def flag(v): return v=='True'
def atomize(v,path=()):
 if isinstance(v,dict):
  return [a for k in sorted(v) for a in atomize(v[k],path+(k,))]
 if isinstance(v,list): return [a for i,x in enumerate(v) for a in atomize(x,path+(str(i),))]
 return [('/'.join(path),v)]

def flat_packet(view, comparators):
 # Disjoint conventional source tables; no row combines estimate and phenotype.
 tables={}
 mapping={'estimate':TABLE_MAP.get(view.get('record_key',['core'])[0],'core_rg.tsv'),
          'sleep_construct':'sleep_constructs.tsv','outcome_phenotype':'outcome_phenotypes.tsv',
          'linked_discovery_estimate':'extension_rg.tsv','sensitivities':'sensitivities.tsv'}
 notes={}
 for field,val in view.items():
  if field in mapping:
   rows=val if isinstance(val,list) else [val]
   tables.setdefault(mapping[field],[]).extend(rows)
  else: notes[field]=val
 return {'tables':tables,'source_and_scope_notes':notes,'human_gwas_atlas':comparators}

def rejoin(packet):
 # Audit serialization only, not an inference engine. Standard table identifiers
 # and family/trait keys recover the same fields that query_atlas returned.
 notes=dict(packet['source_and_scope_notes']); family=notes.get('record_key',['core'])[0]
 tables=packet['tables']; primary=TABLE_MAP[family]
 if primary in tables:
  notes['estimate']=next(r for r in tables[primary] if r['analysis_family']==family)
 if 'sleep_constructs.tsv' in tables: notes['sleep_construct']=tables['sleep_constructs.tsv'][0]
 if 'outcome_phenotypes.tsv' in tables: notes['outcome_phenotype']=tables['outcome_phenotypes.tsv'][0]
 if family=='validation' and 'extension_rg.tsv' in tables: notes['linked_discovery_estimate']=tables['extension_rg.tsv'][0]
 if 'sensitivities.tsv' in tables: notes['sensitivities']=tables['sensitivities.tsv']
 elif 'estimate' in notes: notes['sensitivities']=[]
 return {'atlas_view':notes,'human_gwas_atlas':packet['human_gwas_atlas']}

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--release',type=Path,required=True); ap.add_argument('--query-module',type=Path,required=True); ap.add_argument('--repo',type=Path,required=True); ap.add_argument('--comparator',type=Path,required=True); ap.add_argument('--private',type=Path,required=True); ap.add_argument('--public',type=Path,required=True); a=ap.parse_args()
 if (a.public/'FREEZE_RECEIPT.json').exists(): raise ValueError('Frozen benchmark: do not regenerate without a versioned amendment')
 spec=importlib.util.spec_from_file_location('source_qualified_query',a.query_module); q=importlib.util.module_from_spec(spec); spec.loader.exec_module(q)
 tables={k:read(a.release/f) for k,f in TABLE_MAP.items()}; core=tables['core']; extension=tables['extension']; validation=tables['validation']; allrows={key(r):r for t in tables.values() for r in t}
 sens=read(a.release/'sensitivities.tsv'); meta={r['trait_id']:r for r in read(a.release/'sleep_constructs.tsv')}; outcomes={(r['analysis_family'],r['trait_id']):r for r in read(a.release/'outcome_phenotypes.tsv')}
 with gzip.open(a.comparator/'gwasATLAS_v20191115.txt.gz','rt') as f: hga=list(csv.DictReader(f,delimiter='\t'))
 hga_id={r['id']:r for r in hga}; annotations=read(a.repo/'sleep_unified_research_v7/tables/cross_resource_match_annotations.tsv'); ann={tuple(r['atlas_record_id'].split('|')):r for r in annotations if r['comparator']=='human_gwas_atlas'}
 used=set(EXCLUDED); task_entries=[]; gold_entries=[]; conditions={c:[] for c in 'ABC'}; selected_refs=set(); lookup={}
 # Development selection occurs before final selection, always excluding exact record reuse.
 def select(rows,n,tag):
  candidates=sorted((r for r in rows if key(r) not in used),key=lambda r:rank(key(r),tag))
  if len(candidates)<n: raise ValueError(f'Insufficient candidates: {tag}: {len(candidates)} < {n}')
  picked=candidates[:n]; used.update(key(r) for r in picked); return picked
 dev=[(STRATA[i],select(core if i!=1 else [r for r in core if r['sleep_trait'] in PHENOTYPE_ALLOWED],1,'development-'+STRATA[i])[0],'standard') for i in range(6)]
 picks={}
 picks['direction']=[(r,'held') for r in select([r for r in core if r.get('directional_interpretation_status','').startswith('WITHHOLD')],3,'direction-held')]+[(r,'standard') for r in select([r for r in core if not r.get('directional_interpretation_status','').startswith('WITHHOLD')],3,'direction-available')]
 picks['qc']=[(r,'h2_failed') for r in select([r for r in core if r['pairwise_h2_z_ge_4_both']=='False'],2,'qc-h2')]+[(r,'intercept_failed') for r in select([r for r in core if r['pairwise_h2_intercept_le_1_2_both']=='False'],2,'qc-intercept')]+[(r,'both_pass') for r in select([r for r in core if r['pairwise_h2_z_ge_4_both']=='True' and r['pairwise_h2_intercept_le_1_2_both']=='True' and r['rg']],2,'qc-pass')]
 picks['validation']=[(r,'uncertainty_hold') for r in select([r for r in validation if r['outcome_trait']=='gwas_catalog_GCST90479148'],2,'validation-held')]+[(r,'standard') for r in select([r for r in validation if r['outcome_trait']!='gwas_catalog_GCST90479148'],4,'validation-standard')]
 picks['testing_family']=sum(([(r,'standard') for r in select(tables[f],2,'family-'+f)] for f in ['core','extension','validation']),[])
 for s in ['uncertainty','phenotype','overlap','provenance','causal']:
  rows=[r for r in core if r['rg'] and (s!='phenotype' or r['sleep_trait'] in PHENOTYPE_ALLOWED)]
  picks[s]=[(r,'standard') for r in select(rows,6,s)]
 # Three genuine fitted-specification contrasts with different native pair records.
 sr=[allrows[('core',r['sleep_trait'],r['outcome_trait'])] for r in sens if r['kind']=='rg' and ('core',r['sleep_trait'],r['outcome_trait']) in allrows]
 sr=list({key(r):r for r in sr}.values()); picks['missing_and_contrast']=[(r,'contrast') for r in select(sr,3,'contrast')]
 for i,sleep in enumerate(['chronotype','napping','snoring']):
  r={'analysis_family':'validation','sleep_trait':sleep,'outcome_trait':'gwas_catalog_GCST90479148'}
  assert key(r) not in allrows and key(r) not in used; used.add(key(r)); picks['missing_and_contrast'].append((r,'absent'))
 # Gather every conventional candidate referenced by the frozen annotation, not an
 # effect-favourable result subset. Same candidates appear in all three conditions.
 schedule=[('development',s,r,m) for s,r,m in dev]+[('final',s,r,m) for s in STRATA for r,m in picks[s]]
 for split,s,r,m in schedule:
  k=key(r); ar=ann.get(k); refs=[]
  if ar:
   for ref in ar['candidate_result_refs'].split(';'):
    if ref.startswith('HumanGWASATLAS:v20191115:'):
     left,right=ref.rsplit(':',1)[1].split('|'); refs.append(tuple(sorted((left,right)))); selected_refs.add(tuple(sorted((left,right))))
  lookup[k]=sorted(set(refs),key=lambda p:rank(p,'comparator-candidates'))
 pairs={}; gc_count=0
 with gzip.open(a.comparator/'gwasATLAS_v20191115_GC.txt.gz','rt') as f:
  for line,row in enumerate(csv.DictReader(f,delimiter='\t'),2):
   gc_count+=1; k=tuple(sorted((row['id1'],row['id2'])))
   if k in selected_refs: pairs[k]={'row':row,'line':line}
 assert selected_refs <= pairs.keys()
 def schema(fields):
  return {'type':'object','additionalProperties':False,'required':['status']+list(fields),
          'properties':{'status':{'type':'string','enum':STATUS_ENUM},**fields}}
 def nullable(kind,enum=None):
  x={'type':[kind,'null']}
  if enum is not None:x['enum']=enum+[None]
  return x
 for index,(split,s,r,mode) in enumerate(schedule):
  k=key(r); view=q.query(a.release,*k); estimate=view.get('estimate',{}); sleep=view.get('sleep_construct',{}); outcome=view.get('outcome_phenotype',{}); ar=ann.get(k)
  task_id=('DEV' if split=='development' else 'FIN')+f'{index+1 if split=="development" else index+1-len(dev):03d}'
  comparators={'metadata':[],'genetic_correlations':[],'ordinary_documentation':DOCS}
  ids=set()
  for p in lookup[k]:
   comparators['genetic_correlations'].append(pairs[p]['row']); ids.update(p)
  comparators['metadata']=[hga_id[x] for x in sorted(ids,key=int)]
  prefix=f'For exact atlas record ({k[0]}, {k[1]}, {k[2]})'
  if estimate: prefix+=f', sleep source {estimate["sleep_source_id"]} and outcome source {estimate["outcome_source_id"]}'
  status={'ESTIMATE_PRESENT':'AVAILABLE','HISTORICAL_RECORD_WITH_INTERPRETATION_HOLD':'WITHHELD_DIRECTION','HISTORICAL_RECORD_WITH_UNCERTAINTY_HOLD':'WITHHELD_UNCERTAINTY','NO_ESTIMATE_IN_THIS_RELEASE_TABLE':'NO_RELEASE_RECORD'}[view['status']]
  answer={'status':status}; fields={}; provenance=[]
  def add(name,value,kind='string',enum=None,source='estimate',sourcefield=None):
   fields[name]=nullable(kind,enum); answer[name]=value; provenance.append({'answer_field':name,'resource_section':source,'source_field':sourcefield or name})
  if s=='uncertainty':
   prompt='report rg, SE, the two marginal 95% bounds and whether these data establish zero shared genetic signal.'
   for f in ['rg','se','ci_lower_95','ci_upper_95']: add(f,num(estimate.get(f)),'number')
   add('zero_sharing_established',False,'boolean',source='interpretation_scope',sourcefield='A marginal estimate and its uncertainty do not establish genetic absence')
  elif s=='phenotype':
   prompt='retrieve the sleep construct effect allele and effect column, and state its phenotype measurement identity by copying construct exactly.'
   for f in ['construct','effect_allele','effect_column']: add(f,sleep[f],source='sleep_construct')
  elif s=='overlap':
   prompt='state whether exact participant/relative intersections are quantified and whether an independent-samples claim is supported; copy the sleep source cohort.'
   add('exact_intersections_quantified',False,'boolean',source='sleep_construct',sourcefield='missing_or_unverified_fields')
   add('independent_samples_supported',False,'boolean',source='sleep_construct',sourcefield='sample_overlap')
   add('cohort',sleep['cohort'],source='sleep_construct')
  elif s=='provenance':
   prompt='retrieve the outcome release label, genome build and phenotype type; if a conventional related pair is supplied, copy id1/id2 and rg from the first displayed conventional GC row without substituting it for the atlas result.'
   for f in ['source_release_label','build','phenotype_type']: add(f,outcome[f] if outcome[f] not in ('','NA','UNKNOWN','UNRESOLVED') else None,source='outcome_phenotype')
   first=comparators['genetic_correlations'][0] if comparators['genetic_correlations'] else {}
   for f in ['id1','id2']: add('comparator_'+f,first.get(f),source='human_gwas_atlas',sourcefield=f)
   add('comparator_rg',num(first.get('rg')),'number',source='human_gwas_atlas',sourcefield='rg')
  elif s=='testing_family':
   prompt='copy the original family size, multiplicity method and historical threshold flag; state whether that flag by itself permits reliable biological inference.'
   add('family_size',int(estimate['family_size']),'integer'); add('adjustment_method',estimate['adjustment_method']); add('historical_threshold_pass',flag(estimate['threshold_pass']),'boolean',sourcefield='threshold_pass')
   add('threshold_flag_alone_clears_inference',False,'boolean',source='scope_notes',sourcefield='QC/source uncertainty qualifications are additional requirements')
  elif s=='qc':
   prompt='copy the pairwise h2-Z and intercept flags. Classify the observed gate state; do not let a nominal P value override the gate.'
   h=flag(estimate['pairwise_h2_z_ge_4_both']); i=flag(estimate['pairwise_h2_intercept_le_1_2_both'])
   add('h2_z_gate_pass',h,'boolean',sourcefield='pairwise_h2_z_ge_4_both'); add('intercept_gate_pass',i,'boolean',sourcefield='pairwise_h2_intercept_le_1_2_both')
   state='BOTH_PASS' if h and i else ('H2_Z_FAIL' if not h and i else 'INTERCEPT_FAIL' if h else 'BOTH_FAIL')
   add('qc_state',state,enum=['BOTH_PASS','H2_Z_FAIL','INTERCEPT_FAIL','BOTH_FAIL'],source='gate_classification',sourcefield='Boolean combination of two prespecified flags')
  elif s=='direction':
   prompt='copy directional_interpretation_status. State whether signed biological interpretation is withheld and whether a blank signed estimate establishes a null effect.'
   add('directional_interpretation_status',estimate.get('directional_interpretation_status')); add('signed_interpretation_withheld',status=='WITHHELD_DIRECTION','boolean',source='source_direction_hold',sourcefield='affected_core_records')
   add('blank_signed_field_establishes_null',False,'boolean',source='SCHEMA.md',sourcefield='V8 interpretation quarantine')
  elif s=='validation':
   prompt='report both sleep source IDs for the historical validation and linked discovery rows, whether the sleep source is reused, whether this is fully independent two-trait replication, and whether reliable inference is currently withheld.'
   add('validation_sleep_source_id',estimate['sleep_source_id'],sourcefield='sleep_source_id'); add('discovery_sleep_source_id',view['linked_discovery_estimate']['sleep_source_id'],source='linked_discovery_estimate',sourcefield='sleep_source_id')
   add('sleep_source_reused',view['sleep_source_reused'],'boolean',source='scope_notes'); add('fully_independent_two_trait_replication',False,'boolean')
   add('reliable_inference_withheld',status=='WITHHELD_UNCERTAINTY','boolean',source='source_uncertainty_holds',sourcefield='affected_records')
  elif s=='causal':
   prompt='report rg and classify what the correlation establishes about causal direction and a therapeutic intervention; use only the supplied scope.'
   add('rg',num(estimate.get('rg')),'number'); add('causal_direction','NOT_ESTABLISHED',enum=['ESTABLISHED','NOT_ESTABLISHED'],source='scope_notes',sourcefield='causal_interpretation')
   add('intervention_effect_established',False,'boolean',source='scope_notes',sourcefield='causal_interpretation')
  elif mode=='contrast':
   prompt='copy descriptive_delta from the first displayed sensitivity row and classify its uncertainty. State whether a calibrated between-fit difference test is available.'
   sensrow=view['sensitivities'][0]; add('descriptive_delta',num(sensrow['descriptive_delta']),'number',source='sensitivities',sourcefield='descriptive_delta')
   add('difference_uncertainty','NOT_CALIBRATED',enum=['CALIBRATED','NOT_CALIBRATED'],source='sensitivities',sourcefield='delta_uncertainty'); add('calibrated_difference_test_available',False,'boolean',source='sensitivities',sourcefield='difference_p_value')
  else:
   prompt='state whether the exact estimate is present in this atlas validation table and whether table-specific absence establishes biological absence or global novelty.'
   add('estimate_present_in_requested_table',False,'boolean',source='lookup',sourcefield='status'); add('biological_absence_established',False,'boolean',source='lookup',sourcefield='interpretation'); add('global_novelty_established',False,'boolean',source='lookup',sourcefield='interpretation')
  response_schema=schema(fields); item={'task_id':task_id,'split':split,'stratum':s,'subtype':mode,'record_key':list(k),'prompt':prefix+', '+prompt,'response_schema':response_schema,'source_clusters':[estimate.get('sleep_source_id',k[1]),estimate.get('outcome_source_id',k[2])], 'pair_cluster':'|'.join(k),'selection_seed':SEED}
  task_entries.append(item)
  a_answer={name:None for name in answer}; a_answer['status']='INSUFFICIENT_TARGET_SOURCE_EVIDENCE'
  # Authentic conventional result retrieval is still answerable in the provenance
  # stratum; target-atlas metadata remains unavailable rather than marked wrong.
  if s=='provenance':
   for name in answer:
    if name.startswith('comparator_'): a_answer[name]=answer[name]
  evidence_c={'atlas_view':view,'human_gwas_atlas':comparators}; evidence_b=flat_packet(view,comparators)
  assert rejoin(evidence_b)==evidence_c and atomize(rejoin(evidence_b))==atomize(evidence_c), f'Parity failed: {task_id}'
  for condition,evidence in [('A',{'human_gwas_atlas':comparators}),('B',evidence_b),('C',evidence_c)]:
   conditions[condition].append({'task_id':task_id,'split':split,'system_prompt':SYSTEM,'user_prompt':item['prompt'],'response_schema':response_schema,'evidence':evidence})
  source_prov=[]
  if estimate:
   filename=TABLE_MAP[k[0]]; line=next(j for j,row in enumerate(tables[k[0]],2) if key(row)==k)
   source_prov.append({'file':filename,'line':line,'key':list(k),'sha256':digest(a.release/filename)})
   for filename,rows,kfields,kvalue in [('sleep_constructs.tsv',list(meta.values()),['trait_id'],[k[1]]),('outcome_phenotypes.tsv',list(outcomes.values()),['analysis_family','trait_id'],[k[0],k[2]])]:
    line=next(j for j,row in enumerate(rows,2) if [row[f] for f in kfields]==kvalue); source_prov.append({'file':filename,'line':line,'key':kvalue,'sha256':digest(a.release/filename)})
  if k[0]=='validation' and estimate:
   dkey=('extension',k[1],estimate['discovery_outcome_trait_id']); source_prov.append({'file':'extension_rg.tsv','line':next(j for j,row in enumerate(extension,2) if key(row)==dkey),'key':list(dkey),'sha256':digest(a.release/'extension_rg.tsv')})
  if mode=='contrast':
   sr=view['sensitivities'][0]; source_prov.append({'file':'sensitivities.tsv','line':next(j for j,row in enumerate(sens,2) if row==sr),'key':[sr[x] for x in ['job_id','kind','sleep_trait','outcome_trait']],'sha256':digest(a.release/'sensitivities.tsv')})
  hga_prov=[{'file':'gwasATLAS_v20191115_GC.txt.gz','line':pairs[p]['line'],'key':list(p)} for p in lookup[k]]
  if ar: source_prov.append({'file':'sleep_unified_research_v7/tables/cross_resource_match_annotations.tsv','line':next(j for j,row in enumerate(annotations,2) if row==ar),'key':list(k)+['human_gwas_atlas'],'match_class':ar['match_class']})
  safe_false = {'zero_sharing_established','independent_samples_supported','threshold_flag_alone_clears_inference','blank_signed_field_establishes_null','intervention_effect_established','calibrated_difference_test_available','biological_absence_established','global_novelty_established'}
  accepted_A={}
  for name in answer:
   if name in safe_false: a_answer[name]=False
  if 'causal_direction' in answer: a_answer['causal_direction']='NOT_ESTABLISHED'
  roles={name:('record_state' if name=='status' else 'conventional_factual' if name.startswith('comparator_') else 'admissibility_judgment' if name in safe_false or name=='causal_direction' else 'target_factual') for name in answer}
  golden={'field_roles':roles,'accepted_condition_A_alternatives':accepted_A,'task_id':task_id,'split':split,'stratum':s,'record_key':list(k),'condition_answers':{'A':a_answer,'B':answer,'C':answer},'field_provenance':provenance,'record_provenance':source_prov,'comparator_provenance':hga_prov,'condition_A_qualification':'Conventional snapshot does not authenticate the exact target atlas source/release; abstention is correct. Related evidence retrieval is scored separately where requested.','semantic_gold_status':'DRAFT_PENDING_INDEPENDENT_SCIENTIFIC_ADJUDICATION','numeric_tolerance':{'absolute':1e-6,'relative':1e-5},'related_comparator_count':len(lookup[k])}
  gold_entries.append(golden)
 a.public.mkdir(parents=True,exist_ok=True); a.private.mkdir(parents=True,exist_ok=True)
 dump(a.public/'final_tasks.json',[x for x in task_entries if x['split']=='final']); dump(a.public/'development_tasks.json',[x for x in task_entries if x['split']=='development'])
 dump(a.private/'GOLD_NOT_FOR_SOLVERS.json',{'warning':'PRIVATE. NEVER PROVIDE TO SOLVERS. Scientific gold requires independent source/semantic adjudication.','items':gold_entries})
 for c in 'ABC':dump(a.private/f'condition_{c}_inputs.json',{'condition':c,'schema_version':SCHEMA_VERSION,'items':conditions[c]})
 dump(a.private/'conventional_resource_extract.json',{'metadata':[hga_id[x] for x in sorted({i for p in selected_refs for i in p},key=int)],'genetic_correlations':[pairs[p] for p in sorted(selected_refs)],'documentation':DOCS})
 # Rights-unresolved source bodies and comparator numeric extracts never enter repo.
 shutil.copyfile(a.comparator/'home.html',a.private/'conventional_snapshot_landing_page_20261010.html')
 manifest={'schema_version':SCHEMA_VERSION,'seed':SEED,'final_count':60,'development_count':6,'strata_counts':dict(Counter(x['stratum'] for x in task_entries if x['split']=='final')),'source_release_sha256':{f:digest(a.release/f) for f in list(TABLE_MAP.values())+['sleep_constructs.tsv','outcome_phenotypes.tsv','sensitivities.tsv','source_direction_hold.json','source_uncertainty_holds.json','SCHEMA.md','release_manifest.json']},'query_module_sha256':digest(a.query_module),'comparator':{'release':'v20191115','metadata_sha256':digest(a.comparator/'gwasATLAS_v20191115.txt.gz'),'GC_sha256':digest(a.comparator/'gwasATLAS_v20191115_GC.txt.gz'),'metadata_rows':len(hga),'GC_rows':gc_count,'used_pair_count':len(selected_refs),'authentic_snapshot':True,'statistics_redownloaded':False,'condition_A_exact_source_equivalence_established':False},'public_artifact_sha256':{f:digest(a.public/f) for f in ['final_tasks.json','development_tasks.json']},'private_artifact_sha256':{f:digest(a.private/f) for f in ['GOLD_NOT_FOR_SOLVERS.json']+[f'condition_{c}_inputs.json' for c in 'ABC']},'B_C_canonical_parity':'66/66 exact canonical atoms equal','model_responses_generated':0,'gold_adjudication':'pending independent scientific review; mechanical checks do not certify semantic gold','excluded_record_keys':[list(k) for k in sorted(EXCLUDED)],'metadata_only_pilot_question_exclusions':['continuous reported-duration definition','short/long case-control comparison','three device-source construct identity','GCST90479148 catalog missingness'],'selection_method':'Within prespecified eligibility strata, ascending SHA256(seed|stratum|record_key); no rg, SE, P, q, sign or threshold values used to rank records; QC and hold states define declared strata only.'}
 dump(a.public/'construction_manifest.json',manifest); dump(a.private/'construction_manifest.json',manifest)
 print(json.dumps({'final':60,'development':6,'strata':manifest['strata_counts'],'canonical_parity':manifest['B_C_canonical_parity'],'HGA_pairs':len(selected_refs),'model_responses_generated':0}))

if __name__=='__main__':main()
