#!/usr/bin/env python3
"""Run official AlphaGenome gene/splice/regulatory scorers; resumable per allele."""
import argparse,csv,hashlib,importlib.metadata,json,os,pathlib,sys,time
from alphagenome.data import genome
from alphagenome.models.v1 import dna_client,dna_model,variant_scorers
SCORERS=['RNA_SEQ','SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS','ATAC','DNASE','CHIP_HISTONE','CHIP_TF','CAGE','PROCAP']
p=argparse.ArgumentParser();p.add_argument('--queue',default='outputs/stage5_score_queue.csv');p.add_argument('--limit',type=int,default=0);p.add_argument('--validate-only',action='store_true');p.add_argument('--key-env',default='ALPHAGENOME_API_KEY');p.add_argument('--store',default='/Volumes/Extreme SSD/codex-sleep-frailty-work/stage5_alphagenome_scores');args=p.parse_args()
rows=list(csv.DictReader(open(args.queue)));assert rows and len({r['variant_key'] for r in rows})==len(rows)
for x in rows:
 assert x['chr'].startswith('chr') and int(x['position_GRCh38'])>0 and set(x['reference_allele']+x['alternate_allele'])<=set('ACGT')
 v=genome.Variant(chromosome=x['chr'],position=int(x['position_GRCh38']),reference_bases=x['reference_allele'],alternate_bases=x['alternate_allele']);i=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB);assert v.start==int(x['position_GRCh38'])-1 and i.width==1048576
for k in SCORERS:assert k in variant_scorers.RECOMMENDED_VARIANT_SCORERS
if args.validate_only:
 print(json.dumps({'queue_rows':len(rows),'unique_variant_keys':len({r['variant_key'] for r in rows}),'all_variants_hg38_1based_validated':True,'sequence_length':1048576,'scorers':SCORERS,'client_version':importlib.metadata.version('alphagenome'),'live_inference_performed':False},indent=2));sys.exit(0)
key=os.environ.get(args.key_env) or (os.environ.get('ALPHA_GENOME_API_KEY') if args.key_env=='ALPHAGENOME_API_KEY' else None)
if not key:sys.exit('AlphaGenome API key is not configured; no remote scoring was attempted.')
store=pathlib.Path(args.store);store.mkdir(parents=True,exist_ok=True);scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[x] for x in SCORERS];model=dna_client.create(key,model_version=dna_model.ModelVersion.ALL_FOLDS,timeout=600)
manifest={'client_version':importlib.metadata.version('alphagenome'),'model_version':'ALL_FOLDS','organism':'HOMO_SAPIENS','genome_build':'GRCh38.p13','sequence_length':1048576,'scorers':SCORERS,'api':'official AlphaGenome model API','queue_sha256':hashlib.sha256(pathlib.Path(args.queue).read_bytes()).hexdigest()}
(store/'run_config.json').write_text(json.dumps(manifest,indent=2)+'\n')
metadata=model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS); md_dir=store/'metadata';md_dir.mkdir(exist_ok=True)
for label in ['rna_seq','splice_sites','splice_site_usage','splice_junctions','atac','dnase','chip_histone','chip_tf','cage','procap']:
 table=getattr(metadata,label,None)
 if table is not None:table.to_csv(md_dir/(label+'.csv'),index=False)
selected=rows[:args.limit] if args.limit>0 else rows;ok=0;failed=0
for x in selected:
 digest=hashlib.sha256(x['variant_key'].encode()).hexdigest()[:20];out=store/(digest+'.csv.gz');done=store/(digest+'.json')
 if out.exists() and done.exists():
  old=json.loads(done.read_text())
  previous=old.get('config',{})
  if old.get('variant_key')==x['variant_key'] and old.get('status')=='complete' and old.get('sha256')==hashlib.sha256(out.read_bytes()).hexdigest() and all(previous.get(k)==manifest[k] for k in ('client_version','model_version','organism','genome_build','sequence_length','scorers')):ok+=1;continue
 v=genome.Variant(chromosome=x['chr'],position=int(x['position_GRCh38']),reference_bases=x['reference_allele'],alternate_bases=x['alternate_allele'],name=x['variant_key']);interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
 try:
  scores=model.score_variant(interval=interval,variant=v,variant_scorers=scorers,organism=dna_client.Organism.HOMO_SAPIENS);tidy=variant_scorers.tidy_scores([scores],match_gene_strand=True,include_extended_metadata=True)
  if tidy is None:raise ValueError('score_variant returned no tidy data')
  tmp=store/(digest+'.tmp.csv.gz');tidy.to_csv(tmp,index=False,compression='gzip');tmp.replace(out);receipt={'variant_key':x['variant_key'],'status':'complete','rows':len(tidy),'columns':list(tidy.columns),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'input_labels':x['input_labels'],'file':out.name,'scores':'all returned tracks; gene-matched strand; original raw and quantile scores','config':manifest};done.write_text(json.dumps(receipt,indent=2,default=str)+'\n');ok+=1;print('complete',x['variant_key'],'rows',len(tidy),flush=True)
 except Exception as e:
  err=str(e).replace(key,'[redacted]');done.write_text(json.dumps({'variant_key':x['variant_key'],'status':'failed','error':err[:1200],'config':manifest},indent=2)+'\n');failed+=1;print('failed',x['variant_key'],type(e).__name__,flush=True)
  if 'quota' in err.lower() or 'resource_exhausted' in err.lower() or 'permission_denied' in err.lower() or 'unauthenticated' in err.lower():print('Stopping to avoid repeating an access or quota error.',flush=True);break
print('completed',ok,'failed',failed,'requested',len(selected),flush=True)
if failed:sys.exit(1)
