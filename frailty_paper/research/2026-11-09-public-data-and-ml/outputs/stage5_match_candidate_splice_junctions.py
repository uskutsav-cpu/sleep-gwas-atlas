#!/usr/bin/env python3
"""Find nearby-coordinate candidate event matches; do not infer concordance."""
import argparse,csv,gzip,hashlib,json,re
from pathlib import Path
from collections import defaultdict,Counter
p=argparse.ArgumentParser();p.add_argument('--store',required=True);p.add_argument('--output',default='outputs/stage5_candidate_splice_junction_matches.csv');a=p.parse_args();store=Path(a.store)
def norm(s):return re.sub('[^a-z0-9]+','',s.lower())
with open('outputs/stage5_variant_gene_tissue_evidence_atlas.csv',newline='') as f:qtl=list(csv.DictReader(f))
index=defaultdict(list)
for q in qtl:
 event=q['molecular_event'];m=re.match(r'^(?:chr)?(\d+|X|Y):(\d+):(\d+):',event,re.I)
 if not m or not q['gene_id'] or not q['qtl_alt_effect']:continue
 if 'splice' not in (q['molecular_measure']+' '+event).lower() and q['molecular_measure'].lower()!='sqtl':continue
 index[(q['variant_GRCh38_ref_alt'],q['gene_id'].split('.')[0])].append((q,int(m[2]),int(m[3])))
fields=['variant_GRCh38_ref_alt','gene_id','gene_name','biosample_name','gtex_tissue','track_name','assay','predicted_junction_start','predicted_junction_end','qtl_event','qtl_source','qtl_dataset','qtl_tissue','qtl_measure','qtl_alt_effect','qtl_p_value','qtl_start','qtl_end','start_offset_pred_minus_qtl','end_offset_pred_minus_qtl','model_raw_score','model_quantile_score','match_type','interpretation']
n=Counter()
with open(a.output,'w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for rp in store.glob('*.json'):
  if rp.name.startswith('._') or rp.name=='run_config.json':continue
  r=json.loads(rp.read_text());key=r.get('variant_key');f=store/r.get('file','')
  if r.get('status')!='complete' or not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=r.get('sha256'):continue
  n['files']+=1
  with gzip.open(f,'rt',newline='') as h:
   for s in csv.DictReader(h):
    if s['output_type']!='SPLICE_JUNCTIONS' or not s['gene_id'] or not s['junction_Start'] or not s['junction_End']:continue
    try:start=int(float(s['junction_Start']));end=int(float(s['junction_End']))
    except ValueError:continue
    for q,qs,qe in index.get((key,s['gene_id'].split('.')[0]),[]):
     if abs(start-qs)>1 or abs(end-qe)>1:continue
     if s.get('gtex_tissue') and norm(q['tissue_or_cell'])==norm(s['gtex_tissue']):match='measured tissue label normalized to model GTEx tissue'
     elif q['tissue_or_cell']=='CD4_Naive' and s.get('biosample_name')=='naive thymus-derived CD4-positive, alpha-beta T cell':match='naive CD4 label; samples differ'
     else:continue
     w.writerow({'variant_GRCh38_ref_alt':key,'gene_id':s['gene_id'],'gene_name':s['gene_name'],'biosample_name':s['biosample_name'],'gtex_tissue':s.get('gtex_tissue',''),'track_name':s['track_name'],'assay':s.get('Assay title',''),'predicted_junction_start':start,'predicted_junction_end':end,'qtl_event':q['molecular_event'],'qtl_source':q['qtl_source'],'qtl_dataset':q['qtl_dataset'],'qtl_tissue':q['tissue_or_cell'],'qtl_measure':q['molecular_measure'],'qtl_alt_effect':q['qtl_alt_effect'],'qtl_p_value':q['qtl_p_value'],'qtl_start':qs,'qtl_end':qe,'start_offset_pred_minus_qtl':start-qs,'end_offset_pred_minus_qtl':end-qe,'model_raw_score':s['raw_score'],'model_quantile_score':s.get('quantile_score',''),'match_type':match,'interpretation':'Candidate coordinate match only; verify coordinate conventions, junction orientation, samples, effect definition and scale before any comparison; no concordance assigned'})
     n['rows']+=1
print(dict(n))
