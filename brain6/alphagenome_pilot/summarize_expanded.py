#!/usr/bin/env python3
"""Build audit-friendly six-region summaries from completed track windows."""
import collections,csv,gzip,hashlib,json,math,struct
from pathlib import Path
P=Path(__file__).resolve().parent;E=P/'results/avi_mirror/variant_evidence_avi.tsv';D=P/'results/jbrowse_expanded';N=2001

def sha(b):return hashlib.sha256(b).hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,cols):
 tmp=p.with_name(p.name+'.tmp')
 with tmp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 tmp.replace(p)
rows=read(E);binding=json.loads((D/'input_binding.json').read_text());assert sha(E.read_bytes())==binding['evidence_sha256'];out=[]; counts=collections.Counter()
for row in rows:
 key=sha(row['atlas_variant'].encode())[:20];folder=D/'variants'/key;meta=folder/'provenance.json';r=dict(row);r['preparation_alphagenome_status']=r.pop('alphagenome_status');r['tissue_prediction_status']='NOT_YET_EXTRACTED'
 for name in ('dnase_log2_ratio_501bp','adult_rna_log2_ratio_501bp','embryonic_plus_rna_log2_ratio_501bp','embryonic_minus_rna_log2_ratio_501bp','prediction_metadata_sha256','prediction_window_sha256'):r[name]=''
 if meta.exists():
  m=json.loads(meta.read_text());v=m['variant'];assert f"{v['chromosome']}:{v['position']}:{v['reference_bases']}>{v['alternate_bases']}"==row['atlas_variant']
  raw=(folder/'track_window.f32.gz').read_bytes();assert sha(raw)==m['track_window_sha256'];arrays=struct.unpack('<16008f',gzip.decompress(raw));assert len(m['tracks'])==4
  for i,t in enumerate(m['tracks']):
   ref=arrays[i*4002:i*4002+2001];alt=arrays[i*4002+2001:(i+1)*4002];rs=sum(ref[750:1251]);as_=sum(alt[750:1251]);assert abs(rs-t['ref_sum_501bp'])<1e-7 and abs(as_-t['alt_sum_501bp'])<1e-7
   ratio=math.log2((as_+1)/(rs+1));assert abs(ratio-t['local_log2_ratio_plus_one'])<1e-10
   tm=t['metadata'];typ=t['output_type']; life=tm.get('biosample_life_stage');strand=tm.get('strand');source=tm.get('data_source')
   if typ=='dnase' and life=='adult':name='dnase_log2_ratio_501bp'
   elif typ=='rna_seq' and life=='adult' and source=='gtex':name='adult_rna_log2_ratio_501bp'
   elif typ=='rna_seq' and life=='embryonic' and strand=='+':name='embryonic_plus_rna_log2_ratio_501bp'
   elif typ=='rna_seq' and life=='embryonic' and strand=='-':name='embryonic_minus_rna_log2_ratio_501bp'
   else:raise ValueError(f'unrecognized track metadata {row["snp"]}')
   if r[name]!='':raise ValueError('duplicate track kind')
   r[name]=ratio
  if any(r[n]=='' for n in ('dnase_log2_ratio_501bp','adult_rna_log2_ratio_501bp','embryonic_plus_rna_log2_ratio_501bp','embryonic_minus_rna_log2_ratio_501bp')):raise ValueError('missing track')
  r['tissue_prediction_status']='PREDICTED';r['prediction_metadata_sha256']=sha(meta.read_bytes());r['prediction_window_sha256']=sha(raw)
 out.append(r);counts[(row['candidate_locus_id'],r['tissue_prediction_status'])]+=1
cols=list(out[0]);write(D/'expanded_variant_evidence.tsv',out,cols)
regions=[]
for cid in sorted({r['candidate_locus_id'] for r in rows}):
 total=sum(1 for r in rows if r['candidate_locus_id']==cid);n=counts[(cid,'PREDICTED')]
 regions.append({'candidate_locus_id':cid,'n_candidate_variants':total,'n_predicted':n,'n_remaining':total-n,'status':'COMPLETE' if n==total else 'INCOMPLETE'})
write(D/'region_coverage.tsv',regions,list(regions[0]))
n=sum(r['tissue_prediction_status']=='PREDICTED' for r in out)
status_text=f'**Status: {n} of {len(out)} variants have complete extracted track windows.**'
if n<len(out):status_text+=' Remaining variants are explicitly marked `NOT_YET_EXTRACTED`; this is not a negative prediction.'
lines=['# Expanded frontal-cortex AlphaGenome analysis','',status_text,'', 'The same four frontal-cortex tracks, 1,048,576-bp input and 2,001-bp extracted window are used for every completed variant. The adult DNase and adult GTEx RNA columns in the evidence table are log₂((alternate 501-bp signal sum + 1)/(reference signal sum + 1)). RNA numbers are local-window summaries, not gene-level expression scores. Two embryonic strand-specific RNA summaries are also retained. These predictions are not measured effects, hypothesis-test p-values, causal-gene assignments, or proof of insomnia–ADHD interactions.','', '| Candidate region | Predicted | Total | Status |','|---|---:|---:|---|']
for r in regions:lines.append(f"| {r['candidate_locus_id']} | {r['n_predicted']} | {r['n_candidate_variants']} | {r['status']} |")
if n==len(out):
 lines += ['', '![Predicted DNase changes across all six regions](dnase_six_regions.svg)', '', 'Each dot is a variant; its vertical value is the predicted alternate-versus-reference DNase ratio in a 501 bp window. Panels use their own vertical scales. The chromosome and coordinate range in each title use GRCh37 candidate IDs; plotted positions use GRCh38.', '']
lines += ['', '## Exploratory largest local DNase changes', '', 'Within each currently completed subset, the row below has the largest absolute predicted adult DNase log₂ ratio. The sign compares the GRCh38 alternate allele with the reference allele; it is not a direction of sleep or disease risk. These are descriptive extremes among variants already extracted, not calibrated significance tests; when a region is incomplete its maximum can change. They do not establish causal variants or trait effects.', '', '| Candidate region | Variant | Signed local DNase log₂ ratio | AVI PHRED interval | Coverage |', '|---|---|---:|---:|---|']
for region in regions:
 subset=[x for x in out if x['candidate_locus_id']==region['candidate_locus_id'] and x['tissue_prediction_status']=='PREDICTED']
 if not subset:lines.append(f"| {region['candidate_locus_id']} | not yet available | — | — | {region['n_predicted']}/{region['n_candidate_variants']} |");continue
 top=sorted(subset,key=lambda x:(-abs(float(x['dnase_log2_ratio_501bp'])),x['atlas_variant']))[0]
 low=float(top['avi_phred_lower']);high=float(top['avi_phred_upper']);ph=f'{low:.3f}' if low==high else f'{low:.3f}–{high:.3f}'
 lines.append(f"| {region['candidate_locus_id']} | [{top['snp']}]({top['atlas_url']}) | {float(top['dnase_log2_ratio_501bp']):+.6f} | {ph} | {region['n_predicted']}/{region['n_candidate_variants']} |")
lines += ['', 'The output table joins predictions to the frozen GWAS/AVI evidence without changing Brain6 evidence classes. The preparation-status column records the earlier pre-API-key snapshot; `tissue_prediction_status` is the current result. Each per-variant gzip file contains the 2,001-bp reference and alternate signal vectors; metadata and source byte-range hashes are retained in its companion provenance file. Model release was not exposed by the proxy. The public service is used for non-commercial research and its outputs are not used to train other ML models.','']
(D/'EXPANDED_REPORT.md').write_text('\n'.join(lines))
receipt={'input_sha256':sha(E.read_bytes()),'expected':len(out),'predicted':n,'status':'COMPLETE' if n==len(out) else 'INCOMPLETE','output_sha256':{f:sha((D/f).read_bytes()) for f in ['expanded_variant_evidence.tsv','region_coverage.tsv','EXPANDED_REPORT.md']},'method':'one variant prediction per unique candidate variant; 1 Mb sequence context; 4 frontal-cortex tracks; local 501-bp summaries; not gene-level inference'}
if n==len(out) and (D/'dnase_six_regions.svg').exists():receipt['output_sha256']['dnase_six_regions.svg']=sha((D/'dnase_six_regions.svg').read_bytes())
(D/'summary_provenance.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps({'predicted':n,'expected':len(out),'status':receipt['status']}))
