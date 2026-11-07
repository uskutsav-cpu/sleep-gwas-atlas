import csv,gzip,json,hashlib,math,collections
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'retrieval_receipt.json').read_text());candidates={117920728:('rs62264764',frozenset(('G','A')),'A'),117923158:('rs62264767',frozenset(('A','C')),'C')};rows=[];summary={}
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(2**20),b''):h.update(b)
 return h.hexdigest()
for src in r['members']:
 cell=src['cell'];file=src['file'];assert sha(file)==src['inner_gzip_sha256'];counts=collections.Counter();peak_near=set();locs=collections.Counter()
 with gzip.open(file,'rt') as f:
  assert f.readline().split()==['snp','beta_coef','pvalue','peak','EffectAllele']
  for line in f:
   counts['rows']+=1;x=line.split();assert len(x)==5
   var,beta,pv,peak,effect=x;v=var.split(':');assert len(v)==4 and v[0]=='chr3';pos=int(v[1])
   if pos not in candidates:continue
   name,alleles,alt=candidates[pos];counts['same_coordinate']+=1
   if frozenset(v[2:])!=alleles:counts['allele_mismatch']+=1;continue
   b=float(beta);q=float(pv);assert math.isfinite(b) and 0<=q<=1 and effect in alleles
   parts=peak.split('_');assert len(parts)==3 and parts[0]=='chr3';start,end=int(parts[1]),int(parts[2]);assert start<end
   altbeta=b if effect==alt else -b
   rows.append({'cell':cell,'rsid':name,'variant':var,'beta_alt':altbeta,'source_beta':b,'source_effect_allele':effect,'p':q,'peak':peak,'peak_start':start,'peak_end':end,'variant_inside_peak_bed':start<=pos-1<end,'peak_overlaps_501bp_window_bed':start<pos+250 and end>pos-251})
   counts['allele_compatible']+=1;locs[name]+=1
 summary[cell]={'counts':dict(counts),'candidate_rows':dict(locs),'full_gzip_read':True,'inner_sha256_verified':True}
for x in rows:x['bonferroni_observed_candidate_rows']=min(1,x['p']*len(rows))
fields=['cell','rsid','variant','beta_alt','source_beta','source_effect_allele','p','bonferroni_observed_candidate_rows','peak','peak_start','peak_end','variant_inside_peak_bed','peak_overlaps_501bp_window_bed']
with (p/'candidate_rows.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'datasets':summary,'candidate_positions_and_alleles_grch38':{str(k):[a,sorted(b),c] for k,(a,b,c) in candidates.items()},'total_candidate_rows':len(rows),'bonferroni_family':'all matching rows for two markers across both files, post hoc and not correction for all atlas or research searches','coordinate_convention':'peak overlap assumes 0-based BED and marker positions 1-based; source peak semantics and marker build require independent verification','candidate_rows_sha256':sha(p/'candidate_rows.tsv'),'source_selection':'Model1A prespecified candidate mediation design; absence is not proof of no effect'}
(p/'scan_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print('MIN rows',sorted(rows,key=lambda x:x['p'])[:5])
