import csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent;atlas=p.parents[1];idx=atlas/'replication/finngen_adhd_raw_lookup/source.tbi';pair=Path('/Volumes/Extreme SSD/brain6-mvp-finngen-hfrs-results-20261028.tsv');pred=atlas.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv';url='https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
q=read(pair);c=read(pred);by={x['snp']:x for x in c};assert len(q)==len(by)==425 and len({x['rsid'] for x in q})==425
assert hashlib.md5(idx.read_bytes()).hexdigest()=='dd3757119e402f1184174653496f1599'
groups=defaultdict(list)
for x in q:
 y=by[x['rsid']];assert y['atlas_variant']==x['grch38'] and y['candidate_locus_id']==x['region']; groups[x['region']].append(x)
assert len(groups)==6
observed={}; receipts={}
with pysam.TabixFile(url,index=str(idx)) as tb:
 for region,items in groups.items():
  c,positions=[],[]
  for x in items:
   chrom,pos,_=x['grch38'].split(':');c.append(chrom[3:]);positions.append(int(pos))
  assert len(set(c))==1; chrom=c[0]; start=min(positions);end=max(positions); records=list(tb.fetch(chrom,start-1,end)); wanted=set(positions); h=hashlib.sha256()
  for line in records:
   h.update((line+'\n').encode());s=line.split('\t');assert len(s)==13
   if int(s[1]) in wanted:observed.setdefault((s[0],int(s[1])),[]).append(s)
  receipts[region]={'chromosome':chrom,'start_1based':start,'end_1based':end,'candidate_count':len(items),'fetched_records':len(records),'text_sha256':h.hexdigest()}
out=[]
for x in q:
 chrom,pos,alleles=x['grch38'].split(':');ref,alt=alleles.split('>');matches=[z for z in observed.get((chrom[3:],int(pos)),[]) if x['rsid'] in z[4].split(',') and z[2]==ref and z[3]==alt]; row={'rsid':x['rsid'],'region':x['region'],'grch38':x['grch38'],'adhd_status':'not_found' if not matches else 'multiple' if len(matches)>1 else 'matched','mvp_status':x['mvp_status'],'hfrs_status':x['hfrs_status'],'mvp_p':x['mvp_p'],'hfrs_p':x['hfrs_p'],'adhd_p':'','adhd_beta_alt':'','adhd_se':'','adhd_af':'','max_p':'','max_p_times_425':'','same_direction_all':''}
 if len(matches)==1:
  z=matches[0]
  try:
   pp,beta,se,af=map(float,[z[6],z[8],z[9],z[10]]);assert all(map(math.isfinite,[pp,beta,se,af])) and 0<=pp<=1 and se>0 and 0<=af<=1;row.update(adhd_p=pp,adhd_beta_alt=beta,adhd_se=se,adhd_af=af)
  except (ValueError,AssertionError):row['adhd_status']='invalid'
 if row['adhd_status']=='matched' and x['mvp_status']=='matched' and x['hfrs_status']=='matched':
  ps=[float(x['mvp_p']),float(x['hfrs_p']),float(row['adhd_p'])];assert all(math.isfinite(z) and 0<=z<=1 for z in ps);m=max(ps);signs=[float(x['mvp_beta_alt']),float(x['hfrs_beta_alt']),float(row['adhd_beta_alt'])];row.update(max_p=m,max_p_times_425=min(1,425*m),same_direction_all=(all(z>0 for z in signs) or all(z<0 for z in signs)))
 out.append(row)
assert len(out)==425
with (p/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys(),delimiter='\t');w.writeheader();w.writerows(out)
matched=[x for x in out if x['max_p']!='']; hit=[x for x in matched if x['max_p_times_425']<.05]
(p/'receipt.json').write_text(json.dumps({'url':url,'source_index_md5':hashlib.md5(idx.read_bytes()).hexdigest(),'input_sha256':{str(z):hashlib.sha256(z.read_bytes()).hexdigest() for z in (pair,pred,idx)},'query_intervals':receipts,'adhd_matched':sum(x['adhd_status']=='matched' for x in out),'three_matched':len(matched),'crossings':len(hit),'partial_random_access_only':True},indent=2)+'\n')
print(json.dumps({'adhd_matched':sum(x['adhd_status']=='matched' for x in out),'three_matched':len(matched),'crossings':len(hit),'best':sorted(matched,key=lambda x:x['max_p'])[:5]},indent=2))
