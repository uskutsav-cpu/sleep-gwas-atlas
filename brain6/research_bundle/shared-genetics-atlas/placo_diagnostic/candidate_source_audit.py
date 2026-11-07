#!/usr/bin/env python3
"""Source-match completed PLACO+ candidate rows without silently dropping failures."""
import argparse,collections,csv,gzip,hashlib,json,math,pathlib,statistics
ap=argparse.ArgumentParser();ap.add_argument('--audit',type=pathlib.Path,required=True);ap.add_argument('--s',type=pathlib.Path,required=True);ap.add_argument('--fi',type=pathlib.Path,required=True);ap.add_argument('--ffs',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);args=ap.parse_args()
summary=json.loads((args.audit/'completed_audit.json').read_text());hitsfile=args.audit/'audited_threshold_crossings.tsv';assert hashlib.sha256(hitsfile.read_bytes()).hexdigest()==summary['hits_sha256']
with hitsfile.open() as f:hits=list(csv.DictReader(f,delimiter='\t'))
wanted={}
for h in hits:
 key=(h['SNP'],h['A1'],h['A2']);assert h['SNP'] not in wanted or wanted[h['SNP']]==key;wanted[h['SNP']]=key
files={'S':args.s,'FI':args.fi,'FFS':args.ffs};cols={'S':('SNP','A1','A2','BETA','SE','P'),'FI':('variant_id','effect_allele','other_allele','beta','standard_error','p_value'),'FFS':('SNP','A1','A2','BETA','SE','P')};records={};counts={}
for name,path in files.items():
 found=collections.defaultdict(list);cnt=collections.Counter(); keys=cols[name]
 if wanted:
  with (gzip.open(path,'rt') if str(path).endswith('.gz') else path.open()) as f:
   for row in csv.DictReader(f,delimiter='\t'):
    cnt['source_rows']+=1;snp=row[keys[0]]
    if snp not in wanted:continue
    cnt['candidate_rows']+=1;a,b=row[keys[1]].upper(),row[keys[2]].upper(); _,refa,refb=wanted[snp]
    if (a,b)==(refa,refb):orientation=1
    elif (a,b)==(refb,refa):orientation=-1
    else:cnt['allele_mismatch']+=1;continue
    try:
     beta,se,p=(float(row[k]) for k in keys[3:]); assert all(math.isfinite(x) for x in (beta,se,p)) and se>0 and 0<p<=1 and (beta!=0 or p==1)
     z_p=math.copysign(-statistics.NormalDist().inv_cdf(p/2),beta*orientation) if p<1 else 0.;z_wald=orientation*beta/se
    except (ValueError,OverflowError,AssertionError):cnt['invalid']+=1;continue
    found[snp].append((beta,se,p,orientation,z_p,z_wald))
 for snp in wanted:
  z=found[snp];status='UNIQUE_MATCH' if len(z)==1 else ('NO_MATCH' if not z else 'MULTIPLE_VALID_MATCHES');records[(name,snp)]=(status,z[0] if len(z)==1 else None);cnt[status]+=1
 counts[name]=dict(cnt)
rows=[]
for h in hits:
 row={k:h[k] for k in ('pair','SNP','A1','A2','Z_insomnia','Z_frailty','p_placo_plus')};ok=True
 for tag,source,zcol in [('sleep','S','Z_insomnia'),('frailty',h['pair'],'Z_frailty')]:
  state,vals=records[(source,h['SNP'])];row[tag+'_source_status']=state
  for field in ('beta','se','p','orientation','z_p','z_wald','z_p_minus_rounded'):row[tag+'_'+field]=''
  if vals:
   beta,se,p,orient,zp,zw=vals;diff=zp-float(h[zcol]);row.update({tag+'_beta':beta,tag+'_se':se,tag+'_p':p,tag+'_orientation':orient,tag+'_z_p':zp,tag+'_z_wald':zw,tag+'_z_p_minus_rounded':diff})
   if abs(diff)>.000501:row[tag+'_source_status']='P_DERIVED_Z_MISMATCH';ok=False
  else:ok=False
 row['sensitivity_eligible']='YES' if ok else 'NO';rows.append(row)
fields=['pair','SNP','A1','A2','Z_insomnia','Z_frailty','p_placo_plus']+[tag+'_'+field for tag in ('sleep','frailty') for field in ('source_status','beta','se','p','orientation','z_p','z_wald','z_p_minus_rounded')]+['sensitivity_eligible']
args.out.mkdir(parents=True,exist_ok=True)
with (args.out/'candidate_source_matches.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
(args.out/'candidate_source_audit.json').write_text(json.dumps({'n_threshold_crossings':len(hits),'n_unique_snps':len(wanted),'input_hit_sha256':summary['hits_sha256'],'source_paths':{k:str(v) for k,v in files.items()},'source_sizes_bytes':{k:v.stat().st_size for k,v in files.items()},'counts':counts,'eligible':sum(r['sensitivity_eligible']=='YES' for r in rows)},indent=2)+'\n')
print('crossings',len(hits),'eligible',sum(r['sensitivity_eligible']=='YES' for r in rows),'source_counts',dict(counts))
