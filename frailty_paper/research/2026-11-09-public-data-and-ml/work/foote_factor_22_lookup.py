import sys,csv,subprocess,hashlib,json,pathlib,collections,math
acc=sys.argv[1];digests={'GCST90624046':'37f7f24f76a11a99e8840687ed847034','GCST90624049':'8d7842206a2d604f1611a9cd9cdf3471','GCST90624052':'0e29a67c2128dfa76780e2e01181c1f4'};assert acc in digests
url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/{acc}.tsv'; candidate=list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv')));ids={r['variant'] for r in candidate};found=collections.defaultdict(list);md5=hashlib.md5(); count=0
proc=subprocess.Popen(['curl','-fsSL','--retry','2','--max-time','240',url],stdout=subprocess.PIPE)
header=None
for raw in proc.stdout:
 md5.update(raw)
 if header is None:header=raw.decode().rstrip('\r\n').split('\t');index=header.index('rs_id');continue
 count+=1;v=raw.decode().rstrip('\r\n').split('\t');rs=v[index]
 if rs in ids:found[rs].append(dict(zip(header,v)))
rc=proc.wait();assert rc==0 and md5.hexdigest()==digests[acc],(rc,md5.hexdigest())
res=[]
for r in candidate:
 a=r['counted_allele'];b=r['other_allele'];hits=[v for v in found[r['variant']] if v['chromosome']==r['chr_GRCh37_assumed'] and v['base_pair_location']==r['position_GRCh37_assumed'] and {v['effect_allele'],v['other_allele']}=={a,b}];o={'variant':r['variant'],'allele':a,'other_allele':b,'rsid_rows':len(found[r['variant']]),'exact_rows':len(hits)}
 if len(hits)==1:
  v=hits[0];sign=1 if v['effect_allele']==a else -1; beta=sign*float(v['beta']);se=float(v['standard_error']);p=float(v['p_value']);assert math.isfinite(beta) and math.isfinite(se) and se>0 and 0<=p<=1
  o.update(beta=beta,se=se,p=p,q_metric_p=v['Q_metric_pvalue'],expected_n=v['n'])
 res.append(o)
out={'accession':acc,'url':url,'source_md5':md5.hexdigest(),'rows_scanned':count,'unique_matches':sum(r['exact_rows']==1 for r in res),'candidates':res,'limitations':'Post hoc selected candidates; factor GWAS Q-filtered, missing rows are not null; expected N is Nhat not count; overlapping contributing studies, factor semantics and allele coding/build require review. No individual data.'}
p=pathlib.Path('outputs/foote_'+acc+'_candidate_lookup.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='candidates'},indent=2))
