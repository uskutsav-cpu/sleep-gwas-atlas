import collections,csv,hashlib,json,pathlib,subprocess,sys
p=pathlib.Path(__file__).resolve().parent;acc=sys.argv[1];assert acc in [f'GCST{x}' for x in range(90624046,90624053)]
with (p/'selected.tsv').open() as f: selected={x['rsid']:x for x in csv.DictReader(f,delimiter='\t')}
with (p.parent/'catalog_manifest.tsv').open() as f:man={x['accession']:x for x in csv.DictReader(f,delimiter='\t')}
url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/{acc}.tsv'
proc=subprocess.Popen(['curl','-fsSL','--retry','0','--connect-timeout','30','--max-time','900',url],stdout=subprocess.PIPE)
md5=hashlib.md5();sha=hashlib.sha256();count=0;matches=collections.defaultdict(list);rawbytes=0
header=None
for line in proc.stdout:
 md5.update(line);sha.update(line);rawbytes+=len(line)
 if header is None:
  header=line.decode().rstrip('\r\n').split('\t');assert header==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rs_id','n','Q_metric','Q_metric_pvalue'],header
  continue
 count+=1; fields=line.decode().rstrip('\r\n').split('\t');assert len(fields)==len(header),(count,len(fields))
 if fields[8] in selected:matches[fields[8]].append(dict(zip(header,fields)))
code=proc.wait();assert code==0,('curl exit',code);assert md5.hexdigest()==man[acc]['publisher_summary_md5'],('md5 mismatch',md5.hexdigest());assert count>0
out=[]
for rs,s in selected.items():
 candidates=matches[rs];exact=[r for r in candidates if r['chromosome'].removeprefix('chr')==s['chrom'] and r['base_pair_location']==s['pos37'] and {r['effect_allele'],r['other_allele']}=={s['allele_1'],s['allele_2']}]
 status='UNIQUE_MATCH' if len(exact)==1 else 'MULTIPLE_EXACT' if len(exact)>1 else 'ID_PRESENT_MISMATCH' if candidates else 'ABSENT'
 x={'accession':acc,'rsid':rs,'selection':s['selection'],'status':status,'n_id_rows':len(candidates),'n_exact':len(exact)}
 x.update({col:exact[0][col] if status=='UNIQUE_MATCH' else '' for col in header});out.append(x)
with (p/f'{acc}_crosswalk.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'url':url,'md5':md5.hexdigest(),'publisher_md5':man[acc]['publisher_summary_md5'],'sha256':sha.hexdigest(),'bytes':rawbytes,'source_rows':count,'selection_sha256':hashlib.sha256((p/'selected.tsv').read_bytes()).hexdigest(),'status_counts':dict(collections.Counter(x['status'] for x in out)),'note':'Missing does not distinguish Q filtering, input support, ID mismatch or other exclusion; source output is not independent replication.'}
(p/f'{acc}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(acc,rawbytes,count,receipt['status_counts'],flush=True)
