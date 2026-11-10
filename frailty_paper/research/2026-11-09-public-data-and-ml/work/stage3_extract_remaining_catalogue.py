import pysam,pathlib,csv,json,hashlib,collections,concurrent.futures,time
w=list(csv.DictReader(open('outputs/stage3_remaining_molecular_worklist.csv')));w=[x for x in w if x['candidate']];pathlib.Path('work/stage3_remaining_qtl').mkdir(exist_ok=True)
cols='molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split()
def run(x):
 d=x['dataset_r7'];c=x['candidate'];ch,pos=x['marker_grch38'].split(':');ch=ch.removeprefix('chr');pos=int(pos);lo=max(1,pos-1000000);hi=pos+1000000;idx=pathlib.Path('work')/(d+'.all.tsv.gz.tbi');url=f'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{d}/{d}.all.tsv.gz';out=pathlib.Path('work/stage3_remaining_qtl')/(c+'.tsv');rec={'candidate':c,'dataset':d,'gene_id':x['gene_id'],'gene_symbol':x['gene_symbol'],'tissue':x['tissue_v10'],'chromosome':ch,'region':[lo,hi],'url':url,'index':str(idx)}
 if not idx.exists():rec['status']='index_not_available';return rec
 try:
  f=pysam.TabixFile(url,index=str(idx));rows=0;kept=[];t=time.time()
  for line in f.fetch(ch,lo-1,hi):
   rows+=1;v=line.split('\t')
   if len(v)!=len(cols):raise ValueError('unexpected columns')
   if v[16].split('.')[0]==x['gene_id']:kept.append(dict(zip(cols,v)))
  f.close()
  with out.open('w',newline='') as o:
   writer=csv.DictWriter(o,fieldnames=cols,delimiter='\t');writer.writeheader();writer.writerows(kept)
  rec.update(status='extracted',region_rows_read=rows,gene_rows=len(kept),candidate_rows=sum(v['rsid']==c and v['position']==str(pos) and v['ref']==x['candidate_ref'] and v['alt']==x['candidate_alt'] for v in kept),local_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),elapsed_seconds=round(time.time()-t,1))
 except Exception as e:rec.update(status='error',error=str(e)[:500])
 print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:res=list(ex.map(run,w))
pathlib.Path('outputs/stage3_remaining_catalogue_extract_receipt.json').write_text(json.dumps({'note':'2-Mb candidate-centered windows; release-7 GTEx v8 reprocessed; exact marker allele verification required; remote source checksum not provider-verified','items':res},indent=2)+'\n')
