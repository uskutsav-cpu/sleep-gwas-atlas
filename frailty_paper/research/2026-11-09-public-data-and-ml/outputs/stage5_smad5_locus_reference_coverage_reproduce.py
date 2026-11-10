import csv,collections,pathlib,json,hashlib
from pyliftover import LiftOver
f=LiftOver('work/hg38ToHg19.over.chain.gz');b=LiftOver('work/hg19ToHg38.over.chain.gz');d=list(csv.DictReader(open('outputs/stage5_smad5_splice_Mak_FI_with_credible_indels_aligned.csv')));ds=[x for x in d if len(x['ref'])==len(x['alt'])==1];lookup=collections.defaultdict(list);out=[]
for x in ds:
 k=f"chr5:{x['position_GRCh38']}:{x['ref']}>{x['alt']}";v=f.convert_coordinate('chr5',int(x['position_GRCh38'])-1)
 if len(v)!=1 or v[0][2]!='+':out.append({'variant_GRCh38':k,'status':'nonunique_or_nonpositive_liftover'});continue
 z=v[0];w=b.convert_coordinate(z[0],z[1]);
 if len(w)!=1 or w[0][:3]!=('chr5',int(x['position_GRCh38'])-1,'+'):out.append({'variant_GRCh38':k,'status':'reciprocal_failure'});continue
 key=(z[0].removeprefix('chr'),str(z[1]+1),frozenset((x['ref'],x['alt'])));lookup[key].append({'variant_GRCh38':k,'chr37':key[0],'pos37':key[1],'alt38':x['alt'],'variant_id':x['variant']})
found=collections.defaultdict(list)
with open('work/stage4_magma_reference/g1000_eur.bim') as fp:
 for idx,line in enumerate(fp):
  v=line.split();key=(v[0],v[3],frozenset(v[4:6]))
  if key in lookup:found[key].append((idx,v))
for key,items in lookup.items():
 for x in items:
  hits=found[key];r={'variant_GRCh38':x['variant_GRCh38'],'chr37':x['chr37'],'pos37':x['pos37'],'source_id':x['variant_id'],'status':'unique allele match' if len(hits)==1 and len(items)==1 else 'missing' if not hits else 'ambiguous', 'bim_index_zero_based':hits[0][0] if len(hits)==1 else '', 'bim_rsid':hits[0][1][1] if len(hits)==1 else '', 'bim_a1':hits[0][1][4] if len(hits)==1 else '', 'bim_a2':hits[0][1][5] if len(hits)==1 else '', 'alt38':x['alt38']};out.append(r)
fields=['variant_GRCh38','status','chr37','pos37','source_id','bim_index_zero_based','bim_rsid','bim_a1','bim_a2','alt38'];
with open('outputs/stage5_smad5_splice_snp_ld_reference_coverage.csv','w',newline='') as fp:w=csv.DictWriter(fp,fieldnames=fields);w.writeheader();w.writerows(out)
print(json.dumps({'input_snvs':len(ds),'status':dict(collections.Counter(x['status'] for x in out)),'reference':'MAGMA 1000G EUR phase 3 GRCh37; 503 individuals','qtl_gwas_source_file':'stage5_smad5_splice_Mak_FI_with_credible_indels_aligned.csv'},indent=2))
