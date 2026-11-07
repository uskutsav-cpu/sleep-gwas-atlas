import csv,pathlib,collections,json,hashlib,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;src=p/'selected_rows.tsv';chain=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
with src.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==20 and all(x['CHROM']=='1' for x in rows)
url='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
vcf=pysam.VariantFile(url);lo=LiftOver(str(chain));result=[]
def query(pos,rs,alleles):
 records=[v for v in vcf.fetch('1',pos-1,pos) if v.pos==pos]
 matches=[v for v in records if v.id and rs in v.id.split(';') and alleles.issubset(set(v.alleles))]
 return len(records),len(matches)
for i,x in enumerate(rows):
 rs=x['ID'];pos=int(x['GENPOS']);alleles={x['ALLELE0'],x['ALLELE1']};d_records,d_match=query(pos,rs,alleles);mapping=lo.convert_coordinate('chr1',pos-1)
 item={'rsid':rs,'reported_pos':pos,'allele0':x['ALLELE0'],'allele1':x['ALLELE1'],'direct_reference_records':d_records,'direct_exact_matches':d_match,'lift_status':'UNIQUE_POSITIVE' if len(mapping)==1 and mapping[0][0]=='chr1' and mapping[0][2]=='+' else 'NOT_UNIQUE_POSITIVE'}
 if item['lift_status']=='UNIQUE_POSITIVE':
  target=mapping[0][1]+1;records,matches=query(target,rs,alleles);item.update(lifted_hg19_pos=target,lifted_reference_records=records,lifted_exact_matches=matches)
  item['classification']='BOTH' if d_match and matches else 'DIRECT_ONLY' if d_match else 'LIFT_ONLY' if matches else 'NEITHER'
 else:item['classification']='UNVERIFIABLE_LIFT'
 result.append(item);print(i+1,rs,item['classification'],flush=True)
fields=['rsid','reported_pos','allele0','allele1','direct_reference_records','direct_exact_matches','lift_status','lifted_hg19_pos','lifted_reference_records','lifted_exact_matches','classification']
with (p/'build_matches.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(result)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
(p/'reference_receipt.json').write_text(json.dumps({'selected_rows_sha256':sha(src),'chain_sha256':sha(chain),'vcf_url':url,'reference':'1000 Genomes phase 3 GRCh37; remotely indexed, whole VCF not hashed','counts':dict(collections.Counter(x['classification'] for x in result)),'output_sha256':sha(p/'build_matches.tsv')},indent=2)+'\n')
