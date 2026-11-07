import csv,json,pathlib,hashlib,subprocess,sys,pysam
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;acc=sys.argv[1];assert acc in [f'GCST{x}' for x in range(90624046,90624053)]
chain=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/hg19ToHg38.over.chain.gz');selector=p.parent/'selected.tsv';chosen={'rs62264764','rs817509','rs2431108','rs7703746','rs1476535'}
with selector.open() as f: rows=[x for x in csv.DictReader(f,delimiter='\t') if x['rsid'] in chosen]
assert len(rows)==5
base=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/harmonised/{acc}.h.tsv.gz'
metadata=subprocess.check_output(['curl','-fsSL','--max-time','45',base+'-meta.yaml']);(p/f'{acc}.h.tsv.gz-meta.yaml').write_bytes(metadata)
lift=LiftOver(str(chain));tab=pysam.TabixFile(base);fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','rs_id','n','Q_metric','Q_metric_pvalue','hm_coordinate_conversion','hm_code','variant_id'];out=[]
for x in rows:
 hits=lift.convert_coordinate('chr'+x['chrom'],int(x['pos37'])-1);assert len(hits)==1 and hits[0][0]=='chr'+x['chrom'] and hits[0][2]=='+',(x['rsid'],hits)
 pos=hits[0][1]+1;records=[]
 for raw in tab.fetch(x['chrom'],pos-1,pos):
  vals=raw.split('\t');assert len(vals)==len(fields);r=dict(zip(fields,vals))
  if r['base_pair_location']==str(pos) and x['rsid'] in (r['rsid'],r['rs_id']):records.append(r)
 exact=[r for r in records if {r['effect_allele'],r['other_allele']}=={x['allele_1'],x['allele_2']}]
 status='UNIQUE_MATCH' if len(exact)==1 else 'MULTIPLE_EXACT' if len(exact)>1 else 'ID_PRESENT_MISMATCH' if records else 'ABSENT'
 a={'accession':acc,'rsid_requested':x['rsid'],'pos37':x['pos37'],'pos38':pos,'n_id_records':len(records),'n_exact':len(exact),'status':status};a.update({k:exact[0][k] if len(exact)==1 else '' for k in fields});out.append(a)
with (p/f'{acc}_anchors.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
def sh(x):return hashlib.sha256(x.read_bytes()).hexdigest()
(p/f'{acc}_receipt.json').write_text(json.dumps({'source':base,'index':base+'.tbi','full_file_checksum_verified':False,'metadata_sha256':hashlib.sha256(metadata).hexdigest(),'chain_sha256':sh(chain),'selector_sha256':sh(selector),'script_sha256':sh(pathlib.Path(__file__)),'selection':sorted(chosen),'records':len(out),'limitation':'Remote tabix access is not a full-file checksum; absent record is not evidence of no association or no published locus.'},indent=2)+'\n')
print(acc,[(x['rsid_requested'],x['status']) for x in out],flush=True)
