import csv,json,pathlib,hashlib,sys,pysam,collections
from pyliftover import LiftOver
p=pathlib.Path(__file__).resolve().parent;acc=sys.argv[1];assert acc in [f'GCST{x}' for x in range(90624046,90624053)]
chain=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/hg19ToHg38.over.chain.gz');selector=p.parent/'selected.tsv'
with selector.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==87
lift=LiftOver(str(chain));todo=[]
for x in rows:
 y=lift.convert_coordinate('chr'+x['chrom'],int(x['pos37'])-1)
 if len(y)!=1 or not y[0][0].startswith('chr') or y[0][0][3:] not in [str(n) for n in range(1,23)]:
  x.update(lift_status='FAILED_OR_AMBIGUOUS',chr38='',pos38='',alleles38='');continue
 ch=y[0][0][3:];pos=y[0][1]+1;strand=y[0][2];assert strand in ('+','-')
 conv=lambda a:a if strand=='+' else a.translate(str.maketrans('ACGT','TGCA'))
 x.update(lift_status='UNIQUE_'+strand,chr38=ch,pos38=str(pos),alleles38=','.join(sorted((conv(x['allele_1']),conv(x['allele_2'])))));todo.append(x)
groups=[]
for x in sorted(todo,key=lambda x:(int(x['chr38']),int(x['pos38']))):
 pos=int(x['pos38']);ch=x['chr38']
 if not groups or groups[-1]['chrom']!=ch or pos-groups[-1]['end']>500000:groups.append({'chrom':ch,'start':pos,'end':pos,'targets':[x]})
 else:groups[-1]['end']=pos;groups[-1]['targets'].append(x)
base=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/harmonised/{acc}.h.tsv.gz';tab=pysam.TabixFile(base);fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','rs_id','n','Q_metric','Q_metric_pvalue','hm_coordinate_conversion','hm_code','variant_id'];found=collections.defaultdict(list)
for g in groups:
 ids={x['rsid'] for x in g['targets']}
 for raw in tab.fetch(g['chrom'],g['start']-1,g['end']):
  vals=raw.split('\t');assert len(vals)==len(fields);v=dict(zip(fields,vals))
  for rs in ids & {v['rsid'],v['rs_id']}:found[rs].append(v)
out=[]
for x in rows:
 cand=found[x['rsid']];exact=[r for r in cand if r['chromosome']==x['chr38'] and r['base_pair_location']==x['pos38'] and {r['effect_allele'],r['other_allele']}==set(x['alleles38'].split(','))]
 status='LIFT_FAILED' if not x['chr38'] else 'UNIQUE_MATCH' if len(exact)==1 else 'MULTIPLE_EXACT' if len(exact)>1 else 'ID_PRESENT_MISMATCH' if cand else 'ABSENT'
 row={'accession':acc,'selected_rsid':x['rsid'],'selection':x['selection'],'lift_status':x['lift_status'],'chr38':x['chr38'],'pos38':x['pos38'],'alleles38':x['alleles38'],'status':status,'n_id_rows':len(cand),'n_exact':len(exact)};row.update({k:exact[0][k] if len(exact)==1 else '' for k in fields});out.append(row)
with (p/f'{acc}_extended.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
def sh(x):return hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'source':base,'index':base+'.tbi','chain_sha256':sh(chain),'selector_sha256':sh(selector),'script_sha256':sh(pathlib.Path(__file__)),'harmonized_metadata_sha256':sh(p/f'{acc}.h.tsv.gz-meta.yaml'),'full_file_checksum_verified':False,'unique_lifts':len(todo),'regions':[{k:v for k,v in g.items() if k!='targets'} for g in groups],'status_counts':dict(collections.Counter(x['status'] for x in out))}
(p/f'{acc}_extended_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(acc,receipt['status_counts'],'regions',len(groups),flush=True)
