import csv,pathlib,sys,re
trait=sys.argv[1]
base=pathlib.Path('work/stage4_magma')
log=(base/(trait+'.log')).read_text();assert 'End time is' in log
expected=re.findall(r'found (\d+) genes containing valid SNPs',log);assert len(expected)==1
symbols={r['GENE']:r['gene_symbol_v92'] for r in csv.DictReader(open('outputs/stage4_mak_fi_brown_magma_gene_results.csv'))}
with (base/(trait+'.genes.out')).open() as f:
 head=f.readline().split(); rows=[dict(zip(head,line.split())) for line in f if line.strip()]
assert len(rows)==int(expected[0]) and all(len(r)==len(head) and 0<=float(r['P'])<=1 for r in rows)
order=sorted(range(len(rows)),key=lambda i:float(rows[i]['P'])); q=1.
for rank in range(len(order),0,-1):
 i=order[rank-1];q=min(q,float(rows[i]['P'])*len(rows)/rank);rows[i]['bh_q_all_genes']=q
for row in rows:row['gene_symbol_v92']=symbols.get(row['GENE'],'')
out=pathlib.Path('outputs')/('stage4_'+trait+'_magma_gene_results.csv')
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=head+['gene_symbol_v92','bh_q_all_genes']);w.writeheader();w.writerows(rows)
print(trait,len(rows), 'BH<.05',sum(r['bh_q_all_genes']<.05 for r in rows),'symbols',sum(bool(r['gene_symbol_v92']) for r in rows))
