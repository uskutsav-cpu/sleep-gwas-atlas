import csv,pathlib
p=pathlib.Path(__file__).resolve().parent;source=p.parents[2]/'placo_diagnostic/completed_audit/candidate_positions.tsv'
# p.parents[2] is the atlas root.
with source.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
variants={}
def add(rs,ch,pos,a,b,origin):
 key=(str(ch),str(pos),frozenset((a,b)))
 if rs in variants:assert variants[rs]['key']==key,(rs,variants[rs]['key'],key)
 else:variants[rs]={'key':key,'origin':set()}
 variants[rs]['origin'].add(origin)
for x in rows:add(x['SNP'],x['CHR'],x['BP_GRCH37'],x['A1'],x['A2'],x['pair']+'_PLACO_selected')
for v in [('rs62264764',3,117639575,'G','A'),('rs817509',3,117750861,'G','A'),('rs2431108',5,103947968,'T','C'),('rs7703746',5,103907337,'G','A'),('rs1476535',7,114071035,'C','T')]:add(*v,'prior_anchor')
with (p/'selected.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['rsid','chrom','pos37','allele_1','allele_2','selection']);
 for rs,x in sorted(variants.items()):
  ch,pos,alleles=x['key'];a,b=sorted(alleles);w.writerow([rs,ch,pos,a,b,';'.join(sorted(x['origin']))])
print('unique selected variants',len(variants))
