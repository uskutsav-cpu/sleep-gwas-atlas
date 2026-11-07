import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent;root=p.parents[1];prev=root/'replication/childhood_symptoms/selected_variants.tsv';ctx=root/'variant_context/selected_variant_gene_context.tsv';files=[root/'frailty_coloc/R13_extracted.tsv',root/'frailty_coloc_chr5/R13_extracted.tsv'];receipts=[root/'frailty_coloc/R13_extraction.json',root/'frailty_coloc_chr5/R13_extraction.json']
def read(q):
 with q.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def sha(q):return hashlib.sha256(q.read_bytes()).hexdigest()
for f,r in zip(files,receipts):assert sha(f)==json.loads(r.read_text())['extraction_sha256']
sym=read(prev);assert [r['variant'] for r in sym]==['rs62264764','rs77960','rs1476535']
context={r['snp']:r for r in read(ctx) if r['snp'] in {x['variant'] for x in sym}};raw=[r for f in files for r in read(f)];res=[]
for z in sym:
 snp=z['variant'];c=context[snp];ch,pos,variant=c['grch38_ref_alt'].removeprefix('chr').split(':');ref,alt=variant.split('>');q=[r for r in raw if r['#chrom']==ch and r['pos']==pos and r['rsids']==snp and r['ref']==ref and r['alt']==alt]
 out={'variant':snp,'effect_allele':c['source_A1'],'status':'NOT_ESTIMABLE','sleep_beta_aligned':'','sleep_se':'','sleep_p':'','symptom_beta_aligned':z['symptom_beta_aligned'],'symptom_p':z['symptom_p'],'sleep_direction_agreement':'','joint_max_p':'','joint_times3':''}
 if len(q)==1 and z['status']=='MATCHED':
  r=q[0];b,se,pv=map(float,(r['beta'],r['sebeta'],r['pval']));assert all(map(math.isfinite,(b,se,pv))) and se>0 and 0<=pv<=1
  aligned=b*(1 if c['source_A1']==alt else -1);discovery=float(c['insomnia_beta_A1']);symp=float(z['symptom_p']);joint=max(pv,symp)
  out.update(status='MATCHED',sleep_beta_aligned=aligned,sleep_se=se,sleep_p=pv,sleep_direction_agreement='YES' if aligned*discovery>0 else ('NO' if aligned*discovery<0 else 'ZERO'),joint_max_p=joint,joint_times3=min(1,3*joint))
 elif len(q)>1:out['status']='AMBIGUOUS'
 res.append(out)
with (p/'comparison.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=res[0].keys(),delimiter='\t');w.writeheader();w.writerows(res)
record={'inputs_sha256':{str(q):sha(q) for q in [prev,ctx,*files,*receipts]},'n_requested':3,'n_matched':sum(r['status']=='MATCHED' for r in res),'source_identity':'FinnGen raw source previously matched official R13 object MD5; extracted table hashes match frozen extraction receipts.','note':'Post-selected retrospective comparison; hashes and max-p arithmetic do not prove independence or control selection error.'};(p/'receipt.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(res,indent=2))
