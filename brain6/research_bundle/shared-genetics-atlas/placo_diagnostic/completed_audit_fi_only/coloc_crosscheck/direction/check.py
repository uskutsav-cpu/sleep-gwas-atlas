import csv,json,hashlib,math,pathlib
p=pathlib.Path(__file__).resolve().parent;ld=p.parent/'pairs.tsv';assoc=pathlib.Path('outputs/shared-genetics-atlas/frailty_coloc/matched.tsv')
with ld.open() as f: pairs=[x for x in csv.DictReader(f,delimiter='\t') if x['variant_2'] in ('rs10953765','rs1456031')]
with assoc.open() as f: rows=[x for x in csv.DictReader(f,delimiter='\t') if x['rsid'] in {z['variant_1'] for z in pairs}|{z['variant_2'] for z in pairs} and x['chr']=='7']
lookup={x['rsid']:x for x in rows};reference=json.loads((p.parent/'receipt.json').read_text())['variants'];assert len(pairs)==4 and len(rows)==4 and len(lookup)==4
out=[]
for r in pairs:
 a,b=(lookup[r['variant_1']],lookup[r['variant_2']]);assert r['status']=='OK';x={**r}
 for i,(row,count) in enumerate(((a,r['counted_allele_1']),(b,r['counted_allele_2'])),1):
  alleles=(row['effect_allele'],row['other_allele']);ref=reference[row['rsid']];assert count==ref['counted'] and set(alleles)=={ref['counted'],ref['other']} and len(set(alleles))==2
  orient=1 if alleles[0]==count else -1
  x.update({f'original_effect_allele_{i}':alleles[0],f'original_other_allele_{i}':alleles[1],f'orientation_to_counted_{i}':orient})
  for trait in ('fi','r13'):
   x.update({f'{trait}_beta_counted_{i}':orient*float(row[trait+'_beta']),f'{trait}_se_{i}':float(row[trait+'_se']),f'{trait}_p_{i}':float(row[trait+'_p'])})
 for trait in ('fi','r13'):
  product=x[f'{trait}_beta_counted_1']*x[f'{trait}_beta_counted_2'];rr=float(r['r']);x[trait+'_ld_sign_consistency']='UNDEFINED_ZERO' if product==0 or rr==0 else ('CONSISTENT' if product*rr>0 else 'OPPOSITE')
 out.append(x)
with (p/'aligned_signs.tsv').open('w',newline='') as f:w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
receipt={'ld_sha256':hashlib.sha256(ld.read_bytes()).hexdigest(),'association_sha256':hashlib.sha256(assoc.read_bytes()).hexdigest(),'reference_receipt_sha256':hashlib.sha256((p.parent/'receipt.json').read_bytes()).hexdigest(),'pairs':len(out),'fi_consistent':sum(x['fi_ld_sign_consistency']=='CONSISTENT' for x in out),'r13_consistent':sum(x['r13_ld_sign_consistency']=='CONSISTENT' for x in out),'caveat':'heuristic sign consistency only; source GWAS reused and LD reference is not study sample'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
