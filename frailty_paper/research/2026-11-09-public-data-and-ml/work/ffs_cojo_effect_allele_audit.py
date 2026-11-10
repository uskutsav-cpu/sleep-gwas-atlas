import csv,collections,hashlib,json,pathlib,math
file=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv')
rows=list(csv.DictReader(open('outputs/published_ffs_37_locus_inventory.csv')));wanted={r['ffs_lead_variant'] for r in rows};found=collections.defaultdict(list);n=0
with file.open() as f:
 header=next(f).rstrip('\n').split('\t'); ix={a:header.index(a) for a in ('CHROM','POS','A1','A2','BETA','SNP')}
 for line in f:
  n+=1;a=line.rstrip('\n').split('\t');id=a[ix['SNP']]
  if id in wanted:found[id].append({k:a[v] for k,v in ix.items()})
result=[]
for r in rows:
 hits=[v for v in found[r['ffs_lead_variant']] if v['CHROM']==r['chromosome'] and v['POS']==r['position_reported']]; d={'locus':r['source_locus_id'],'variant':r['ffs_lead_variant'],'source_id_matches':len(found[r['ffs_lead_variant']]),'coordinate_matches':len(hits),'cojo_refA':r['source_refA_not_yet_author_confirmed_effect_allele'],'cojo_marginal_beta':r['ffs_marginal_beta_reported']}
 if len(hits)==1:
  h=hits[0]; ref=d['cojo_refA'];raw=float(h['BETA']);b=float(d['cojo_marginal_beta']);direction=1 if ref==h['A1'] else -1 if ref==h['A2'] else 0
  d.update(source_A1=h['A1'],source_A2=h['A2'],source_beta=h['BETA'],refA_match='A1' if direction==1 else 'A2' if direction==-1 else 'neither',beta_aligned_diff=(direction*raw-b) if direction else '',beta_diff_if_A1_effect=raw-b,beta_diff_if_A2_effect=-raw-b)
 result.append(d)
out=pathlib.Path('outputs/ffs_cojo_effect_allele_audit.csv')
with out.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in result for k in r)));w.writeheader();w.writerows(result)
matched=[x for x in result if x['coordinate_matches']==1];res={'source_rows':n,'published_loci':len(result),'unique_id_coordinate_matches':len(matched),'reference_allele_matches':dict(collections.Counter(x['refA_match'] for x in matched)),'beta_refA_aligned_within_1e-6':sum(x['beta_aligned_diff']!='' and abs(x['beta_aligned_diff'])<=1e-6 for x in matched),'beta_if_a1_effect_within_1e-6':sum(abs(x['beta_diff_if_A1_effect'])<=1e-6 for x in matched),'beta_if_a2_effect_within_1e-6':sum(abs(x['beta_diff_if_A2_effect'])<=1e-6 for x in matched),'unmatched':[x['variant'] for x in result if x['coordinate_matches']!=1],'source_sha256_from_prior_audit':'912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55','output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()};pathlib.Path('outputs/ffs_cojo_effect_allele_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
