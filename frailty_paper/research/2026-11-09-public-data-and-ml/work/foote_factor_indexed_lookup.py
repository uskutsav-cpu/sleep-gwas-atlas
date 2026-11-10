import sys,csv,json,hashlib,pathlib,pysam,math
acc=sys.argv[1];u=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/{acc}/harmonised/{acc}.h.tsv.gz';idx=pathlib.Path('work')/(acc+'.h.tsv.gz.tbi');f=pysam.TabixFile(u,index=str(idx));rows=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))); assert len(rows)==22
fields='chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rsid rs_id n Q_metric Q_metric_pvalue hm_coordinate_conversion hm_code variant_id'.split();out=[]
for r in rows:
 ch=r['chr_GRCh38']; pos=int(r['position_GRCh38']);a=r['score_counted_allele'];b=r['other_allele'];candidates=[]
 for line in f.fetch(ch,pos-1,pos):
  vals=line.split('\t');assert len(vals)==len(fields),(acc,len(vals));z=dict(zip(fields,vals))
  if z['base_pair_location']==str(pos) and z['rsid']==r['variant'] and {z['effect_allele'],z['other_allele']}=={a,b}:candidates.append(z)
 row={'variant':r['variant'],'position_GRCh38':pos,'score_allele':a,'other_allele':b,'matched_rows':len(candidates)}
 if len(candidates)==1:
  z=candidates[0];beta=float(z['beta'])*(1 if z['effect_allele']==a else -1);row.update(beta=beta,se=float(z['standard_error']),p=float(z['p_value']),expected_n=z['n'],q_metric_p=z['Q_metric_pvalue'],harmonisation_code=z['hm_code'])
 out.append(row);print(acc,r['variant'],len(candidates),flush=True)
result={'accession':acc,'url':u,'index_md5':hashlib.md5(idx.read_bytes()).hexdigest(),'factor_candidate_rows':out,'limitations':'Remote indexed subset only; full data checksum NOT verified; Q-filtered, absence not null; positive strand presumed from earlier candidate liftover; harmonised release GRCh38, no individual data; factor N is Nhat.'}
pathlib.Path('outputs/foote_'+acc+'_indexed_candidate_lookup.json').write_text(json.dumps(result,indent=2)+'\n');print('complete',acc,flush=True)
