import csv,json,math,collections
raw=list(csv.DictReader(open('work/raw_external_prototype_matches.csv'))); markers={x['variant']:x for x in csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv'))}; coords={x['variant']:x for x in csv.DictReader(open('work/grch37_to_grch38_mapping.csv'))}; ref={x['variant']:x for x in csv.DictReader(open('work/variant_reference_grch38.csv'))}; out=[]
for x in raw:
 m=markers[x['variant']];co=coords[x['variant']]; re=ref[x['variant']];mapping=json.loads(re['Ensembl_GRCh38_mappings']) if re['status']=='success' else []
 a=(m['counted_allele'],m['other_allele']);b=(x['A1'],x['A2']);mode='direct' if a==b else 'swap' if a==b[::-1] else 'incompatible'
 same_coord=(co['status']=='unique_mapping' and co['mapping_strand']=='1' and x['chr']==co['grch38_chr'] and x['position']==co['grch38_pos']);position_ref=any(str(q['chr'])==x['chr'] and str(q['pos'])==x['position'] for q in mapping);pair_ref=any(all(a0 in (q.get('alleles') or '').split('/') for a0 in a) for q in mapping)
 if not (same_coord and position_ref and pair_ref and mode in ('direct','swap') and set(a) not in ({'A','T'},{'C','G'})):status='do_not_interpret';aligned='';concordance='not_assessed'
 else:
  value=float(x['effect_value']);value=math.log(value) if x['source']=='MVP' else value
  aligned=value*(1 if mode=='direct' else -1); prod=aligned*float(m['insomnia_beta_weight']);concordance='same' if prod>0 else 'opposite' if prod<0 else 'zero';status='coordinate_allele_match_provisional'
  if x['source']=='MVP' and float(x['effect_value'])==1.0:concordance='zero_from_rounded_OR'
 out.append({'variant':x['variant'],'source':x['source'],'source_trait':'insomnia PheCode327.4' if x['source']=='MVP' else 'F5_INSOMNIA','chr_GRCh37':m['chr'],'pos_GRCh37':m['position_GRCh37_inferred'],'chr_GRCh38':x['chr'],'pos_GRCh38':x['position'],'GRCh37_to_38_map_strand':co['mapping_strand'],'candidate_effect_allele':a[0],'candidate_other_allele':a[1],'source_effect_allele':b[0],'source_other_allele':b[1],'allele_relation':mode,'Ensembl_mapping_and_alleles_match':str(same_coord and position_ref and pair_ref),'source_effect':x['effect_value'],'source_effect_type':x['effect_type'],'aligned_beta_to_candidate_allele':aligned,'Jansen_insomnia_beta':m['insomnia_beta_weight'],'sign_concordance':concordance,'source_p':x['p'],'source_SE':x['SE'],'source_EAF':x['frequency'],'variant_has_more_than_two_reference_alleles':m['GRCh37_reference_multiallelic'],'status':status,'limitations':'Descriptive post-selection lookup; sample independence unknown; source QC not rerun; multiallelic study contrast unconfirmed; not independent frailty validation or predictive risk score'})
assert len(out)==30 and len({(x['variant'],x['source']) for x in out})==30
with open('outputs/prototype_15_variants_full_external_insomnia_lookup.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
for s in ('MVP','R13'):
 a=[x for x in out if x['source']==s];print(s,'matched',sum(x['status']=='coordinate_allele_match_provisional' for x in a),'signs',dict(collections.Counter(x['sign_concordance'] for x in a)), 'raw p<0.05',sum(float(x['source_p'])<.05 for x in a),'raw p<5e-8',sum(float(x['source_p'])<5e-8 for x in a))
