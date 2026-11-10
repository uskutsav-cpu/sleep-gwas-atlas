import csv, json, subprocess, hashlib, pathlib
base=pathlib.Path('outputs'); table=list(csv.DictReader((base/'shared22_eight_source_evidence_matrix.csv').open())); rows=[r for r in table if r['all_eight_nominal_same_sign']=='True']; assert len(rows)==5
res=[]; rawpath=pathlib.Path('work/vep_shared_five_responses');rawpath.mkdir(exist_ok=True)
for r in rows:
 rs=r['variant']; assert rs.startswith('rs') and rs[2:].isdigit();url='https://grch37.rest.ensembl.org/vep/human/id/'+rs+'?content-type=application/json'
 d=subprocess.run(['curl','-fsSL','--max-time','30',url],capture_output=True,check=True).stdout; (rawpath/(rs+'.json')).write_bytes(d); hits=json.loads(d)
 matches=[x for x in hits if str(x['seq_region_name'])==r['chr_GRCh37_assumed'] and str(x['start'])==r['position_GRCh37_assumed']]; assert len(matches)==1,(rs,len(matches))
 h=matches[0]; alls=h['allele_string'].split('/'); ref=alls[0]; pair={r['score_counted_allele'],r['other_allele']};assert ref in pair and len(pair)==2; alt=next(iter(pair-{ref})); assert alt in alls
 tc=[t for t in h.get('transcript_consequences',[]) if t.get('variant_allele')==alt]
 terms=sorted({s for t in tc for s in t.get('consequence_terms',[])})
 genes=sorted({t.get('gene_symbol') or t.get('gene_id') for t in tc if t.get('gene_symbol') or t.get('gene_id')})
 regulatory=sorted({s for t in h.get('regulatory_feature_consequences',[]) if t.get('variant_allele')==alt for s in t.get('consequence_terms',[])})
 res.append(dict(variant=rs,GRCh37_chr_assumed=r['chr_GRCh37_assumed'],GRCh37_pos_assumed=r['position_GRCh37_assumed'],shared_count_allele=r['score_counted_allele'],other_allele=r['other_allele'],nearest_gene_body=r['nearest_gene_body_positional'],reference_allele=ref,annotated_alternate_allele=alt,vep_alleles=h['allele_string'],alternate_transcript_terms=';'.join(terms),alternate_transcript_genes=';'.join(genes),alternate_regulatory_terms=';'.join(regulatory),n_transcript_records=len(tc),source_url=url,response_sha256=hashlib.sha256(d).hexdigest()))
 print(rs,alt,terms,genes,flush=True)
out=base/'shared5_allele_specific_vep.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=res[0].keys());w.writeheader();w.writerows(res)
print('output SHA256',hashlib.sha256(out.read_bytes()).hexdigest())
