import csv,json,pathlib,math,hashlib,collections,pysam
p=pathlib.Path(__file__).resolve().parent;root=p.parent;panel=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/1kg_phase3_panel.tsv');eur={r['sample'] for r in csv.DictReader(panel.open(),delimiter='\t') if r['super_pop']=='EUR'};assert len(eur)==503
source=root.parent/'brain6-alphagenome-pilot/results/jbrowse_expanded/expanded_variant_evidence.tsv';allrows=list(csv.DictReader(source.open(),delimiter='\t'));regs={'placo_insomnia__adhd_chr3_117142005_118142005':('rs817509','G','A'),'placo_insomnia__adhd_chr7_113571035_114571035':('rs12536335','G','A')};base='https://hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz';out=[];receipt=[]
for reg,(anchor,anchor_a,anchor_b) in regs.items():
 rows=[x for x in allrows if x['candidate_locus_id']==reg];selected={};
 for col,label in [('insomnia_conditional_abf_pip','top4_insomnia_pip'),('adhd_conditional_abf_pip','top4_adhd_pip')]:
  for r in sorted(rows,key=lambda x:float(x[col]),reverse=True)[:4]:selected.setdefault(r['snp'],{'row':r,'why':[]})['why'].append(label)
 for col,label,fun in [('dnase_log2_ratio_501bp','max_abs_dnase',lambda r:abs(float(r['dnase_log2_ratio_501bp']))),('avi_phred_lower','max_avi',lambda r:float(r['avi_phred_lower']))]:
  r=max(rows,key=fun);selected.setdefault(r['snp'],{'row':r,'why':[]})['why'].append(label)
 r=next(x for x in rows if x['snp']==anchor);selected.setdefault(anchor,{'row':r,'why':[]})['why'].append('frailty_anchor')
 chrom=next(iter(selected.values()))['row']['chromosome'].removeprefix('chr');vcf=pysam.VariantFile(base.format(chrom=chrom));calls={}
 for rs,entry in selected.items():
  r=entry['row'];pos=int(r['pos_grch37_1based']);counted=r['source_a1'];other=r['source_a2']
  if rs==anchor:counted,other=anchor_a,anchor_b
  note={'snp':rs,'GRCh37':chrom+':'+str(pos),'allele_counted':counted,'other_allele':other,'status':'UNKNOWN'}
  try:
   records=[x for x in vcf.fetch(chrom,pos-1,pos) if x.pos==pos and x.id==rs and counted in x.alleles and other in x.alleles];assert len(records)==1,(rs,len(records));rec=records[0];a=rec.alleles.index(counted);b=rec.alleles.index(other);v={}
   for sample in sorted(eur.intersection(vcf.header.samples)):
    c=rec.samples[sample];g=c.get('GT')
    if g and len(g)==2 and c.phased and all(t in (a,b) for t in g):v[sample]=tuple(int(t==a) for t in g)
   calls[rs]=v;note.update({'status':'OK','matched_reference_alleles':list(rec.alleles),'n_eligible':len(v)})
  except Exception as e:note['status']='FAILED '+str(e)
  receipt.append(note)
 if anchor not in calls:continue
 for rs,e in selected.items():
  r=e['row'];z={'region':reg,'frailty_anchor':anchor,'anchor_counted_allele':anchor_a,'brain_variant':rs,'brain_counted_allele':r['source_a1'] if rs!=anchor else anchor_a,'grch37_position':r['pos_grch37_1based'],'selection':';'.join(e['why']),'insomnia_PIP':r['insomnia_conditional_abf_pip'],'ADHD_PIP':r['adhd_conditional_abf_pip'],'insomnia_credible_member':r['insomnia_credible_set_member'],'ADHD_credible_member':r['adhd_credible_set_member'],'dnase_alt_ref_log2ratio':r['dnase_log2_ratio_501bp'],'avi_phred':r['avi_phred_lower']}
  if rs in calls:
   samples=sorted(set(calls[anchor])&set(calls[rs]));x=[];y=[]
   for s in samples:x.extend(calls[anchor][s]);y.extend(calls[rs][s])
   n=len(x);px=sum(x)/n;py=sum(y)/n;pxy=sum(a*b for a,b in zip(x,y))/n;den=math.sqrt(px*(1-px)*py*(1-py));rr=(pxy-px*py)/den if den else None;z.update({'status':'OK','n_people':len(samples),'r':rr,'r2':rr*rr if rr is not None else None})
  else:z['status']='FAILED_REFERENCE_MATCH'
  out.append(z)
with (p/'credible_representative_ld.tsv').open('w',newline='') as f:
 cols=list(dict.fromkeys(k for r in out for k in r));w=csv.DictWriter(f,cols,delimiter='\t');w.writeheader();w.writerows(out)
(p/'credible_representative_receipt.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(panel.read_bytes()).hexdigest(),'n_source_rows':len(allrows),'selection':'top 4 marginal ABF PIPs per trait, top absolute DNase, top AVI and frailty anchor in each region; deduplicated','reference_records':receipt,'vcf_urls':[base.format(chrom='3'),base.format(chrom='7')],'note':'post-selected exact allele matches; remote indexed VCF not fully hashed; no individual calls saved'},indent=2)+'\n')
print('rows',len(out),'OK',sum(x['status']=='OK' for x in out),'failed',sum(x['status']!='OK' for x in out),flush=True)
