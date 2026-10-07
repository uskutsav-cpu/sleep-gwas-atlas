#!/usr/bin/env python3
import csv,json,pathlib,sys
import numpy as np
from reconstruct_signed_ld import regional_rows,header_and_samples,scalar_pearson,write,sha
HERE=pathlib.Path(__file__).resolve().parent;RAW=HERE.parents[3]/'work/ld_genotypes_research_v1'
def main():
 m=json.loads((HERE/'additional_prior_lead_acquisition_manifest.json').read_text())
 if m['protocol_sha256']!=sha(HERE/'ADDITIONAL_PRIOR_LEADS_PROTOCOL.md'):raise ValueError('Protocol changed')
 plans=json.loads((HERE/'genotype_acquisition_manifest.json').read_text())['regions'];panel={r['sample']:r for r in csv.DictReader((RAW/'integrated_call_samples_v3.20130502.ALL.panel').open(),delimiter='\t')}
 lead_specs={'rs171697':(5,103956516,{'C','G'}),'rs30266':(5,103972357,{'A','G'}),'rs7105462':(11,112912048,None),'rs9485410':(6,101261007,None)}
 results=[];extra={}
 for obj in plans:
  c=obj['chr'];head,samples=header_and_samples(next(RAW.glob(f'chr{c}_header_*.bgzf')));eur=[s for s in samples if panel[s]['super_pop']=='EUR'];ix=[samples.index(s)+9 for s in eur]
  paths=sorted(RAW.glob(f'chr{c}_genotypes_*.bgzf'),key=lambda p:int(p.name.split('_')[2]))
  for row in regional_rows(paths,c,obj['start'],obj['stop']):
   pos=int(row[1]);wanted=[(rs,spec) for rs,spec in lead_specs.items() if spec[:2]==(c,pos)]
   if not wanted:continue
   ref,alt=row[3:5]
   for rs,spec in wanted:
    expected=spec[2];mapping=json.loads((RAW/(rs+'_grch37.json')).read_text());maps=[x for x in mapping['mappings'] if x['seq_region_name']==str(c) and x['assembly_name']=='GRCh37' and x['strand']==1 and x['start']==pos]
    if len(maps)!=1:raise ValueError('Mapping mismatch')
    if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':continue
    if not {ref,alt}.issubset(set(maps[0]['allele_string'].split('/'))):continue
    if expected and {ref,alt}!=expected:continue
    tokens=[row[i] for i in ix]
    if not all(x in ['0|0','0|1','1|0','1|1','0/0','0/1','1/0','1/1'] for x in tokens):raise ValueError('Invalid diploidGT')
    gt=np.array([int(x[0])+int(x[2]) for x in tokens]);key=f'{c}:{pos}:{ref}:{alt}'
    results.append({'published_lead':rs,'variant_key':key,'n_reference':len(eur),'counted_allele':alt,'ALT_AF':float(gt.mean()/2),'primary_source_allele_pair_available':bool(expected),'genotype_status':'EXACT_G37_COORDINATE_AND_KNOWN_SOURCE_ALLELES' if expected else 'G37_COORDINATE_AND_DBSNP_ALLELES_SOURCE_EFFECT_ALLELES_UNRESOLVED','primary_finemap_inclusion':'PALINDROMIC_EXCLUDED' if {ref,alt} in [{'A','T'},{'C','G'}] else 'NONPAL_REFERENCE_AVAILABLE'})
    if expected:extra[rs]=(gt,key,alt,eur)
 write(HERE/'additional_prior_lead_availability.tsv',results)
 factors=np.load(HERE/'chr5_genotype_factors.npz');v=list(csv.DictReader((HERE/'chr5_variants.tsv').open(),delimiter='\t'));by={r['mapped_lead_rsid']:i for i,r in enumerate(v) if r['mapped_lead_rsid']};comparisons=[]
 for rs,(g,key,alt,eur) in extra.items():
  for brain6 in ['rs2431108','rs77960']:
   b=by[brain6];g2=factors['ALT_dosage'][:,b];assert list(factors['samples'])==eur
   for pop in ['EUR','CEU','FIN','GBR','IBS','TSI']:
    ix=[i for i,s in enumerate(eur) if pop=='EUR' or panel[s]['pop']==pop];r=scalar_pearson(g[ix],g2[ix]);comparisons.append({'published_lead':rs,'published_variant_key':key,'published_counted_allele':alt,'comparison_lead':brain6,'comparison_key':v[b]['variant_key'],'comparison_counted_allele':v[b]['ALT'],'population':pop,'n_reference':len(ix),'signed_r':r,'r_squared':r*r,'scope':'RETROSPECTIVE_REFERENCE_LD_NOT_CONDITIONAL_OR_REPLICATION_GWAS','published_palin_omitted_primary_fit':rs=='rs171697'})
 write(HERE/'additional_prior_vs_brain6_LD.tsv',comparisons)
 print(json.dumps({'available_markers':len(results),'LD_population_comparisons':len(comparisons),'EUR':[x for x in comparisons if x['population']=='EUR']}),flush=True)
if __name__=='__main__':main()
