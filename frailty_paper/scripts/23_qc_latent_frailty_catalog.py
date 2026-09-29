#!/usr/bin/env python3
"""Bounded-memory ordered-coordinate and row-N QC of latent GWAS Catalog TSVs.

The existing acquired-resource manifest supplies file hashes. Numeric coordinate order and adjacent duplicate keys are diagnostic only; duplicate detection is not global when a file is unsorted. This scan reports findings only and does not alter source data.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
FIELDS=['accession','file','rows','coordinate_order','coordinate_order_descents','adjacent_exact_variant_key_duplicates','invalid_chromosomes','invalid_positions','invalid_alleles','n_missing_or_invalid','n_distinct_values','n_min','n_max','n_values_outside_metadata','interpretation','examples']
CHR={str(i):i for i in range(1,23)}|{'X':23,'Y':24,'M':25,'MT':25}
ALLELES=set('ACGTN')
def scan(path, accession, expected):
 c={k:0 for k in ('coordinate_order_descents','adjacent_exact_variant_key_duplicates','invalid_chromosomes','invalid_positions','invalid_alleles','n_missing_or_invalid','n_values_outside_metadata')}; c['rows']=0
 prev_coord=prev_variant=None; nvals=set(); nmin=nmax=None; ex=[]
 with path.open(newline='',encoding='utf-8-sig') as f:
  r=csv.reader(f,delimiter='\t'); header=next(r); ix={x:i for i,x in enumerate(header)}
  needed=['chromosome','base_pair_location','effect_allele','other_allele','n']
  if any(x not in ix for x in needed): raise ValueError(f'{accession}: missing columns')
  ci,pi,ei,oi,ni=(ix[x] for x in needed); ri=ix.get('rs_id',-1)
  for row in r:
   c['rows']+=1; ch=row[ci].strip().upper().removeprefix('CHR'); ps=row[pi].strip(); ea=row[ei].strip().upper(); oa=row[oi].strip().upper(); issues=[]
   rank=CHR.get(ch,1000)
   if rank==1000: c['invalid_chromosomes']+=1;issues.append('chromosome')
   try: pos=int(ps); assert pos>0
   except (ValueError,AssertionError): pos=-1;c['invalid_positions']+=1;issues.append('position')
   coord=(rank,pos); variant=(ch,ps,ea,oa)
   if prev_coord is not None and coord<prev_coord:c['coordinate_order_descents']+=1;issues.append('order')
   if variant==prev_variant:c['adjacent_exact_variant_key_duplicates']+=1;issues.append('adjacent_exact_duplicate')
   prev_coord,prev_variant=coord,variant
   if not ea or not oa or (set(ea)|set(oa))-ALLELES:c['invalid_alleles']+=1;issues.append('allele')
   try:
    n=int(row[ni]); assert n>0
    nvals.add(n);nmin=n if nmin is None else min(nmin,n);nmax=n if nmax is None else max(nmax,n)
    if n!=expected:c['n_values_outside_metadata']+=1
   except (ValueError,AssertionError):c['n_missing_or_invalid']+=1;issues.append('n')
   if issues and len(ex)<10:ex.append({'row':c['rows'],'chromosome':ch,'position':ps,'rs_id':row[ri] if ri>=0 else '', 'issues':issues})
 return {'accession':accession,'file':str(path),'rows':c['rows'],'coordinate_order':'NO_NUMERIC_DESCENTS' if not c['coordinate_order_descents'] else 'NUMERIC_DESCENTS_OBSERVED','coordinate_order_descents':c['coordinate_order_descents'],'adjacent_exact_variant_key_duplicates':c['adjacent_exact_variant_key_duplicates'],'invalid_chromosomes':c['invalid_chromosomes'],'invalid_positions':c['invalid_positions'],'invalid_alleles':c['invalid_alleles'],'n_missing_or_invalid':c['n_missing_or_invalid'],'n_distinct_values':len(nvals),'n_min':nmin or '','n_max':nmax or '','n_values_outside_metadata':c['n_values_outside_metadata'],'interpretation':'STRUCTURAL_FIELDS_VALID; DUPLICATE_ASSESSMENT_ADJACENT_ONLY' if not any(c[k] for k in ('invalid_chromosomes','invalid_positions','invalid_alleles')) else 'STRUCTURAL_FIELD_ANOMALIES_REPORTED_NO_FILTERING','examples':json.dumps(ex,separators=(',',':'))}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--accession',action='append');ap.add_argument('--output',type=Path,default=Path('frailty_paper/manifests/gwas_latent_deep_qc.tsv'));a=ap.parse_args();repo=a.repo.resolve()
 idx=list(csv.DictReader((repo/'frailty_paper/manifests/latent_frailty_accession_index.tsv').open(),delimiter='\t')); wanted=set(a.accession or [x['accession'] for x in idx]); out=[]
 for m in idx:
  if m['accession'] not in wanted:continue
  p=repo/'frailty_paper/data/gwas/latent_frailty_catalog'/m['accession']/f"{m['accession']}.tsv"; row=scan(p,m['accession'],int(m['sample_size'])) if p.is_file() else {'accession':m['accession'],'file':str(p),'interpretation':'MISSING'};row['file']=p.relative_to(repo).as_posix();out.append(row)
  print(f"{row['accession']} rows={row['rows']} coordinate_order={row['coordinate_order']} adjacent_duplicate_keys={row['adjacent_exact_variant_key_duplicates']} N={row.get('n_min')}..{row.get('n_max')} metadata_N={m['sample_size']} rows_N_mismatch={row['n_values_outside_metadata']}",flush=True)
 dest=a.output if a.output.is_absolute() else repo/a.output; dest.parent.mkdir(parents=True,exist_ok=True)
 write_header=not a.accession or not dest.exists() or not dest.stat().st_size
 with dest.open('a' if a.accession else 'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n',extrasaction='ignore')
  if write_header:w.writeheader()
  w.writerows(out)
if __name__=='__main__':main()
