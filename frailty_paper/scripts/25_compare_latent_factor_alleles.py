#!/usr/bin/env python3
"""Compare allele orientation across the seven latent frailty factor GWAS.

Inputs are required to be coordinate-sorted (checked separately); this is a
read-only source QC report and does not harmonize or change effects.
"""
from __future__ import annotations
import argparse,csv,heapq,itertools
from collections import defaultdict
from itertools import groupby
from pathlib import Path
from typing import Iterator
FIELDS=['accession_a','accession_b','shared_coordinates','single_variant_coordinate_comparisons','same_effect_allele','swapped_effect_allele','complemented_same','complemented_swapped','palindromic_coordinate_comparisons','allele_mismatches','multi_variant_coordinates','interpretation']
CHROM={str(i):i for i in range(1,23)}|{'X':23,'Y':24,'M':25,'MT':25}
COMP={'A':'T','T':'A','C':'G','G':'C'}
def rows(path:Path, index:int)->Iterator[tuple]:
 with path.open(newline='',encoding='utf-8-sig') as f:
  r=csv.DictReader(f,delimiter='\t')
  for row in r:
   ch=(row.get('chromosome') or '').strip().upper().removeprefix('CHR'); pos=int(row['base_pair_location'])
   yield (CHROM.get(ch,1000),pos,index,ch,(row.get('rs_id') or '').strip(),(row.get('effect_allele') or '').strip().upper(),(row.get('other_allele') or '').strip().upper())
def orient(x,y):
 _,_,_,_,_,ea1,oa1=x; _,_,_,_,_,ea2,oa2=y
 if (ea1,oa1)==(ea2,oa2): return 'same_effect_allele'
 if (ea1,oa1)==(oa2,ea2): return 'swapped_effect_allele'
 if len(ea1)==len(oa1)==len(ea2)==len(oa2)==1 and all(a in COMP for a in (ea1,oa1,ea2,oa2)):
  c1,c2=COMP[ea1],COMP[oa1]
  if (c1,c2)==(ea2,oa2): return 'complemented_same'
  if (c1,c2)==(oa2,ea2): return 'complemented_swapped'
 return 'allele_mismatches'
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[2]); ap.add_argument('--output',type=Path,default=Path('frailty_paper/manifests/gwas_latent_allele_compatibility.tsv'));a=ap.parse_args();repo=a.repo.resolve()
 idx=[r for r in csv.DictReader((repo/'frailty_paper/manifests/latent_frailty_accession_index.tsv').open(),delimiter='\t') if not r['accession'].endswith('053')]
 streams=[]
 for i,m in enumerate(idx):
  p=repo/'frailty_paper/data/gwas/latent_frailty_catalog'/m['accession']/f"{m['accession']}.tsv"
  streams.append(rows(p,i))
 merged=heapq.merge(*streams,key=lambda r:(r[0],r[1]))
 counts={(i,j):defaultdict(int) for i in range(len(idx)) for j in range(i+1,len(idx))}
 for coord,events in groupby(merged,key=lambda r:(r[0],r[1])):
  by={i:[] for i in range(len(idx))}
  for event in events: by[event[2]].append(event)
  present=[i for i,v in by.items() if v]
  for i,j in itertools.combinations(present,2):
   c=counts[(i,j)]; c['shared_coordinates']+=1
   if len(by[i])==len(by[j])==1:
    x,y=by[i][0],by[j][0]; c['single_variant_coordinate_comparisons']+=1
    c[orient(x,y)]+=1
    if len(x[5])==len(x[6])==len(y[5])==len(y[6])==1 and x[5]+x[6] in {'AT','TA','CG','GC'}: c['palindromic_coordinate_comparisons']+=1
   else:c['multi_variant_coordinates']+=1
 out=[]
 for (i,j),c in counts.items():
  row={'accession_a':idx[i]['accession'],'accession_b':idx[j]['accession'],**{k:c[k] for k in ('shared_coordinates','single_variant_coordinate_comparisons','same_effect_allele','swapped_effect_allele','complemented_same','complemented_swapped','palindromic_coordinate_comparisons','allele_mismatches','multi_variant_coordinates')}}
  row['interpretation']='NO_SINGLETON_COORDINATE_ALLELE_MISMATCH' if not c['allele_mismatches'] else 'ALLELE_MISMATCH_REVIEW_REQUIRED'
  out.append(row);print(f"{row['accession_a']} vs {row['accession_b']} overlap={row['shared_coordinates']} single={row['single_variant_coordinate_comparisons']} mismatch={row['allele_mismatches']} palindromic={row['palindromic_coordinate_comparisons']}",flush=True)
 dest=a.output if a.output.is_absolute() else repo/a.output;dest.parent.mkdir(parents=True,exist_ok=True)
 with dest.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
 print(f'LATENT_ALLELE_COMPARISON_WRITTEN rows={len(out)} output={dest}')
if __name__=='__main__':main()
