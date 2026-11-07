#!/usr/bin/env python3
import csv,re,sys
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'sample.tsv').open() as f:ids={x['SNP'] for x in csv.DictReader(f,delimiter='\t')}
head=sys.stdin.readline().rstrip('\n').split('\t'); assert head[:5]==['#chrom','pos','ref','alt','rsids']
count=0;hits=[]
for line in sys.stdin:
 count+=1;fields=line.rstrip('\n').split('\t');assert len(fields)==len(head)
 if any(x in ids for x in re.split(r'[,;]',fields[4])):hits.append(fields)
with (p/'raw_matches.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(head);w.writerows(hits)
print('raw_rows',count,'matched_rows',len(hits),file=sys.stderr,flush=True)
