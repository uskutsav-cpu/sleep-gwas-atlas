#!/usr/bin/env python3
"""Provisional day/hour-selected GSE56931 profiles; phase metadata conflict; no genotype or clinical outcome."""
import csv, gzip, hashlib, json, math, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
root=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
matrix=root/'GSE56931_series_matrix.txt.gz'; platform=root/'GPL10379.txt'
genes={'CLOCK','ARNTL','PER1','PER2','PER3','CRY1','CRY2','DBP'}
mode='candidates' if '--candidates' in sys.argv else 'clock'
if mode=='candidates':
    genes=set()
    for filename,column in [('outputs/candidate_risk_gene_review_15.csv','nearest_gene_body'),('outputs/unvalidated_22_shared_direction_score.csv','nearest_gene_body_positional')]:
        with open(filename) as source:
            for row in csv.DictReader(source):
                genes.update(x.strip() for x in re.split(r'[,;/]',row[column]) if x.strip())
with platform.open() as f:
    for line in f:
        if line.startswith('!platform_table_begin'): break
    probes={r['ID']:r['GeneSymbol'] for r in csv.DictReader(f,delimiter='\t') if r.get('GeneSymbol') in genes}
with gzip.open(matrix,'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles=next(csv.reader([line],delimiter='\t'))[1:]
        if line.startswith('!series_matrix_table_begin'):break
    reader=csv.reader(f,delimiter='\t'); header=next(reader)[1:]
    assert len(header)==len(titles)
    selected=[]
    for i,title in enumerate(titles):
        match=re.search(r'patient (BH\d+).*day\s*(\d+),\s*hour\s*(\d+)',title,re.I)
        if not match: raise ValueError(title)
        subject,day,hour=match.group(1),int(match.group(2)),int(match.group(3))
        if (day==1 and hour in (8,12,16,20)) or (day==2 and hour in (0,4)):
            selected.append((i,subject,hour))
    assert len({(s,h) for _,s,h in selected})==len(selected)
    rows=[]
    for row in reader:
        if row[0] not in probes:continue
        values=defaultdict(list)
        for i,subject,hour in selected:
            try:v=float(row[i+1])
            except ValueError:continue
            if math.isfinite(v):values[subject].append((hour,v))
        centered=defaultdict(list)
        for subject,obs in values.items():
            if len(obs)<4:continue
            mean=np.mean([v for _,v in obs])
            for hour,v in obs:centered[hour].append(v-mean)
        if len(centered)<6:continue
        hours=sorted(centered); y=np.array([np.mean(centered[h]) for h in hours]); theta=2*np.pi*np.array(hours)/24
        coef=np.linalg.lstsq(np.column_stack([np.ones(len(hours)),np.cos(theta),np.sin(theta)]),y,rcond=None)[0]
        amp=math.hypot(coef[1],coef[2]); peak=(math.atan2(coef[2],coef[1])%(2*math.pi))*24/(2*math.pi)
        rows.append({'gene':probes[row[0]],'probe':row[0],'subjects':sum(len(v)>=4 for v in values.values()),'samples':sum(len(v) for v in centered.values()),'amplitude_log2_units':round(amp,4),'fitted_peak_clock_hour':round(peak,2),'note':'descriptive; no statistical test'})
rows.sort(key=lambda x:(x['gene'],x['probe']))
out=Path('outputs/stage6_gse56931_'+mode+'_baseline_descriptive.csv')
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps({'samples':len(selected),'subjects':len(set(s for _,s,_ in selected)),'probes':len(rows),'requested_genes':len(genes),'matched_genes':len({r['gene'] for r in rows}),'unmatched_genes':sorted(genes-{r['gene'] for r in rows}),'matrix_sha256':hashlib.file_digest(matrix.open('rb'),'sha256').hexdigest(),'platform_sha256':hashlib.file_digest(platform.open('rb'),'sha256').hexdigest()},indent=2))
