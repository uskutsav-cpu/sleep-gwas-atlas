#!/usr/bin/env python3
"""Inventory and descriptive within-person harmonic fit in public GSE39445."""
import csv,gzip,hashlib,json,math,re
from collections import defaultdict
from pathlib import Path
import numpy as np
root=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
matrix=root/'GSE39445_series_matrix.txt.gz'; platform=root/'GPL15331.txt'
genes={'CLOCK','ARNTL','PER1','PER2','PER3','CRY1','CRY2','DBP'}
for filename,col in [('outputs/candidate_risk_gene_review_15.csv','nearest_gene_body'),('outputs/unvalidated_22_shared_direction_score.csv','nearest_gene_body_positional')]:
    with open(filename) as f:
        for row in csv.DictReader(f):genes.update(x.strip() for x in re.split(r'[,;/]',row[col]) if x.strip())
with platform.open() as f:
    for line in f:
        if line.startswith('!platform_table_begin'):break
    probes={r['ID']:r['GENE_SYMBOL'] for r in csv.DictReader(f,delimiter='\t') if r.get('GENE_SYMBOL') in genes}
with gzip.open(matrix,'rt') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            titles=next(csv.reader([line],delimiter='\t'))[1:]; meta=[{} for _ in titles]
        elif line.startswith('!Sample_characteristics_ch1'):
            for d,v in zip(meta,next(csv.reader([line],delimiter='\t'))[1:]):
                key,sep,val=v.partition(': ')
                if sep:d[key]=val
        elif line.startswith('!series_matrix_table_begin'):break
    reader=csv.reader(f,delimiter='\t');header=next(reader)[1:]
    assert len(header)==len(meta)
    selected=[]
    recovered=0
    for i,d in enumerate(meta):
        title_match=re.search(r'(AF\d+)_([A-Z]+)_T\d+',titles[i])
        if not title_match:raise ValueError(f'Unparseable title: {titles[i]}')
        title_subject=title_match.group(1)
        if 'subject' in d and d['subject']!=title_subject:
            raise ValueError(f'Subject mismatch: {titles[i]} vs {d["subject"]}')
        if d.get('sleepprotocol')!='Sleep Extension' or 'circadianphase' not in d:continue
        if title_match.group(2)!='SE':raise ValueError(f'Condition mismatch: {titles[i]}')
        try:phase=float(d['circadianphase'])
        except ValueError:continue
        if 'subject' not in d: recovered+=1
        selected.append((i,d.get('subject',title_subject),phase))
    rows=[]
    for r in reader:
        if r[0] not in probes:continue
        groups=defaultdict(list)
        for i,subject,phase in selected:
            try: v=float(r[i+1])
            except ValueError: continue
            if math.isfinite(v):groups[subject].append((phase,v))
        x=[];y=[]
        for subject,obs in groups.items():
            if len(obs)<4:continue
            vals=np.array([v for _,v in obs]); sd=vals.std()
            if not sd:continue
            phases=np.array([t for t,_ in obs]);th=2*np.pi*phases/24
            cos=np.cos(th);sin=np.sin(th)
            x.extend(zip(cos-cos.mean(),sin-sin.mean()));y.extend((vals-vals.mean())/sd)
        if not x:continue
        coeff=np.linalg.lstsq(np.array(x),np.array(y),rcond=None)[0]
        rows.append({'gene':probes[r[0]],'probe':r[0],'n_subjects':sum(len(v)>=4 and np.std([x[1] for x in v])>0 for v in groups.values()),'n_samples':len(y),'standardized_harmonic_amplitude':round(float(np.hypot(*coeff)),4),'peak_relative_to_phase_zero_h':round(float((math.atan2(coeff[1],coeff[0])%(2*math.pi))*24/(2*math.pi)),2),'interpretation':'descriptive only; no significance test'})
rows.sort(key=lambda x:(x['gene'],x['probe']))
out=Path('outputs/stage6_gse39445_candidate_clock_probe_inventory.csv')
with out.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print(json.dumps({'eligible_samples':len(selected),'subject_ids_recovered_from_title':recovered,'eligible_subjects':len({s for _,s,_ in selected}),'requested_genes':len(genes),'matched_genes':len({r['gene'] for r in rows}),'probes':len(rows),'matrix_sha256':hashlib.file_digest(matrix.open('rb'),'sha256').hexdigest(),'platform_sha256':hashlib.file_digest(platform.open('rb'),'sha256').hexdigest()},indent=2))
