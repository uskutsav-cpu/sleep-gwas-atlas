#!/usr/bin/env python3
"""Construct pre-coded complete-case raw scores, without performing a GWAS."""
import argparse,csv,json,math,pathlib,sys
ap=argparse.ArgumentParser();ap.add_argument('--input',type=pathlib.Path,required=True);ap.add_argument('--mapping',type=pathlib.Path,required=True);ap.add_argument('--output',type=pathlib.Path,required=True);a=ap.parse_args();m=json.loads(a.mapping.read_text());fi=m['fi_deficits'];ffs=m['ffs_components'];ids=m['id_columns'];sleep=m['fi_sleep'];fatigue=m['fi_fatigue'];exhaust=m['ffs_exhaustion']
assert len(fi)==49 and len(set(fi))==49 and sleep in fi and fatigue in fi and sleep!=fatigue
assert len(ffs)==5 and len(set(ffs))==5 and exhaust in ffs
assert len(ids)>0 and len(set(ids))==len(ids) and not set(ids)&(set(fi)|set(ffs))
assert len(set(fi+ffs))==54, 'Input pre-coded score columns must have distinct names even if derived from the same raw response'
assert a.input.resolve()!=a.output.resolve();a.output.parent.mkdir(parents=True,exist_ok=True)
if a.output.exists():raise FileExistsError('Choose a new output file')
seen=set();counts={'input':0,'complete':0,'missing':0};missing={'','NA','N/A','NULL','.'}
created_output=a.output.open('x',newline='')
try:
 with created_output as g,a.input.open(newline='') as f:
  r=csv.DictReader(f,delimiter='\t');assert set(fi+ffs+ids).issubset(r.fieldnames or []),'Input is missing required columns';assert len(set(r.fieldnames))==len(r.fieldnames),'Duplicate column header'
  w=csv.writer(g,delimiter='\t');w.writerow(ids+['FI49_raw','FI48_no_sleep_raw','FI47_no_sleep_fatigue_raw','FFS5_raw','FFS4_no_exhaustion_raw'])
  for row in r:
   counts['input']+=1;ident=tuple(row[k].strip() for k in ids)
   if any(not x or x.upper() in missing for x in ident):raise ValueError('Missing participant ID')
   if ident in seen:raise ValueError('Duplicate participant ID')
   seen.add(ident)
   if any(row[k].strip().upper() in missing for k in fi+ffs):counts['missing']+=1;continue
   vals={k:float(row[k]) for k in fi+ffs}
   if any(not math.isfinite(vals[k]) or not 0<=vals[k]<=1 for k in fi):raise ValueError('FI deficit outside [0,1] or nonfinite')
   if any(vals[k] not in (0.,1.) for k in ffs):raise ValueError('FFS component must be 0 or 1')
   total=sum(vals[k] for k in fi);ftotal=sum(vals[k] for k in ffs)
   w.writerow([*ident,total/49,(total-vals[sleep])/48,(total-vals[sleep]-vals[fatigue])/47,ftotal,ftotal-vals[exhaust]]);counts['complete']+=1
except Exception:
 a.output.unlink(missing_ok=True)
 raise
print(json.dumps(counts,sort_keys=True));print('Raw research scores only; validate mapping, normalize as planned, and perform approved GWAS separately.',file=sys.stderr)
