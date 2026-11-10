import csv,gzip,glob,os,statistics
G=['FOXP2','SMAD5','RBM6','MEIS1']; groups={g:{'FRA':[],'NOR':[]} for g in G}; libs=[]
for f in sorted(glob.glob('work/handls/*.count.txt.gz')):
 label=os.path.basename(f).split('_')[1]; counts={}; total=0
 with gzip.open(f,'rt') as h:
  next(h)
  for row in csv.reader(h,delimiter='\t'):
   if len(row)<2: continue
   try: n=float(row[1])
   except ValueError: continue
   if row[0] in counts: raise ValueError('duplicate symbol '+row[0])
   counts[row[0]]=n; total+=n
 if not total: raise ValueError('empty library')
 libs.append(total)
 for g in G:
  groups[g][label].append(None if g not in counts else counts[g]/total*1e6)
with open('outputs/stage6_handls_candidate_expression.csv','w',newline='') as h:
 w=csv.writer(h);w.writerow(['gene','group','n_files','n_nonzero','median_cpm','min_cpm','max_cpm'])
 for g in G:
  for label in ['NOR','FRA']:
   vals=groups[g][label]; valid=[x for x in vals if x is not None]
   w.writerow([g,label,len(valid),sum(x>0 for x in valid),statistics.median(valid) if valid else '',min(valid) if valid else '',max(valid) if valid else ''])
print('libraries',len(libs),'library total range',min(libs),max(libs))
