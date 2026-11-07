import csv, collections, math, pathlib, json
from pyliftover import LiftOver
p=pathlib.Path('outputs/shared-genetics-atlas/frailty_coloc')
lo=LiftOver('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/reference/hg38ToHg19.over.chain.gz')
regions={'3':(117142005,118142005),'7':(113571035,114571035)}
for chrom,pos in [('3',117920728),('7',114430980)]:print('lift_check',chrom,pos,lo.convert_coordinate('chr'+chrom,pos-1))
def read(name):
 with (p/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
a=read('FI_extracted.tsv');b=read('R13_extracted.tsv'); counts=collections.Counter(); index=collections.defaultdict(list)
for r in b:
 chrom=r['#chrom']; loc=lo.convert_coordinate('chr'+chrom,int(r['pos'])-1)
 if len(loc)!=1 or loc[0][0]!='chr'+chrom or loc[0][2]!='+':counts['R13_bad_liftover']+=1;continue
 bp=loc[0][1]+1
 if not regions[chrom][0]<=bp<=regions[chrom][1]:continue
 counts[f'chr{chrom}_R13_in_region']+=1
 index[(chrom,bp,r['rsids'])].append(r)
fiindex=collections.Counter((r['CHR'],int(r['BP']),r['SNP']) for r in a)
out=[]
for r in a:
 c=r['CHR']; counts[f'chr{c}_FI_rows']+=1
 try:
  a1,a2=r['A1'].upper(),r['A2'].upper();f=float(r['FRQ']);beta=float(r['BETA']);se=float(r['SE'])
  if a1 not in 'ACGT' or a2 not in 'ACGT' or len(a1)!=1 or len(a2)!=1 or a1==a2 or frozenset((a1,a2)) in (frozenset('AT'),frozenset('CG')):raise ValueError('allele')
  if not all(map(math.isfinite,[f,beta,se])) or not .01<=f<=.99 or se<=0:raise ValueError('numeric')
 except (ValueError,KeyError):counts[f'chr{c}_FI_ineligible']+=1;continue
 counts[f'chr{c}_FI_eligible']+=1
 k=(c,int(r['BP']),r['SNP'])
 if fiindex[k]!=1:counts[f'chr{c}_FI_duplicate']+=1;continue
 matches=index.get(k,[])
 if len(matches)!=1:counts[f'chr{c}_no_unique_R13_match']+=1;continue
 s=matches[0];ref,alt=s['ref'].upper(),s['alt'].upper()
 if {a1,a2}!={ref,alt}:counts[f'chr{c}_allele_mismatch']+=1;continue
 try:
  freq=float(s['af_alt']);bet=float(s['beta']);sse=float(s['sebeta'])
  if not all(map(math.isfinite,[freq,bet,sse])) or not .01<=freq<=.99 or sse<=0:raise ValueError()
 except ValueError:counts[f'chr{c}_R13_ineligible']+=1;continue
 if a1!=alt:freq=1-freq;bet=-bet
 if abs(f-freq)>.15:counts[f'chr{c}_frequency_discordant']+=1;continue
 out.append({'chr':c,'bp37':r['BP'],'rsid':r['SNP'],'effect_allele':a1,'other_allele':a2,'fi_freq':f,'r13_effect_freq':freq,'fi_beta':beta,'fi_se':se,'r13_beta':bet,'r13_se':sse,'fi_p':r['P'],'r13_p':s['pval']});counts[f'chr{c}_matched']+=1
with (p/'matched.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]) if out else ['chr'],delimiter='\t');w.writeheader();w.writerows(out)
(p/'matching_counts.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts,indent=2))
