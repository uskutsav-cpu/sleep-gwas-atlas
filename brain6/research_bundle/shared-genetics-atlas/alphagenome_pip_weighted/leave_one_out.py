"""Post-hoc winner-removal sensitivity; descriptive, not a significance test."""
import csv, collections, math, pathlib
p=pathlib.Path(__file__).resolve().parent
with (p/'variant_contributions.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
groups=collections.defaultdict(list)
for row in rows: groups[(row['region'],row['trait'])].append(row)
results=[]
for (region,trait), group in sorted(groups.items()):
    pp=lambda x:float(x['pip'])
    effect=lambda x:float(x['abs_dnase'])
    winner=sorted(group,key=lambda x:(-pp(x)*effect(x),x['snp']))[0]
    rest=[x for x in group if x is not winner]
    mass=math.fsum(map(pp,group)); rest_mass=math.fsum(map(pp,rest))
    assert len(rest)>0 and 0<rest_mass<=mass
    original=math.fsum(pp(x)*effect(x) for x in group)/mass
    removed=math.fsum(pp(x)*effect(x) for x in rest)/rest_mass
    unweighted_rest=math.fsum(map(effect,rest))/len(rest)
    results.append({'region':region,'trait':trait,'removed_variant':winner['snp'],
      'original_weighted_mean':original,'removed_pip':pp(winner),
      'remaining_pip_mass':rest_mass,'weighted_mean_without_winner':removed,
      'unweighted_mean_without_winner':unweighted_rest,'weighted_change':removed-original})
with (p/'leave_one_out.tsv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=results[0],delimiter='\t');w.writeheader();w.writerows(results)
print('Checked',len(results),'eligible region-trait pairs')
