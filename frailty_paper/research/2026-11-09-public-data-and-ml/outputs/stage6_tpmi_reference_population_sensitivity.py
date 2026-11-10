#!/usr/bin/env python3
"""Sensitivity of 17-marker conditional statistic to EAS reference subset."""
import csv,json,math,runpy
from pathlib import Path
import numpy as np

out=Path(__file__).resolve().parent
# This also rebuilds the base 17-marker result and supplies exactly its input vectors.
base=runpy.run_path(str(out/'stage6_tpmi_17_marker_exploratory_z.py'))
x=base['x'];markers=base['markers'];names=base['samples'];chroms=base['chroms'];orient=base['orientation']
panel_path=out.parent/'work/1kg_phase3_panel.tsv'
with panel_path.open(newline='') as f: pop={r['sample']:r['pop'] for r in csv.DictReader(f,delimiter='\t')}
groups={'all_EAS':{'CHB','CHS','JPT','KHV','CDX'},'CHB_CHS':{'CHB','CHS'},'JPT':{'JPT'},'KHV_CDX':{'KHV','CDX'}}
records=[]
for label, group in groups.items():
    chosen=np.array([pop[name] in group for name in names]);g=x[chosen,:];sd=g.std(axis=0,ddof=1)
    zero=[markers[i]['variant'] for i in range(len(markers)) if sd[i]==0]
    # Nonvariable sites contribute zero in this reference. Correlations are
    # defined only for variable sites and denominator is computed as covariance.
    cov=np.cov(g,rowvar=False,ddof=1);cov=np.where(np.equal.outer(chroms,chroms),cov,0.)
    denom=float(np.sqrt(np.ones(len(markers))@cov@np.ones(len(markers))))
    for phenotype in ('insomnia','persistent'):
        z=np.array([float(r[phenotype+'_effect'])/float(r[phenotype+'_se']) for r in markers])*orient
        stat=float(sd@z/denom)
        records.append({'reference_group':label,'reference_n':int(chosen.sum()),'phenotype':phenotype,'weighted_z':stat,'nominal_p':math.erfc(abs(stat)/math.sqrt(2)),'reference_nonvariable_markers':zero})
(out/'stage6_tpmi_reference_population_sensitivity.json').write_text(json.dumps({'records':records,'limitations':'Postselected and conditional; groups overlap all_EAS, small pooled reference subsets are not TPMI target LD; no outcome-level uncertainty or prediction evaluated.'},indent=2)+'\n')
print(json.dumps(records,indent=2))
