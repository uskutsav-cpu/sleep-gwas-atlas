#!/usr/bin/env python3
"""Check the exact local score definitions and post hoc result artifacts."""
from pathlib import Path
import pandas as pd,json,hashlib
O=Path('outputs');p15=O/'unvalidated_15_variant_score_prototype.csv';p22=O/'unvalidated_22_shared_direction_score.csv';s15=pd.read_csv(p15);s22=pd.read_csv(p22);assert len(s15)==15 and s15.variant.is_unique and len(s22)==22 and s22.variant.is_unique and (s22.weight==1).all();w=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'];assert s15[w].notna().all().all()
h=pd.read_csv(O/'stage6_hfrs_sparse_score_summary_test.csv');h=h[h.omitted_variant.eq('none')].set_index('weight').loc[w];i=pd.read_csv(O/'stage6_insomnia_sparse_score_summary.csv');rows={}
for source,pop in [('MVP','EUR503'),('R13','FIN99')]:
 x=i[(i.source==source)&(i.reference==pop)&(i.omitted=='none')].set_index('weight').loc[w];assert (x.n==15).all();rows[source]=dict(zip(w,x.weighted_z.astype(float)))
receipt=json.loads((O/'stage6_insomnia_sparse_score_summary_receipt.json').read_text());assert receipt['weight_sha256']==hashlib.sha256(p15.read_bytes()).hexdigest();m=sorted(set(s15.variant)&set(s22.variant));out={'prototype_15':{'n':len(s15),'weight_columns':w,'sha256':hashlib.sha256(p15.read_bytes()).hexdigest()},'prototype_22_equal_weight':{'n':len(s22),'sha256':hashlib.sha256(p22.read_bytes()).hexdigest()},'shared_rsid_count':len(m),'posthoc_weighted_z':{'HFRS_FIN99':dict(zip(w,h.weighted_z.astype(float))),**rows},'scope':'Source score definitions use different allele orientations and selection; shared IDs do not prove same counted allele. Summary statistics do not validate predictive performance.'};(O/'stage6_score_evidence_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
