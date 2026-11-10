#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,json
root=Path(__file__).resolve().parents[1];o=root/'outputs';w=root/'work';d=pd.read_csv(o/'stage6_foxp2_gnomad_narrow_variants.csv');n=len(d);ids=d.rs.tolist();a=ids.index('rs12536335');b=ids.index('rs35624276');results=[]
for pop in ['fin','nfe']:
 R=np.fromfile(w/f'stage6_foxp2_gnomad_{pop}_narrow_ld.bin',dtype='<f8').reshape(n,n);eig=np.linalg.eigvalsh(R);results.append({'population':pop,'n_variants':n,'min_eigen':float(eig.min()),'n_negative_over_1e_8':int(sum(eig < -1e-8)),'anchor_r':float(R[a,b]),'anchor_r2':float(R[a,b]**2),'source_note':'full genome file checksum not verified; source r sign cross-check via 1000G/SISu'})
(o/'stage6_foxp2_gnomad_narrow_qc.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
