import csv,json,hashlib
import numpy as np
from pathlib import Path
root=Path(__file__).resolve().parents[1]; archive=np.load(root/'outputs/stage5_smad5_muscle_eur_reference_ld.npz');R=archive['r'];names=archive['variant'].tolist();assert len(names)==R.shape[0]
def convert(k):c,p,r,a=k.split(':');return f'chr{c}:{p}:{r}>{a}'
rows={convert(r['variant_GRCh38']):r for r in csv.DictReader((root/'outputs/stage6_tenk10k_peak_pal_freq_Mak_FI_aligned.csv').open())}
idx=[i for i,k in enumerate(names) if k in rows];assert len(idx)>0
out=root/'work/stage6_tenk10k_susie_diagnostic';out.mkdir(exist_ok=True)
selected=[rows[names[i]] for i in idx]
with (out/'aligned.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=selected[0].keys());w.writeheader();w.writerows(selected)
a=np.asarray(R[np.ix_(idx,idx)],dtype='<f8',order='F');(out/'ld_f64.bin').write_bytes(a.tobytes(order='F'))
print('variants',len(idx),'sample ref',archive['n_samples'],'ld bytes',a.nbytes)
