#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent.parent
queries=[('work/niagads_tad_finnish_ld_anchor_4m.json','5:137615744:G:A',['5:136202653:A:C','5:136206558:C:T'],4000000),
         ('work/niagads_smad5_finnish_ld_anchor.json','5:136202653:A:C',['5:136206558:C:T'],100000)]
out=[]
for path,anchor,targets,window in queries:
 p=root/path; j=json.loads(p.read_text()); rows=j['ld']; assert rows and all(r['variation1']==anchor for r in rows)
 for target in targets:
  found=[r for r in rows if r['variation2']==target];assert len(found)==1
  out.append({'anchor':anchor,'target':target,'window_parameter':window,'r2':found[0]['r2'],
              'response_rows':len(rows),'response_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
receipt={'panel':'sisu42','source':'https://api.finngen.fi/api/ld','results':out,
         'limits':'Public reference-panel pairwise LD; rounded r2; not study-specific signed LD, fine-mapping or score validation.'}
(root/'outputs/stage6_niagads_sisu_pair_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(out,indent=2))
