#!/usr/bin/env python3
"""Read public VETSA concordance and capture relevant variable descriptions."""
from pathlib import Path
import openpyxl,re,json,hashlib
O=Path('outputs');p=Path('work/stage6_vetsa_concordance.xlsx');w=openpyxl.load_workbook(p,read_only=True,data_only=True);terms={'sleep':r'\b(insomnia|sleep|sleeping)\b','grip':r'\bgrip\b','walking':r'\b(walking|walk|gait)\b','fatigue':r'\b(fatigue|tired|exhaustion|exhausted)\b','weight_loss':r'\b(weight\s+loss|lost\s+weight|lose\s+weight)\b','physical_activity':r'\b(physical\s+activity|exercise)\b','frailty':r'\bfrail(?:ty)?\b','genotypes':r'\b(genotype|genotypes|genotyping|polygenic|genome-wide|single.nucleotide|snp)\b'};matches={k:[] for k in terms}
for sheet in w:
 for index,row in enumerate(sheet.iter_rows(values_only=True),1):
  cells=[str(v).strip() for v in row if v is not None];text=' | '.join(cells)
  for key,pattern in terms.items():
   if re.search(pattern,text,re.I):matches[key].append({'sheet':sheet.title,'row':index,'source_text':text[:700]})
r={'source_url':'https://psychiatry.ucsd.edu/research/programs-centers/vetsa/_files/VETSAComposite-concordance_2023_12_01.xlsx','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'publisher_checksum_verified':False,'sheets':w.sheetnames,'matches':matches,'scope':'Term-based inventory only. Concordance is not a complete inventory of public-use files; neither genotype availability nor clinical outcome construction is verified.'};(O/'stage6_vetsa_variable_screen.json').write_text(json.dumps(r,indent=2)+'\n');print({k:len(v) for k,v in matches.items()})
