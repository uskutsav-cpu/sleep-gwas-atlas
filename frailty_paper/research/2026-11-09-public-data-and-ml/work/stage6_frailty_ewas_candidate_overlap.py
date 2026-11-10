import openpyxl,csv,json,hashlib
from pathlib import Path
p=Path('work/mak2024_ewas/ACE-23-0795_supplement_2.xlsx');w=openpyxl.load_workbook(p,read_only=True,data_only=True)
candidates=list(csv.DictReader(open('outputs/stage6_gsa_reference_tag_screen.csv')));out=[];counts={}
for sheet in ['ST3','ST4']:
 rows=list(w[sheet].values)[3:];counts[sheet]=len(rows)
 for r in rows:
  if not isinstance(r[2],(int,float)) or not isinstance(r[3],(int,float)):continue
  for x in candidates:
   if str(int(r[2]))!=x['chr_GRCh37']:continue
   distance=int(r[3])-int(x['position_GRCh37'])
   if abs(distance)<=100000:out.append(dict(sheet=sheet,candidate=x['variant'],score_lists=x['score_lists'],candidate_pos=x['position_GRCh37'],chromosome=int(r[2]),cpg=r[0],cpg_position=int(r[3]),gene_annotation=r[1],cpg_region=r[4],distance_bp=distance,meta_beta=r[26],meta_se=r[27],meta_p=r[28],meta_fdr=r[29],heterogeneity_I2=r[34]))
assert len(candidates)==33
fields=['sheet','candidate','score_lists','candidate_pos','chromosome','cpg','cpg_position','gene_annotation','cpg_region','distance_bp','meta_beta','meta_se','meta_p','meta_fdr','heterogeneity_I2']
with open('outputs/stage6_frailty_ewas_candidate_overlap.csv','w',newline='') as f:
 z=csv.DictWriter(f,fieldnames=fields);z.writeheader();z.writerows(out)
receipt={'selected_ewas_only':True,'window_bp':100000,'n_candidates':len(candidates),'rows_in_source':counts,'overlap_rows':len(out),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'limitations':'Position only; no mQTL, colocalization, gene regulatory relationship, causality or independent score prediction inferred.'}
Path('outputs/stage6_frailty_ewas_candidate_overlap.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print(out)
