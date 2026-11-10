#!/usr/bin/env python3
"""Extract relevant AlphaGenome paper track metadata without changing source labels."""
import csv,hashlib,json,openpyxl,pathlib
p=pathlib.Path('work/alphagenome_2026_supp_tables.xlsx');h=hashlib.sha256(p.read_bytes()).hexdigest()
w=openpyxl.load_workbook(p,read_only=True,data_only=True);s=w['Suppl Table 2 Track metadata (f'];it=s.iter_rows(values_only=True);header=next(it)
fields=['output_type','name','strand','track_index','biosample_name','ontology_curie','data_source','Experiment accession','File accession','provenance_flag'];out=[]
contexts={'SK-N-SH','Jurkat','Jurkat, Clone E6-1','CD4-positive, alpha-beta T cell','naive thymus-derived CD4-positive, alpha-beta T cell'}
for row in it:
 d=dict(zip(header,row))
 if d.get('organism')!='human' or d.get('output_type') not in ('ATAC','DNASE','CAGE') or d.get('biosample_name') not in contexts:continue
 item={k:d.get(k) for k in fields if k!='provenance_flag'};label=str(d.get('File accession') or '').lower()
 item['provenance_flag']='Model ontology label says naive CD4; source file name says memory conventional T cells; confirm version and original metadata' if d.get('biosample_name')=='naive thymus-derived CD4-positive, alpha-beta T cell' and d.get('output_type')=='CAGE' and 'memory' in label else 'No specific sample-label discrepancy assessed'
 out.append(item)
with open('outputs/stage5_alphagenome_priority_track_provenance.csv','w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(out)
pathlib.Path('outputs/stage5_alphagenome_paper_metadata_receipt.json').write_text(json.dumps({'source':'Avsec et al. 2026 supplementary table 2','url':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-025-10014-0/MediaObjects/41586_2025_10014_MOESM3_ESM.xlsx','sha256':h,'sheet':s.title,'selected_rows':len(out),'note':'Published paper metadata; exact equivalence to the live API model version and track provenance has not been independently established.'},indent=2)+'\n')
print('selected rows',len(out),'label flags',sum('source file name says memory' in x['provenance_flag'] for x in out))
