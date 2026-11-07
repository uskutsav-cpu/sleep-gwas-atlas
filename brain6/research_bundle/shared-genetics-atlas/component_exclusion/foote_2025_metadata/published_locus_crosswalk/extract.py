import csv,hashlib,json,pathlib,openpyxl
p=pathlib.Path(__file__).resolve().parent;src=pathlib.Path('work/foote2025_supplement_springer.xlsx');w=openpyxl.load_workbook(src,read_only=True,data_only=True)
wanted={'Supplementary Table 13':{'rs1363103','rs2396766','rs8089807'},'Supplementary Table 14':{'rs7703746'}}
out=[]
for sheet,ids in wanted.items():
 ws=w[sheet]; header=[str(c or '') for c in next(ws.iter_rows(min_row=3,max_row=3,values_only=True))[:11]];assert len(set(header))==11
 for rownum,r in enumerate(ws.iter_rows(min_row=4,values_only=True),4):
  if r[1] in ids:
   row={'source_sheet':sheet,'source_row':rownum};row.update(dict(zip(header,r[:11])));out.append(row)
 assert {x['Lead SNP'] for x in out if x['source_sheet']==sheet}==ids
with (p/'published_locus_rows.tsv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');writer.writeheader();writer.writerows(out)
rows=[]
for i,r in enumerate(w['Supplementary Table 1'].iter_rows(values_only=True),1):
 if r[0] and str(r[0]) in ['Freqency of tiredness/lethargy in the last 2 weeks','Insomnia']:
  rows.append({'source_sheet':'Supplementary Table 1','source_row':i,'phenotype':r[0],'domain':r[1],'original_gwas':r[2],'sample_size_as_reported':r[3],'measure':r[4],'source':r[5],'description':r[6],'removal_reason':r[13]})
assert len(rows)==2
with (p/'sleep_input_rows.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
ws=w['Supplementary Table 6'];hdr=[str(v or '') for v in next(ws.iter_rows(min_row=3,max_row=3,values_only=True))[:14]];assert len(set(hdr))==14
loci=[]
for i,r in enumerate(ws.iter_rows(min_row=4,values_only=True),4):
 if r[2]=='rs2431108':
  row={'source_sheet':'Supplementary Table 6','source_row':i};row.update(dict(zip(hdr,r[:14])));loci.append(row)
assert len(loci)==1
with (p/'published_general_locus_row.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(loci[0]),delimiter='\t');wr.writeheader();wr.writerows(loci)
sha=hashlib.sha256(src.read_bytes()).hexdigest();receipt={'source_url':'https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02269-0/MediaObjects/41588_2025_2269_MOESM4_ESM.xlsx','downloaded_workbook_sha256':sha,'workbook_path':str(src),'sheets':'Supplementary Tables 1, 6, 13 and 14','selected_locus_rows':len(out),'general_locus_rows':len(loci),'selected_indicator_rows':len(rows),'note':'Read-only extraction of selected author-reported cells; their replication label is distinct from independent two-trait replication.'};(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('wrote',len(out),'locus rows and',len(rows),'indicator rows; sha256',sha)
