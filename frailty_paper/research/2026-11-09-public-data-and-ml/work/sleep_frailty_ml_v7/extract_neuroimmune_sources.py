"""Read-only extraction of verified publisher worksheets into raw CSV intermediates."""
import csv,json,hashlib,posixpath,zipfile,xml.etree.ElementTree as ET,re
from pathlib import Path
BASE=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work/sleep_frailty_ml_v7');NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read_tables(path,names):
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None;strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for event,n in ET.iterparse(z.open('xl/sharedStrings.xml'),events=('end',)):
                if n.tag.endswith('}si'):strings.append(''.join(n.itertext()));n.clear()
        rel={n.attrib['Id']:n.attrib['Target'] for n in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};book=ET.fromstring(z.read('xl/workbook.xml'));result={}
        for sheet in book.find('m:sheets',NS):
            if sheet.attrib['name'] not in names:continue
            target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);rows=[];formula_count=0
            for event,node in ET.iterparse(z.open(target),events=('end',)):
                if not node.tag.endswith('}row'):continue
                cells={}
                for c in node:
                    if not c.tag.endswith('}c'):continue
                    col=re.match('[A-Z]+',c.attrib['r']).group();index=0
                    for char in col:index=index*26+ord(char)-64
                    v=c.find('m:v',NS);value=v.text if v is not None and v.text is not None else '';kind=c.attrib.get('t','n')
                    if kind=='s':value=strings[int(value)] if value else ''
                    elif kind=='inlineStr':value=''.join(c.find('m:is',NS).itertext())
                    if c.find('m:f',NS) is not None:formula_count+=1
                    cells[index-1]=value
                if cells and any(v!='' for v in cells.values()):rows.append([cells.get(k,'') for k in range(max(cells)+1)])
                node.clear()
            result[sheet.attrib['name']]={'rows':rows,'source_formula_cells':formula_count}
        assert set(result)==set(names);return result
def main():
    source=BASE/'sources/neuroimmune_publisher_supplement.xlsx';mirror=BASE/'sources/S25D8RVYD_41467_2026_77114_MOESM5_ESM.xlsx';assert sha(source)==sha(mirror);names={'SI Data 2.9':'design','SI_Data 2.12':'baseline','SI_Data 2.13':'IL6','SI_Data 2.14':'IFNa'};tables=read_tables(source,names);folder=BASE/'neuroimmune_source_tables';folder.mkdir(exist_ok=True);receipt={}
    for name,label in names.items():
        content=tables[name];rows=content['rows'];header=rows[0];assert len(set(header))==len(header) and len(header)>1;path=folder/(label+'.csv')
        with path.open('w') as f:w=csv.writer(f);w.writerows(rows)
        receipt[name]={'raw_table':str(path),'sha256':sha(path),'rows_excluding_header':len(rows)-1,'columns':header,'formula_cells':content['source_formula_cells']};print('SOURCE TABLE',label,len(rows)-1,'rows',flush=True)
    (BASE/'outputs/neuroimmune_raw_extract_receipt.json').write_text(json.dumps({'publisher_workbook_sha256':sha(source),'archive_identical_to_publisher':True,'sheets':receipt,'source_workbook_modified':False,'training_or_model_performance_use':False,'purpose':'Additional candidate source provenance and data availability, separate from current fixed V7 experiment.'},indent=2)+'\n')
if __name__=='__main__':main()
