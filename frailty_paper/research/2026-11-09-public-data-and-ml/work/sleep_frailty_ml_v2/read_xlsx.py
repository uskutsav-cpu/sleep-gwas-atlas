"""Read-only XLSX value extraction with the Python standard library.

No workbook writing, formulas, styling or arbitrary code execution.
"""
import zipfile,xml.etree.ElementTree as E,re
N={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def rows(path,sheet):
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=[''.join(x.itertext()) for x in E.fromstring(z.read('xl/sharedStrings.xml'))]
        rels={r.attrib['Id']:r.attrib['Target'] for r in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        wb=E.fromstring(z.read('xl/workbook.xml'));s=next(s for s in wb.find('m:sheets',N) if s.attrib['name']==sheet)
        target=rels[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
        with z.open(target) as f:
            for event,row in E.iterparse(f,events=['end']):
                if row.tag!='{'+N['m']+'}row':continue
                values={}
                for c in row:
                    ref=re.match(r'[A-Z]+',c.attrib['r'])[0];col=0
                    for a in ref:col=26*col+ord(a)-64
                    v=c.find('m:v',N);text=v.text if v is not None else ''
                    if c.attrib.get('t')=='s':text=strings[int(text)]
                    if c.attrib.get('t')=='inlineStr':text=''.join(c.find('m:is',N).itertext())
                    values[col-1]=text
                yield [values.get(i,'') for i in range(max(values,default=-1)+1)]
                row.clear()
def dictionaries(path,sheet):
    it=rows(path,sheet);header=next(it)
    for row in it:yield dict(zip(header,row))
