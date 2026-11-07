import pathlib,zipfile,xml.etree.ElementTree as ET,csv,hashlib,math,json
src=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo/work/song_supplement.docx');out=pathlib.Path('/Users/swethasunilkumar/Documents/Codex/2026-10-27/brain6/outputs/shared-genetics-atlas/local_evidence');ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
with zipfile.ZipFile(src) as f:root=ET.fromstring(f.read('word/document.xml'))
tables=root.findall('.//w:tbl',ns);rows=[]
for tr in tables[1].findall('./w:tr',ns):rows.append([''.join((t.text or '') for t in c.findall('.//w:t',ns)) for c in tr.findall('./w:tc',ns)])
assert rows[0]==['Chr','Start','End','Number of SNPs','local_rhog','SE','P'] and rows[1][:4]==['3','47727212','49316972','967']
with (out/'song_local_reported.tsv').open('w',newline='') as f:w=csv.writer(f,delimiter='\t');w.writerows(rows[:2])
r=rows[1];z=float(r[4])/float(r[5]);approx=math.erfc(abs(z)/math.sqrt(2));receipt={'source_path':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'supplement_table_index_zero_based':1,'test_statistic_from_rounded_estimate_and_se':z,'normal_approx_p_from_rounded_values':approx,'reported_p':float(r[6]),'coordinate_build':'NOT VERIFIED in located source','note':'Reported local_rhog is consistent with HESS local covariance output; exact paper version and transformations not verified. Not an independent rerun.'}
(out/'song_local_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
