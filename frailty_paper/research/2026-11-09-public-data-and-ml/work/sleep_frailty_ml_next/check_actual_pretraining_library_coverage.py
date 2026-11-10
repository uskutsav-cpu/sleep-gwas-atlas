"""Check that the conservative published-library scan covers real author examples."""
import csv,gzip,json,zipfile,xml.etree.ElementTree as ET,posixpath,re,hashlib,time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def canonical(s):return min(s,s.translate(str.maketrans('ACGT','TGCA'))[::-1])


def main():
    records=json.loads((OUT/'all_pretraining_fold_manifest_receipt.json').read_text())['records'];bycell={r['cell']:r for r in records}
    for r in records:assert hashlib.sha256(Path(r['manifest_path']).read_bytes()).hexdigest()==r['manifest_sha256']
    path=OUT/'author_large_scale_library_design.xlsx';summary={}
    with zipfile.ZipFile(path) as z:
        strings=[]
        for event,n in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
            if n.tag.endswith('}si'):strings.append(''.join(n.itertext()));n.clear()
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};wb=ET.fromstring(z.read('xl/workbook.xml'))
        for sheet in wb.find('m:sheets',NS):
            cell=sheet.attrib['name'].split()[0];target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);sequences=set();header=None
            for event,row in ET.iterparse(z.open(target),events=['end']):
                if not row.tag.endswith('}row'):continue
                values={}
                for c in row:
                    if not c.tag.endswith('}c'):continue
                    v=c.find('m:v',NS);text=v.text if v is not None and v.text is not None else ''
                    if c.attrib.get('t')=='s':text=strings[int(text)]
                    elif c.attrib.get('t')=='inlineStr':text=''.join(c.find('m:is',NS).itertext())
                    values[re.sub('[0-9]','',c.attrib['r'])]=text
                if header is None:
                    if 'name' in values.values() and 'chr.hg38' in values.values():header=values;seqcol=next(k for k,v in values.items() if v.startswith('230nt sequence'))
                    row.clear();continue
                seq=values.get(seqcol,'').upper()
                if len(seq)==230 and not set(seq)-set('ACGT'):sequences.add(canonical(seq))
                row.clear()
            total=0;missing=0;folds={str(i):0 for i in range(1,11)};missing_ids=[]
            with gzip.open(bycell[cell]['manifest_path'],'rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    total+=1;folds[r['fold']]+=1
                    if canonical(r['seq']) not in sequences:
                        missing+=1
                        if len(missing_ids)<10:missing_ids.append(r['seq_id'])
            summary[cell]={'actual_author_rows':total,'author_fold_counts':folds,'rows_whose_full_230bp_DNA_not_in_published_library_in_either_orientation':missing,'example_unmapped_ids':missing_ids,
                           'all_actual_training_and_validation_sequences_covered_by_published_library_scan':missing==0}
            print('ACTUAL PRETRAINING COVERAGE',cell,total,'uncovered_rows',missing,flush=True)
    complete=all(r['all_actual_training_and_validation_sequences_covered_by_published_library_scan'] for r in summary.values())
    (OUT/'actual_pretraining_library_coverage_audit.json').write_text(json.dumps({'by_cell':summary,'all_actual_sequences_covered':complete,
        'strict_source_subset_verified_against_actual_pretraining_examples':complete,
        'interpretation':'Published-library profiling is conservative only if every actual example is covered. Uncovered examples require an additional exposure scan; do not infer absence of leakage from the paper table alone.',
        'source_coordinates_for_noncanonical_or_unmapped_library_rows_still_limited':True,'new_experimental_model_performance_read':False},indent=2)+'\n')


if __name__=='__main__':main()
