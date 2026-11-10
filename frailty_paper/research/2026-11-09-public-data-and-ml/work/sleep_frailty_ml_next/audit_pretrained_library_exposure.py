"""Full published library exposure profile, with no new experimental outcomes."""
import csv,json,hashlib,zipfile,xml.etree.ElementTree as ET,posixpath,re,bisect,time,collections
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'
OUT=SOURCE/'candidate_pretrained_readers'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def main():
    started=time.monotonic();path=OUT/'author_large_scale_library_design.xlsx';receipt=json.loads((OUT/'author_pretraining_library_retrieval_receipt.json').read_text());assert hashlib.sha256(path.read_bytes()).hexdigest()==receipt['sha256']
    with (SOURCE/'Tardaguila_reserved_no_historical_overlap.csv').open() as f:cases=list(csv.DictReader(f))
    core={}
    for i,r in enumerate(cases):
        for role in ['reference','alternate']:
            sequence=r[role]
            for s in [sequence,rc(sequence)]:core.setdefault(s,set()).add((i,role))
    prefixes={s[:12] for s in core};exposure={};profiles={}
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None;strings=[]
        for event,n in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
            if n.tag.endswith('}si'):strings.append(''.join(n.itertext()));n.clear()
        rel={n.attrib['Id']:n.attrib['Target'] for n in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};wb=ET.fromstring(z.read('xl/workbook.xml'))
        for sheet in wb.find('m:sheets',NS):
            cell=sheet.attrib['name'].split()[0];assert cell in ['HepG2','K562','WTC11'];target=rel[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
            header=None;intervals=collections.defaultdict(list);matches=set();rows=0;valid_sequences=0;windows=0;unmapped=0
            for event,row in ET.iterparse(z.open(target),events=['end']):
                if not row.tag.endswith('}row'):continue
                values={}
                for c in row:
                    if not c.tag.endswith('}c'):continue
                    column=re.sub('[0-9]','',c.attrib['r']);node=c.find('m:v',NS);value=node.text if node is not None and node.text is not None else ''
                    if c.attrib.get('t')=='s':value=strings[int(value)]
                    elif c.attrib.get('t')=='inlineStr':value=''.join(c.find('m:is',NS).itertext())
                    values[column]=value
                if header is None:
                    if 'name' in values.values() and 'chr.hg38' in values.values():header={v:k for k,v in values.items() if v};seqcol=next(k for v,k in header.items() if v.startswith('230nt sequence'))
                    row.clear();continue
                if not values.get(header['name']):row.clear();continue
                rows+=1;sequence=values.get(seqcol,'').upper()
                if len(sequence)==230 and not set(sequence)-set('ACGT'):
                    valid_sequences+=1
                    for j in range(len(sequence)-148+1):
                        windows+=1
                        if sequence[j:j+12] in prefixes:
                            hit=core.get(sequence[j:j+148])
                            if hit:matches.update(hit)
                chrom=values.get(header['chr.hg38'],'');chrom='chr'+chrom.removeprefix('chr') if chrom else ''
                try:start=int(float(values.get(header['start.hg38'],'')));stop=int(float(values.get(header['stop.hg38'],'')))
                except ValueError:unmapped+=1;row.clear();continue
                if chrom not in ['chr'+str(j) for j in range(1,23)]+['chrX','chrY','chrM'] or not 0<=start<stop:unmapped+=1
                else:intervals[chrom].append((start,stop))
                row.clear()
            index={}
            for chrom,v in intervals.items():
                v.sort();starts=[x[0] for x in v];ends=[];maximum=-1
                for _,end in v:maximum=max(maximum,end);ends.append(maximum)
                index[chrom]=(starts,ends)
            per=[]
            for i,r in enumerate(cases):
                pos=int(r['position_hg38'])-1;starts,ends=index.get(r['chromosome'],([],[]));j=bisect.bisect_right(starts,pos+2000)-1;near=j>=0 and ends[j]>=pos-2000
                per.append({'reference_exact_148bp_training_library_match':(i,'reference') in matches,'alternate_exact_148bp_training_library_match':(i,'alternate') in matches,'within_2kb_of_published_library_interval':near})
            exposure[cell]=per;profiles[cell]={'published_library_rows':rows,'unambiguous_actual_230bp_sequences':valid_sequences,'all_148bp_substrings_examined':windows,'unmapped_or_noncanonical_interval_rows':unmapped,
                'candidate_reference_exact_match_designs':sum(r['reference_exact_148bp_training_library_match'] for r in per),'candidate_alternate_exact_match_designs':sum(r['alternate_exact_148bp_training_library_match'] for r in per),'candidate_near_library_designs':sum(r['within_2kb_of_published_library_interval'] for r in per)}
            print('PUBLISHED LIBRARY EXPOSURE SCAN',cell,profiles[cell],flush=True)
    output=[]
    for i,r in enumerate(cases):
        result={k:r[k] for k in ['source_SNP','comparison','window','variant_GRCh38','component_root']}
        for cell in exposure:
            for k,v in exposure[cell][i].items():result[cell+'__'+k]=v
        output.append(result)
    with (SOURCE/'candidate_pretraining_exposure.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    # Propagate potential exposure across the complete candidate component.
    exposed_roots={int(cases[i]['component_root']) for cell in exposure for i,r in enumerate(exposure[cell]) if any(r.values())}
    strict=[r for r in cases if int(r['component_root']) not in exposed_roots]
    if strict:
        with (SOURCE/'candidate_strict_no_published_pretraining_neighborhood.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(strict[0]));w.writeheader();w.writerows(strict)
    summary={'passed':True,'source_library_sha256':receipt['sha256'],'by_reader_library':profiles,'candidate_designs':len(cases),'strict_no_exposure_designs':len(strict),'strict_no_exposure_distinct_variants':len({r['variant_GRCh38'] for r in strict}),'strict_no_exposure_regions':len({r['component_root'] for r in strict}),
        'exposure_scope':'Whole published large-scale library, regardless of whether a particular row belongs to the selected reader training, validation or test folds. Potential exposure is conservative and is not proof that a selected reader trained on the row.',
        'new_source_training_and_test_eligibility_not_changed_by_this_profile':True,'new_outcomes_used_to_choose_model_or_calculate_performance':False,
        'full_biological_or_donor_independence_established':False,'elapsed_seconds':time.monotonic()-started}
    (SOURCE/'full_published_pretraining_exposure_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
