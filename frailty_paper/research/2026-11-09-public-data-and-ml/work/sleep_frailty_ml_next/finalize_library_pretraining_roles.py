"""Pretraining neighborhood and exact-DNA guards before the one coverage expansion."""
import sys,zipfile,xml.etree.ElementTree as ET,posixpath,re,bisect,time,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *

SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'
PRE=SCREEN/'candidate_pretrained_readers'
DEST=SCREEN/'Siraj_unused_libraries'
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def main():
    started=time.monotonic();assert json.loads((PRE/'actual_pretraining_library_coverage_audit.json').read_text())['all_actual_sequences_covered']
    path=PRE/'author_large_scale_library_design.xlsx';receipt=json.loads((PRE/'author_pretraining_library_retrieval_receipt.json').read_text());assert sha(path)==receipt['sha256']
    rows=list(csv_rows(DEST/'full_reference_and_shared_core_roles.csv'));reserved=[r for r in rows if r['phase_after_supervised_history_and_shared_core']=='reserved_candidate'];near=set();exact=set();profiles={}
    root_of=lambda r:int(r['joint_component_after_exact_shared_DNA'])
    with zipfile.ZipFile(path) as z:
        strings=[]
        for event,node in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
            if node.tag.endswith('}si'):strings.append(''.join(node.itertext()));node.clear()
        relations={n.attrib['Id']:n.attrib['Target'] for n in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};wb=ET.fromstring(z.read('xl/workbook.xml'));sheets={}
        for sheet in wb.find('m:sheets',NS):
            target=relations[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];sheets[sheet.attrib['name'].split()[0]]=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        def records(cell):
            header=None
            for event,node in ET.iterparse(z.open(sheets[cell]),events=['end']):
                if not node.tag.endswith('}row'):continue
                values={}
                for c in node:
                    if not c.tag.endswith('}c'):continue
                    n=c.find('m:v',NS);v=n.text if n is not None and n.text is not None else ''
                    if c.attrib.get('t')=='s':v=strings[int(v)]
                    elif c.attrib.get('t')=='inlineStr':v=''.join(c.find('m:is',NS).itertext())
                    values[re.sub('[0-9]','',c.attrib['r'])]=v
                if header is None:
                    if 'name' in values.values() and 'chr.hg38' in values.values():header={v:k for k,v in values.items() if v};seqcol=next(k for v,k in header.items() if v.startswith('230nt sequence'))
                    node.clear();continue
                if values.get(header['name']):yield values.get(header['chr.hg38'],''),values.get(header['start.hg38'],''),values.get(header['stop.hg38'],''),values.get(seqcol,'').upper()
                node.clear()
        for cell in sheets:
            intervals=collections.defaultdict(list);unmapped=0;count=0
            for c,s,e,seq in records(cell):
                count+=1;c='chr'+c.removeprefix('chr') if c else ''
                try:s,e=int(float(s)),int(float(e))
                except ValueError:unmapped+=1;continue
                if c not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY','chrM'] or not 0<=s<e:unmapped+=1;continue
                intervals[c].append((s,e))
            index={}
            for c,v in intervals.items():
                v.sort();starts=[s for s,e in v];ends=[];maximum=-1
                for s,e in v:maximum=max(maximum,e);ends.append(maximum)
                index[c]=(starts,ends)
            roots=set()
            for r in reserved:
                p=int(r['position_hg38'])-1;starts,ends=index.get(r['chromosome'],([],[]));j=bisect.bisect_right(starts,p+2000)-1
                if j>=0 and ends[j]>p-2000:roots.add(root_of(r))
            near.update(roots);profiles[cell]={'published_rows':count,'unmapped_or_noncanonical_intervals':unmapped,'reserved_components_near_mapped_library':len(roots)}
            print('PRETRAINED NEIGHBORHOOD GUARD',cell,len(roots),flush=True)
        # Exact comparison is needed only for reserved components that survive
        # the broad mapped-neighborhood screen, keeping the working set small.
        cases=[r for r in reserved if root_of(r) not in near];n=sum(4*(min(int(r['changed_position']),len(r['reference'])-148)-max(0,int(r['changed_position'])-147)+1) for r in cases)
        dtype=[('hash','<u8'),('case','<u4'),('start','u1'),('allele','u1'),('reverse','u1')];bank=np.empty(n,dtype=dtype);prefixes=set();cursor=0
        for i,r in enumerate(cases):
            p=int(r['changed_position']);lo=max(0,p-147);hi=min(p,len(r['reference'])-148)
            for allele,field in enumerate(['reference','alternate']):
                for start in range(lo,hi+1):
                    fragment=r[field][start:start+148]
                    for reverse,s in enumerate([fragment,rc(fragment)]):
                        bank[cursor]=(int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),'big'),i,start,allele,reverse);prefixes.add(s[:12]);cursor+=1
        assert cursor==n
        if n:
            order=np.argsort(bank['hash'],kind='stable');bank=bank[order];del order;keys=np.ascontiguousarray(bank['hash']);cached={}
            def matching_roots(h,sequence):
                key=(int(h),sequence)
                if key not in cached:
                    lo=np.searchsorted(keys,h,'left');hi=np.searchsorted(keys,h,'right');found=set()
                    for record in bank[lo:hi]:
                        r=cases[int(record['case'])];field='alternate' if record['allele'] else 'reference';start=int(record['start']);s=r[field][start:start+148];s=rc(s) if record['reverse'] else s
                        if s==sequence:found.add(root_of(r))
                    cached[key]=found
                return cached[key]
            def flush(batch):
                if not batch:return set()
                hh=np.array([int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),'big') for s in batch],dtype='uint64');pp=np.searchsorted(keys,hh);ids=np.flatnonzero(pp<len(keys));ids=ids[keys[pp[ids]]==hh[ids]];found=set()
                for i in ids:found.update(matching_roots(hh[i],batch[i]))
                return found
            for cell in sheets:
                batch=[];roots=set();windows=0
                for c,s,e,seq in records(cell):
                    if len(seq)!=230 or set(seq)-set('ACGT'):continue
                    for j in range(83):
                        windows+=1
                        if seq[j:j+12] in prefixes:batch.append(seq[j:j+148])
                        if len(batch)>=8192:roots.update(flush(batch));batch=[]
                roots.update(flush(batch));exact.update(roots);profiles[cell].update(exact_scan_windows=windows,additional_literal_shared_DNA_components=len(roots));print('PRETRAINED LITERAL DNA GUARD',cell,len(roots),flush=True)
    output=[]
    for r in rows:
        root=root_of(r);old=r['phase_after_supervised_history_and_shared_core'];role=old
        if old=='reserved_candidate':role='reserved_secondary_pretraining_exposure' if root in near or root in exact else 'reserved_primary'
        output.append({'original_variant_GRCh37':r['original_variant_GRCh37'],'variant_GRCh38':r['variant_GRCh38'],'library':r['library'],'component_root':root,'role':role,
                       'mapped_pretraining_neighborhood_exposed':root in near,'literal_pretraining_DNA_exposed':root in exact,'source_reference_design_row':rows.index(r) if False else '',
                       'signed_source_coefficient_convention_established':False})
    write_csv(DEST/'pretraining_guarded_source_roles.csv',output)
    summary={'passed':True,'profiles':profiles,'by_library':{lib:dict(collections.Counter(r['role'] for r in output if r['library']==lib)) for lib in ['OL27','OL33']},
             'all_actual_author_pretraining_sequences_covered_by_published_library':True,'exact_hash_hits_verified_by_full_literal_DNA':True,
             'whole_component_exposure_propagated':True,'pretraining_exposure_in_development_allowed_not_a_new_test':True,
             'scope':'Conservative whole-published-library screening, including author held-out folds. Unmapped/noncanonical source intervals limit claims about their wider neighborhoods; all source sequences still enter exact comparison.',
             'new_source_effects_or_model_performance_used':False,'assay_and_target_availability_audit_pending':True,'elapsed_seconds':time.monotonic()-started}
    write_json(DEST/'final_pretraining_guard_audit.json',summary);assert_previous_sources();print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
