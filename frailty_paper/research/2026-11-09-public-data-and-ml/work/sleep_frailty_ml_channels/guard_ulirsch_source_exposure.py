"""Complete-source literal DNA components, history and pretrained-library guards."""
import re,bisect,zipfile,xml.etree.ElementTree as ET,posixpath,time
from channel_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
SOURCE=DEST/'source_candidates';PREP=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/candidate_pretrained_readers';NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

def main():
    assert_frozen_previous();start=time.monotonic();audit=json.loads((SOURCE/'Ulirsch_physical_and_source_count_audit.json').read_text());assert audit['passed'];assert sha(SOURCE/'Ulirsch_full_reference_and_spatial_audit.csv')==audit['physical_table_sha256'];assert json.loads((PREP/'actual_pretraining_library_coverage_audit.json').read_text())['all_actual_sequences_covered']
    rows=list(csv_rows(SOURCE/'Ulirsch_full_reference_and_spatial_audit.csv'));designs=list(csv.reader(gzip.open(SOURCE/'GSE87711_design.txt.gz','rt'),delimiter='\t'));variants=sorted({identifier.split('_')[0] for identifier,sequence in designs});lift=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));points=collections.defaultdict(set);mapped={}
    for variant in variants:
        c,p,ra,aa=variant.split(':');m=lift.convert_coordinate('chr'+c.removeprefix('chr'),int(p)-1)
        if len(m)==1:chrom,q,strand,_=m[0];points[chrom].add(int(q)+1);mapped[variant]=(chrom,int(q)+1)
    spatial={};root=0
    for chrom,positions in sorted(points.items()):
        previous=None
        for p in sorted(positions):
            if previous is None or p-previous>2000:root+=1
            spatial[(chrom,p)]=root;previous=p
    node={variant:spatial[mapped[variant]] if variant in mapped else -i-1 for i,variant in enumerate(variants)};parent={v:v for v in node.values()}
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    lookup={};joins=0
    for identifier,sequence in designs:
        v=node[identifier.split('_')[0]];canonical=min(sequence,rc(sequence))
        if canonical in lookup:
            if find(v)!=find(lookup[canonical]):joins+=1
            union(v,lookup[canonical])
        else:lookup[canonical]=v
    core={}
    for sequence,n in lookup.items():core[sequence]=find(n);core[rc(sequence)]=find(n)
    prefixes={sequence[:12] for sequence in core};old_positions=collections.defaultdict(set);history_hits=set();history_profiles=[]
    for path,expected in audit['historical_files_sha256'].items():
        assert sha(path)==expected;count=0;windows=0;hits=set()
        for row in csv_rows(path):
            count+=1;old_positions[row['chromosome']].add(int(row['position_hg38']));j=int(row['changed_position'])
            for field in ['reference','alternate']:
                sequence=row[field]
                for offset in range(max(0,j-144),min(j,len(sequence)-145)+1):
                    windows+=1
                    if sequence[offset:offset+12] in prefixes:
                        hit=core.get(sequence[offset:offset+145])
                        if hit is not None:hits.add(hit)
        history_hits.update(hits);history_profiles.append({'path':path,'sha256':expected,'physical_rows':count,'variant_containing_145bp_windows':windows,'shared_components':len(hits)});print('ULIRSCH HISTORY GUARD',Path(path).parent.name,count,windows,flush=True)
    old_positions={c:sorted(v) for c,v in old_positions.items()};near_history=set()
    for point,n in spatial.items():
        c,p=point;values=old_positions.get(c,[]);i=bisect.bisect_left(values,p)
        if any(0<=j<len(values) and abs(values[j]-p)<=2000 for j in [i-1,i]):near_history.add(find(n))
    library=PREP/'author_large_scale_library_design.xlsx';library_hash=json.loads((PREP/'author_pretraining_library_retrieval_receipt.json').read_text())['sha256'];assert sha(library)==library_hash;near_pretrained=set();literal_pretrained=set();profiles={}
    with zipfile.ZipFile(library) as z:
        strings=[]
        for event,n in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
            if n.tag.endswith('}si'):strings.append(''.join(n.itertext()));n.clear()
        relations={n.attrib['Id']:n.attrib['Target'] for n in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))};book=ET.fromstring(z.read('xl/workbook.xml'))
        for sheet in book.find('m:sheets',NS):
            cell=sheet.attrib['name'].split()[0];target=relations[sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target);header=None;intervals=collections.defaultdict(list);matches=set();count=0;windows=0;unmapped=0
            for event,row in ET.iterparse(z.open(target),events=['end']):
                if not row.tag.endswith('}row'):continue
                values={}
                for c in row:
                    if not c.tag.endswith('}c'):continue
                    n=c.find('m:v',NS);v=n.text if n is not None and n.text is not None else ''
                    if c.attrib.get('t')=='s':v=strings[int(v)]
                    elif c.attrib.get('t')=='inlineStr':v=''.join(c.find('m:is',NS).itertext())
                    values[re.sub('[0-9]','',c.attrib['r'])]=v
                if header is None:
                    if 'name' in values.values() and 'chr.hg38' in values.values():header={v:k for k,v in values.items() if v};seqcol=next(k for v,k in header.items() if v.startswith('230nt sequence'))
                    row.clear();continue
                if not values.get(header['name']):row.clear();continue
                count+=1;sequence=values.get(seqcol,'').upper()
                if len(sequence)==230 and not set(sequence)-set('ACGT'):
                    for offset in range(86):
                        windows+=1
                        if sequence[offset:offset+12] in prefixes:
                            hit=core.get(sequence[offset:offset+145])
                            if hit is not None:matches.add(hit)
                chrom=values.get(header['chr.hg38'],'');chrom='chr'+chrom.removeprefix('chr') if chrom else ''
                try:begin,end=int(float(values.get(header['start.hg38'],''))),int(float(values.get(header['stop.hg38'],'')))
                except ValueError:unmapped+=1;row.clear();continue
                if chrom not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY','chrM'] or not 0<=begin<end:unmapped+=1
                else:intervals[chrom].append((begin,end))
                row.clear()
            index={}
            for chrom,values in intervals.items():
                values.sort();starts=[s for s,e in values];ends=[];maximum=-1
                for s,e in values:maximum=max(maximum,e);ends.append(maximum)
                index[chrom]=(starts,ends)
            near=set()
            for (chrom,p),n in spatial.items():
                p-=1;starts,ends=index.get(chrom,([],[]));j=bisect.bisect_right(starts,p+2000)-1
                if j>=0 and ends[j]>p-2000:near.add(find(n))
            near_pretrained.update(near);literal_pretrained.update(matches);profiles[cell]={'published_library_rows':count,'literal_145bp_windows_examined':windows,'components_near_library':len(near),'components_with_literal_shared_DNA':len(matches),'unmapped_or_noncanonical_library_intervals':unmapped};print('ULIRSCH PRETRAIN GUARD',cell,profiles[cell],flush=True)
    result=[];strict=[]
    for row in rows:
        original=row['original_GRCh37_variant'].removeprefix('chr');original_root=node[original];assert original_root==int(row['component_root']);r=find(original_root);reasons=[]
        if r in near_history:reasons.append('whole_component_within_2kb_of_revealed_history')
        if r in history_hits:reasons.append('whole_component_shares_literal_145bp_with_history')
        if r in near_pretrained:reasons.append('whole_component_near_published_pretraining_library')
        if r in literal_pretrained:reasons.append('whole_component_shares_literal_145bp_with_pretraining_library')
        output={**row,'joint_component_after_all_source_DNA':r,'strict_history_and_pretraining_eligible':not reasons,'exclusion_reasons':'|'.join(reasons),'exact_shared_DNA_and_pretraining_audits_pending':False};result.append(output)
        if not reasons:strict.append(output)
    write_csv(SOURCE/'Ulirsch_complete_exposure_profile.csv',result)
    if strict:write_csv(SOURCE/'Ulirsch_candidate_strict_no_history_or_published_pretraining.csv',strict)
    final={'passed':True,'complete_source_allele_DNA_constructs_joined':len(designs),'additional_shared_DNA_component_joins':joins,'history_profiles':history_profiles,'pretraining_profiles':profiles,
           'strict_candidate_designs':len(strict),'strict_distinct_variants':len({r['variant_GRCh38'] for r in strict}),'strict_components':len({r['joint_component_after_all_source_DNA'] for r in strict}),
           'strict_1Mb_windows':len({(r['chromosome'],(int(r['position_hg38'])-1)//1000000) for r in strict}),'actual_author_pretraining_sequences_all_covered_by_scanned_published_library':True,
           'all_literal_matches_are_full_string_comparisons_both_orientations':True,'outcome_coefficients_or_model_performance_calculated':False,'eligible_for_final_validation_yet':False,
           'pretraining_scope':'Entire published library screened conservatively, including author-held-out folds. Unmapped/noncanonical library coordinates limit wider-neighborhood claims; literal sequences are still scanned.',
           'remaining':['Predeclare and audit coefficient reconstruction from released original counts and replicate/QC rules','Fix a scientifically justified candidate and confidence rules before outcome comparison'],
           'source_library_sha256':library_hash,'profile_sha256':sha(SOURCE/'Ulirsch_complete_exposure_profile.csv'),'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
    write_json(SOURCE/'Ulirsch_complete_history_and_pretraining_guard.json',final);assert_frozen_previous();print(json.dumps(final),flush=True)

if __name__=='__main__':main()
