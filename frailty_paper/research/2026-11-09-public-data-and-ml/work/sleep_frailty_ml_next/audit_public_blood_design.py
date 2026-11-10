"""Outcome-blind physical-reference and historical-overlap audit of new blood tiles."""
import sys,importlib.util,collections,bisect,re,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
sys.path.append(str(ROOT/'work/sleep_frailty_ml_v7'))
spec=importlib.util.spec_from_file_location('read_only_genome_reader',ROOT/'work/sleep_frailty_ml_v7/reference.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)

DEST=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'


def main():
    start_time=time.monotonic();assert_previous_sources()
    receipt=json.loads((DEST/'targeted_archive_retrieval_receipt.json').read_text())
    for f in receipt['files']:assert sha(f['path'])==f['sha256']
    genome_path=V6/'sources/official_hg38.2bit';genome_receipt=json.loads((V6/'outputs/additional_reference_audit.json').read_text());assert sha(genome_path)==genome_receipt['sha256']
    genome=ref.TwoBit(genome_path);lift=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));back=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'))
    with (DEST/'Tardaguila2026_sequence_information.tsv').open() as f:raw=list(csv.DictReader(f,delimiter='\t'))
    verified=[];excluded=collections.Counter();source_keys=set();all_source_positions=collections.defaultdict(set)
    for r in raw:
        key=(r['SNP'],r['comparison'],r['window']);assert key not in source_keys;source_keys.add(key)
        original=re.findall(r'(chr(?:\d+|X|Y))_(\d+)_([ACGT]+)_([ACGT]+)',r['SNP'])
        for chrom,position,_,_ in original:
            mapped=lift.convert_coordinate(chrom,int(position)-1)
            if len(mapped)==1:all_source_positions[mapped[0][0]].add(int(mapped[0][1])+1)
        a,b=r['sequence_REF'].upper(),r['sequence_ALT'].upper()
        if len(a)!=148 or len(b)!=148:excluded['not_equal_actual_148bp']+=1;continue
        if set(a+b)-set('ACGT'):excluded['ambiguous_DNA']+=1;continue
        changes=[j for j,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(changes)!=1:excluded[f'actual_substitution_count_{len(changes)}']+=1;continue
        j=changes[0]
        if r['SNP_rel_pos']!=str(j):excluded['source_changed_position_annotation_disagrees_with_actual_DNA']+=1;continue
        # These author start/end coordinates are 1-based inclusive, verified
        # against the separately labeled SNP position and the actual change.
        if int(r['end'])-int(r['start'])+1!=148:excluded['source_interval_length_disagrees_with_actual_DNA']+=1;continue
        candidates=[(chrom,int(position),ra,aa) for chrom,position,ra,aa in original if chrom==r['chr'] and int(position)==int(r['start'])+j and ra==a[j] and aa==b[j]]
        if len(candidates)!=1:excluded['source_position_or_allele_not_unique_match']+=1;continue
        chrom,pos,ra,aa=candidates[0];maps=lift.convert_coordinate(chrom,pos-1)
        if len(maps)!=1:excluded['nonunique_build_conversion']+=1;continue
        c,q,strand,_=maps[0];reverse=back.convert_coordinate(c,q)
        if len(reverse)!=1 or reverse[0][:3]!=(chrom,pos-1,strand):excluded['nonreciprocal_build_conversion']+=1;continue
        offset=j if strand=='+' else 147-j;start=int(q)-offset
        expected=a if strand=='+' else rc(a)
        if start<0 or genome.sequence(c,start,start+148)!=expected:excluded['full_actual_reference_mismatch_to_GRCh38']+=1;continue
        gra,gaa=(ra,aa) if strand=='+' else (rc(ra),rc(aa))
        verified.append({'source_SNP':r['SNP'],'comparison':r['comparison'],'window':r['window'],'source_ID':r['ID'],'source_KEY':r['KEY'],'source_Project':r['Project'],
                         'original_GRCh37_variant':f'{chrom}:{pos}:{ra}:{aa}','variant_GRCh38':f'{c}:{int(q)+1}:{gra}:{gaa}','chromosome':c,'position_hg38':int(q)+1,
                         'reference':a,'alternate':b,'changed_position':j,'construct_strand_in_GRCh38':strand,'physical_reference_length':148,'full_physical_reference_verified':True})
    assert len(verified)>100,(len(verified),excluded)
    # Connect the complete candidate-source point universe before exclusions.
    # This includes positions belonging to unsupported multi-change constructs.
    point_root={};root_members={};root=0
    for chrom,positions in sorted(all_source_positions.items()):
        previous=None
        for pos in sorted(positions):
            if previous is None or pos-previous>2000:root+=1
            point_root[(chrom,pos)]=root;root_members.setdefault(root,[]).append((chrom,pos));previous=pos
    for r in verified:r['component_root']=point_root[(r['chromosome'],r['position_hg38'])]
    histories=[V6/'sequence_bank.csv.gz',V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',V8/'matoba_full_verified_physical_design.csv']
    historical_positions=collections.defaultdict(set);candidate_sequences=collections.defaultdict(set)
    for r in verified:
        for s in [r['reference'],r['alternate']]:
            candidate_sequences[s].add(r['component_root']);candidate_sequences[rc(s)].add(r['component_root'])
    prefixes={s[:12] for s in candidate_sequences};shared_roots=set();shared_matches=collections.Counter();history_stats=[]
    for path in histories:
        windows=0;rows=0
        for r in csv_rows(path):
            rows+=1;historical_positions[r['chromosome']].add(int(r['position_hg38']));changed=int(r['changed_position']);length=len(r['reference'])
            if length<148:continue
            low=max(0,changed-147);high=min(changed,length-148)
            for field in ['reference','alternate']:
                sequence=r[field]
                for k in range(low,high+1):
                    windows+=1
                    if sequence[k:k+12] in prefixes:
                        hit=candidate_sequences.get(sequence[k:k+148])
                        if hit:shared_roots.update(hit);shared_matches[str(path)]+=1
        history_stats.append({'path':str(path),'sha256':sha(path),'physical_rows':rows,'all_variant_containing_148bp_windows_examined':windows})
        print('HISTORICAL 148bp SCAN',path.parent.name,rows,windows,flush=True)
    historical_positions={c:sorted(ps) for c,ps in historical_positions.items()};near_roots=set()
    for root,points in root_members.items():
        for chrom,pos in points:
            ps=historical_positions.get(chrom,[]);i=bisect.bisect_left(ps,pos)
            if any(0<=j<len(ps) and abs(ps[j]-pos)<=2000 for j in [i-1,i]):near_roots.add(root);break
    accepted=[];decisions=[]
    for r in verified:
        root=r['component_root'];why=[]
        if root in near_roots:why.append('whole_component_within_2kb_of_revealed_history')
        if root in shared_roots:why.append('whole_component_shares_actual_148bp_sequence_with_history')
        row={**r,'eligible_after_historical_exclusion':not why,'exclusion_reason':'|'.join(why)};decisions.append(row)
        if not why:accepted.append(row)
    write_csv(DEST/'Tardaguila_full_reference_verified_design.csv',decisions)
    if accepted:write_csv(DEST/'Tardaguila_reserved_no_historical_overlap.csv',accepted)
    duplicates=collections.Counter((r['variant_GRCh38'],r['reference'],r['alternate']) for r in accepted)
    summary={'passed':True,'source_design_rows':len(raw),'full_reference_verified_one_substitution_rows':len(verified),'physical_or_build_exclusions':dict(excluded),
             'historical_exclusion_counts':dict(collections.Counter(r['exclusion_reason'] for r in decisions if r['exclusion_reason'])),
             'accepted_tile_comparisons':len(accepted),'accepted_distinct_physical_pairs':len(duplicates),'accepted_distinct_GRCh38_variants':len({r['variant_GRCh38'] for r in accepted}),
             'accepted_connected_2kb_regions':len({r['component_root'] for r in accepted}),'accepted_1Mb_windows':len({(r['chromosome'],(r['position_hg38']-1)//1000000) for r in accepted}),
             'duplicate_physical_pair_aliases':sum(v-1 for v in duplicates.values()),'shared_historical_sequence_matches':dict(shared_matches),'historical_scans':history_stats,
             'source_start_end_convention':'1-based inclusive, reconciled to labeled original SNP and actual changed position; no guessed coordinate offset.',
             'source_build':'GRCh37 to GRCh38 with unique reciprocal conversion and full 148bp physical-reference match.',
             'target_outcomes_or_significance_values_used':False,'new_model_fitted_or_performance_calculated':False,
             'pretraining_overlap_excluded_for_future_model':False,'remaining':['Reconcile original source model coefficients/comparators and raw assay QC','Audit any future pretrained-model sequence exposure','Freeze next-model and validation/calibration roles before using outcomes'],
             'scope':'Candidate laboratory validation resource, not independent human genetic-score validation. Historical near/shared-region exclusion is conservative and does not prove donor or experimental independence.',
             'elapsed_seconds':time.monotonic()-start_time}
    write_json(DEST/'new_source_physical_and_historical_overlap_audit.json',summary);assert_previous_sources();print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
