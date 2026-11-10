"""Original construct/count joins and reciprocal full-genome reference checks.

Counts are used for identities and nonnegative-integer availability only. No
effect estimate, p-value or model prediction determines physical eligibility.
"""
import re,importlib.util,bisect,time
from channel_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
sys.path.append(str(ROOT/'work/sleep_frailty_ml_v7'))
spec=importlib.util.spec_from_file_location('ulirsch_read_only_genome',ROOT/'work/sleep_frailty_ml_v7/reference.py');reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
SOURCE=DEST/'source_candidates'

def main():
    assert_frozen_previous();start=time.monotonic();receipt=json.loads((SOURCE/'Ulirsch_original_retrieval_receipt.json').read_text())
    for entry in receipt['files']:assert sha(entry['path'])==entry['sha256']
    pattern=re.compile(r'^([^:]+):(\d+):([ACGT]+):([ACGT]+)_(Ref|Mut)_Bot_([^_]+)_Top_([^_]+)(_mod)?$');constructs={};by_variant=collections.defaultdict(dict);metadata={}
    with gzip.open(SOURCE/'GSE87711_design.txt.gz','rt') as handle:
        for identifier,sequence in csv.reader(handle,delimiter='\t'):
            assert identifier not in constructs;match=pattern.fullmatch(identifier);assert match,identifier;c,p,ra,aa,allele,bot,top,mod=match.groups();assert bot in ['1/3','1/2','2/3'];assert top in ['1/3','1/2','2/3'];assert len(sequence)==145 and not set(sequence)-set('ACGT')
            constructs[identifier]=sequence;variant=f'chr{c.removeprefix("chr")}:{p}:{ra}:{aa}';key=(variant,bot,top);assert allele not in by_variant[key];by_variant[key][allele]=identifier;metadata[identifier]=(c,p,ra,aa,allele,bot,top,'mod' if mod else 'var')
    barcodes={}
    with gzip.open(SOURCE/'GSE87711_design_barcodes.txt.gz','rt') as handle:
        for identifier,values in csv.reader(handle,delimiter='\t'):
            assert identifier not in barcodes;codes=values.split(',');assert len(codes)==len(set(codes))==14;barcodes[identifier]=codes
    assert constructs.keys()==barcodes.keys();barcode_lengths=collections.Counter(len(b) for values in barcodes.values() for b in values);keys_to_id={v:k for k,v in metadata.items()};assert len(keys_to_id)==len(constructs);count_rows=collections.Counter();readout_names=None;total=0
    with gzip.open(SOURCE/'GSE87711_RBC_MPRA_minP_raw.txt.gz','rt') as handle:
        reader=csv.DictReader(handle,delimiter='\t');readout_names=reader.fieldnames[11:]
        for row in reader:
            key=tuple(row[field] for field in ['chr','pos','ref','alt','type','bot','top','clean']);assert key in keys_to_id,key;count_rows[keys_to_id[key]]+=1;total+=1
            for field in readout_names:value=float(row[field]);assert np.isfinite(value) and value>=0 and value==int(value)
    assert set(count_rows)==set(constructs) and set(count_rows.values())=={14}
    genome_path=V6/'sources/official_hg38.2bit';assert sha(genome_path)==json.loads((V6/'outputs/additional_reference_audit.json').read_text())['sha256'];genome=reference.TwoBit(genome_path);lift=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));back=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));verified=[];excluded=collections.Counter();points=collections.defaultdict(set)
    # Build candidate components from every mapped original variant before exclusions.
    for (variant,bot,top),alleles in sorted(by_variant.items()):
        chrom,p,ra,aa=variant.split(':');p=int(p);mapped=lift.convert_coordinate(chrom,p-1)
        if len(mapped)!=1:excluded['nonunique_build_conversion']+=1;continue
        c,q,strand,_=mapped[0];points[c].add(int(q)+1);reciprocal=back.convert_coordinate(c,q)
        if len(reciprocal)!=1 or reciprocal[0][:3]!=(chrom,p-1,strand):excluded['nonreciprocal_build_conversion']+=1;continue
        if set(alleles)!=set(['Ref','Mut']):excluded['no_unique_actual_Ref_and_Mut_design']+=1;continue
        if len(ra)!=1 or len(aa)!=1:excluded['not_labeled_single_base_variant']+=1;continue
        a,b=[constructs[alleles[name]] for name in ['Ref','Mut']];changes=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(changes)!=1:excluded['actual_pair_not_one_substitution']+=1;continue
        j=changes[0]
        if (a[j],b[j])!=(ra,aa):excluded['labeled_alleles_do_not_match_actual_DNA']+=1;continue
        begin=int(q)-(j if strand=='+' else 144-j);truth=genome.sequence(c,begin,begin+145) if begin>=0 else '';amatch=truth==(a if strand=='+' else rc(a));bmatch=truth==(b if strand=='+' else rc(b))
        if amatch==bmatch:excluded['no_unique_full_145bp_genomic_reference_match']+=1;continue
        rr,alt=(a,b) if amatch else (b,a);gr,ga=(rr[j],alt[j]) if strand=='+' else (rc(rr[j]),rc(alt[j]))
        verified.append({'original_GRCh37_variant':variant,'variant_GRCh38':f'{c}:{int(q)+1}:{gr}:{ga}','window_bot':bot,'window_top':top,'source_reference_construct':alleles['Ref'],'source_alternate_construct':alleles['Mut'],
                         'reference':rr,'alternate':alt,'changed_position':j,'chromosome':c,'position_hg38':int(q)+1,'construct_strand_in_GRCh38':strand,'GRCh38_reference_is_source_Ref':amatch,'full_actual_reference_verified':True})
    roots={};root=0
    for chrom,positions in sorted(points.items()):
        previous=None
        for p in sorted(positions):
            if previous is None or p-previous>2000:root+=1
            roots[(chrom,p)]=root;previous=p
    screen=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen';history=[V6/'sequence_bank.csv.gz',V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',V8/'matoba_full_verified_physical_design.csv',screen/'Tardaguila_full_reference_verified_design.csv',screen/'Siraj_unused_libraries/full_reference_and_shared_core_roles.csv'];past=collections.defaultdict(set)
    for path in history:
        for row in csv_rows(path):past[row['chromosome']].add(int(row['position_hg38']))
    past={c:sorted(values) for c,values in past.items()};bad=set()
    for (chrom,p),root in roots.items():
        values=past.get(chrom,[]);i=bisect.bisect_left(values,p)
        if any(0<=j<len(values) and abs(values[j]-p)<=2000 for j in [i-1,i]):bad.add(root)
    for row in verified:
        component=roots[(row['chromosome'],row['position_hg38'])];row.update(component_root=component,spatial_history_eligible=component not in bad,exact_shared_DNA_and_pretraining_audits_pending=True)
    write_csv(SOURCE/'Ulirsch_full_reference_and_spatial_audit.csv',verified);eligible=[r for r in verified if r['spatial_history_eligible']]
    report={'passed':True,'original_count_rows':total,'original_allele_constructs':len(constructs),'all_constructs_exactly_14_count_rows':True,'count_joins_use_full_alleles_not_ambiguous_byallele_column':True,
            'source_barcode_lengths':dict(barcode_lengths),'metadata_discrepancy':'GEO prose says eleven nucleotide barcodes; actual supplied barcode strings have the lengths reported here. Genomic inserts are actual 145 bp.',
            'original_mapped_variant_window_pairs':len(by_variant),'full_reference_one_SNV_pairs':len(verified),'physical_exclusions':dict(excluded),'after_whole_component_2kb_history_pairs':len(eligible),
            'after_spatial_history_distinct_variants':len({r['variant_GRCh38'] for r in eligible}),'after_spatial_history_regions':len({r['component_root'] for r in eligible}),
            'after_spatial_history_1Mb_windows':len({(r['chromosome'],(r['position_hg38']-1)//1000000) for r in eligible}),
            'sample_columns':readout_names,'historical_files_sha256':{str(path):sha(path) for path in history},'genome_sha256':sha(genome_path),'physical_table_sha256':sha(SOURCE/'Ulirsch_full_reference_and_spatial_audit.csv'),
            'effect_coefficients_or_p_values_selected':False,'model_performance_calculated':False,'eligible_as_final_validation':False,
            'remaining':['Variant-containing shared DNA core/history audit including complete source components','Whole author pretraining library and actual sequence coverage audit','Count-to-coefficient/replicate/noise reconstruction and protocol freeze','Frozen candidate confidence and predictions before any outcome comparison'],
            'not_human_insomnia_frailty_or_genetic_score_validation':True,'code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
    write_json(SOURCE/'Ulirsch_physical_and_source_count_audit.json',report);assert_frozen_previous();print(json.dumps(report),flush=True)

if __name__=='__main__':main()
