"""Physical-DNA and spatial audit of unused multi-cell libraries, no model scoring."""
import sys,importlib.util,re,bisect,collections,zipfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'work/sleep_frailty_ml_v9'))
from v9_common import *
sys.path.append(str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
sys.path.append(str(ROOT/'work/sleep_frailty_ml_v7'))
spec=importlib.util.spec_from_file_location('additional_library_reference',ROOT/'work/sleep_frailty_ml_v7/reference.py');ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)

SCREEN=ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen'
DEST=SCREEN/'Siraj_unused_libraries'


def main():
    start_time=time.monotonic();assert_previous_sources();receipt=json.loads((DEST/'targeted_retrieval_receipt.json').read_text())
    for f in receipt['files']:assert sha(f['path'])==f['sha256']
    bylibrary={};wanted=set()
    for library in ['OL27','OL33']:
        per=[]
        for cell in ['HepG2','K562','SKNSH']:
            counts=collections.Counter()
            with gzip.open(DEST/f'{library}_{cell}_emVAR_glm_paired_20221220.out.gz','rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    if r['comb'].endswith('_center_fwd_ref'):counts[r['SNP']]+=1
            per.append({s for s,n in counts.items() if n==1})
        bylibrary[library]=set.intersection(*per);wanted.update('chr'+s for s in bylibrary[library])
    table=BASE/'siraj2026_tables.zip';sequences={};multiplicity=collections.Counter()
    with zipfile.ZipFile(table) as z:
        with z.open('tables/stable2.txt.gz') as handle:
            with gzip.open(handle,'rt') as f:
                for r in csv.DictReader(f,delimiter='\t'):
                    if r['variant'] in wanted:sequences[r['variant']]=r;multiplicity[r['variant']]+=1
    genome_path=V6/'sources/official_hg38.2bit';assert sha(genome_path)==json.loads((V6/'outputs/additional_reference_audit.json').read_text())['sha256'];genome=ref.TwoBit(genome_path)
    lift=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));back=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));verified={};excluded=collections.Counter();all_mapped=collections.defaultdict(set)
    for key in sorted(wanted):
        parts=key.split(':')
        if len(parts)!=4 or len(parts[2])!=1 or len(parts[3])!=1 or set(parts[2]+parts[3])-set('ACGT'):excluded['non_single_nucleotide_variant']+=1;continue
        chrom,p,ra,aa=parts;p=int(p);maps=lift.convert_coordinate(chrom,p-1)
        if len(maps)!=1:excluded['nonunique_build_conversion']+=1;continue
        c,q,strand,_=maps[0];all_mapped[c].add(int(q)+1);reverse=back.convert_coordinate(c,q)
        if len(reverse)!=1 or reverse[0][:3]!=(chrom,p-1,strand):excluded['nonreciprocal_build_conversion']+=1;continue
        if key not in sequences or multiplicity[key]!=1:excluded['no_unique_published_physical_design']+=1;continue
        design=sequences[key]
        if design['duplicate']!='no':excluded['published_duplicate_construct_flag']+=1;continue
        a,b=design['allele1_oligo'].upper(),design['allele2_oligo'].upper()
        if len(a)!=200 or len(b)!=200 or set(a+b)-set('ACGT'):excluded['not_unambiguous_actual_200bp_DNA']+=1;continue
        changes=[j for j,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(changes)!=1 or changes[0] not in [99,100]:excluded['not_actual_centered_single_substitution']+=1;continue
        j=changes[0]
        if (a[j],b[j])!=(ra,aa):excluded['original_labeled_allele_pair_disagrees_with_DNA']+=1;continue
        offset=j if strand=='+' else 199-j;begin=int(q)-offset;truth=genome.sequence(c,begin,begin+200) if begin>=0 else '';amatch=truth==(a if strand=='+' else rc(a));bmatch=truth==(b if strand=='+' else rc(b))
        if amatch==bmatch:excluded['full_actual_reference_not_unique_GRCh38_match']+=1;continue
        rref,ralt=(a,b) if amatch else (b,a);gra,gaa=(rref[j],ralt[j]) if strand=='+' else (rc(rref[j]),rc(ralt[j]))
        verified[key]={'original_variant_GRCh37':key,'variant_GRCh38':f'{c}:{int(q)+1}:{gra}:{gaa}','chromosome':c,'position_hg38':int(q)+1,'reference':rref,'alternate':ralt,
                       'changed_position':j,'construct_strand_in_GRCh38':strand,'GRCh38_reference_is_source_allele1':amatch,'full_physical_reference_verified':True}
    # Components use the entire mapped candidate-source universe before any
    # physical or historical exclusions, preventing a removed bridge leak.
    roots={};root=0
    for chrom,positions in sorted(all_mapped.items()):
        previous=None
        for p in sorted(positions):
            if previous is None or p-previous>2000:root+=1
            roots[(chrom,p)]=root;previous=p
    histories=[V6/'sequence_bank.csv.gz',V7/'sequence_bank.csv.gz',DATA/'extra_sequence_bank.csv.gz',V8/'matoba_full_verified_physical_design.csv',SCREEN/'Tardaguila_full_reference_verified_design.csv'];old=collections.defaultdict(set)
    for path in histories:
        for r in csv_rows(path):old[r['chromosome']].add(int(r['position_hg38']))
    old={c:sorted(ps) for c,ps in old.items()};bad=set()
    for (chrom,p),root in roots.items():
        points=old.get(chrom,[]);i=bisect.bisect_left(points,p)
        if any(0<=j<len(points) and abs(points[j]-p)<=2000 for j in [i-1,i]):bad.add(root)
    rows=[];reserved_roots={roots[(verified[key]['chromosome'],verified[key]['position_hg38'])] for key in verified if key.removeprefix('chr') in bylibrary['OL33']}
    for library in ['OL27','OL33']:
        for key in sorted(verified):
            if key.removeprefix('chr') not in bylibrary[library]:continue
            r=verified[key];component=roots[(r['chromosome'],r['position_hg38'])];spatial_ok=component not in bad
            role='historical_overlap_excluded' if not spatial_ok else 'reserved_candidate' if component in reserved_roots else 'development_candidate'
            rows.append({**r,'library':library,'component_root':component,'after_spatial_exclusion':spatial_ok,'provisional_role':role,'shared_sequence_exclusion_still_required':True})
    write_csv(DEST/'full_reference_verified_with_spatial_roles.csv',rows)
    summary={'physical_audit_passed':True,'primary_new_test_eligibility_finalized':False,'source_union_variant_labels':len(wanted),'unique_full_reference_verified_variants':len(verified),'physical_or_reference_exclusions':dict(excluded),
             'by_library':{lib:dict(collections.Counter(r['provisional_role'] for r in rows if r['library']==lib)) for lib in ['OL27','OL33']},
             'new_source_measurements_in_three_exact_previously_modeled_cell_lines':'HepG2, K562 and SK-N-SH; each physical design has an original centered-forward measurement in all three.',
             'source_sequence_table_sha256':sha(table),'historical_source_sha256':{str(p):sha(p) for p in histories},
             'no_favorable_effect_or_significance_selection':True,'new_model_predictions_or_performance_calculated':False,
             'role_policy':'All spatial components touching OL33 are reserved, including overlapping OL27 designs. No OL33 component enters training. These roles remain provisional until shared-DNA and pretraining exposure checks.',
             'not_independent_publications_or_donors':True,'original_signed_coefficient_convention_for_these_source_files_not_independently_established':True,
             'remaining':['Connect shared variant-containing DNA cores across distant genomic components before final roles','Historical shared-core exclusion','Published/actual pretraining sequence and neighborhood exposure','Source-coefficient/assay QC and protocol lock'],
             'elapsed_seconds':time.monotonic()-start_time}
    write_json(DEST/'physical_and_spatial_audit.json',summary);assert_previous_sources();print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
