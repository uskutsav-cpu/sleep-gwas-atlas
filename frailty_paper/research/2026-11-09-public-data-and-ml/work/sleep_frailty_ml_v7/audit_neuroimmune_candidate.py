"""Audit new original paired sequences/summary effects; no model predictions."""
import sys,gzip,collections,bisect,shutil
from common import *
from reference import get_reference
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def numeric(value):
    try:v=float(value);return v if np.isfinite(v) else None
    except (ValueError,TypeError):return None
def main():
    receipt=json.loads((OUT/'neuroimmune_raw_extract_receipt.json').read_text());folder=DATA/'neuroimmune_source_tables'
    for r in receipt['sheets'].values():assert sha(r['raw_table'])==r['sha256']
    old_audit=OUT/'neuroimmune_candidate_source_audit.json'
    if old_audit.exists() and 'exact_duplicate_design_rows_removed' not in json.loads(old_audit.read_text()):
        archive=DATA/'neuroimmune_initial_duplicate_diagnostic';archive.mkdir(exist_ok=True)
        for path in [old_audit,DATA/'neuroimmune_physical_design_verified.csv',DATA/'neuroimmune_verified_candidate_effects.csv']:
            if path.exists():shutil.copyfile(path,archive/path.name)
        write_json(archive/'archive_receipt.json',{'files':{p.name:sha(p) for p in archive.iterdir() if p.is_file()},'reason':'Initial diagnostic excluded192exactly duplicated biallelic groups. Source inspection established384identical duplicate design rows and no conflicting repeated target values. Correct by deduplicating exact constructs before eligibility, preserving true multiallelic/missing exclusions. No model performance used.'})
    genome=get_reference();chain=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));groups=collections.defaultdict(list);excluded=collections.Counter();accepted=[];mapping=collections.Counter();duplicates=0;unique_constructs=0
    for r in csv.DictReader((folder/'design.csv').open()):
        parts=r['SNP'].split(':')
        if len(parts)!=4 or not parts[0].startswith('chr') or not parts[2].isdigit() or len(parts[3])!=1 or parts[3] not in 'ACGT':excluded['control_or_unsupported_ID']+=1;continue
        groups[':'.join(parts[:3])].append((r,parts))
    for identity,items in groups.items():
        before=len(items);items=list({(r['SNP'],r['Sequence']):(r,parts) for r,parts in items}.values());duplicates+=before-len(items);unique_constructs+=len(items)
        if len(items)!=2:excluded['not_exactly_two_constructs']+=1;continue
        a,b=[r['Sequence'].upper() for r,_ in items]
        if len(a)!=len(b) or set(a+b)-set('ACGT') or not 51<=len(a)<=1024:excluded['unsupported_sequence_pair']+=1;continue
        changes=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(changes)!=1:excluded['not_a_single_substitution']+=1;continue
        pi=changes[0];c,rsid,p,_=items[0][1];p=int(p);c={'chr23':'chrX','chr24':'chrY'}.get(c,c);candidates={}
        for offset in [0,1]:
            source_pos=p+offset;locations=[(c,source_pos,'direct_hg38',offset)];lift=chain.convert_coordinate(c,source_pos-1)
            if lift and len(lift)==1:locations.append((lift[0][0],int(lift[0][1])+1,'hg19_lift_to_hg38',offset))
            for chrom,pos,build,off in locations:
                for strand in ['+','-']:
                    aa,bb=(a,b) if strand=='+' else (rc(a),rc(b));index=pi if strand=='+' else len(a)-1-pi;start=pos-1-index
                    try:reference=genome.sequence(chrom,start,start+len(aa))
                    except (KeyError,AssertionError):continue
                    if reference==aa:refseq,altseq,ra,al=a,b,aa[index],bb[index];refid,altid=items[0][0]['SNP'],items[1][0]['SNP']
                    elif reference==bb:refseq,altseq,ra,al=b,a,bb[index],aa[index];refid,altid=items[1][0]['SNP'],items[0][0]['SNP']
                    else:continue
                    key=f'{chrom}:{pos}:{ra}:{al}';record={'source_ID':identity,'rsid':rsid,'variant_GRCh38':key,'chromosome':chrom,'position_hg38':pos,'reference':refseq,'alternate':altseq,'changed_position':pi,'construct_strand':strand,'source_ref_construct':refid,'source_alt_construct':altid,'source_build_interpretation':build,'source_position_offset':off};candidates.setdefault(key,[]).append(record)
        if len(candidates)!=1:excluded['no_unique_full_sequence_genomic_mapping']+=1;continue
        choices=next(iter(candidates.values()));choices.sort(key=lambda r:(r['source_position_offset'],r['source_build_interpretation']!='direct_hg38'));r=choices[0];r['build_alternatives_same_exact_hg38_variant']=';'.join(sorted({x['source_build_interpretation'] for x in choices}));mapping[r['source_build_interpretation']+'|offset'+str(r['source_position_offset'])]+=1;accepted.append(r)
    assert accepted
    # Audit the complete source against all current development inputs and
    # the primary reserved study before choosing any additional test subset.
    points=collections.defaultdict(set);old_cores=set()
    def add_old(r):
        points[r['chromosome']].add(int(r['position_hg38']));p=int(r['changed_position'])
        for s in [r['reference'],r['alternate']]:
            start=max(0,min(p-75,len(s)-150));core=s[start:start+150];old_cores.add(min(core,rc(core)))
    with gzip.open(DATA/'sequence_bank.csv.gz','rt') as f:
        for r in csv.DictReader(f):add_old(r)
    for r in csv.DictReader((DATA/'reserved_physical_design.csv').open()):add_old(r)
    points={c:sorted(v) for c,v in points.items()};chroms=collections.defaultdict(list)
    for i,r in enumerate(accepted):chroms[r['chromosome']].append(i)
    for c,ids in chroms.items():
        ids.sort(key=lambda i:accepted[i]['position_hg38']);chunks=[]
        for i in ids:
            if not chunks or accepted[i]['position_hg38']-accepted[chunks[-1][-1]]['position_hg38']>2000:chunks.append([])
            chunks[-1].append(i)
        for chunk in chunks:
            reasons=set();ps=points.get(c,[])
            for i in chunk:
                r=accepted[i];p=r['position_hg38'];j=bisect.bisect_left(ps,p)
                if any(0<=k<len(ps) and abs(ps[k]-p)<=2000 for k in [j-1,j]):reasons.add('within2kb_of_development_or_primary_test')
                pi=r['changed_position']
                for seq in [r['reference'],r['alternate']]:
                    if any(min(seq[j:j+150],rc(seq[j:j+150])) in old_cores for j in range(max(0,pi-149),min(pi,len(seq)-150)+1)):reasons.add('variant_containing150bp_matches_a_prior_central_context')
            for i in chunk:accepted[i]['additional_test_overlap_exclusion']=';'.join(sorted(reasons))
    for i,r in enumerate(accepted):r['design_index']=i
    path=DATA/'neuroimmune_physical_design_verified.csv'
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(accepted[0]));w.writeheader();w.writerows(accepted)
    eligible={r['source_ID']:r for r in accepted if not r['additional_test_overlap_exclusion']};effect_rows=[];counts={};effect_exclusions={}
    for label in ['baseline','IL6','IFNa']:
        byid=collections.defaultdict(list)
        for r in csv.DictReader((folder/(label+'.csv')).open()):byid[r['ID']].append(r)
        bad=collections.Counter();all_verified=0;before_overlap=0
        for design in accepted:
            values=byid.get(design['source_ID'],[]);source_duplicates=collections.defaultdict(list)
            for row in values:source_duplicates[row['CRS']].append(row)
            if any(len({r.get('emVar.logFC','') for r in rows})>1 for rows in source_duplicates.values()):bad['repeated_construct_target_values_disagree']+=1;continue
            byconstruct={key:rows[0] for key,rows in source_duplicates.items()}
            if not all(design[k] in byconstruct for k in ['source_ref_construct','source_alt_construct']):bad['both_physical_alleles_not_in_postQC_summary']+=1;continue
            rr=byconstruct[design['source_ref_construct']];aa=byconstruct[design['source_alt_construct']];a=numeric(rr.get('emVar.logFC'));b=numeric(aa.get('emVar.logFC'))
            if a is None or b is None:bad['missing_measured_allelic_effect_not_a_negative']+=1;continue
            if abs(abs(a)-abs(b))>1e-8:bad['source_allele_rows_disagree_on_effect_magnitude']+=1;continue
            all_verified+=1;before_overlap+=1
            if design['source_ID'] not in eligible:continue
            effect_rows.append({'design_index':design['design_index'],'source_ID':design['source_ID'],'rsid':design['rsid'],'variant_GRCh38':design['variant_GRCh38'],'condition':label,'absolute_log2_effect':abs(a),'source_logFC_ref_row':a,'source_logFC_alt_row':b,'source_q':rr.get('emVar.qvalue',''),'source_significance_label':rr.get('emVar.group',''),'source_table':str(folder/(label+'.csv'))})
        counts[label]={'source_postQC_rows':sum(len(v) for v in byid.values()),'source_postQC_unique_constructs':len({r['CRS'] for v in byid.values() for r in v}),'verified_biallelic_measured_effects_before_overlap':all_verified,'retained_after_genomic_overlap_exclusion':sum(r['condition']==label for r in effect_rows)};effect_exclusions[label]=dict(bad)
    effects=DATA/'neuroimmune_verified_candidate_effects.csv'
    if effect_rows:
        with effects.open('w') as f:w=csv.DictWriter(f,fieldnames=list(effect_rows[0]));w.writeheader();w.writerows(effect_rows)
    study=json.loads((SOURCES/'OpenMPRA_S25D8RVYD_metadata.json').read_text());report={'study':'Retallick-Townsley et al.2026,doi10.1038/s41467-026-77114-x','source_variant_groups_before_genome_checks':len(groups),'exact_duplicate_design_rows_removed':duplicates,'unique_allele_constructs_in_source':unique_constructs,'full_reference_genome_verified_pairs':len(accepted),'design_exclusions':dict(excluded),'coordinate_mappings':dict(mapping),'physical_pairs_outside_all_current_development_and_primary_test_neighborhoods':len(eligible),'published_effect_measurement_inventory':counts,'effect_exclusions':effect_exclusions,'source_design_sha256':sha(path),'verified_effect_file_sha256':sha(effects) if effects.exists() else None,'publisher_and_archive_workbook_identical':True,'raw_barcode_QC_independently_reconstructed':False,'postQC_source_selection':'Deduplicate identical construct rows; preserve complete biallelic physical pairs and finite matching absolute original fitted effects. Do not filter by activity, significance or positive calls; missing effects never become zero.','source_limitations':['9244annotated allele-design rows reconcile to8860unique constructs after removing384duplicates. Some source IDs have three alleles or one available construct and remain unsupported by this biallelic endpoint.','The experiment has two donors and repeated cultures, not thousands of independent people.','IL6 and IFNa are stimulus contexts absent from the present model inputs; no condition-specific model accuracy is established.','Sequence/build and neighborhood checks do not prove participant, ancestry-LD or gene-wide independence.','This supplies laboratory effects only, not independent human sleep/frailty genotype-outcome score validation.'],'status':'Audited source availability; not used in current V7 training, selection or calibration. No candidate model performance calculated.','current_primary_test_evaluation_has_started':(OUT/'reserved_evaluation_started.json').exists()};write_json(OUT/'neuroimmune_candidate_source_audit.json',report);print('ADDITIONAL NEURAL SOURCE AUDITED',json.dumps(report),flush=True)
if __name__=='__main__':main()
