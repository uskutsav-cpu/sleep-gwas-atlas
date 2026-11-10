"""Freeze physical SCZ-library designs; original-effect provenance is not model selection."""
import gzip,collections,sys
from common import *
from reference import get_reference
sys.path.insert(0,str(ROOT/'work/stage5_env/lib/python3.13/site-packages'))
from pyliftover import LiftOver
COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
def json_list(value):
    v=json.loads(value)
    if isinstance(v,str):v=json.loads(v)
    assert isinstance(v,list)
    return v
def main():
    assert not (OUT/'model_search_started.json').exists()
    assert not (OUT/'reserved_evaluation_started.json').exists()
    previous=assert_previous_intact();genome=get_reference();source=SOURCES/'IGVFFI6542CYCI.gz'
    with gzip.open(source,'rt') as f:raw=list(csv.DictReader(f,delimiter='\t'))
    refs=collections.defaultdict(list);alts=[];excluded=collections.Counter()
    for r in raw:
        if r['category']!='variant' or r['class']!='test' or r['ref']!='GRCh38':excluded['not_test_variant_GRCh38']+=1;continue
        spdi,allele,kind,positions=map(json_list,[r['SPDI'],r['allele'],r['variant_class'],r['variant_pos']])
        if len(spdi)!=1 or kind!=['SNV'] or len(positions)!=1 or allele not in [['ref'],['alt']]:excluded['not_supported_single_SNV']+=1;continue
        r={**r,'SPDI':spdi[0],'changed_position':int(positions[0])}
        if allele==['ref']:refs[r['SPDI']].append(r)
        else:alts.append(r)
    alt_counts=collections.Counter(r['SPDI'] for r in alts);kept=[]
    for alt in alts:
        key=alt['SPDI']
        if len(refs.get(key,[]))!=1 or alt_counts[key]!=1:excluded['nonunique_or_missing_physical_pair']+=1;continue
        ref=refs[key][0];a,b=ref['sequence'].upper(),alt['sequence'].upper()
        if len(a)!=len(b) or len(a)<51 or len(a)>1024 or set(a+b)-set('ACGT'):excluded['invalid_sequence_length_or_bases']+=1;continue
        changes=[i for i,(x,y) in enumerate(zip(a,b)) if x!=y]
        if len(changes)!=1:excluded['not_exactly_one_substitution']+=1;continue
        j=changes[0];acc,pos,ra,aa=key.split(':');pos=int(pos);c=ref['chr'];start,end=int(ref['start']),int(ref['end']);strand=ref['strand']
        if any(ref[k]!=alt[k] for k in ['chr','start','end','strand']) or j!=ref['changed_position'] or j!=alt['changed_position'] or end-start!=len(a):excluded['physical_pair_annotations_disagree']+=1;continue
        genomic_position=start+j if strand=='+' else end-j-1
        observed=(a[j],b[j]) if strand=='+' else (rc(a[j]),rc(b[j]))
        if genomic_position!=pos or observed!=(ra,aa):excluded['SPDI_position_or_forward_allele_disagreement']+=1;continue
        actual=a if strand=='+' else rc(a)
        if genome.sequence(c,start,end)!=actual:excluded['entire_reference_construct_does_not_match_GRCh38']+=1;continue
        kept.append({'design_index':len(kept),'SPDI':key,'variant_GRCh38':f'{c}:{pos+1}:{ra}:{aa}','chromosome':c,'position_hg38':pos+1,'start':start,'end':end,'construct_strand':strand,'reference_name':ref['name'],'alternate_name':alt['name'],'reference':a,'alternate':b,'changed_position':j})
    assert len(kept)>1000 and len({r['SPDI'] for r in kept})==len(kept)
    design=DATA/'reserved_physical_design.csv'
    with design.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(kept[0]));w.writeheader();w.writerows(kept)
    by_variant={(r['chromosome'],int(r['position_hg38'])-1,r['SPDI'].split(':')[-2],r['SPDI'].split(':')[-1]):r for r in kept}
    lift=LiftOver(str(ROOT/'work/hg19ToHg38.over.chain.gz'));back=LiftOver(str(ROOT/'work/hg38ToHg19.over.chain.gz'));source_audits={};sealed=[];annotations=[]
    specifications=[('neural_progenitor','HNP_primary','IGVFFI7229XIGT','IGVFFI6735WZEI','IGVFFI5341BYYN'),('kidney_derived','HEK293','IGVFFI3812LRGS','IGVFFI7790OLFT','IGVFFI0365SBZK')]
    for cell,cell_line,published_id,canonical_id,bed_id in specifications:
        with gzip.open(SOURCES/(canonical_id+'.gz'),'rt') as f:canonical={r['variant_id']:r for r in csv.DictReader(f,delimiter='\t')}
        with gzip.open(SOURCES/(published_id+'.gz'),'rt') as f:published=list(csv.DictReader(f,delimiter='\t'))
        matched=[];counts=collections.Counter();originals=collections.defaultdict(list)
        for r in published:
            c='chr'+r['CHR'];p=int(r['BP'])-1;maps=lift.convert_coordinate(c,p)
            if len(maps)!=1:counts['original_nonunique_liftover']+=1;continue
            cc,pp,strand,_=maps[0];reverse=back.convert_coordinate(cc,pp)
            if len(reverse)!=1 or reverse[0][:3]!=(c,p,strand):counts['original_nonreciprocal_liftover']+=1;continue
            a,b=r['A1'],r['A2']
            if len(a)!=1 or len(b)!=1 or set(a+b)-set('ACGT'):counts['original_non_SNV']+=1;continue
            if strand=='-':a,b=rc(a),rc(b)
            found=by_variant.get((cc,pp,a,b)) or by_variant.get((cc,pp,b,a))
            if found is None:counts['not_in_verified_physical_design']+=1;continue
            originals[found['SPDI']].append(r)
        for d in kept:
            key=d['SPDI'];original=originals.get(key,[]);current=canonical.get(key)
            if len(original)!=1 or current is None:counts['physical_pair_without_unique_original_and_current_result']+=1;continue
            r=original[0]
            try:
                value,other,p,q=(float(r['MPRA_logFC']),float(current['log2FoldChange']),float(r['MPRA_P']),float(r['MPRA_FDR']))
                assert np.isfinite([value,other,p,q]).all() and 0<=p<=1 and 0<=q<=1
            except (ValueError,AssertionError):counts['missing_or_unassessed_original_effect']+=1;continue
            error=abs(abs(value)-abs(other))
            if error>1e-10:counts['original_vs_current_effect_magnitude_disagrees']+=1;continue
            matched.append(error);sealed.append({'design_index':d['design_index'],'SPDI':key,'cell_group':cell,'cell_line':cell_line,'absolute_log2_effect':abs(value),'original_logFC':value,'original_p':p,'original_q':q,'rsID':r['rsID'],'source':'McAfee2023','original_effect_file':published_id})
            annotations.append({'design_index':d['design_index'],'SPDI':key,'cell_group':cell,'cell_line':cell_line,'original_rsid':r['rsID'],'original_build':'GRCh37,reciprocal mapping verified','original_published_file':published_id,'canonical_replacement_file':canonical_id,'current_genomic_BED_file':bed_id,'original_vs_current_absolute_effect_agrees':True})
        source_audits[cell]={'published_rows':len(published),'canonical_rows':len(canonical),'unique_eligible_original_effects':len(matched),'maximum_original_vs_current_absolute_effect_error':max(matched) if matched else None,'exclusions':dict(counts),'source_files':{x:sha(SOURCES/(x+'.gz')) for x in [published_id,canonical_id,bed_id]}}
    target=DATA/'reserved_original_effects.csv.gz'
    with gzip.open(target,'wt',newline='') as f:w=csv.DictWriter(f,fieldnames=list(sealed[0]));w.writeheader();w.writerows(sealed)
    with (DATA/'reserved_design_annotations.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(annotations[0]));w.writeheader();w.writerows(annotations)
    write_json(OUT/'reserved_design_and_original_source_audit.json',{'study':'McAfee2023,doi10.1016/j.xgen.2023.100404','source_library':'IGVFDS8315BLIH','physical_designs':len(kept),'physical_design_exclusions':dict(excluded),'physical_sequence_and_entire_reference_genome_verified':True,'source_build':'GRCh38 SPDI0based; original publication coordinates GRCh37 reciprocally mapped','original_sources':source_audits,'sealed_effect_rows':len(sealed),'physical_design_sha256':sha(design),'original_effects_file_sha256':sha(target),'design_annotations_sha256':sha(DATA/'reserved_design_annotations.csv'),'source_role_plan_sha256':sha(WORK/'SOURCE_ASSIGNMENT_AMENDMENT.txt'),'V6_immutable_lock_sha256':previous,'models_fitted_or_compared':False,'scope':'Actual original published fitted MPRA log2effect magnitudes, independently matched to current released summaries and physical sequence pairs. This is not a reanalysis of original barcodes. Reprocessed long-count files have no barcode counts and were not independently reconciled to the fitted original coefficients; do not claim barcode-QC reconstruction or reinterpret those raw ratios as the published estimates. Direction/risk-allele claims are outside this magnitude endpoint.'})
    print('RESERVED DESIGN AND ORIGINAL EFFECT PROVENANCE CHECKED',len(kept),'physical pairs;',len(sealed),'cell-specific published effects; NO V7 MODEL PERFORMANCE COMPUTED',flush=True)
if __name__=='__main__':main()
