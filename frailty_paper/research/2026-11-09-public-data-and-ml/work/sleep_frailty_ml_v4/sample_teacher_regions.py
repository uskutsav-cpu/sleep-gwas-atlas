"""Select score-independent novel genomic regions before fetching teacher values."""
import gzip,collections,bisect
from common import *
def main():
    assert not (OUT/'teacher_region_plan.json').exists();sizes=dict(line.split() for line in (DATA/'hg38.chrom.sizes').read_text().splitlines());sizes={c:int(n) for c,n in sizes.items()};rng=np.random.default_rng(41207);blocked=collections.defaultdict(list);sources=[]
    for p in [ROOT/'work/sleep_frailty_ml_v3/salomon_manifest_local.csv',ROOT/'work/sleep_frailty_ml_v3/siraj_manifest_local.csv',ROOT/'work/sleep_frailty_ml_v2/external_verified_all_local.csv']:
        sources.append({'file':str(p),'sha256':sha(p)})
        for r in csv.DictReader(p.open()):blocked['chr'+r['chromosome'].removeprefix('chr')].append(int(r['position_hg38']))
    project=V3/'project_raw_inputs.json';sources.append({'file':str(project),'sha256':sha(project)})
    for r in json.loads(project.read_text()):blocked['chr'+str(r['chromosome']).removeprefix('chr')].append(int(r['position_hg38']))
    for c in blocked:blocked[c].sort()
    reference=json.loads((V3/'reference_chromatin_sources.json').read_text());peaks=[];peak_receipts=[]
    for r in reference['files']:
        if r['accession'] not in ['ENCFF752OZB','GSE113480']:continue
        if r['accession']=='GSE113480':continue # These arehg19; do not silently usehg38 coordinates.
        path=ROOT/r['local_file'];assert sha(path)==r['sha256'];peak_receipts.append({'file':str(path),'sha256':sha(path),'assembly':'GRCh38'})
        with gzip.open(path,'rt') as f:
            for line in f:
                x=line.split('\t')
                if x[0] in ['chr'+str(c) for c in range(1,23)]:peaks.append((x[0],(int(x[1])+int(x[2]))//2+1))
    assert peaks;chosen=collections.defaultdict(list);rows=[]
    def allowed(c,p):
        if not 500<p<sizes[c]-500:return False
        for array,distance in [(blocked[c],2150),(chosen[c],4096)]:
            j=bisect.bisect_left(array,p-distance)
            if j<len(array) and array[j]<=p+distance:return False
        return True
    # Pilot subsets are fixed before any score. The full reserved teacher set
    # remains unused by pilot training/checkpoint selection.
    for split,chroms,count,pilot in [('train',list(range(6,23)),1120,280),('validation',[4,5],240,60),('reserved_teacher_test',[1,2,3],240,0)]:
        names=['chr'+str(c) for c in chroms];pool=[(c,p) for c,p in peaks if c in names];lengths=np.array([sizes[c] for c in names],dtype='float64');probs=lengths/lengths.sum();n=0
        while n<count:
            stratum='SKNSH_accessible' if n%2==0 else 'genome_background'
            if n%2==0:c,p=pool[int(rng.integers(len(pool)))]
            else:c=str(rng.choice(names,p=probs));p=int(rng.integers(500,sizes[c]-500))
            if not allowed(c,p):continue
            bisect.insort(chosen[c],p);rows.append({'anchor_id':f'{split}_{n:04d}','chromosome':c,'position_hg38':p,'split':split,'sampling_stratum':stratum,'pilot':n<pilot,'window_offsets':[-15,15]});n+=1
    path=DATA/'teacher_regions.json';write_json(path,rows);write_json(OUT/'teacher_region_plan.json',{'seed':41207,'anchor_manifest_sha256':sha(path),'total_anchors':len(rows),'pilot_training_anchors':280,'pilot_validation_anchors':60,'pilot_test_anchors':0,'expanded_train_anchors':1120,'expanded_validation_anchors':240,'reserved_teacher_test_anchors':240,'synthetic_variants_per_region':93,'allele_rule':'Every genomic-reference-to-other-base SNV at31consecutive central positions; native270base paired contexts','distance_from_any_previous_assay_or_project_variant_bp':2150,'anchor_minimum_separation_bp':4096,'sampling':'Equal SK-N-SH called-accessible peak centers and whole-genome background; backgrounds are not guaranteed inaccessible. Regions chosen without teacher scores or lab outcomes. OnlyGRCh38peak file used;hg19profiles not silently reused.','source_exclusions':sources,'peak_source_receipts':peak_receipts,'independence':'Reserved teacher predictions are not measured biological truth. Prior V3 lab tests have already been observed and cannot be claimed fresh for V4.','training_use':'Only designated permissive AVI SNV raw/PHRED scores; detailed API outputs excluded.','old_V3_registry_sha256':sha(V3/'final_model_registry.json'),'old_V3_prediction_sha256':sha(V3/'reserved_test_predictions.npz')});print('Selected1600novel teacher regions; pilot340anchors; reserved240untouched',flush=True)
if __name__=='__main__':main()
