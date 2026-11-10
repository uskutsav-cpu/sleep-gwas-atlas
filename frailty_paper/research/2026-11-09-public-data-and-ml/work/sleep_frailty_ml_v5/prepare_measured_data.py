"""Unify measured examples; reserve the final study before any model fitting."""
import gzip,collections,bisect,re,math
from common import *
from prepare_blind_test_design import fasta,rc
V3W=ROOT/'work/sleep_frailty_ml_v3';V3O=ROOT/'outputs/sleep_frailty_ai_v3'
LUT=np.zeros(256,dtype='uint8');LUT[[ord(c) for c in 'ACGT']]=[7,8,9,10]
def tokens(s):return LUT[np.frombuffer(s.encode(),dtype='uint8')]
def decode(x):return np.frombuffer(b'ACGT',dtype='uint8')[x-7].tobytes().decode()
def core_hashes(a,b,p,width=150):
    lo=max(0,p-width+1);hi=min(p,len(a)-width)
    for seq in [a,b]:
        for i in range(lo,hi+1):
            part=seq[i:i+width];yield hashlib.sha256(min(part,rc(part)).encode()).digest()
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists()
    assert not (OUT/'model_search_started.json').exists()
    test=list(csv.DictReader((DATA/'kosicki_blind_design.csv').open()))
    assert sha(DATA/'kosicki_blind_design.csv')==json.loads((OUT/'kosicki_blind_design_audit.json').read_text())['design_sha256']
    tf=fasta(SOURCES/'kosicki_ENCFF752NFF.fasta.gz');test_pos=collections.defaultdict(set);test_alleles=set();test_cores=set()
    for r in test:
        test_pos[r['chromosome']].add(int(r['position_hg38']));test_alleles.add(r['variant_GRCh38'])
        test_cores.update(core_hashes(tf[r['reference_name']],tf[r['name']],int(r['changed_position'])))
    test_pos={k:sorted(v) for k,v in test_pos.items()};print('Reserved design exact150core hashes',len(test_cores),flush=True)
    exclusions=collections.Counter();bank=[];observations=[];siraj_map={};orig_receipt=json.loads((V3O/'feature_receipt.json').read_text())
    def eligible(row,a,b,p):
        c=row['chromosome'];pos=int(row['position_hg38']);points=test_pos.get(c,[]);j=bisect.bisect_left(points,pos)
        if row['variant_GRCh38'] in test_alleles:return 'exact_reserved_allele'
        if any(0<=q<len(points) and abs(points[q]-pos)<=1000 for q in [j-1,j]):return 'within1kb_reserved_mutation'
        if any(h in test_cores for h in core_hashes(a,b,p)):return 'reserved_exact150bp_variant_core'
        return None
    for stem,cell,group in [('salomon','NGN2_induced_neurons','neuronal_induced'),('siraj','SK-N-SH','neuronal_tumor')]:
        rows=list(csv.DictReader((V3W/(stem+'_manifest_local.csv')).open()));sp=V3W/(stem+'_sequences_local.npz');fp=V3W/(stem+'_features_local.npz');prior=next(r for r in orig_receipt['datasets'] if r['study']==stem)
        assert sha(sp)==prior['input_sha256'] and sha(fp)==prior['features_sha256']
        with np.load(sp) as z:arrays={k:z[k] for k in ['ref','alt','position','label','activity_mean','absolute_effect']}
        for i,r in enumerate(rows):
            a=decode(arrays['ref'][i]);b=decode(arrays['alt'][i]);p=int(arrays['position'][i]);r={**r,'chromosome':'chr'+r['chromosome']};reason=eligible(r,a,b,p)
            if reason:exclusions[stem+'_'+reason]+=1;continue
            idx=len(bank);bank.append({'bank_id':idx,'variant_GRCh38':r['variant_GRCh38'],'chromosome':r['chromosome'],'position_hg38':int(r['position_hg38']),'reference':a,'alternate':b,'changed_position':p,'sequence_source':stem,'original_row_index':i,'parent_locus':''})
            if stem=='siraj':siraj_map[r['variant']]=idx
            observations.append({'bank_id':idx,'study':'Salomon' if stem=='salomon' else 'Siraj','experiment':r['study'],'cell':cell,'cell_group':group,'label':int(arrays['label'][i]),'absolute_effect':float(arrays['absolute_effect'][i]),'activity_mean':float(arrays['activity_mean'][i]),'active_screen':int(r['active_screen']),'pvalue':'','tags':'','parent_locus':''})
        print('Existing measured source',stem,'accepted bank',len(bank),'observations',len(observations),flush=True)
    source_cell_receipt=json.loads((SOURCES/'siraj_additional_cell_receipt.json').read_text());cell_exclusions=collections.Counter()
    for item in source_cell_receipt['files']:
        path=Path(item['file']);assert sha(path)==item['sha256'];cell='HepG2' if '_HepG2_' in path.name else 'HCT116' if '_HCT116_' in path.name else None;assert cell
        group='liver' if cell=='HepG2' else 'colon'
        with gzip.open(path,'rt') as f:raw=[r for r in csv.DictReader(f,delimiter='\t') if r['comb'].endswith('_center_fwd_ref')]
        counts=collections.Counter(r['SNP'] for r in raw)
        for r in raw:
            if counts[r['SNP']]!=1:cell_exclusions[cell+'_duplicate_variant']+=1;continue
            key='chr'+r['SNP']
            if key not in siraj_map:cell_exclusions[cell+'_not_in_eligible_reference_checked_neural_construct_universe']+=1;continue
            try:
                aa,bb,ap,bp,skew=(float(r[k]) for k in ['A_log2FC','B_log2FC','A_logPadj_BF','B_logPadj_BF','Log2Skew'])
                assert np.isfinite([aa,bb,skew]).all() and all(not math.isnan(v) and v>=0 for v in [ap,bp])
                active=(abs(aa)>=1 and ap>=2) or (abs(bb)>=1 and bp>=2)
                if active:q=float(r['Skew_logFDR_act']);assert not math.isnan(q) and q>=0;label=int(q>=1)
                else:label=0
            except (ValueError,AssertionError):cell_exclusions[cell+'_unassessable_activity_or_skew']+=1;continue
            observations.append({'bank_id':siraj_map[key],'study':'Siraj','experiment':'OL42_'+cell,'cell':cell,'cell_group':group,'label':label,'absolute_effect':abs(skew),'activity_mean':.5*(aa+bb),'active_screen':int(active),'pvalue':'','tags':'','parent_locus':''})
        print('Added measured cell',cell,'observations',len(observations),flush=True)
    kp={r['variant_GRCh38']:r for r in csv.DictReader((DATA/'kircher_sequence_pairs.csv').open())};km={}
    for key,r in kp.items():
        a,b,p=r['reference'],r['alternate'],int(r['changed_position']);reason=eligible(r,a,b,p)
        if reason:exclusions['kircher_'+reason]+=1;continue
        idx=len(bank);km[key]=idx;bank.append({'bank_id':idx,'variant_GRCh38':key,'chromosome':r['chromosome'],'position_hg38':int(r['position_hg38']),'reference':a,'alternate':b,'changed_position':p,'sequence_source':'kircher','original_row_index':-1,'parent_locus':r['parent_locus']})
    for r in csv.DictReader((DATA/'kircher_observations.csv').open()):
        if r['variant_GRCh38'] not in km:continue
        observations.append({'bank_id':km[r['variant_GRCh38']],**{k:r[k] for k in ['study','experiment','cell','cell_group','label','absolute_effect','activity_mean','active_screen','pvalue','tags','parent_locus']}})
    del test_cores,tf
    # Loci are connected by any <=2kb adjacency across all studies/cells;
    # common reference cores also join distant duplicate contexts.
    parent=list(range(len(bank)))
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    bychr=collections.defaultdict(list);seen_core={}
    for r in bank:
        i=r['bank_id'];bychr[r['chromosome']].append(i);p=r['changed_position'];a=r['reference'];start=max(0,min(p-75,len(a)-151));s=a[start:start+151];h=hashlib.sha256(min(s,rc(s)).encode()).digest()
        if h in seen_core:union(i,seen_core[h])
        else:seen_core[h]=i
    for c,ids in bychr.items():
        ids.sort(key=lambda i:bank[i]['position_hg38'])
        for a,b in zip(ids,ids[1:]):
            if bank[b]['position_hg38']-bank[a]['position_hg38']<=2000:union(a,b)
    roots=collections.defaultdict(list)
    for r in bank:roots[find(r['bank_id'])].append(r['bank_id'])
    for ids in roots.values():
        coords=sorted((bank[i]['chromosome'],bank[i]['position_hg38']) for i in ids);name='|'.join(f'{c}:{p}' for c,p in coords);h=hashlib.sha256(('V5_locus_split_20261107|'+name).encode()).hexdigest();fold=int(h[:8],16)%10;split='calibration' if fold<2 else 'validation' if fold<4 else 'train';gid=h[:20]
        for i in ids:bank[i]['locus_group']=gid;bank[i]['split']=split
    for r in observations:r['locus_group']=bank[int(r['bank_id'])]['locus_group'];r['split']=bank[int(r['bank_id'])]['split']
    local_ref=np.empty((len(bank),151),dtype='uint8');local_alt=np.empty_like(local_ref);local_pos=np.empty(len(bank),dtype='int32')
    for r in bank:
        i=r['bank_id'];p=r['changed_position'];a,b=r['reference'],r['alternate'];start=max(0,min(p-75,len(a)-151));local_ref[i]=tokens(a[start:start+151]);local_alt[i]=tokens(b[start:start+151]);local_pos[i]=p-start;assert np.sum(local_ref[i]!=local_alt[i])==1
    np.savez_compressed(DATA/'local151_sequences.npz',ref=local_ref,alt=local_alt,position=local_pos)
    with gzip.open(DATA/'sequence_bank.csv.gz','wt',newline='') as f:w=csv.DictWriter(f,fieldnames=list(bank[0]));w.writeheader();w.writerows(bank)
    with (DATA/'measured_observations.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(observations[0]));w.writeheader();w.writerows(observations)
    counts=collections.defaultdict(lambda:collections.Counter())
    for r in observations:counts[r['cell']][r['split']]+=1;counts[r['cell']][r['split']+'_detected']+=int(r['label'])
    assert all(len({r['split'] for r in observations if int(r['bank_id'])==idx})==1 for idx in list(range(min(200,len(bank)))))
    report={'sequence_bank_rows':len(bank),'observations':len(observations),'locus_groups':len(roots),'cells':{k:dict(v) for k,v in counts.items()},'overlap_and_source_exclusions':dict(exclusions),'additional_cell_eligibility_exclusions':dict(cell_exclusions),'partition':'SHA256 seed V5_locus_split_20261107, connected <=2kb loci, common canonical151reference cores joined; folds0–1calibration,2–3validation,4–9train. Shared variants/conditions stay together.','reserved_study_labels_used':False,'old_V3_test_role':'Former V3 tests are now previously observed development data; only the new whole-study Kosicki outcome test is reserved. No V3 artifact changed.','sequence_scope':'Actual published oligos for Salomon/Siraj; native genomic proxies for Kircher. Local151window is clipped within each source context without padding.','outputs':{n:sha(DATA/n) for n in ['sequence_bank.csv.gz','local151_sequences.npz','measured_observations.csv']},'V3_registry_sha256':sha(V3O/'final_model_registry.json')}
    write_json(OUT/'measured_data_audit.json',report);print('MEASURED DATA AUDIT',json.dumps({k:report[k] for k in ['sequence_bank_rows','observations','locus_groups','cells','overlap_and_source_exclusions']}),flush=True)
if __name__=='__main__':main()
