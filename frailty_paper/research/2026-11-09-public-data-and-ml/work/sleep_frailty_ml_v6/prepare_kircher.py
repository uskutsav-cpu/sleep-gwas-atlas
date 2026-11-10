"""Measured saturation data with human-cell and reference checks.

The input sequence is a GRCh38 genomic-context proxy, not an assertion that
the Sanger-confirmed physical amplicon exactly matches the reference genome.
"""
import gzip,collections,concurrent.futures,math
from common import *
CELL={
 'F9':('HepG2','liver'),'FOXE1':('HeLa','cervical'),'GP1BA':('HEL92.1.7','blood'),
 'HBB':('HEL92.1.7','blood'),'HBG1':('HEL92.1.7','blood'),'HNF4A':('HEK293T','kidney'),
 'LDLR':('HepG2','liver'),'LDLR.2':('HepG2','liver'),'MSMB':('HEK293T','kidney'),
 'PKLR-24h':('K562','blood'),'PKLR-48h':('K562','blood'),
 'TERT-GBM':('SF7996','glioma'),'TERT-HEK':('HEK293T','kidney'),
 'BCL11A':('HEL92.1.7','blood'),'IRF4':('SK-MEL-28','melanoma'),'IRF6':('HaCaT','skin'),
 'MYCrs6983267':('HEK293T','kidney'),'MYCrs11986220':('LNCaP','prostate'),
 'SORT1':('HepG2','liver'),'SORT1.2':('HepG2','liver'),'SORT1-flip':('HepG2','liver')}
def parent(name):
    if name.startswith('TERT-'):return 'TERT'
    if name.startswith('PKLR-'):return 'PKLR'
    if name.startswith('SORT1'):return 'SORT1'
    return name.removesuffix('.2')
def main():
    source=SOURCES/'kircher_elements.tsv.gz';rows=[];excluded=collections.Counter()
    with gzip.open(source,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['Release']!='GRCh38':continue
            e=r['Element']
            if e not in CELL:excluded['nonhuman_or_unresolved_cell_condition']+=1;continue
            if len(r['Ref'])!=1 or len(r['Alt'])!=1 or set(r['Ref']+r['Alt'])-set('ACGT') or r['Ref']==r['Alt']:excluded['not_SNV']+=1;continue
            try:
                tags=int(r['Barcodes']);coef=float(r['Coefficient']);p=float(r['pValue']);dna=int(r['DNA']);rna=int(r['RNA']);assert np.isfinite([coef,p]).all() and 0<=p<=1
            except (ValueError,AssertionError):excluded['unassessable_effect_or_pvalue']+=1;continue
            if tags<10 or dna<=0 or rna<0:excluded['fewer_than_10_tags_or_invalid_count']+=1;continue
            rows.append({**r,'position':int(r['Pos']),'tags':tags,'coefficient':coef,'pvalue':p,'parent_locus':parent(e)})
    regions={}
    for r in rows:
        key=(r['Chrom'],r['parent_locus']);p=r['position']
        if key not in regions:regions[key]=[p,p]
        else:regions[key]=[min(regions[key][0],p),max(regions[key][1],p)]
    cache=DATA/'kircher_reference_cache';cache.mkdir(exist_ok=True)
    def one(item):
        (chrom,locus),(lo,hi)=item;path=cache/(locus+'.json');url=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom=chr{chrom}&start={lo-1}&end={hi}'
        if not path.exists():fetch(url,path,100000)
        d=json.loads(path.read_text());assert d['start']==lo-1 and d['end']==hi and d['chrom']=='chr'+chrom
        dna=d['dna'].upper();assert len(dna)==hi-lo+1 and len(dna)<=1024
        return (chrom,locus),(lo,dna,sha(path))
    with concurrent.futures.ThreadPoolExecutor(4) as pool:reference=dict(pool.map(one,regions.items()))
    out=[];pairs={}
    for r in rows:
        lo,ref,h=reference[(r['Chrom'],r['parent_locus'])];i=r['position']-lo
        if set(ref)-set('ACGT'):excluded['ambiguous_genomic_context']+=1;continue
        if ref[i]!=r['Ref']:excluded['genomic_reference_base_mismatch']+=1;continue
        alt=ref[:i]+r['Alt']+ref[i+1:];key=f"chr{r['Chrom']}:{r['position']}:{r['Ref']}:{r['Alt']}"
        # Same genomic contrast across full replicates/cell conditions shares
        # one sequence key and one eventual locus partition.
        if key not in pairs:pairs[key]={'variant_GRCh38':key,'reference':ref,'alternate':alt,'changed_position':i,'chromosome':'chr'+r['Chrom'],'position_hg38':r['position'],'parent_locus':r['parent_locus']}
        cell,group=CELL[r['Element']]
        out.append({'variant_GRCh38':key,'study':'Kircher2019','experiment':r['Element'],'cell':cell,'cell_group':group,'label':int(r['pvalue']<1e-5),'absolute_effect':abs(r['coefficient']),'signed_effect':r['coefficient'],'pvalue':r['pvalue'],'tags':r['tags'],'activity_mean':'','active_screen':'','parent_locus':r['parent_locus'],'label_definition':'Source coefficient p<1e-5; other assessed rows are non-detections, not proven inactive.'})
    for name,data in [('kircher_observations.csv',out),('kircher_sequence_pairs.csv',list(pairs.values()))]:
        with (DATA/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    report={'source_sha256':sha(source),'source_revision':json.loads((SOURCES/'kircher_download_receipt.json').read_text())['revision'],'observations':len(out),'unique_genomic_contrasts':len(pairs),'human_cell_lines':dict(collections.Counter(r['cell'] for r in out)),'experiments':dict(collections.Counter(r['experiment'] for r in out)),'exclusions':dict(excluded),'called_effects':sum(r['label'] for r in out),'source_name_caution':'The portal uses GP1BA in result rows, but its table labels this regulatory locus GP1BB. Preserve the source identifier and genomic coordinates; do not infer a gene assignment.','sequence_scope':'Native GRCh38 genomic-context proxy with every changed reference base checked, not recovered physical amplicons.','replicate_scope':'LDLR/SORT1 full replicates and PKLR timepoints retained as distinct measured observations with shared sequence/locus partition. Not independent variants.','excluded_cell_conditions':'All mouse-cell assays and TERT-GAa/GSc unresolved conditions excluded from primary training.','outputs':{name:sha(DATA/name) for name in ['kircher_observations.csv','kircher_sequence_pairs.csv']},'reference_payloads':{c+'_'+l:v[2] for (c,l),v in reference.items()}}
    write_json(OUT/'kircher_preparation_audit.json',report);print(json.dumps({k:report[k] for k in ['observations','unique_genomic_contrasts','human_cell_lines','exclusions','called_effects']}),flush=True)
if __name__=='__main__':main()
