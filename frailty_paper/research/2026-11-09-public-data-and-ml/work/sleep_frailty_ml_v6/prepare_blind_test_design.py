"""Audit final-study sequences and coordinates without opening outcomes."""
import gzip,re,collections,concurrent.futures
from common import *
COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
def fasta(path):
    result={};name=None;parts=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    assert name not in result;result[name]=''.join(parts).upper()
                name=line[1:].strip();parts=[]
            else:parts.append(line.strip())
    if name is not None:assert name not in result;result[name]=''.join(parts).upper()
    return result
def main():
    assert not (OUT/'final_test_evaluation_started.json').exists()
    fp=SOURCES/'kosicki_ENCFF752NFF.fasta.gz';bp=SOURCES/'kosicki_ENCFF538QRO.bed.gz'
    sequences=fasta(fp)
    with gzip.open(bp,'rt') as f:bed={r[3]:r for r in (line.rstrip().split('\t') for line in f)}
    assert set(bed)==set(sequences)
    exclusion=collections.Counter();conventions=collections.Counter();candidates=[];refs={}
    for name,alt in sequences.items():
        m=re.search(r'\|([ACGT])-([ACGT])-(\d+)$',name)
        if not m:exclusion['not_named_SNV_design']+=1;continue
        prefix=name.rsplit('|',1)[0];rn=prefix+'|reference'
        if rn not in sequences:exclusion['missing_reference_construct']+=1;continue
        ref=sequences[rn]
        if len(ref)!=270 or len(alt)!=270 or set(ref+alt)-set('ACGT'):exclusion['non270_or_ambiguous_sequence']+=1;continue
        changes=[i for i,(a,b) in enumerate(zip(ref,alt)) if a!=b]
        if len(changes)!=1:exclusion['not_one_base_difference']+=1;continue
        i=changes[0]
        coord=re.search(r'\|(chr[^-]+)-(\d+)([+-])(\d+)$',prefix)
        if not coord:exclusion['unparsed_coordinates']+=1;continue
        chrom,start,strand,end=coord[1],int(coord[2]),coord[3],int(coord[4])
        if start>end or end-start+1!=270 or bed[name][0]!=chrom or int(bed[name][1])!=start-1 or int(bed[name][2])!=end:exclusion['coordinate_or_bed_disagreement']+=1;continue
        if (ref[i],alt[i],i+1)==(m[1],m[2],int(m[3])):convention='construct_orientation_descriptor'
        elif strand=='-' and (ref[i],alt[i],270-i)==(rc(m[1]),rc(m[2]),int(m[3])):convention='forward_genomic_descriptor_on_reverse_construct'
        else:exclusion['descriptor_disagrees_with_sequence_in_both_supported_conventions']+=1;continue
        conventions[convention]+=1
        pos=start+i if strand=='+' else end-i
        a,b=(ref[i],alt[i]) if strand=='+' else (rc(ref[i]),rc(alt[i]))
        candidates.append({'name':name,'reference_name':rn,'chromosome':chrom,'position_hg38':pos,'ref':a,'alt':b,'construct_strand':strand,'descriptor_convention':convention,'changed_position':i,'variant_GRCh38':f'{chrom}:{pos}:{a}:{b}'})
        refs[rn]={'chromosome':chrom,'start':start,'end':end,'strand':strand,'reference_sequence':ref}
    cache=DATA/'kosicki_reference_cache';cache.mkdir(exist_ok=True)
    def verify(item):
        name,r=item;key=hashlib.sha256(name.encode()).hexdigest()[:24];p=cache/(key+'.json');url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom={r['chromosome']}&start={r['start']-1}&end={r['end']}"
        for attempt in range(3):
            try:
                if not p.exists():fetch(url,p,100_000)
                d=json.loads(p.read_text());assert d['chrom']==r['chromosome'] and d['start']==r['start']-1 and d['end']==r['end'];seq=d['dna'].upper();assert len(seq)==270
                expected=r['reference_sequence'] if r['strand']=='+' else rc(r['reference_sequence'])
                return name,seq,expected==seq,sha(p)
            except Exception:
                if attempt==2:raise
                time.sleep(attempt+1)
    checked={}
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        for n,result in enumerate(pool.map(verify,refs.items()),1):
            checked[result[0]]=result[1:]
            if n%50==0:print('Blind construct reference checks',n,'/',len(refs),flush=True)
    kept=[];reference_mismatches=0
    for r in candidates:
        g,full,h=checked[r['reference_name']];info=refs[r['reference_name']]
        if g[r['position_hg38']-info['start']]!=r['ref']:exclusion['genomic_reference_base_mismatch']+=1;continue
        r['full_reference_construct_matches_hg38']=full
        if not full:reference_mismatches+=1
        kept.append(r)
    # Store design only. No measured outcomes, significance, or assay activity.
    with (DATA/'kosicki_blind_design.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(kept[0]));writer.writeheader();writer.writerows(kept)
    write_json(OUT/'kosicki_blind_design_audit.json',{'study':'Kosicki2025','fasta_sha256':sha(fp),'bed_sha256':sha(bp),'design_sha256':sha(DATA/'kosicki_blind_design.csv'),'designed_named_SNVs':len(candidates),'eligible_reference_checked_pairs':len(kept),'unique_genomic_contrasts':len(set(r['variant_GRCh38'] for r in kept)),'unique_reference_tiles':len(refs),'nonmatching_full_context_pairs':reference_mismatches,'exclusions':dict(exclusion),'descriptor_conventions':dict(conventions),'full_reference_checks':{k:{'full_construct_match':v[1],'payload_sha256':v[2]} for k,v in checked.items()},'measured_outcomes_opened':False,'evaluation_or_model_selection_performed':False,'scope':'Designed pairs; actual experimental QC and assessability will be applied by the frozen final evaluator.'})
    print('BLIND DESIGN COMPLETE',len(kept),'pairs;',len(set(r['variant_GRCh38'] for r in kept)),'unique contrasts;',dict(exclusion),flush=True)
if __name__=='__main__':main()
