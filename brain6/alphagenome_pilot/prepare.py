#!/usr/bin/env python3
"""Prepare Brain6's six eligible insomnia–ADHD ABF credible-set variants.

Coordinate conventions: input POS is 1-based GRCh37, pyliftover and UCSC
sequence windows are 0-based, output positions are 1-based GRCh38.
"""
import argparse, collections, csv, hashlib, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path
from pyliftover import LiftOver

ROOT = Path(__file__).resolve().parents[2]
UP = ROOT / 'brain6/results/brain6_exploratory_finemap_coloc_v2'
DEFAULT_EXTERNAL = Path('/Volumes/Extreme SSD/brain6-work/exploratory-finemap-coloc-v2')
COMP = str.maketrans('ACGT', 'TGCA')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def read(path):
    with path.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows, fields):
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--external',type=Path,default=DEFAULT_EXTERNAL); p.add_argument('--chain',type=Path,required=True); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True); cache=a.out/'sequence_cache'; cache.mkdir(exist_ok=True)
    prov=json.loads((UP/'integrity_provenance.json').read_text())
    for name,expected in prov['result_sha256'].items():
        if digest(UP/name)!=expected: raise ValueError(f'Upstream hash mismatch: {name}')
    pip=a.external/'variant_pip.tsv'
    if digest(pip)!=prov['full_variant_pip_sha256']: raise ValueError('Full PIP source hash mismatch')
    coloc=[x for x in read(UP/'trait_coloc_25.tsv') if x['pair_id']=='insomnia__adhd' and x['status']=='EXPLORATORY_TRAIT_COLOC_ABF']
    if len(coloc)!=6 or len({x['candidate_locus_id'] for x in coloc})!=6: raise ValueError('Expected six distinct eligible candidates')
    ids={x['candidate_locus_id'] for x in coloc}
    credible=[x for x in read(UP/'fine_mapping_credible_sets.tsv') if x['candidate_locus_id'] in ids]
    members=collections.defaultdict(dict)
    for x in credible:
        key=(x['candidate_locus_id'],x['snp'].lower()); trait=x['trait_id']
        if trait in members[key]: raise ValueError('Duplicate credible set membership')
        members[key][trait]=x
    manifests={(x['candidate_locus_id'],x['trait_id']):x for x in read(UP/'fine_mapping_input_manifest.tsv') if x['candidate_locus_id'] in ids}
    source={}; source_hash={}
    for key,m in manifests.items():
        path=a.external/'inputs'/Path(m['input_tsv']).name
        if m['status']!='READY' or not path.exists() or digest(path)!=m['input_sha256']: raise ValueError(f'Input unavailable or hash mismatch: {key} {path}')
        source_hash[path.name]=digest(path)
        for row in read(path):
            k=(key[0],row['SNP'].lower(),key[1])
            if k in source: raise ValueError(f'Duplicate source SNP {k}')
            source[k]=row
    pips={}
    for row in read(pip):
        key=(row['candidate_locus_id'],row['snp'].lower(),row['trait_id'])
        if key[:2] in members:
            if key in pips: raise ValueError('Duplicate SNP PIP')
            pips[key]=row
    lift=LiftOver(str(a.chain)); sequence_hashes={}
    def base(build,chrom,pos):
        # Cache a reference sequence chunk with its coordinates in filename.
        start=(pos//100000)*100000; end=start+100000
        path=cache/f'{build}_{chrom}_{start}_{end}.json'
        if not path.exists():
            url='https://api.genome.ucsc.edu/getData/sequence?'+urllib.parse.urlencode({'genome':build,'chrom':chrom,'start':start,'end':end},safe=':')
            last=None
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url,timeout=60) as response: raw=response.read()
                    obj=json.loads(raw)
                    if obj.get('genome')!=build or obj.get('chrom')!=chrom or obj.get('start')!=start or obj.get('end')!=end or len(obj.get('dna',''))!=end-start: raise ValueError('Unexpected UCSC sequence response')
                    path.write_bytes(raw); break
                except Exception as e:
                    last=e; time.sleep(2**attempt)
            else: raise RuntimeError(f'UCSC sequence fetch failed: {build} {chrom}:{start}: {last}')
        obj=json.loads(path.read_text()); sequence_hashes[path.name]=digest(path)
        return obj['dna'][pos-start].upper()
    rows=[]
    for (cid,snp), traits in sorted(members.items()):
        available=[source.get((cid,snp,t)) for t in ('insomnia','adhd')]; avail=[s for s in available if s]
        if not avail: raise ValueError(f'Credible variant lacks source {cid} {snp}')
        orig=avail[0]; chrom='chr'+orig['CHR']; pos37=int(orig['POS']); alleles={orig['A1'].upper(),orig['A2'].upper()}
        if any((s['CHR'],s['POS'],s['A1'],s['A2'])!=(orig['CHR'],orig['POS'],orig['A1'],orig['A2']) for s in avail): raise ValueError(f'Inconsistent source alleles {cid} {snp}')
        row={'candidate_locus_id':cid,'snp':snp,'chromosome':chrom,'pos_grch37_1based':pos37,'source_a1':orig['A1'],'source_a2':orig['A2']}
        for t,s in zip(('insomnia','adhd'),available):
            pp=pips.get((cid,snp,t)); member=traits.get(t)
            row.update({f'{t}_beta_a1':s['BETA'] if s else '',f'{t}_se':s['SE'] if s else '',f'{t}_conditional_abf_pip':pp['conditional_pip'] if pp else '', f'{t}_credible_set_member':'YES' if member else 'NO',f'{t}_credible_set_rank':member['rank'] if member else ''})
        row.update({'mapping_status':'','pos_grch38_1based':'','ref_grch37':'','ref_grch38':'','alt_grch38':'','atlas_variant':'','atlas_url':'','alphagenome_status':'NOT_RUN_NO_API_KEY'})
        if len(alleles)!=2 or not alleles <= set('ACGT') or any(len(x)!=1 for x in alleles): row['mapping_status']='UNSUPPORTED_NON_SNV_OR_ALLELES'
        elif alleles in ({'A','T'},{'C','G'}): row['mapping_status']='HOLD_PALINDROMIC_ALLELES'
        else:
            hits=lift.convert_coordinate(chrom,pos37-1)
            if len(hits)!=1 or hits[0][0]!=chrom: row['mapping_status']='HOLD_NONUNIQUE_OR_CROSS_CHROMOSOME'
            else:
                chr38,p38,strand,_=hits[0]; ref37=base('hg19',chrom,pos37-1); ref38=base('hg38',chr38,p38)
                row.update(ref_grch37=ref37,ref_grch38=ref38,pos_grch38_1based=p38+1)
                # Resolve LD-reference strand against hg19 before lift; reject ambiguous cases.
                forward=alleles if ref37 in alleles else {x.translate(COMP) for x in alleles}
                if ref37 not in forward: row['mapping_status']='HOLD_GRCH37_ALLELE_MISMATCH'
                else:
                    alt37=next(x for x in forward if x!=ref37)
                    projected_ref=ref37.translate(COMP) if strand=='-' else ref37
                    alt38=alt37.translate(COMP) if strand=='-' else alt37
                    if ref38!=projected_ref: row['mapping_status']='HOLD_REFERENCE_CHANGED_BETWEEN_ASSEMBLIES'
                    elif alt38==ref38: row['mapping_status']='HOLD_ALT_EQUALS_REF'
                    else:
                        v=f'{chr38}:{p38+1}:{ref38}>{alt38}'
                        row.update(mapping_status='READY_UNIQUE_REFERENCE_VERIFIED',alt_grch38=alt38,atlas_variant=v,atlas_url='https://deepmind.google.com/science/alphagenome/atlas?'+urllib.parse.urlencode({'q':v,'m':'variant'}))
        rows.append(row)
    fields=['candidate_locus_id','snp','chromosome','pos_grch37_1based','source_a1','source_a2','insomnia_beta_a1','insomnia_se','insomnia_conditional_abf_pip','insomnia_credible_set_member','insomnia_credible_set_rank','adhd_beta_a1','adhd_se','adhd_conditional_abf_pip','adhd_credible_set_member','adhd_credible_set_rank','mapping_status','pos_grch38_1based','ref_grch37','ref_grch38','alt_grch38','atlas_variant','atlas_url','alphagenome_status']
    write(a.out/'variant_evidence.tsv',rows,fields)
    sens=read(UP/'trait_coloc_prior_sensitivity.tsv'); fm=read(UP/'fine_mapping_50.tsv')
    regional=[]
    for c in sorted(coloc,key=lambda x: (not x['candidate_locus_id'].startswith('placo_insomnia__adhd_chr5_103447968'),x['candidate_locus_id'])):
        cid=c['candidate_locus_id']; rr=[r for r in rows if r['candidate_locus_id']==cid]
        lo=next(s for s in sens if s['candidate_locus_id']==cid and float(s['p12'])==1e-6)
        results={f['trait_id']:f for f in fm if f['candidate_locus_id']==cid}
        regional.append({'candidate_locus_id':cid,'n_union_credible_variants':len(rr),'n_ready':sum(x['mapping_status']=='READY_UNIQUE_REFERENCE_VERIFIED' for x in rr),'n_held':sum(x['mapping_status']!='READY_UNIQUE_REFERENCE_VERIFIED' for x in rr),'insomnia_set_status':results['insomnia']['status'],'adhd_set_status':results['adhd']['status'],'PP_H4_p12_1e-5':c['PP.H4'],'PP_H4_p12_1e-6':lo['PP.H4'],'label':'EXPLORATORY_SINGLE_SIGNAL_ABF_NO_INDEPENDENT_LOCUS_REPLICATION'})
    write(a.out/'regions.tsv',regional,list(regional[0]))
    final=ROOT/'brain6/paper/final_package_v1'
    final_prov=json.loads((final/'BRAIN6_FINAL_PROVENANCE.json').read_text())
    context_names=['BRAIN6_FINAL_CANDIDATES.tsv','BRAIN6_FINAL_REGULATORY.tsv','BRAIN6_FINAL_EQTL_COLOC.tsv','BRAIN6_FINAL_SQTL_COLOC.tsv']
    for name in context_names:
        path=final/name
        if digest(path)!=final_prov['output_sha256'][name]: raise ValueError(f'Final package hash mismatch: {name}')
        all_rows=read(path); selected=[x for x in all_rows if x['candidate_locus_id'] in ids]
        if not selected: raise ValueError(f'No context rows: {name}')
        write(a.out/name.lower(),selected,list(all_rows[0]))
    receipt={'source_git_commit':'303c905c229104b7bf623df6de7431abb87e0727','chain_sha256':digest(a.chain),'upstream_integrity_sha256':digest(UP/'integrity_provenance.json'),'full_variant_pip_sha256':digest(pip),'input_sha256':source_hash,'sequence_sha256':sequence_hashes,'sequence_source':'UCSC Genome Browser API hg19 and hg38','n_regions':len(regional),'n_variant_rows':len(rows),'status_counts':dict(collections.Counter(x['mapping_status'] for x in rows)),'output_sha256':{x:digest(a.out/x) for x in ['variant_evidence.tsv','regions.tsv']},'alphagenome':'NOT_RUN_NO_API_KEY','method':'exploratory only; single signal credible sets; missing or held is not negative'}
    receipt['final_package_provenance_sha256']=digest(final/'BRAIN6_FINAL_PROVENANCE.json')
    receipt['final_context_source_sha256']={name:digest(final/name) for name in context_names}
    receipt['output_sha256'].update({name.lower():digest(a.out/name.lower()) for name in context_names})
    (a.out/'provenance.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'n_regions':len(regional),'n_variant_rows':len(rows),'mapping':receipt['status_counts']},indent=2))
if __name__=='__main__': main()
