#!/usr/bin/env python3
"""Query AlphaGenome Atlas after explicit selection of scorers and tissue ontology IDs."""
import argparse, csv, gzip, hashlib, importlib.metadata, json, math, os
from datetime import datetime, timezone
from pathlib import Path


def hash_file(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def clean(v):
    if hasattr(v,'item'):
        try: v=v.item()
        except (ValueError,AttributeError): pass
    if v is None or isinstance(v,float) and not math.isfinite(v): return None
    if isinstance(v,(str,int,float,bool)): return v
    return str(v)

def flatten(results):
    """Retain all returned track and observation metadata with each score."""
    out=[]
    for scorer, adata in results.items():
        for i in range(adata.n_obs):
            obs={str(k):clean(v) for k,v in adata.obs.iloc[i].items()}
            for j in range(adata.n_vars):
                track={str(k):clean(v) for k,v in adata.var.iloc[j].items()}
                quantiles=adata.layers.get('quantiles')
                out.append({'scorer':scorer,'observation':obs,'track':track,'raw_score':clean(adata.X[i,j]),'quantile_score':clean(quantiles[i,j]) if quantiles is not None else None})
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    meta=sub.add_parser('metadata'); meta.add_argument('--out',type=Path,required=True)
    run=sub.add_parser('score'); run.add_argument('--evidence',type=Path,required=True); run.add_argument('--out',type=Path,required=True); run.add_argument('--scorer',action='append',required=True); run.add_argument('--ontology',action='append',required=True); run.add_argument('--region',default='placo_insomnia__adhd_chr5_103447968_104447968'); run.add_argument('--limit',type=int)
    args=parser.parse_args(); key=os.environ.get('ALPHAGENOME_API_KEY')
    if not key: parser.error('Set ALPHAGENOME_API_KEY in your environment; do not pass it on the command line.')
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    client=atlas.create(key,timeout=60); metadata=client.scorer_metadata()
    if args.command=='metadata':
        args.out.mkdir(parents=True,exist_ok=False)
        rows=[]
        for name,m in metadata.items():
            for _,r in m.track_metadata.iterrows():
                rows.append({'scorer':name,'is_signed':m.is_signed,**{str(k):clean(v) for k,v in r.items()}})
        (args.out/'scorer_metadata.json').write_text(json.dumps(rows,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
        print(f'Saved {len(rows)} scorer-track records to {args.out}.'); return
    if args.limit is not None and args.limit<1: parser.error('--limit must be positive')
    unknown=set(args.scorer)-set(metadata)
    if unknown: parser.error(f'Unknown scorer(s): {sorted(unknown)}; inspect metadata first')
    known=set()
    for scorer in args.scorer:
        tracks=metadata[scorer].track_metadata
        if 'ontology_curie' in tracks: known.update(str(x) for x in tracks['ontology_curie'].dropna())
    unknown_ont=set(args.ontology)-known
    if unknown_ont: parser.error(f'Ontology ID(s) absent from chosen scorer metadata: {sorted(unknown_ont)}')
    with args.evidence.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    regions={r['candidate_locus_id'] for r in rows}
    if args.region!='all' and args.region not in regions: parser.error('Region not in evidence')
    eligible=[r for r in rows if r['mapping_status']=='READY_UNIQUE_REFERENCE_VERIFIED' and (args.region=='all' or r['candidate_locus_id']==args.region)]
    variants=sorted({r['atlas_variant'] for r in eligible})
    if args.limit: variants=variants[:args.limit]
    if not variants: parser.error('No eligible variants selected')
    args.out.mkdir(parents=True,exist_ok=False); (args.out/'scores').mkdir()
    config={'evidence_sha256':hash_file(args.evidence),'sdk_version':importlib.metadata.version('alphagenome'),'scorers':args.scorer,'ontology_terms':args.ontology,'region':args.region,'variant_count':len(variants),'service_model_release':'not exposed by this client; unknown','retrieved_utc':datetime.now(timezone.utc).isoformat(),'disclaimer':'Atlas outputs are molecular predictions, not validation or disease-risk predictions; API terms restrict training other ML models.'}
    (args.out/'run_config.json').write_text(json.dumps(config,indent=2)+'\n')
    status=[]
    for v in variants:
        try:
            result=client.query_variant(genome.Variant.from_str(v),requested_scorers=args.scorer,ontology_terms=args.ontology)
            scores=flatten(result)
            name=hashlib.sha256(v.encode()).hexdigest()+'.json.gz'
            with gzip.open(args.out/'scores'/name,'wt',encoding='utf-8') as f: json.dump({'variant':v,'scores':scores},f,allow_nan=False)
            status.append({'variant':v,'status':'SCORED' if scores else 'NO_RETURNED_TRACKS','n_score_rows':len(scores),'file':f'scores/{name}','error_type':''})
        except Exception as e:
            status.append({'variant':v,'status':'QUERY_FAILED','n_score_rows':0,'file':'','error_type':type(e).__name__})
        # Preserve progress even if a later request fails or the process is interrupted.
        with (args.out/'query_status.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['variant','status','n_score_rows','file','error_type'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(status)
        print(f'{len(status)}/{len(variants)} {status[-1]["status"]}',flush=True)
    print('Completed queries. Inspect query_status.tsv; failed or empty queries are not negative evidence.')
if __name__=='__main__': main()
