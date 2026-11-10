"""Stream pinned author benchmarks; exclude entire linked candidate components.

Only identity/location metadata are used. New-study outcome files are unread.
"""
import bisect,csv,hashlib,json,ssl,time,urllib.request
from pathlib import Path
import certifi

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/sleep_frailty_ml_native_stable/benchmark_exposure'
SOURCE=ROOT/'outputs/sleep_frailty_ml_channels/source_candidates/Ulirsch_complete_exposure_profile.csv'
COMMIT='a20dbc7d0b08daae3af1aebcfb8a71d9f198ab5d'
FILES=[
 ('datasets/asb/HepG2__hepatoblastoma_.asb.tsv',29194322,'460ead596917111c509dc953b3936a0da4286287'),
 ('datasets/asb/K562__myelogenous_leukemia_.asb.tsv',26055028,'ba68aebdc00cf37531b49eb430bb05781452d922'),
 ('datasets/coverage/HepG2__hepatoblastoma_.atac.tsv',3808746,'fe1a69e83dad82880d9e92cf25e05046941ab70c'),
 ('datasets/coverage/HepG2__hepatoblastoma_.dnase.tsv',3391645,'bec5cb2f76fc0519a98e93423bdd1785b4b11f89'),
 ('datasets/coverage/K562__myelogenous_leukemia_.atac.tsv',2138034,'12ed62c32db69226afd092222ccd719794e132e5'),
 ('datasets/coverage/K562__myelogenous_leukemia_.dnase.tsv',7683928,'ed05d116b0f5eb41c1b3adbb6640ea3048af74e1')]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    start=time.monotonic();OUT.mkdir(exist_ok=True)
    rows=list(csv.DictReader(SOURCE.open()));fields=list(rows[0]);intervals={}
    for r in rows:
        component=int(r['joint_component_after_all_source_DNA'])
        for field,build in [('original_GRCh37_variant','GRCh37_conservative'),('variant_GRCh38','GRCh38')]:
            chrom,pos,*_=r[field].split(':');pos=int(pos)
            # Source construct strand can differ from genome. The symmetric
            # 145bp envelope covers either orientation, plus a 2kb neighborhood.
            intervals.setdefault(chrom,set()).add((pos-144-2000,pos+144+2000,component,build))
    lookup={c:sorted(a) for c,a in intervals.items()};starts={c:[a[0] for a in v] for c,v in lookup.items()}
    touched=set();receipts=[];hits=set()
    ctx=ssl.create_default_context(cafile=certifi.where())
    for path,size,expected in FILES:
        url=f'https://raw.githubusercontent.com/autosome-ru/human_legnet/{COMMIT}/{path}'
        git=hashlib.sha1(f'blob {size}\0'.encode());h=hashlib.sha256();nbytes=0;nrows=0;local=set()
        def consume(line):
            nonlocal nbytes
            git.update(line);h.update(line);nbytes+=len(line)
        with urllib.request.urlopen(url,context=ctx,timeout=60) as stream:
            header=stream.readline();consume(header);names=header.decode().rstrip('\r\n').split('\t')
            ci=names.index('#chr');pi=names.index('pos' if '/asb/' in path else 'start')
            ri=names.index('ref');ai=names.index('alt')
            for line in stream:
                consume(line)
                if not line.strip():continue
                parts=line.decode().rstrip('\r\n').split('\t');chrom=parts[ci]
                # Original asb_predict CLI passes one_indexed=True;
                # coverage_predict CLI explicitly passes one_indexed=False.
                position=int(parts[pi])+(0 if '/asb/' in path else 1);nrows+=1
                if chrom not in lookup:continue
                end=bisect.bisect_right(starts[chrom],position)
                # Fixed 4288bp span permits a bounded bisect search.
                begin=bisect.bisect_left(starts[chrom],position-4288)
                for lo,hi,component,build in lookup[chrom][begin:end]:
                    if lo<=position<=hi:
                        local.add(component);hits.add((path,chrom,position,parts[ri],parts[ai],component,build))
        assert nbytes==size and git.hexdigest()==expected,(path,nbytes,git.hexdigest())
        touched.update(local);receipts.append({'path':path,'url':url,'rows':nrows,'bytes':nbytes,'verified_Git_blob_SHA1':expected,
            'sha256':h.hexdigest(),'position_conversion':'1-based pos' if '/asb/' in path else '0-based start + 1',
            'source_components_near_benchmark':len(local)})
        print('AUTHOR BENCHMARK SCREEN',path,nrows,'components',len(local),flush=True)
    hit_path=OUT/'author_benchmark_component_hits.csv'
    with hit_path.open('w') as handle:
        writer=csv.writer(handle);writer.writerow(['benchmark','chromosome','position_1based','reference','alternate','whole_source_component','candidate_coordinate_build']);writer.writerows(sorted(hits))
    candidates=[r for r in rows if r['strict_history_and_pretraining_eligible']=='True']
    kept=[r for r in candidates if int(r['joint_component_after_all_source_DNA']) not in touched]
    target=OUT/'Ulirsch_candidates_after_author_benchmark_guard.csv'
    with target.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(kept)
    result={'passed_identity_and_component_screen':True,'original_reader_repository_commit':COMMIT,
      'source_profile_sha256':sha(SOURCE),'source_designs_screened_before_candidate_filter':len(rows),
      'benchmark_sources':receipts,'candidate_designs_before':len(candidates),'candidate_designs_after':len(kept),
      'candidate_components_before':len({r['joint_component_after_all_source_DNA'] for r in candidates}),
      'candidate_components_after':len({r['joint_component_after_all_source_DNA'] for r in kept}),
      'candidate_variants_after':len({r['variant_GRCh38'] for r in kept}),
      'exclusion_scope':'Entire connected component excluded if any full-source member is within 2kb plus a symmetric 145bp source-envelope of any published author benchmark position; all benchmark rows included regardless of association strength.',
      'build_policy':'ADASTRA and UDACHA publications use hg38. Conservatively screen both verified candidate hg38 and original hg19 positions; the extra hg19 screen does not assert benchmark data are hg19.',
      'primary_build_sources':['https://pmc.ncbi.nlm.nih.gov/articles/PMC9252736/','https://www.nature.com/articles/s41467-024-55513-2'],
      'outcome_files_read':False,'model_performance_calculated':False,'outside_validation_ready':False,
      'remaining':'Complete the four Kircher saturation-locus benchmark inventory and supplementary-methods exposure review, and freeze a justified model and confidence rules before any new-study comparison.',
      'candidate_file_sha256':sha(target),'hits_sha256':sha(hit_path),'code_sha256':sha(Path(__file__)),
      'elapsed_seconds':time.monotonic()-start}
    (OUT/'author_benchmark_exposure_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k.startswith('candidate_')}),flush=True)


if __name__=='__main__':main()
