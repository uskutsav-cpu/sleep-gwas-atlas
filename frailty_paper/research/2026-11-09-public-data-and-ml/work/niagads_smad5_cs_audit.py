#!/usr/bin/env python3
"""Inspect chromosome 5 SMAD5 credible-set SNPs without extracting the archive."""
import csv, gzip, hashlib, io, json, tarfile
from pathlib import Path

root = Path(__file__).resolve().parent.parent
archive = root / 'work/niagads_snuc_cs95.tar'
digest = hashlib.md5(archive.read_bytes()).hexdigest()
assert digest == 'a278bf209cf7156cf90edabf6f7c0a4f', 'archive MD5 mismatch'
results = {}
with tarfile.open(archive) as tar:
    for member in tar:
        name = member.name
        if 'ROSMAP_CUIMC1_2_MIT.snuc-eQTL.' not in name or '_v1_5_snp_scfmCs95.' not in name or not name.endswith('.bed.gz'):
            continue
        context = name.split('.snuc-eQTL.')[1].split('.')[0]
        if context not in ('exc', 'inh', 'oli'):
            continue
        assert context not in results
        decoded = gzip.decompress(tar.extractfile(member).read()).decode()
        records = [r for r in csv.DictReader(io.StringIO(decoded), delimiter='\t') if r['target_ensembl_id'] == 'ENSG00000113658']
        records.sort(key=lambda r: float(r['PIP']), reverse=True)
        target = [r for r in records if r['variant_id'] == 'chr5:136202653:A:C']
        results[context] = {'archive_member': name, 'n_snp_rows': len(records),
                            'credible_sets': sorted(set(r['cs_95'] for r in records)),
                            'highest_pip_variant': records[0]['variant_id'], 'highest_pip': float(records[0]['PIP']),
                            'candidate': [{'variant_id': r['variant_id'], 'pip': float(r['PIP']), 'cs_95': r['cs_95']} for r in target]}
assert set(results) == {'exc', 'inh', 'oli'}
receipt = {'archive_md5': digest, 'source': 'NIAGADS NG00184.v1 snuc-eQTL single-context 95% credible sets',
           'target': 'ENSG00000113658', 'candidate_GRCh38': 'chr5:136202653:A:C', 'results': results}
(root / 'outputs/stage6_niagads_smad5_cs_audit.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(results, indent=2))
