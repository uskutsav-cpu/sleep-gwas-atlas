"""Add a separate development-only sidecar of original paired GLM contrasts.

Direction is traced to published software and physical allele identities, not
guessed from agreement with activity means. Old targets remain immutable.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sleep_frailty_ml_measurement'))
from measurement_common import *
OUT=ROOT/'outputs/sleep_frailty_ml_direction'


def main():
    assert_prior_versions()
    code=(OUT/'sources/MPRAmodel_e853506.R').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(code)).encode()+b'\0'+code).hexdigest()=='ce908dacaf0a4ff580e8200bc6cc59fdb0288f0a'
    text=code.decode(); assert 'contrast=list("materialRNA.allelealt","materialDNA.allelealt")' in text
    assert 'allele=factor(rep(c("ref","alt")' in text
    assert (OUT/'paired_contrast_independent_R_receipt.tsv').exists()
    with np.load(OLD/'coverage_training_arrays.npz') as z: d={k:z[k] for k in z.files}
    sourcecode=list(d['study_names']).index('Siraj')
    eligible=(d['study']==sourcecode)&(d['split']!=2)
    by_bank_cell={}
    for i in np.flatnonzero(eligible):
        key=(int(d['bank_index'][i]),str(d['cell_line_names'][d['line'][i]]))
        assert key not in by_bank_cell; by_bank_cell[key]=int(i)
    old_source_rows={int(r['bank_id']):int(r['original_row_index']) for r in csv_rows(V6/'sequence_bank.csv.gz') if r['sequence_source']=='siraj'}
    manifest=list(csv_rows(ROOT/'work/sleep_frailty_ml_v3/siraj_manifest_local.csv'))
    identity={};physical={}
    for r in csv_rows(V7/'sequence_bank.csv.gz'):
        if r['context_feature_source']!='v6' or int(r['context_feature_index']) not in old_source_rows: continue
        m=manifest[old_source_rows[int(r['context_feature_index'])]]; bank=int(r['bank_id'])
        assert m['variant'] not in identity
        identity[m['variant']]=bank;physical[bank]=r
    ol27={}
    for r in csv_rows(OLD/'OL27_development_sequence_bank.csv.gz'):
        assert r['original_variant_GRCh37'] not in ol27
        bank=int(r['bank_id']);ol27[r['original_variant_GRCh37']]=bank;physical[bank]=r
    cr=json.loads((V6/'sources/siraj_additional_cell_receipt.json').read_text())
    br=json.loads((V6/'sources/siraj_pooled_K562_receipt.json').read_text())
    files=[(BASE/'siraj2026_neural/OL42_SKNSH_emVAR_glm_paired_20221220.out.gz','SK-N-SH',identity)]
    files += [(Path(r['file']),'HepG2' if '_HepG2_' in r['file'] else 'HCT116',identity) for r in cr['files']]
    files += [(Path(br['file']),'K562',identity)]
    retrieval=json.loads((ROOT/'outputs/sleep_frailty_ai_v9_validation/next_source_screen/Siraj_unused_libraries/targeted_retrieval_receipt.json').read_text())
    line_map={'SKNSH':'SK-N-SH','HepG2':'HepG2','K562':'K562'}
    files += [(ROOT/r['path'],line_map[r['cell']],ol27) for r in retrieval['files'] if r['library']=='OL27']
    expected={r['path']:r['sha256'] for r in json.loads((V8/'outputs/measured_training_data_audit.json').read_text())['source_files']}
    expected.update({str(ROOT/r['path']):r['sha256'] for r in retrieval['files'] if r['library']=='OL27'})
    values=np.full(len(d['effect']),np.nan,dtype='float32');rows=[];sources=[];counts=collections.Counter();seen=set()
    for path,line,mapping in files:
        digest=sha(path);assert digest==expected[str(path)];sources.append({'path':str(path),'sha256':digest,'line':line})
        with gzip.open(path,'rt') as handle:
            for r in csv.DictReader(handle,delimiter='\t'):
                if r['window']!='center' or r['strand']!='fwd' or r['haplotype']!='ref': continue
                bank=mapping.get('chr'+r['SNP'])
                if bank is None: continue
                i=by_bank_cell.get((bank,line))
                if i is None: continue
                assert i not in seen,(path,i);seen.add(i)
                assert r['ID']==r['SNP']+':R:wC' and r['allele']=='ref'
                p=physical[bank];at=int(p['changed_position']);a,b=p['reference'],p['alternate']
                assert len(a)==len(b) and [j for j,(x,y) in enumerate(zip(a,b)) if x!=y]==[at]
                parts=p['variant_GRCh38'].split(':')
                if [a[at],b[at]]==parts[2:]: orientation='+'
                elif [rc(a[at]),rc(b[at])]==parts[2:]: orientation='-'
                else: raise AssertionError((bank,p['variant_GRCh38'],a[at],b[at],r['SNP']))
                if [r['ref_allele'],r['alt_allele']]!=[a[at],b[at]]:
                    counts['excluded_unresolved_source_to_physical_allele_orientation']+=1;continue
                value=float(r['Log2Skew']);assert np.isfinite(value)
                assert abs(abs(value)-float(d['effect'][i]))<1e-5
                values[i]=value;counts[line]+=1
                rows.append({'observation_index':i,'bank_index':bank,'cell_line':line,
                  'original_source_path':str(path),'source_row_ID':r['ID'],
                  'source_variant_GRCh37':r['SNP'],'variant_GRCh38':p['variant_GRCh38'],
                  'construct_genomic_orientation':orientation,
                  'signed_ALT_minus_REF_log2_effect':value,'split':int(d['split'][i])})
    assert np.isfinite(values).sum()==len(rows)
    assert not np.isfinite(values[d['split']==2]).any()
    assert np.allclose(abs(values[np.isfinite(values)]),d['effect'][np.isfinite(values)],rtol=1e-5,atol=1e-6)
    old=np.isfinite(d['signed_effect']);assert not (old&np.isfinite(values)).any()
    target=OUT/'original_paired_signed_development.npy';np.save(target,values)
    with gzip.open(OUT/'original_paired_signed_development.csv.gz','wt',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write_json(OUT/'original_paired_signed_recovery_receipt.json',{'passed':True,
       'new_original_signed_development_targets':len(rows),'new_signed_training_rows':int(np.isfinite(values[d['split']==0]).sum()),
       'new_signed_validation_rows':int(np.isfinite(values[d['split']==1]).sum()),'calibration_rows_added':0,
       'by_cell_line':dict(counts),'source_files':sources,'target_sha256':sha(target),
       'CSV_sha256':sha(OUT/'original_paired_signed_development.csv.gz'),
       'original_training_arrays_sha256':sha(OLD/'coverage_training_arrays.npz'),
       'source_software_Git_blob_SHA1':'ce908dacaf0a4ff580e8200bc6cc59fdb0288f0a',
       'source_software_sha256':sha(OUT/'sources/MPRAmodel_e853506.R'),
       'contrast_independent_R_receipt_sha256':sha(OUT/'paired_contrast_independent_R_receipt.tsv'),
       'target_is_original_paired_GLM_Log2Skew':True,'inferred_from_RNA_mean_correlation':False,
       'original_coefficients_independently_refitted':False,'exact_2022_executed_binary_recovered':False,
       'source_scope':'Published paired contrast and unambiguous actual allele mapping, not proof of the exact historical execution environment.',
       'no_significance_or_activity_or_effect_size_or_SE_filter':True,'magnitude_targets_and_partitions_unchanged':True,
       'new_model_performance_or_reserved_Ulirsch_outcomes_read':False})
    assert_prior_versions();print('ORIGINAL SIGNED DEVELOPMENT TARGETS RECOVERED',len(rows),dict(counts),flush=True)


if __name__=='__main__': main()
