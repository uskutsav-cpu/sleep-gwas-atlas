"""Unchanged maintained PlinkBEDFile on a <=288-SNP byte-exact source window."""
import argparse
import importlib.util
import io
from pathlib import Path
import sys
from signed_reference_source_common_v1 import load_frozen,regular_sha,binary_new,save

PIN='df98f324ca775ece05a6e8f93e0584c88bde9535c360fda6ff96bb7f475e726b'


def export(metadata,expected,np,pd):
    spec=load_frozen(metadata,expected)
    if spec['samples']!=503 or not 0<len(spec['rows'])<=288:raise RuntimeError('DECODER_B_BOUNDED503_WINDOW_REQUIRED')
    if regular_sha(spec['bed_path'])!=spec['bed_sha256']:raise RuntimeError('DECODER_B_SOURCE_WINDOW_CHANGED_BEFORE_READ')
    path=Path(spec['maintained_reader_path'])
    if regular_sha(path)!=PIN:raise RuntimeError('MAINTAINED_DECODER_B_CODE_CHANGED')
    loader=importlib.util.spec_from_file_location('maintained_plink_B',path);module=importlib.util.module_from_spec(loader);loader.loader.exec_module(module)
    if regular_sha(path)!=PIN:raise RuntimeError('MAINTAINED_DECODER_B_CHANGED_AFTER_IMPORT')
    class SNPList:
        IDList=[row['SNP'] for row in spec['rows']]
        df=pd.DataFrame([[r['CHR'],r['SNP'],r['BP'],0] for r in spec['rows']],columns=['CHR','SNP','BP','CM'])
    # Invoke unchanged maintained code, including its original filters. They
    # must remove no already-QC-eligible window row, never define the new J.
    obj=module.PlinkBEDFile(spec['bed_path'],503,SNPList(),mafMin=None)
    if obj.m!=len(spec['rows']) or list(obj.kept_snps)!=list(range(len(spec['rows']))):raise RuntimeError('DECODER_B_UNEXPECTED_ORIGINAL_FILTER_REMOVAL')
    raw=np.array(obj.geno.decode(obj._bedcode),dtype=np.int8).reshape((obj.m,obj.nru)).T[:503,:]
    original_y=obj.nextSNPs(obj.m,minorRef=None)
    if not np.all(np.isin(raw,[0,1,2,9])) or not np.all(np.isfinite(original_y)):raise RuntimeError('DECODER_B_NONFINITE_OR_DISCRETE_FAILED')
    if regular_sha(spec['bed_path'])!=spec['bed_sha256']:raise RuntimeError('DECODER_B_WINDOW_CHANGED_AFTER_READ')
    # Source counts A2. Preserve literal raw/standardized B outputs; the A
    # control independently applies A1=2-A2 and target=-Y_B*recorded_sign.
    payload=io.BytesIO();np.savez(payload,A2_dose=raw,standardized_A2=original_y)
    out=Path(spec['output_path']);out_sha=binary_new(out,payload.getvalue())
    receipt={'schema':'maintained_plink_decoderB_export_v1','source_window_sha256':spec['bed_sha256'],
      'metadata_sha256':expected,'maintained_reader_sha256':PIN,'kept_source_rows':list(obj.kept_snps),
      'counts_BIM_A2':True,'minorRef':None,'original_mafMin_default':0,'rows':spec['rows'],
      'samples':503,'output_path':str(out),'output_sha256':out_sha,
      'original_decoder_and_standardization_unchanged':True,'zero_variance_fallback_not_eligible_contract':True,
      'padding_rejection_not_provided_by_original_B':True}
    save(spec['receipt_path'],receipt)
    return receipt


def main():
    p=argparse.ArgumentParser();p.add_argument('--metadata',required=True);p.add_argument('--metadata-sha256',required=True);a=p.parse_args()
    import numpy as np
    import pandas as pd
    import bitarray
    if sys.version_info[:3]!=(3,9,23) or np.__version__!='1.21.5' or pd.__version__!='1.3.3' or bitarray.__version__!='2.8.3':raise RuntimeError('DECODER_B_DECLARED_RUNTIME_VERSION_DIFFERS')
    export(a.metadata,a.metadata_sha256,np,pd)


if __name__=='__main__':main()
