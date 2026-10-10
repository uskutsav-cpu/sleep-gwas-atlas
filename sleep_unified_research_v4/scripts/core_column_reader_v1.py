"""Preserve native whole-column inference before bounded core harmonization.

This is a preparation candidate. It supplies no scientific admission, new
filter, dtype emulator, or executable campaign. Its caller must hold the
reviewed owned-worker guard and private SSD namespace for every operation.
"""
import csv
import gc
import gzip
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

ROWS=50_000
NUMERIC={'CHR','BP','FRQ','BETA','SE','P','N','N_EFF','N_EFF_HALF','NCASE','NCONTROL','INFO'}


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def dump_new(path,value):
    with Path(path).open('xb') as f:pickle.dump(value,f,protocol=5)


def load(path):
    # All pickle files are private, produced by this candidate and hash-bound.
    with Path(path).open('rb') as f:return pickle.load(f)


def c_shape_gate(path,read_kwargs,expected_rows):
    """Reject unsupported physical shapes; do not silently alter parser rows."""
    sep=read_kwargs['sep']
    if sep not in ['\t',',']:raise RuntimeError('ONLY_ORIGINAL_TAB_OR_COMMA_C_SOURCE_SUPPORTED')
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='') as f:
        for _ in range(read_kwargs['skiprows']):next(f)
        reader=csv.reader(f,delimiter=sep)
        header=next(reader);width=len(header);rows=0;previous=reader.line_num
        if not header or len(set(header))!=width:raise RuntimeError('UNIQUE_COMPLETE_C_HEADER_REQUIRED')
        for row in reader:
            current=reader.line_num
            if current-previous!=1:raise RuntimeError('MULTILINE_SOURCE_RECORD_NOT_QUALIFIED')
            previous=current
            if not row:continue
            if len(row)!=width:raise RuntimeError('C_SOURCE_IMPLICIT_INDEX_OR_RAGGED_RECORD_NOT_QUALIFIED')
            rows+=1
        if rows!=expected_rows:raise RuntimeError('FULL_C_SOURCE_SHAPE_ROW_COUNT_DIFFERS')
    return {'header':header,'full_header_width':width,'rows':rows,'raw_gzip_EOF_and_CRC_reached':True}


def python_tokens(path,read_kwargs,spool,expected_rows):
    """Spool original PythonParser tokens before any dtype/NA conversion."""
    kwargs=dict(read_kwargs);kwargs['iterator']=True
    reader=pd.read_csv(path,**kwargs);engine=reader._engine
    if type(engine).__name__!='PythonParser' or engine._implicit_index or engine.index_col:
        reader.close();raise RuntimeError('ORIGINAL_PYTHON_PARSER_WITHOUT_IMPLICIT_INDEX_REQUIRED')
    pieces={column:[] for column in read_kwargs['usecols']};rows=0
    try:
        while True:
            try:content=engine._get_lines(ROWS)
            except StopIteration:break
            if not content:continue
            data,names=engine._exclude_implicit_index(engine._rows_to_cols(content))
            if set(data)!=set(pieces):raise RuntimeError('PYTHON_PREINFERENCE_COLUMN_SET_DIFFERS')
            lengths={len(v) for v in data.values()}
            if lengths!={len(content)}:raise RuntimeError('PYTHON_PREINFERENCE_ROW_ALIGNMENT_DIFFERS')
            for i,column in enumerate(pieces):
                dest=spool/('tokens_%02d_%08d.pkl'%(i,rows))
                dump_new(dest,data[column]);pieces[column].append({'path':str(dest),'sha256':digest(dest),'rows':len(content)})
            rows+=len(content)
        if rows!=expected_rows:raise RuntimeError('PYTHON_PREINFERENCE_FULL_SOURCE_COUNT_DIFFERS')
        # Conversion receives one complete original token array, preserving the
        # Python engine's global fallback, bool/numeric interaction and NA rules.
        for column,parts in pieces.items():
            values=np.empty(rows,dtype=object);start=0
            for part in parts:
                if digest(part['path'])!=part['sha256']:raise RuntimeError('PRIVATE_PREINFERENCE_SPOOL_CHANGED')
                chunk=load(part['path']);values[start:start+len(chunk)]=chunk;start+=len(chunk)
                del chunk
            converted=engine._convert_data({column:values})
            if set(converted)!={column}:raise RuntimeError('ORIGINAL_SINGLETON_PYTHON_CONVERSION_DIFFERS')
            yield column,pd.Series(converted[column]),parts
            del converted,values;gc.collect()
    finally:reader.close()


def normalized(series,canonical):
    """The original whole-column conversions, before any QC decision."""
    if canonical=='OR':
        odds_ratio=pd.to_numeric(series,errors='coerce')
        return np.log(odds_ratio.where(odds_ratio>0))
    if canonical in NUMERIC:series=pd.to_numeric(series,errors='coerce')
    if canonical in ['A1','A2','SNP']:
        series=series.astype('string').str.strip()
        series=series.str.lower() if canonical=='SNP' else series.str.upper()
    if canonical=='CHR':
        series=series.astype('string').str.replace('chr','',case=False,regex=False)
        series=pd.to_numeric(series,errors='coerce')
    return series


def freeze_columns(path,read_kwargs,columns,spool,expected_rows):
    """Freeze typed 50k pieces after unchanged complete inference per column."""
    if pd.__version__!='2.2.3' or np.__version__!='1.26.4':raise RuntimeError('DECLARED_CORE_PARSE_RUNTIME_REQUIRED')
    spool=Path(spool);spool.mkdir(parents=True,exist_ok=False)
    if list(read_kwargs['usecols'])!=list(dict.fromkeys(columns.values())):raise RuntimeError('EXACT_ORIGINAL_SELECTED_COLUMN_ORDER_REQUIRED')
    source_to_standard={original:standard for standard,original in columns.items()}
    if len(source_to_standard)!=len(columns):raise RuntimeError('ONE_TO_ONE_ORIGINAL_COLUMN_MAP_REQUIRED')
    if 'BETA' in columns and 'OR' in columns:raise RuntimeError('ORIGINAL_BETA_PREFERENCE_NOT_APPLIED')
    python=read_kwargs.get('engine')=='python'
    shape=None if python else c_shape_gate(path,read_kwargs,expected_rows)
    def c_columns():
        for column in read_kwargs['usecols']:
            kwargs=dict(read_kwargs);kwargs['usecols']=[column]
            data=pd.read_csv(path,**kwargs)
            if not data.index.equals(pd.RangeIndex(expected_rows)) or list(data.columns)!=[column]:raise RuntimeError('FULL_NATIVE_C_COLUMN_SHAPE_DIFFERS')
            # The C parser's low_memory inference block size uses full physical
            # table width. Complete native read avoids converted-chunk merging.
            yield column,data[column],[]
            del data;gc.collect()
    result={'source':str(path),'source_sha256_before':digest(path),'expected_rows':expected_rows,'read_kwargs':read_kwargs,'columns':columns,'physical_shape_gate':shape,'typed_pieces':{},'global_dtype':{},'original_BUILD_values':[],'preinference_token_pieces':{},'status':'PREPARATION_CANDIDATE_NOT_ADMITTED'}
    for index,(column,series,token_parts) in enumerate(python_tokens(path,read_kwargs,spool,expected_rows) if python else c_columns()):
        standard=source_to_standard[column]
        if standard=='BUILD':result['original_BUILD_values']=[str(x) for x in series.dropna().unique() if str(x).strip()]
        series=normalized(series,standard)
        target='BETA' if standard=='OR' else standard
        result['global_dtype'][target]=str(series.dtype);parts=[]
        for start in range(0,expected_rows,ROWS):
            piece=spool/('typed_%02d_%08d.pkl'%(index,start));chunk=series.iloc[start:start+ROWS]
            dump_new(piece,chunk)
            parts.append({'path':str(piece),'sha256':digest(piece),'start':start,'rows':len(chunk),'dtype':str(chunk.dtype)})
        result['typed_pieces'][target]=parts;result['preinference_token_pieces'][target]=token_parts
        del series;gc.collect()
    result['source_sha256_after']=digest(path)
    if result['source_sha256_before']!=result['source_sha256_after']:raise RuntimeError('CORE_RAW_SOURCE_CHANGED_DURING_COLUMN_PREPARATION')
    out=spool/'column_preparation_manifest.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    return result
