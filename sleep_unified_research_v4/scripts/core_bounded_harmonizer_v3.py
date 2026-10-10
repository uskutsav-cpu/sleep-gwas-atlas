#!/usr/bin/env python3
"""Candidate bounded orchestration of the frozen original harmonizer AST.

No production admission is supplied by this module. A separately reviewed
controller must own the complete SSD namespace, source identities, runtime,
resource supervision, process groups, and final template/QC comparison.
"""
import argparse
import ast
from collections import OrderedDict
import copy
import gc
import gzip
import importlib.util
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pandas as pd

from core_column_reader_v4 import digest, dump_new, freeze_columns, load, ROWS

ORIGINAL_SHA = '487074007f8361e7de7f8f771dd0a912f926d4b7a0427cfa3c47b02be1421f7a'
N_FIELDS = ['N', 'N_EFF', 'N_EFF_HALF', 'NCASE', 'NCONTROL']


class DiskSet:
    """Exact membership for original rsID strings or coordinate-key tuples."""
    def __init__(self, path):
        if Path(path).exists():raise RuntimeError('PRIVATE_DISK_SET_ALREADY_EXISTS')
        self.db=sqlite3.connect(path)
        self.db.execute('PRAGMA cache_size=-4096')
        self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('CREATE TABLE keys (key BLOB PRIMARY KEY) WITHOUT ROWID')
    @staticmethod
    def encoded(value):
        import pickle
        return pickle.dumps(value,protocol=5)
    def update(self, values):
        self.db.executemany('INSERT OR IGNORE INTO keys VALUES (?)',((self.encoded(v),) for v in values))
        self.db.commit()
    def __contains__(self, value):
        return self.db.execute('SELECT 1 FROM keys WHERE key=?',(self.encoded(value),)).fetchone() is not None
    def __len__(self):return self.db.execute('SELECT COUNT(*) FROM keys').fetchone()[0]
    def close(self):self.db.close()


class FirstSurvivingSNP:
    """Original first duplicate is chosen after ordinary QC and before N QC."""
    def __init__(self,path):
        if Path(path).exists():raise RuntimeError('PRIVATE_DUPLICATE_INDEX_ALREADY_EXISTS')
        self.db=sqlite3.connect(path)
        self.db.execute('PRAGMA cache_size=-4096')
        self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('CREATE TABLE snps (snp TEXT COLLATE BINARY PRIMARY KEY) WITHOUT ROWID')
    def mask(self,series):
        duplicate=[]
        for value in series:
            if not isinstance(value,str):raise RuntimeError('PRE_N_QC_RSID_MUST_BE_NORMALIZED_STRING')
            cursor=self.db.execute('INSERT OR IGNORE INTO snps VALUES (?)',(value,))
            duplicate.append(cursor.rowcount==0)
        self.db.commit()
        return pd.Series(duplicate,index=series.index,dtype=bool)
    def close(self):self.db.close()


def checked_piece(item):
    path=Path(item['path'])
    if path.is_symlink() or not path.is_file() or digest(path)!=item['sha256']:
        raise RuntimeError('PRIVATE_TYPED_FRAME_OR_COLUMN_CHANGED')
    return load(path)


def frames(columns):
    pieces=columns['typed_pieces'];counts={len(v) for v in pieces.values()}
    if len(counts)!=1:raise RuntimeError('WHOLE_COLUMN_TYPED_PIECE_COUNTS_DIFFER')
    for index in range(next(iter(counts))):
        data={name:checked_piece(items[index]) for name,items in pieces.items()}
        frame=pd.DataFrame(data)
        start=index*ROWS;end=min(start+ROWS,columns['expected_rows'])
        if not frame.index.equals(pd.RangeIndex(start,end)):
            raise RuntimeError('ORIGINAL_SOURCE_ROW_ORDER_OR_INDEX_DIFFERS')
        for name,dtype in columns['global_dtype'].items():
            if str(frame[name].dtype)!=dtype:raise RuntimeError('GLOBAL_TYPED_COLUMN_DTYPE_CHANGED')
        if 'SNP' not in frame:frame['SNP']=pd.Series(pd.NA,index=frame.index,dtype='string')
        for name in ['CHR','BP']:
            if name not in frame:frame[name]=np.nan
        yield frame


def call_name(node,name):
    return isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id==name


def load_original(path):
    path=Path(path)
    if path.is_symlink() or digest(path)!=ORIGINAL_SHA:raise RuntimeError('FROZEN_ORIGINAL_HARMONIZER_CHANGED')
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location('_original_harmonizer_science_v1',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    tree=ast.parse(path.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    return module,main.body


def compiled_stage(namespace,name,parameters,body,drop_node=None,tail='return out, steps'):
    prefix=ast.parse('steps=[]').body
    if drop_node is not None:prefix.append(copy.deepcopy(drop_node))
    function=ast.FunctionDef(name=name,args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=p) for p in parameters],kwonlyargs=[],kw_defaults=[],defaults=[]),body=prefix+copy.deepcopy(body)+ast.parse(tail).body,decorator_list=[])
    tree=ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    exec(compile(tree,'<frozen-original-'+name+'>','exec'),namespace)
    return namespace[name]


def original_stages(module,body):
    namespace=dict(vars(module))
    raw_at=next(i for i,n in enumerate(body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='raw' for t in n.targets))
    drop_at=next(i for i,n in enumerate(body) if isinstance(n,ast.FunctionDef) and n.name=='drop')
    life_at=next(i for i,n in enumerate(body) if i>drop_at and isinstance(n,ast.If) and ast.unparse(n.test)=='args.liftover_chain')
    map_at=next(i for i,n in enumerate(body) if i>life_at and isinstance(n,ast.If) and ast.unparse(n.test)=='args.variant_map')
    duplicate_at=next(i for i,n in enumerate(body) if isinstance(n,ast.Expr) and call_name(n.value,'drop') and len(n.value.args)>1 and isinstance(n.value.args[1],ast.Constant) and n.value.args[1].value=='duplicate SNP ID')
    binary_at=next(i for i,n in enumerate(body) if i>duplicate_at and isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='binary' for t in n.targets))
    mkdir_at=next(i for i,n in enumerate(body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='os.makedirs')
    serial_at=next(i for i,n in enumerate(body) if i>binary_at and isinstance(n,ast.For) and ast.unparse(n.target)=='column')
    prefix=ast.fix_missing_locations(ast.Module(body=copy.deepcopy(body[:raw_at]),type_ignores=[]))
    drop_node=body[drop_at]
    life=body[life_at]
    if not isinstance(life.body[0],ast.Try):raise RuntimeError('ORIGINAL_LIFTOVER_LOAD_AST_LAYOUT_CHANGED')
    stage_life=compiled_stage(namespace,'stage_life',['out','args','chain_index','liftover_provenance'],life.body[1:],drop_node,'return out, steps, reverse_strand_count')
    mapping=body[map_at]
    if not isinstance(mapping.body[0],ast.If) or not isinstance(mapping.body[1],ast.Try):raise RuntimeError('ORIGINAL_MAP_LOAD_AST_LAYOUT_CHANGED')
    required=compiled_stage(namespace,'stage_required_keys',['out','args'],[mapping.body[0]],tail='return required_map_keys')
    stage_map=compiled_stage(namespace,'stage_map',['out','args','variant_index','variant_map_provenance'],mapping.body[2:],drop_node)
    ordinary=copy.deepcopy(body[map_at+1:duplicate_at+1])
    original_duplicate=ast.dump(ordinary[-1],include_attributes=False)
    ordinary[-1].value.args[0]=ast.parse('duplicate_mask(out["SNP"])',mode='eval').body
    stage_qc=compiled_stage(namespace,'stage_qc',['out','duplicate_mask'],ordinary,drop_node)
    stage_n=compiled_stage(namespace,'stage_n',['out','metadata','columns'],body[binary_at:serial_at],drop_node)
    stage_serialize=compiled_stage(namespace,'stage_serialize',['out'],body[serial_at:mkdir_at],tail='return out')
    report=[copy.deepcopy(n) for n in body[mkdir_at:] if not (isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='out.to_csv')]
    return dict(namespace=namespace,prefix=compile(prefix,'<frozen-original-schema>','exec'),life=stage_life,required=required,mapping=stage_map,qc=stage_qc,sample=stage_n,serialize=stage_serialize,report=compile(ast.fix_missing_locations(ast.Module(body=report,type_ignores=[])),'<frozen-original-QC-report>','exec'),original_duplicate_AST=original_duplicate)


def aggregate_steps(all_steps):
    """Add identical ordered per-frame losses; cumulative totals are global."""
    totals=None
    for steps in all_steps:
        if totals is None:totals=[[reason,0] for reason,_,_ in steps]
        if [s[0] for s in steps]!=[t[0] for t in totals]:raise RuntimeError('ORIGINAL_ORDERED_QC_STEP_LABELS_DIFFER_ACROSS_FRAMES')
        for target,(_,removed,_) in zip(totals,steps):target[1]+=removed
    return totals or []


def save_frame(folder,index,frame):
    path=folder/('%08d.pkl'%index);dump_new(path,frame)
    return {'path':str(path),'sha256':digest(path),'rows':len(frame)}


class LengthView:
    def __init__(self,n):self.n=n
    def __len__(self):return self.n


def harmonize(original,original_argv,spool,expected_rows,expected_source_sha256):
    """Candidate only: no resource admission, source replacement or filtering fix."""
    # The schema prefix sniffs and reads the raw header. Bind the source before
    # that first read, using the exact --infile token in frozen original commands.
    source_arguments=[]
    for index,token in enumerate(original_argv):
        if token=='--infile':
            if index+1>=len(original_argv):raise RuntimeError('EXACT_ORIGINAL_INFILE_VALUE_REQUIRED')
            source_arguments.append(original_argv[index+1])
        elif token.startswith('--infile='):source_arguments.append(token.split('=',1)[1])
    if len(source_arguments)!=1 or not source_arguments[0]:raise RuntimeError('EXACT_SINGLE_ORIGINAL_INFILE_ARGUMENT_REQUIRED')
    sealed_source=Path(source_arguments[0])
    if sealed_source.is_symlink() or not sealed_source.is_file() or digest(sealed_source)!=expected_source_sha256:raise RuntimeError('SEALED_CORE_RAW_SOURCE_IDENTITY_DIFFERS_BEFORE_SCHEMA_PREFIX')
    module,body=load_original(original);stages=original_stages(module,body);context=stages['namespace']
    old_argv=sys.argv;sys.argv=[str(original),*original_argv]
    try:exec(stages['prefix'],context)
    finally:sys.argv=old_argv
    args=context['args'];metadata=context['metadata'];columns=context['columns']
    if Path(args.infile)!=sealed_source:raise RuntimeError('PARSED_ORIGINAL_INFILE_DIFFERS_FROM_PREPREFIX_SOURCE_BINDING')
    spool=Path(spool);spool.mkdir(parents=True,exist_ok=False)
    parsed=args.variant_map_strategy=='BY_COORD_ALLELES' and ('CHR' not in columns or 'BP' not in columns)
    if parsed and ('CHR' in columns or 'BP' in columns):module.fail('coordinate mapping found only one of CHR/BP; inspect the source schema')
    prepared=freeze_columns(args.infile,context['read_kwargs'],columns,spool/'columns',expected_rows,expected_source_sha256,parse_coordinate_labels=parsed)
    file_builds={module.normalise_build(v) for v in prepared['original_BUILD_values']}
    expected_build='hg38' if context['liftover_sets_build'] else module.TARGET_BUILD
    if file_builds and file_builds!={expected_build}:module.fail('file build values do not equal declared original build')
    life_folder=spool/'post_liftover';life_folder.mkdir();before_n_folder=spool/'post_ordinary_QC';before_n_folder.mkdir()
    chain_index=None;life_provenance=None
    if args.liftover_chain:chain_index,life_provenance=module.load_chain(args.liftover_chain,args.expected_liftover_chain_sha256,args.expected_liftover_chain_bytes)
    required=DiskSet(spool/'required_map_keys.sqlite3') if args.variant_map else None
    life_parts=[];life_steps=[];reverse_count=0
    for index,frame in enumerate(frames(prepared)):
        if args.liftover_chain:
            frame,steps,reversed_here=stages['life'](frame,args,chain_index,life_provenance)
            reverse_count+=reversed_here;life_steps.append(steps[:4])
        if required is not None:required.update(stages['required'](frame,args))
        life_parts.append(save_frame(life_folder,index,frame));del frame;gc.collect()
    variant_index=None;map_provenance=None
    if args.variant_map:variant_index,map_provenance=module.load_variant_map(args.variant_map,args.variant_map_strategy,args.expected_variant_map_sha256,args.expected_variant_map_bytes,required_keys=required)
    if required is not None:required_count=len(required);required.close()
    else:required_count=0
    duplicate=FirstSurvivingSNP(spool/'first_surviving_SNP.sqlite3');map_steps=[];qc_steps=[];before_n=[]
    for index,part in enumerate(life_parts):
        frame=checked_piece(part)
        if args.variant_map:
            frame,steps=stages['mapping'](frame,args,variant_index,map_provenance);map_steps.append(steps)
        frame,steps=stages['qc'](frame,duplicate.mask);qc_steps.append(steps)
        before_n.append(save_frame(before_n_folder,index,frame));del frame;gc.collect()
    duplicate.close();del variant_index,chain_index;gc.collect()
    fields=[name for name in N_FIELDS if name in prepared['typed_pieces']]
    # Native inference/normalization was already global. Only these numerical
    # columns are combined for the exact original global N branch/any/max.
    n_parts=[]
    for part in before_n:
        frame=checked_piece(part);n_parts.append(frame[fields]);del frame
    n_frame=pd.concat(n_parts,axis=0);del n_parts;gc.collect()
    n_frame,n_steps=stages['sample'](n_frame,metadata,columns);n_values=n_frame['N'];del n_frame;gc.collect()
    output_path=Path(args.outdir)/(args.trait+'.harmonized.tsv.gz')
    output_path.parent.mkdir(parents=True,exist_ok=True)
    if output_path.exists() or output_path.is_symlink():raise RuntimeError('PRIOR_HARMONIZED_OUTPUT_PRESERVED')
    rows_out=0
    with gzip.open(output_path,'xt',newline='') as output:
        first=True
        for part in before_n:
            frame=checked_piece(part);frame=frame.loc[frame.index.isin(n_values.index)].copy()
            frame['N']=n_values.loc[frame.index]
            frame=stages['serialize'](frame)
            frame.to_csv(output,sep='\t',index=False,na_rep='NA',header=first)
            first=False;rows_out+=len(frame);del frame;gc.collect()
    totals=[]
    if parsed:totals.append(['parsed CHR/BP from coordinate-based variant label',0])
    totals.extend(aggregate_steps(life_steps))
    if args.liftover_chain:
        totals.extend([[f'reverse-strand liftover mappings with allele complements: {reverse_count}',0],[f"coordinates lifted hg38 to hg19; chain_sha256={life_provenance['chain_sha256']}",0]])
    totals.extend(aggregate_steps(map_steps));totals.extend(aggregate_steps(qc_steps));totals.extend([[reason,removed] for reason,removed,_ in n_steps])
    remaining=expected_rows;steps=[]
    for reason,removed in totals:remaining-=removed;steps.append((reason,removed,remaining))
    if remaining!=rows_out or rows_out!=len(n_values):raise RuntimeError('GLOBAL_ORIGINAL_ORDERED_QC_TOTALS_OR_N_ROW_COUNT_DIFFERS')
    if digest(args.infile)!=expected_source_sha256:raise RuntimeError('ORIGINAL_RAW_CHANGED_AFTER_GLOBAL_HARMONIZATION')
    context.update(n_input=expected_rows,out=LengthView(rows_out),steps=steps,variant_map_provenance=map_provenance,liftover_provenance=life_provenance)
    exec(stages['report'],context)
    receipt={'status':'CANDIDATE_GLOBAL_HARMONIZATION_PRODUCED_NOT_SCIENTIFICALLY_ADMITTED','original_code_sha256':ORIGINAL_SHA,'source_sha256':expected_source_sha256,'expected_rows':expected_rows,'rows_out':rows_out,'coordinate_label_parsed':parsed,'global_required_map_key_count':required_count,'global_reverse_strand_count':reverse_count,'global_sample_size_rows':len(n_values),'original_QC_steps':steps,'output_sha256':digest(output_path),'column_manifest_sha256':digest(spool/'columns/column_preparation_manifest.json'),'post_liftover_frames':life_parts,'post_ordinary_QC_frames':before_n,'duplicate_rule':'FIRST_SURVIVING_AFTER_ORDINARY_QC_BEFORE_GLOBAL_N_QC','full_original_output_template_and_QC_comparison_required':True}
    with (spool/'global_harmonization_candidate_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    return receipt


def main():
    p=argparse.ArgumentParser();p.add_argument('--original',type=Path,required=True);p.add_argument('--spool',type=Path,required=True);p.add_argument('--expected-rows',type=int,required=True);p.add_argument('--expected-source-sha256',required=True);p.add_argument('original_argv',nargs=argparse.REMAINDER)
    a=p.parse_args();argv=a.original_argv[1:] if a.original_argv[:1]==['--'] else a.original_argv
    harmonize(a.original,argv,a.spool,a.expected_rows,a.expected_source_sha256)


if __name__=='__main__':main()
