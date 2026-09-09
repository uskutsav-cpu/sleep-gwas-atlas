"""Isolated launcher for the original MATLAB pleioFDR implementation.

Requires independently prepared, reference-index-aligned MAT inputs. This does
not manufacture MAT alignment or reinterpret an exit code as valid conjFDR.
"""
from __future__ import annotations
import os
import shutil
import subprocess
from pathlib import Path
from .io import require, read_json, check_hash, write_json, safe_write_path


def matlab_string(value):
    require(isinstance(value,str) and '\n' not in value and '\r' not in value and '\x00' not in value,
            'Invalid MATLAB path/string')
    return "'"+value.replace("'","''")+"'"


def config_entries(path):
    entries={}
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line or line.startswith('#'):continue
        require('=' in line,'Malformed pleioFDR config')
        key,value=line.split('=',1);key=key.strip();value=value.strip()
        require(key and key not in entries,'Duplicate pleioFDR config key')
        entries[key]=value
    return entries


def run_pleiofdr(c,out):
    require(c.get('alignment_reviewed') is True,'Reference-index MAT alignment must be reviewed')
    require(c.get('configuration_reviewed') is True,'Review official method settings before execution')
    alignment=read_json(c['alignment_manifest'])
    require(alignment.get('reviewed') is True and alignment.get('ancestry')=='EUR',
            'Missing reviewed alignment provenance')
    require(alignment.get('genome_build')==c.get('genome_build'),'Reference genome-build mismatch')
    require(alignment.get('reference_order_sha256') and alignment.get('n_reference_variants',0)>0,
            'Alignment must identify exact reference variant order')
    for name in ['reference_mat','trait1_mat','trait2_mat']:
        require(Path(alignment[name]['path']).resolve()==Path(c[name]).resolve(),'MAT path mismatch')
        check_hash(c[name],alignment[name]['sha256'])
    for name in ['trait1_mat','trait2_mat']:
        require(alignment[name].get('reference_order_sha256')==alignment['reference_order_sha256'],
                'MAT index order not verified against the reference')
    binary=shutil.which(c.get('matlab','matlab'));require(binary,'BLOCKED_BY_SOFTWARE: MATLAB missing')
    upstream=Path(c['repo_root']).resolve();copy=out/'upstream';copy.mkdir()
    # Copy code only. runme changes directory; it must never write into the source checkout.
    sources=list(upstream.rglob('*.m'));require(sources,'Missing official MATLAB source')
    for source in sources:
        require(not source.is_symlink(),'MATLAB code symlink requires explicit review')
        target=safe_write_path(copy,str(source.relative_to(upstream)))
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
    require((copy/'runme.m').is_file(),'Official runme.m missing')
    entries=config_entries(c['config_file'])
    require(entries.get('stattype')=='conjfdr','This adapter is for reviewed conjunction-FDR')
    # Runtime paths are always overridden to preserve source and reference files.
    entries.update(traitfolder='',traitfile1=str(Path(c['trait1_mat']).resolve()),
        traitfiles='{'+matlab_string(str(Path(c['trait2_mat']).resolve()))+'}',
        reffile=str(Path(c['reference_mat']).resolve()),mlibrary=str(copy),
        outputdir=str(out/'results'),onscreen='false',exit_matlab_upon_completion='false',
        randprune_file='',reset_pruneidx='true')
    require(entries.get('dummy_zscore','false')=='false','Do not fabricate unsigned effect directions')
    # Config-sourced refinfo is disallowed unless separately checksum-bound in settings.
    require(not entries.get('refinfo',''),'Supply already complete reference MAT metadata; unbound refinfo forbidden')
    cfgpath=out/'config.resolved.txt'
    cfgpath.write_text('\n'.join(f'{k}={v}' for k,v in entries.items())+'\n')
    driver=out/'driver.m'
    driver.write_text(f"rng({int(c.get('seed',20260908))}, 'twister');\n"
      f"addpath(genpath({matlab_string(str(copy))}));\n"
      f"config={matlab_string(str(cfgpath))};\n"
      f"run({matlab_string(str(copy/'runme.m'))});\n"
      f"fid=fopen({matlab_string(str(out/'matlab_version.txt'))},'w'); fprintf(fid,'%s',version); fclose(fid);\n")
    result=subprocess.run([binary,'-batch',f'run({matlab_string(str(driver))})'],check=False)
    require(result.returncode==0,'Official MATLAB pleioFDR execution failed')
    require(c.get('expected_outputs'),'Declare expected native output filenames')
    for relative in c['expected_outputs']:
        require(safe_write_path(out/'results',relative).is_file(),f'Missing MATLAB output: {relative}')
    write_json(out/'status.json',{'status':'INSUFFICIENT_EVIDENCE',
        'reason':'NATIVE_OUTPUTS_REQUIRE_SCIENTIFIC_QC',
        'note':'Official code completed; conjunction-FDR calibration and locus tables require review before consumption.'})
