"""Bind native adapters to the data, references, and source files they actually read."""
from __future__ import annotations
import csv
from pathlib import Path
from .io import require, file_record, check_hash

DIRECT = {
 'placo': ('placo_source','pair_file','chunk_file','parameters'),
 'susie': ('locus_file','ld_file'), 'coloc': ('fit1','fit2'),
 'lava': ('input_info','sample_overlap_file','loci_file'),
 'genomicsem': ('covariance_file',), 'mr': ('exposure','outcome'),
 'munge': ('munge_script','sumstats','hapmap3'),
 'h2': ('ldsc_script','sumstats'), 'rg': ('ldsc_script','sumstats1','sumstats2'),
 'clump': ('associations',), 'ld': (),
 'pleiofdr': ('config_file','reference_mat','trait1_mat','trait2_mat','alignment_manifest')
}

def discover_inputs(method, settings):
    """Files are checksummed once while preparing a job, then reverified on execution.

    Prefix-based tools read many files; the three PLINK members and chromosome
    LD-score families must not be represented by a meaningless prefix string.
    """
    require(method in DIRECT, 'Unknown native method for data binding')
    files=set()
    def add(value):
        require(isinstance(value,str) and value and value!='UNRESOLVED','Unresolved native input')
        p=Path(value).resolve();require(p.is_file(),f'Missing native input: {p}');files.add(p)
    for key in DIRECT[method]:
        if settings.get(key) is not None:add(settings[key])
    if method=='genomicsem' and settings.get('mode')=='covariance':
        require(settings.get('traits'),'No GenomicSEM trait inputs')
        for value in settings['traits']:add(value)
        for key in ['ld_directory','weights_directory']:
            p=Path(settings[key]);require(p.is_dir(),f'Missing LD directory: {p}')
            matches=[x for x in p.iterdir() if x.is_file() and not x.name.startswith('.')]
            require(matches,f'Empty LD directory: {p}')
            for x in matches:add(str(x))
    if method in {'lava','clump','ld'}:
        require(settings.get('reference_prefix') not in {None,'UNRESOLVED'},'Reference prefix unresolved')
        if method=='lava' and settings.get('reference_format')=='LAVA_CUSTOM':
            require(settings.get('reference_files'), 'Custom LAVA reference needs exact file inventory')
            for record in settings['reference_files']:
                check_hash(record['path'],record['sha256']);add(record['path'])
        else:
            for suffix in ['.bed','.bim','.fam']:add(settings['reference_prefix']+suffix)
    if method in {'h2','rg'}:
        for key in ['ref_ld_chr','weights_chr']:
            prefix=settings[key];require(isinstance(prefix,str) and prefix!='UNRESOLVED','LD-score prefix unresolved')
            for chromosome in range(1,23):
                name=f'{prefix}{chromosome}.l2.ldscore.gz';add(name)
            # LD-score variant counts are also used by LDSC. Preserve all present variants.
            for chromosome in range(1,23):
                for suffix in ['.l2.M','.l2.M_5_50']:
                    p=Path(f'{prefix}{chromosome}{suffix}')
                    if p.is_file():add(str(p))
    if method=='lava':
        # LAVA's info table contains data-file paths; bind each referenced file.
        with Path(settings['input_info']).open() as f:
            lines=[line.split() for line in f if line.strip()]
        require(lines and 'filename' in lines[0],'LAVA input info needs filename column')
        i=lines[0].index('filename')
        for row in lines[1:]:
            require(len(row)==len(lines[0]),'Malformed LAVA input info');add(row[i])
    if method in {'munge','h2','rg'}:
        base=Path(settings.get('ldsc_script',settings.get('munge_script'))).resolve().parent
        for p in (base/'ldscore').glob('*.py'):add(str(p))
    if method=='pleiofdr':
        base=Path(settings['repo_root']).resolve();require(base.is_dir(),'Missing official pleioFDR checkout')
        require((base/'runme.m').is_file(),'Official pleioFDR runme.m missing')
        for p in base.rglob('*.m'):
            require(not p.is_symlink(),'Upstream MATLAB source symlinks must be reviewed explicitly');add(str(p))
    return {f'native_{i:05d}':file_record(p) for i,p in enumerate(sorted(files))}
