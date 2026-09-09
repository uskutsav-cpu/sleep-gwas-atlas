"""Build a complete data-preparation DAG from the actual reviewed pair lock."""
from pathlib import Path
from .pairs import verify_lock
from .io import require, read_json, file_record, write_json


def assemble_preparation(pair_lock, sources_file, output_root, runtime_out, *, reviewed=False):
    require(reviewed, 'Review source identities and pair decisions before building the DAG')
    lock=verify_lock(pair_lock);sources=read_json(sources_file)
    needed=sorted({p[k] for p in lock['pairs'] for k in ['sleep_trait','disease_trait']})
    require(set(needed)<=set(sources), f'Missing source cards: {set(needed)-set(sources)}')
    tasks=[]
    for trait in needed:
        path=Path(sources[trait]).resolve();card=read_json(path)
        require(card['trait_id']==trait,'Source card trait does not match its registry key')
        require(card.get('synthetic',False)==lock['synthetic'],'Synthetic source contamination')
        tasks.append({'id':'normalize_'+trait,'kind':'normalize','parameters':{'source':str(path)}})
    for pair in lock['pairs']:
        identity=pair['pair_id'];s=pair['sleep_trait'];d=pair['disease_trait']
        name='pair_'+identity
        tasks.append({'id':name,'kind':'join','pair_id':identity,
            'depends_on':['normalize_'+s,'normalize_'+d],
            'parameters':{'left':'@normalize_'+s,'right':'@normalize_'+d}})
        tasks.append({'id':'index_'+identity,'kind':'index','pair_id':identity,
                      'depends_on':[name],'parameters':{'pair_dir':'@'+name}})
        if identity!='insomnia__adhd':
            tasks.append({'id':'chunks_'+identity,'kind':'split','pair_id':identity,
                          'depends_on':[name],'parameters':{'pair_dir':'@'+name,'chunk_rows':20000}})
    value={'schema_version':1,'reviewed':True,'synthetic':lock['synthetic'],
           'pair_lock':str(Path(pair_lock).resolve()),'source_registry':file_record(sources_file),
           'output_root':str(Path(output_root).resolve()),'tasks':tasks,
           'note':'Preparation only. No rerun of legacy ADHD PLACO. Native production requires independently reviewed jobs.'}
    write_json(runtime_out,value);return value
