import csv,hashlib,json,subprocess,pathlib,collections
p=pathlib.Path(__file__).resolve().parent;src=p.parent/'candidate_positions.tsv';rows=list(csv.DictReader(src.open(),delimiter='\t'));assert len(rows)==16
base='https://grch37.rest.ensembl.org/overlap/region/human/'; out=[];receipt=[]
for chrom in ['7','18','22']:
 a=[r for r in rows if r['CHR']==chrom];positions=[int(r['BP_GRCH37']) for r in a];lo=min(positions)-250000;hi=max(positions)+250000;url=f'{base}{chrom}:{lo}-{hi}?feature=gene'
 run=subprocess.run(['curl','-fsSL','--max-time','30','-H','Content-Type: application/json',url],check=True,capture_output=True);raw=run.stdout;genes=json.loads(raw);assert isinstance(genes,list)
 assert all(str(g['seq_region_name'])==chrom and g['assembly_name']=='GRCh37' for g in genes)
 (p/f'chr{chrom}_response.json').write_bytes(raw);receipt.append({'chrom':chrom,'range':f'{lo}-{hi}','url':url,'response_sha256':hashlib.sha256(raw).hexdigest(),'n_genes':len(genes)})
 for r in a:
  bp=int(r['BP_GRCH37']);dist=lambda g:max(g['start']-bp,0,bp-g['end'])
  overlap=[g for g in genes if dist(g)==0]; coding=[g for g in genes if g.get('biotype')=='protein_coding']; mind=min((dist(g) for g in coding),default=None);near=[g for g in coding if dist(g)==mind] if mind is not None else []
  pack=lambda gs:';'.join((g.get('external_name') or g['id'])+'|'+g['id']+'|'+g.get('biotype','') for g in gs)
  out.append({'SNP':r['SNP'],'chr37':chrom,'bp37':bp,'overlap_genes':pack(overlap),'nearest_coding_gene':pack(near),'nearest_coding_distance_bp':mind if mind is not None else 'NOT_FOUND_IN_WINDOW'})
with (p/'marker_context.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'candidate_table_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'queries':receipt,'notes':'GRCh37 Ensembl reference annotation, position only, not effector-gene evidence'},indent=2)+'\n')
for row in out:print(row)
