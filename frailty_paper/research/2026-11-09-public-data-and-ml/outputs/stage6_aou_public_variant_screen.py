#!/usr/bin/env python3
"""Read selected-variant public association endpoint; absence is not a null."""
import csv,hashlib,json,re,subprocess,time
from pathlib import Path
from urllib.parse import quote
p=Path(__file__).resolve().parent
def read(n):
 with (p/n).open(newline='') as f:return list(csv.DictReader(f))
score=read('unvalidated_22_shared_direction_score.csv'); inv={r['variant']:r for r in read('stage6_tpmi_shared22_insomnia_inventory.csv')};atlas=read('stage5_functional_atlas_with_priority_splice.csv')
def query(vid,phen):
 u='https://allbyall.researchallofus.org/api/variants/associations/variant/'+quote(vid,safe='-')+'?analysis_id='+quote(phen,safe='._')
 for attempt in range(3):
  r=subprocess.run(['curl','-fsS','--max-time','20',u],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  if r.returncode==0:
   x=json.loads(r.stdout);assert isinstance(x.get('data'),list) and x.get('count')==len(x['data']);return x['data']
  time.sleep(1+attempt)
 raise RuntimeError('public lookup failed: '+vid+' '+phen)
controls=[]
for v,phen in [('chr3-53494203-G-GA','NS_333.2'),('chr8-40846505-CT-C','SS_807.21')]:
 r=query(v,phen);assert len(r)==1 and r[0]['phenotype']==phen and r[0]['pvalue']<5e-8;controls.append(r[0])
rows=[]
for item in score:
 name=item['variant'];i=inv[name]
 if i.get('pos38') and i.get('queue_identity_match')=='True':
  contig='chr'+i['chr38'];position=i['pos38'];ref=i['ref'];alt=i['alt'];mapping='previous TPMI/queue identity check'
 else:
  hits={a['variant_GRCh38_ref_alt'] for a in atlas if a['input_classes']=='exploratory_marker' and a['input_labels']==name}
  if len(hits)!=1:raise ValueError('ambiguous GRCh38 mapping '+name)
  m=re.fullmatch(r'(chr\d+):(\d+):([ACGT]+)>([ACGT]+)',next(iter(hits)))
  if not m:raise ValueError('invalid mapped allele '+name)
  contig,position,ref,alt=m.groups();mapping='previous functional atlas identity'
 if {ref,alt}!={item['counted_allele'],item['other_allele']}:raise ValueError('pair or strand mismatch '+name)
 vid='-'.join([contig,position,ref,alt])
 for phen in ('NS_333.2','SS_807.21'):
  hits=query(vid,phen)
  if len(hits)>1:raise ValueError('unexpected non-unique association '+vid)
  d=hits[0] if hits else {};row={'variant':name,'variant_grch38':vid,'counted_allele':item['counted_allele'],'other_allele':item['other_allele'],'mapping_source':mapping,'phenotype_id':phen,'status':'selected_association_returned' if hits else 'no_row_in_public_selected_association_lookup','reported_ancestry':d.get('ancestry',''),'sequencing_type':d.get('sequencing_type',''),'pvalue':d.get('pvalue',''),'beta_source':d.get('beta',''),'se':d.get('se',''),'af':d.get('af','')};rows.append(row)
with (p/'stage6_aou_public_variant_screen.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
report={'source':'https://allbyall.researchallofus.org/','endpoint':'/api/variants/associations/variant/{variant_id}?analysis_id={phenotype}','browser_version_observed':'0.1.0-317fe0c-202609021906','server_source_version':'https://github.com/broadinstitute/all-by-all-aou-browser/blob/317fe0c/axaou-server/src/variants/annotations.rs','inputs_sha256':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in ['unvalidated_22_shared_direction_score.csv','stage6_tpmi_shared22_insomnia_inventory.csv','stage5_functional_atlas_with_priority_splice.csv']},'positive_controls':controls,'rows':len(rows),'hits':sum(r['status']=='selected_association_returned' for r in rows),'limits':'Selected-association endpoint only, response absence is not a null or p-value; server implementation uses significant_variants table and LIMIT 1 without ancestry/sequencing filtering, so assess reported source carefully. Manual prior mapping and allele build not independently reverified by this screen. Participant and relative overlap unverified; not score validation.'}
(p/'stage6_aou_public_variant_screen.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'rows':report['rows'],'hits':report['hits'],'positive_controls':len(controls)}))
