#!/usr/bin/env python3
"""Record rounded TPMI browser transcriptions; this script performs no GWAS analysis."""
import csv
from pathlib import Path

root = Path(__file__).resolve().parent
# Browser transcriptions: GRCh38 chromosome, position, REF, ALT, then effect, SE,
# p for phecodes 327.4 and 327.41. Effects are NOT verified as score-allele effects.
observed = {
    'rs28528997': (7,2034487,'A','G',-.0082,.013,.52,-.00081,.017,.96),
    'rs2815764': (1,72288440,'A','G',-.00035,.041,.99,.038,.054,.49),
    'rs2155292': (11,113040061,'A','G',.015,.012,.21,.0085,.016,.59),
    'rs817509': (3,118032014,'A','G',-.0047,.035,.89,.0065,.046,.89),
    'rs2838788': (21,45122459,'T','C',-.011,.0099,.27,.0026,.013,.84),
    'rs12671330': (7,114693988,'T','C',.022,.013,.11,.022,.018,.22),
    'rs6452792': (5,88497027,'G','A',-.0038,.0098,.70,-.0058,.013,.66),
    'rs7691487': (4,89910850,'T','C',.0049,.011,.65,.011,.015,.44),
    'rs2604551': (4,15089577,'T','G',.052,.030,.083,.060,.040,.13),
    'rs6058750': (20,32590184,'T','C',.0059,.024,.80,.041,.032,.20),
    'rs6938026': (6,43217995,'A','G',.011,.0098,.26,.012,.013,.35),
    'rs7731706': (5,136235299,'T','G',.0018,.0093,.85,.0083,.012,.50),
    'rs11636381': (15,73842542,'T','C',-.00064,.011,.95,.0048,.014,.73),
}
with (root / 'stage6_tpmi_insomnia_candidate_inventory.csv').open(newline='') as f:
    prior = {r['variant']:r for r in csv.DictReader(f)}
with (root / 'stage5_score_queue.csv').open(newline='') as f:
    queue = list(csv.DictReader(f))
with (root / 'unvalidated_22_shared_direction_score.csv').open(newline='') as f:
    score = list(csv.DictReader(f))
fields = ['variant','counted_allele','other_allele','status','chr38','pos38','ref','alt',
          'insomnia_effect','insomnia_se','insomnia_p','persistent_effect','persistent_se','persistent_p',
          'queue_identity_match','effect_allele_verified','source_url']
rows=[]
for item in score:
    name=item['variant']
    row=dict.fromkeys(fields,''); row.update(variant=name,counted_allele=item['counted_allele'],
        other_allele=item['other_allele'],status='rsid_search_no_exact_match',effect_allele_verified='False')
    if name in prior and prior[name]['status']=='browser_lookup_rounded':
        p=prior[name]; ref,alt=p['browser_pair'].split('/')
        data=(int(float(p['chr38'])),int(float(p['pos38'])),ref,alt,*[p[k] for k in fields[8:14]])
    else:
        data=observed.get(name)
    if data:
        chrom,pos,ref,alt,*values=data
        matches=[q for q in queue if name in q['aliases'].split('|') and
                 q['chr']==f'chr{chrom}' and q['position_GRCh38']==str(pos) and
                 q['reference_allele']==ref and q['alternate_allele']==alt]
        if len(matches)!=1: raise ValueError(f'{name}: queue variant/allele mismatch')
        if {ref,alt}!={item['counted_allele'],item['other_allele']}:
            raise ValueError(f'{name}: score allele pair mismatch')
        row.update(status='browser_lookup_rounded',chr38=chrom,pos38=pos,ref=ref,alt=alt,
                   queue_identity_match='True',source_url=f'https://pheweb.ibms.sinica.edu.tw/variant/{chrom}-{pos}-{ref}-{alt}')
        row.update(dict(zip(fields[8:14],values)))
    rows.append(row)
assert len(rows)==22 and len({r['variant'] for r in rows})==22
with (root / 'stage6_tpmi_shared22_insomnia_inventory.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} markers; {sum(r['status']=='browser_lookup_rounded' for r in rows)} rounded browser matches")
