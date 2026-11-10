from pathlib import Path
import pandas as pd
O=Path('outputs'); s=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv')
data={
'rs8067974':(17,45155199,'A','G',.013,.0078,.11),
'rs12536335':(7,114403104,'A','G',.00029,.0079,.97),
'rs10865954':(3,49174556,'T','C',.037,.021,.084),
'rs7734045':(5,136202653,'A','C',.015,.0079,.062),
'rs1032524':(1,201821698,'T','C',.016,.0096,.09),
'rs1909650':(13,53217760,'G','A',-.022,.011,.026),
'rs1942262':(18,55206086,'G','A',.012,.0084,.14),
'rs11152363':(18,55389957,'G','A',-.0016,.0079,.83),
'rs4886140':(13,59259385,'A','G',.022,.011,.052),
'rs35881094':(2,58695786,'T','G',.0012,.008,.89),
'rs17550821':(2,58751388,'G','A',-.040,.029,.16),
'rs11790060':(9,93440650,'T','C',-.020,.0087,.011)}
assert len(data)==12 and set(data).issubset(set(s.variant))
rows=[]
for x in s.itertuples(index=False):
 row={'variant':x.variant,'prototype_chr37':x.chr,'prototype_pos37':x.position_GRCh37_inferred,'counted':x.counted_allele,'other':x.other_allele,'status':'not_verified','chr38':None,'pos38':None,'browser_pair':None,'browser_effect_rounded':None,'browser_se_rounded':None,'browser_p_rounded':None,'effect_allele_verified':False}
 if x.variant in data:
  chrom,pos,ref,alt,beta,se,p=data[x.variant]; assert set([ref,alt])==set([x.counted_allele,x.other_allele]);row.update(status='browser_lookup',chr38=chrom,pos38=pos,browser_pair=ref+'/'+alt,browser_effect_rounded=beta,browser_se_rounded=se,browser_p_rounded=p)
 rows.append(row)
d=pd.DataFrame(rows); d.to_csv(O/'stage6_tpmi_sleep_candidate_inventory.csv',index=False);print('Candidates:',len(d),'browser lookups:',(d.status=='browser_lookup').sum(),'unverified:',','.join(d.loc[d.status=='not_verified','variant']))
