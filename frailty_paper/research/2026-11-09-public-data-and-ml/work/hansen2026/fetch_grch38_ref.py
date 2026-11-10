import csv,requests,concurrent.futures,hashlib
c=list(csv.DictReader(open('outputs/expanded_22_candidate_risk_gene_table.csv')))
def get(r):
 u='https://api.genome.ucsc.edu/getData/sequence';pos=int(r['GRCh38_pos']);p={'genome':'hg38','chrom':'chr'+r['GRCh38_chr'],'start':pos-1,'end':pos};x=requests.get(u,params=p,timeout=30);x.raise_for_status();s=x.json()['dna'].upper();assert len(s)==1 and s in 'ACGT',repr(s);return {'variant':r['variant'],'chromosome':r['GRCh38_chr'],'position38':str(pos),'reference_base':s,'source_url':x.url}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(get,c))
with open('outputs/stage6_hansen_candidate_grch38_reference.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('reference bases',len(rows),hashlib.sha256(open('outputs/stage6_hansen_candidate_grch38_reference.csv','rb').read()).hexdigest())
