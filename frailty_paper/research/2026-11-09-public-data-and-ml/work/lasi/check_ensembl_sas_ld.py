#!/usr/bin/env python3
"""Retrieve bounded 1000 Genomes SAS pair LD for two candidate pairs."""
import hashlib, json, subprocess

pairs=[('rs113851554','rs78206187'),('rs35881094','rs17550821')]
data=[]
for a,b in pairs:
    url=f'https://grch37.rest.ensembl.org/ld/human/pairwise/{a}/{b}?content-type=application/json;population_name=1000GENOMES:phase_3:SAS'
    body=subprocess.run(['curl','-fsSL','--max-time','30',url],check=True,stdout=subprocess.PIPE).stdout
    data.append({'pair':[a,b],'url':url,'sha256':hashlib.sha256(body).hexdigest(),'response':json.loads(body)})
with open('outputs/stage6_ensembl_sas_chr2_ld.json','w') as f: json.dump(data,f,indent=2)
print(json.dumps(data,indent=2))
