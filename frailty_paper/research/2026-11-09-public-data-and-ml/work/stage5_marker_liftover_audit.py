import csv,concurrent.futures,subprocess,json,hashlib,urllib.parse
m=list(csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv')))
a=[{'marker':x['variant'],'chr37':x['chr_GRCh37_assumed'],'pos37':x['position_GRCh37_assumed'],'chr38':x['chr_GRCh38'],'pos38':x['position_GRCh38']} for x in m];x=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');a.append({'marker':'rs323509','chr37':x['GRCh37_chr'],'pos37':x['GRCh37_pos'],'chr38':x['GRCh37_chr'],'pos38':x['GRCh38_pos']})
def get(x):
 u='https://grch37.rest.ensembl.org/map/human/GRCh37/'+x['chr37']+':'+x['pos37']+'..'+x['pos37']+'/GRCh38?content-type=application/json'
 try:
  b=subprocess.run(['curl','-fsS','--max-time','35',u],capture_output=True,check=True).stdout;j=json.loads(b);maps=j.get('mappings',[]);matching=[z for z in maps if str(z['mapped']['seq_region_name'])==x['chr38'] and str(z['mapped']['start'])==x['pos38'] and str(z['mapped']['end'])==x['pos38']];return {**x,'n_mappings':len(maps),'n_expected_coordinate_matches':len(matching),'mapping_strand':matching[0]['mapped']['strand'] if len(matching)==1 else '', 'status':'unique_expected_positive_strand' if len(maps)==len(matching)==1 and matching[0]['mapped']['strand']==1 else 'needs_review','url':u,'response_sha256':hashlib.sha256(b).hexdigest()}
 except Exception as e:return {**x,'n_mappings':'','n_expected_coordinate_matches':'','mapping_strand':'','status':'request_failed','url':u,'response_sha256':'','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(get,a))
fields=['marker','chr37','pos37','chr38','pos38','n_mappings','n_expected_coordinate_matches','mapping_strand','status','url','response_sha256','error']
with open('outputs/stage5_marker_grch37_grch38_liftover_audit.csv','w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(r)
print({s:sum(x['status']==s for x in r) for s in set(x['status'] for x in r)})
