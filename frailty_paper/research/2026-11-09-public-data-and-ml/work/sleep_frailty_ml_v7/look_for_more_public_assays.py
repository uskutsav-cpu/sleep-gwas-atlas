"""Additional source leads only; never changes the current fixed experiment."""
import urllib.parse
from common import *
BASE='https://mpra.johnomeara.com/api/v1'
def main():
    sources=[]
    for query in ['Tewhey','Matoba','neuro-immune','circadian']:
        url=BASE+'/studies?'+urllib.parse.urlencode({'q':query,'limit':100});path=SOURCES/('OpenMPRA_'+query.replace('-','_')+'_metadata.json');receipt=fetch(url,path,10_000_000);value=json.loads(path.read_text());sources.append({'query':query,'receipt':receipt,'pagination':value.get('pagination'),'studies':value.get('data',[])});print(query,'studies',value.get('pagination'),flush=True)
        for study in value.get('data',[]):print('SOURCE LEAD',study.get('study_id',study.get('id')),study.get('title'),flush=True)
    write_json(OUT/'additional_public_assay_leads.json',{'access_date_UTC':time.strftime('%Y-%m-%d',time.gmtime()),'database':'OpenMPRA,JohnS.OMeara,https://mpra.johnomeara.com','sources':sources,'curation_limitation':'Database describes AI-assisted curation, not independent scientific review. Its processed tables/metadata require verification against original publication and source bytes. Units and quality-control rules are not standardized.','role':'Availability leads for additional molecular studies only. Not accepted training, final validation or human genetic-score data. Current V7 source roles and final-test rules remain fixed.','current_final_study_model_performance_seen':(OUT/'reserved_evaluation_started.json').exists()})
if __name__=='__main__':main()
