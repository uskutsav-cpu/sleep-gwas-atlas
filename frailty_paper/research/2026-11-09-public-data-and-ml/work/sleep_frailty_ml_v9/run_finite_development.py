"""Finite augmented-data comparison, then advance to locked external testing."""
import subprocess,time
from v9_common import *


def main():
    plan={'recipes':['frozen','tuned'],'held_studies':['none','Salomon','Kosicki2025','Siraj','IGVF240K','Koesterich2023','McAfee2023','Retallick2026'],'epochs':8,'seed':42,'primary_selection':'One common recipe/epoch maximizing mean raw magnitude-error skill over the training-only cell/culture-line median, across whole-study omissions with supported cell groups. Report all sources and unsupported lines; no per-held-study epoch selection presented as independent evidence.','RNA_loss':'Direct per-allele activities and contrasts, precision-weighted auxiliary and paired losses; main effect-size MAE stays unmodified by SE.','independent_test':'Preserve the same eligible Matoba cases and source/audit files.','stopping_rule':'After this declared coverage-amended round, fix the best research candidate and complete the matched AVI and fresh-study assessment without additional architecture or source changes before that test. If no neural model beats the simple baseline, retain and state that adverse finding rather than claiming successful training.'}
    path=OUT/'finite_development_plan.json'
    if path.exists():assert json.loads(path.read_text())==plan
    else:write_json(path,plan)
    start=time.monotonic();completed=[]
    for recipe in plan['recipes']:
        for study in plan['held_studies']:
            name=f'{recipe}_{study}_42';record=OUT/(name+'_record.json')
            if record.exists():
                r=json.loads(record.read_text());assert r['completed_epochs']==8 and sha(r['checkpoint'])==r['checkpoint_sha256'];print('VERIFIED COMPLETED V9 FIT',name,flush=True)
            else:
                print('STARTING DECLARED V9 FIT',name,flush=True)
                with (OUT/(name+'.log')).open('w') as f:subprocess.run([sys.executable,str(WORK/'train_development.py'),'--recipe',recipe,'--held-study',study,'--epochs','8'],env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True)
            completed.append(record.name);print('COMPLETED DECLARED V9 FIT',name,'seconds',round(time.monotonic()-start,1),flush=True)
    write_json(OUT/'finite_development_execution_receipt.json',{'fits':completed,'plan_sha256':sha(path),'fresh_model_performance_read':False,'elapsed_seconds':time.monotonic()-start})
    print('FINITE V9 COMPARISONS COMPLETE; LOCKED EXTERNAL ASSESSMENT NEXT',flush=True)


if __name__=='__main__':main()
