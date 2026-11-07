import pathlib,csv,subprocess
p=pathlib.Path('outputs/shared-genetics-atlas/prediction_rank_comparison');subprocess.run(['python3',str(p/'compare.py')],check=True,stdout=subprocess.DEVNULL)
r=list(csv.DictReader((p/'variant_ranks.tsv').open(),delimiter='\t'));assert len(r)==425
x={z['snp']:z for z in r};assert x['rs62264764']['dnase_absolute_rank']=='1' and x['rs62264764']['insomnia_pip_rank']=='6' and x['rs62264764']['adhd_pip_rank']=='2';assert x['rs2431108']['dnase_absolute_rank']=='69' and x['rs2431108']['avi_raw_rank']=='1';assert x['rs4916723']['insomnia_pip_rank']==''
y=list(csv.DictReader((p/'overlap_summary.tsv').open(),delimiter='\t'));assert len(y)==12 and all(z['dnase_pip_top5_overlap']=='' for z in y if z['trait_rank_gate']=='NOT_ESTIMATED')
print('PASS: full variant coverage, example ranks and withheld trait gates')
