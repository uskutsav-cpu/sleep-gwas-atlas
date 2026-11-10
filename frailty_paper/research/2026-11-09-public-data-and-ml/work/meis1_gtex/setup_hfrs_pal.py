from pathlib import Path
p=Path('work/meis1_gtex');o=Path('outputs')
src=(p/'prepare_hfrs_input.py').read_text().replace("or {ref,alt} in ({'A','T'},{'C','G'})", "or ({ref,alt} in ({'A','T'},{'C','G'}) and r['eqtl_rsid'] not in ('rs12621789','rs6546230'))").replace('meis1_hfrs_input.csv','meis1_hfrspal_input.csv')
(p/'prepare_hfrspal_input.py').write_text(src)
src=(p/'hfrs_frequency_extract.py').read_text().replace('meis1_hfrs_input.csv','meis1_hfrspal_input.csv').replace('meis1_hfrs_frequencies.csv','meis1_hfrspal_frequencies.csv')
(p/'hfrspal_frequency_extract.py').write_text(src)
for pop in ['nfe','fin']:
 for task in ['index','ld']:
  s=(o/f'stage6_meis1_hfrs_{pop}_{task}.py').read_text().replace('meis1_hfrs_', 'meis1_hfrspal_')
  (o/f'stage6_meis1_hfrspal_{pop}_{task}.py').write_text(s)
s=(p/'run_hfrs_susie.R').read_text().replace('meis1_hfrs_', 'meis1_hfrspal_')
(o/'stage6_meis1_hfrspal_susie.R').write_text(s)
