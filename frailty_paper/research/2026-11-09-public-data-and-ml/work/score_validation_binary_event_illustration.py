import json, math
from pathlib import Path
# Local-null approximation: one standardized continuous score, equal
# score variance in cases and controls, independent observations.
scenarios=[('published_MESA_imaging_subgroup_illustration',32,1711),('hypothetical_100_cases',100,3906),('hypothetical_500_cases',500,3506)]
rows=[]
for name,cases,controls in scenarios:
 se=math.sqrt(1/cases+1/controls)
 rows.append({'scenario':name,'cases':cases,'controls':controls,'null_log_or_se':se,'null_or_95ci_lower':math.exp(-1.959964*se),'null_or_95ci_upper':math.exp(1.959964*se),'approx_or_for_80pct_power_two_sided':math.exp((1.959964+0.841621)*se)})
out={'method':'Near-null unadjusted logistic approximation; 1 SD score, equal variance, independent observations; power via normal shift; no finite-sample or nuisance adjustment.','scenarios':rows,'limits':'MESA counts are from a published imaging subset, not known genotype/sleep overlap. Other case counts are hypothetical and not cohort forecasts. Not a sample size calculation for incremental AUC, calibration or joint outcome.'}
Path('outputs/stage6_binary_event_precision_illustration.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(rows,indent=2))
