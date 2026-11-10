# Five-component Fried phenotype audit

The reference is Fried et al., J Gerontol A (2001), https://pubmed.ncbi.nlm.nih.gov/11253156/; Appendix details are available at https://pameladwilson.com/wp-content/uploads/4_3_2020-Fried-Frailty-In-Older-Adults.pdf. The five-column dataset retains each component and distinguishes them from the original validated phenotype. These components, scoring choices and output names are in `derived_variable_dictionary.csv` and `supplement_49_item_audit.csv`.

| Component | CHS original | NHANES implementation | Classification and limitation |
|---|---|---|---|
| Shrinking | Baseline self-report >10 lb unintentional loss in prior year; follow-up measured ≥5% weight loss and intent | Difference of WHD050 and WHD020 >10 lb, WHQ060=2 for unintentional loss; no >10 lb loss coded 0 if both weights valid | Close self-report proxy, not the same single question or longitudinal measured follow-up. WHQ060 is asked after substantial loss; missing intention stays missing. |
| Weakness | Dominant-hand grip, sex/BMI stratified cutoffs | Mean of available dominant-hand NHANES dynamometer trials; CHS sex/BMI cutoffs | Close measurement, different devices/protocol; inferred first-hand mapping; both hands and known dominance required. No claim of validation. |
| Exhaustion | Two CES-D effort/get-going statements over prior week, at least 3 days for either | DPQ040 tiredness/little energy at least more than half the days of preceding two weeks | Proxy only: different questions, timeframe and frequency. |
| Slowness | Timed 15-foot walk and sex/height cutoffs | PFQ061B quarter-mile difficulty; PFQ054 equipment-needed skip treated as limitation | Proxy for mobility limitation; **no timed gait speed** in these selected variables. |
| Activity | Minnesota Leisure Time Activity questionnaire and sex-specific kcal/week cutoffs | GPAQ work, travel and recreation multiplied by MET and body weight; original 383/270 kcal thresholds | Proxy: different question coverage and gross-energy computation, so applying CHS thresholds is not validated. |

The adapted five-item sum requires all five observed; the four-item physical mean omits exhaustion and requires the other four. The ≥3 category is explicitly exploratory. The previous four-component measure uses combined grip with fixed exploratory thresholds, quarter-mile difficulty, no reported activity in five domains, and ≥5% loss regardless of intent. The current package does not silently replace that earlier score.
