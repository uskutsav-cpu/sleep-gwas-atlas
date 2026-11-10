# LBC1936 published frailty index: sleep-item and coding audit

A public [LBC1936 frailty-index study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7456784/) supplies 30 top-level items and published coding in its Table A1. We transcribed the table into [a machine-readable inventory](lbc1936_published_fi_items.csv); [the receipt](lbc1936_published_fi_items_receipt.json) records the public XML source and hash. Run `python3 work/lbc1936_fi_item_audit.py` using the saved XML to reproduce it. This is a study-level description, not a participant-level reconstruction.

None of the 30 **top-level item labels** explicitly mentions sleep, insomnia, fatigue, tiredness, exhaustion or lethargy. This makes the index a useful candidate for a sensitivity analysis of the sleep–frailty question. It does not prove that every multi-item scale is free of sleep- or energy-related content or that the index is independent of sleep. The study's depression, anxiety and disability scales and other deficits need item-level review. This FI is also not necessarily identical to the LBC1936 FI used in every other paper.

The FI includes six-metre walk time and grip strength, while the published Fried-type phenotype also uses walking and weakness. Chair-standing and disability items are also physical-function measures. These shared or related inputs mean any FI-versus-physical-frailty comparison is not a comparison of independent constructs. Test the physical phenotype, the published FI and, only if justified and allowed, a separately defined index excluding relevant overlapping items; do not treat a modified index as automatically validated.

## Coding must be clarified before reconstruction

The public Table A1's wording has several ambiguities that should not be silently resolved in code:

- The row labeled “Able to stand up from a chair” lists “Yes (1) or No (0)”; as written, the direction is counterintuitive for a deficit. The wording or underlying variable coding may differ from the actual algorithm. Check the original instrument and analysis code.
- Six-metre walk time is described as greater than 10 seconds or physically unable versus less than 10 seconds; the table does not specify exactly 10 seconds.
- MMSE ranges as written omit exactly 10 and include 20 in two ranges. Verify the intended cut points and missing-data rules.

These are issues in the **published table as read**, not proof the original authors computed the FI incorrectly. Percentile tie handling, missingness, wave-specific thresholds and the permitted linked variables also require documentation. The cohort's published sleep questionnaire and genetic release cannot be used to reconstruct a linked score without the required access and linkage approvals. No data application, participant analysis or score validation has been performed.

## Wave 3 underlying-scale labels

The official [LBC1936 data dictionary](https://lothian-birth-cohorts.ed.ac.uk/sites/default/files/2023-11/Lothian%20Birth%20Cohort%201936%20Data%20Dictionary%20Waves%201%20to%204.xls) lists all 14 HADS anxiety and depression item labels at Wave 3. A search of those labels found no explicit sleep, insomnia, fatigue, tiredness, exhaustion or lethargy keyword; however, the depression item `HADSD4_w3` asks about feeling “slowed down,” a potentially related symptom rather than a direct fatigue measurement. This strengthens the narrow description of the HADS labels but does **not** prove that an FI using those scales is biologically independent of sleep or fatigue, nor that the published Table A1 actually used every field exactly as labeled.

The same dictionary lists the Wave 3 Townsend activities-of-daily-living total (`adl_w3`) and chair-stand fields. A [cohort methods description](https://link.springer.com/article/10.1186/1471-2318-7-28) describes Townsend as a nine-item activities-of-daily-living scale, but the limited search did not identify its individual component labels in the Wave 3 dictionary; do not conclude its items are free of sleep/fatigue wording without the actual instrument. The dictionary also does not resolve whether the Table A1 chair-standing coding is reversed or how scores are calculated.

[Source-field inventory](lbc1936_fi_underlying_scale_fields.csv), [receipt](lbc1936_fi_underlying_scale_fields.json); reproduce with `work/lbc1936_underlying_scale_audit.py` and the saved public XLS. No participant data were inspected. At authorized analysis time, confirm the instrument versions, components, scoring and data linkage before claiming a sleep- or fatigue-excluded index.
