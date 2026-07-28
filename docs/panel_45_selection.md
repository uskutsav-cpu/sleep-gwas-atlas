# Final 45-trait panel: selection contract

`config/panel_45_selection.tsv` is the shared selection contract for the
final Phase 0/1 panel. It has exactly 45 distinct candidate trait IDs:

- 12 sleep/circadian traits, 6 neuro/psychiatric traits, 6 immune traits,
  9 metabolic/cardiovascular traits, 6 aging traits, and 6 cancer traits.
- Codex owns the 24 sleep/circadian, aging, and cancer rows.
- Claude owns the non-overlapping 21 neuro/psychiatric, immune, metabolic,
  and cardiovascular rows on its own branch.

Selection means a trait is in scope; it does **not** mean it is curated,
downloadable, or eligible for LDSC. The 86-row historic registry remains the
source of truth until the two curation branches are merged. At that point the
active `config/traits.tsv` will be reduced to these 45 rows and the excluded
rows will be archived rather than silently discarded.

The panel deliberately prefers common, interpretable disorders and traits
(for example Alzheimer’s disease, major depression, asthma, type 2 diabetes,
coronary artery disease, stroke, breast cancer, and sleep duration) over
duplicated phenotype definitions or poorly powered placeholders. Every source
still has to meet the existing EUR, GRCh37/hg19, public-access, header, and
metadata gates before `status=CURATED` is allowed.
