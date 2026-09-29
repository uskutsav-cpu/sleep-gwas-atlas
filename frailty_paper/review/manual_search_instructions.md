# Manual database search and export instructions

The PubMed searches do not substitute for licensed database searches. Do not automate scraping or bypass institutional access. Run each search in the database interface or an authorized institutional platform and save the exact query, platform/provider, search date, filters, hit count and exported-file checksum in `manual_exports/search_log.tsv`.

The operational query templates below are recorded as `AMEND-001` in
`protocol_amendments.tsv`. Protocol v1's eligibility and outcome definitions
remain unchanged; no licensed manual search or screening had started when this
search syntax was added.

Use the primary concept structure below, adapting only syntax and controlled vocabulary to the database. Preserve each database's exact executed query. Search from database inception through **2026-09-22** for this frozen version; any update search must be logged as an amendment.

## Core concepts

- Sleep/circadian: sleep, insomnia, sleep duration, short sleep, long sleep, sleep quality, fragmentation, efficiency, daytime sleepiness, napping, sleep apnea/sleep-disordered breathing, snoring, chronotype, circadian.
- Frailty: frail*, frailty index, frailty phenotype, Fried frailty, Clinical Frailty Scale, prefrail/pre-frail, cognitive frailty.
- Secondary evidence map: sarcopenia, grip strength, gait/walking speed, exhaustion, physical activity, weight loss, disability, falls, cognition/dementia, Alzheimer, Parkinson, longevity/lifespan, healthspan, mortality.

## Frozen free-text query blocks

Run two searches in each database: a **primary sleep × frailty** search and a
separate **secondary sleep × aging evidence-map** search. The parenthesized
blocks below are the title/abstract/keyword text equivalents. Do not add an
adult, human, language, study-design, or publication-type filter: those filters
can remove eligible records and the review includes observational and genetic
evidence. Apply the date range through 2026-09-22 in the database interface and
record the interface's inclusive/exclusive date behavior and all selected
indexes. If a platform cannot search a named field or interpret a phrase,
adapt only its syntax, retain the same concepts, and save the exact executed
query and any changes in `manual_exports/search_log.tsv`.

**Sleep/circadian block**

```text
sleep OR insomnia OR "sleep duration" OR "short sleep" OR "long sleep" OR "sleep quality" OR "sleep fragmentation" OR "sleep efficiency" OR "daytime sleepiness" OR hypersomnolence OR napping OR "sleep apnea" OR "sleep apnoea" OR "sleep disordered breathing" OR snoring OR chronotype OR circadian
```

**Frailty block**

```text
frail* OR frailty OR "frailty index" OR "frailty phenotype" OR "frailty syndrome" OR "Fried phenotype" OR "Fried frailty" OR prefrail OR pre-frail OR "clinical frailty scale" OR "cognitive frailty"
```

**Secondary aging block**

```text
sarcopenia OR "grip strength" OR "gait speed" OR "walking speed" OR exhaustion OR "physical activity" OR "weight loss" OR disability OR falls OR cognition OR dementia OR Alzheimer* OR Parkinson* OR longevity OR lifespan OR "life span" OR healthspan OR mortality
```

Construct the primary set as `(sleep/circadian block) AND (frailty block)`.
Construct the secondary set as `(sleep/circadian block) AND (secondary aging
block)`. Keep these searches, result counts, records, and later evidence
classification separate. The secondary search is an evidence map, not evidence
that these adjacent outcomes are frailty.

### Platform syntax templates

Use these in the database's advanced-search interface, replacing each label
with the corresponding parenthesized block above. Do not paste the explanatory
labels literally.

| Platform | Primary template | Secondary template |
|---|---|---|
| Embase.com | `((SLEEP BLOCK):ti,ab,kw) AND ((FRAILTY BLOCK):ti,ab,kw)` | `((SLEEP BLOCK):ti,ab,kw) AND ((AGING BLOCK):ti,ab,kw)` |
| Scopus | `TITLE-ABS-KEY((SLEEP BLOCK)) AND TITLE-ABS-KEY((FRAILTY BLOCK))` | `TITLE-ABS-KEY((SLEEP BLOCK)) AND TITLE-ABS-KEY((AGING BLOCK))` |
| Web of Science Core Collection | `TS=((SLEEP BLOCK)) AND TS=((FRAILTY BLOCK))` | `TS=((SLEEP BLOCK)) AND TS=((AGING BLOCK))` |
| PsycINFO, EBSCOhost | Search the free-text blocks in TI and AB as separate OR sets, then map relevant APA Thesaurus descriptors and OR those headings into their matching concept sets. Combine the sleep set with the frailty set using AND. Save the translated search history. | Same field construction, using the aging block. |
| PsycINFO, Ovid | Search the free-text blocks in `.ti,ab.` and combine with AND; separately explode/map relevant APA Thesaurus headings and OR those lines into the matching concept set. Preserve the full Ovid line history. | Same field construction, using the aging block. |

For Embase, also inspect Emtree and add relevant exploded preferred headings
for each concept as OR terms in the appropriate concept set. Record the actual
preferred headings and complete translated query; Emtree headings can change
and should not be guessed from the free-text template. On Scopus and Web of
Science, retain the full query exactly as submitted, including parentheses,
and note the searched collections/indexes. On PsycINFO, record the vendor and
interface because field codes, thesaurus headings, and export behavior differ.

## Embase

In Embase.com, combine Emtree terms and title/abstract/keyword terms for the sleep/circadian concept with frailty terms. Run a separate secondary evidence-map query for the listed adjacent aging outcomes. Save the full search history and export all records with citation, abstract, DOI, PMID, Emtree terms, accession number and database source in RIS or CSV. Do not export only selected records.

## Scopus

Use `TITLE-ABS-KEY(...)` for the sleep/circadian terms AND the frailty terms; run the adjacent-outcome query separately. Export all results with title, authors, year, source, abstract, DOI, PMID where indexed, keywords, affiliations and citation identifier as CSV or RIS.

## Web of Science Core Collection

Use `TS=(...)` topic searches for sleep/circadian AND frailty; run the adjacent-outcome query separately. Record which collections were searched. Export full records and cited references where institutionally permitted, including abstract, DOI, PMID, accession number and database identifier.

## PsycINFO (if available and appropriate)

Search APA Thesaurus terms and title/abstract fields for sleep/circadian AND frailty. Include adult human records; do not apply a language filter. Export all results with abstracts, DOI, PMID, descriptors and source identifiers. Record platform/provider (e.g., EBSCOhost or Ovid) because syntax and export limits differ.

## Import and audit

Place untouched exports in `frailty_paper/review/manual_exports/` using filenames such as `embase_YYYY-MM-DD.ris`. Do not overwrite source exports. Keep database source IDs and export checksums. Run the deduplication builder only after imports have been converted to the agreed schema; every discarded duplicate must map to a retained record with the match key and source identifiers. Log truncated exports, platform limits, inaccessible full texts and records awaiting classification.

Supported input formats are RIS and CSV. Begin filenames with `embase`, `scopus`, `webofscience` (or `wos`), or `psycinfo`, then run `make -C frailty_paper import-review-exports`. This creates a SHA-256 import manifest and a normalized-record table without changing raw exports. Run the review-record builder afterward, before any screening decisions are entered. The builder preserves every source record and records cross-source duplicate links; uncertain identifier/title conflicts are retained separately rather than silently merged.

## Syntax references

Consulted 2026-09-23: [Embase search fields and release notes](https://service.elsevier.com/app/answers/detail/a_id/29534/supporthub/evolve/), [Embase Boolean/wildcard syntax](https://service.elsevier.com/app/answers/detail/a_id/17875/supporthub/embase/kw/proximity/), [Scopus advanced-search fields](https://service.elsevier.com/app/answers/detail/a_id/11365/supporthub/scopus/~/how-can-i-best-use-the-advanced-search%3F/), [Web of Science Core Collection field tags](https://webofscience.help.clarivate.com/en-us/Content/wos-core-collection/woscc-search-field-tags.htm), [Web of Science search rules](https://webofscience.help.clarivate.com/en-us/Content/search-rules.htm), and [Ovid PsycINFO database guide](https://ospguides.ovid.com/OSPguides/psbkdb.htm). These documents support field/operator syntax; the searcher must still preserve the exact query as translated and executed in the institution's licensed interface.
