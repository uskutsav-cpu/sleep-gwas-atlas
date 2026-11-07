# Author-code audit: identifying the sleep-item-excluded FI files

The authors' [pinned sensitivity-analysis code](https://github.com/jonathanklmak/Metabolomics_frailty/blob/6f370810052d08d69ba465d205b34966fc59be2a/Program/2.2_MR_frailty_sensitivity.R) identifies the modified mental-wellbeing-excluded FI as `nomental`. It constructs local input paths for three corresponding GWAS files:

| File implied by the author code | Role in the code |
|---|---|
| `fi_nomental.fastGWA` | Modified FI full-sample version |
| `fi_nomental_no_metabolomics.fastGWA` | Modified FI version without the metabolomics subset |
| `fi_nomental_random50.fastGWA` | Modified FI random-half version |

These names are constructed by a loop and refer to the authors' local `Data/GWAS_summary_statistics_data/Frailty/` directory. They make a targeted request possible, but **do not establish that the files are publicly downloadable or that they are present in a particular release**. The inspected public repository tree at commit `6f370810052d08d69ba465d205b34966fc59be2a` contains eight tracked files and no `Data/` blobs. The existing [version-1 Figshare inventory](https://doi.org/10.6084/m9.figshare.24591573.v1) includes `fi.fastGWA` and full-FI subset files, but none named `nomental`; this does not rule out another repository, differently named file, or access by request.

The sensitivity script describes removing eight items: general health, fatigue, sleep, depressed feelings, nervous personality, anxiety or panic, loneliness and misery. Thus `nomental` is broader than a sleep-only or sleep-and-fatigue-only deletion. Its MR code reads all three files and uses the full, non-metabolomics or random-half results for different biomarker exposures. The expected fields in the script include SNP, BETA, SE, AF1, A1, A2, P, POS and N; allele semantics, genome build, sample definitions, phenotype construction and model must still be verified with the authors or actual file metadata. The [main analysis code](https://github.com/jonathanklmak/Metabolomics_frailty/blob/6f370810052d08d69ba465d205b34966fc59be2a/Program/2.1_MR_frailty_selected_metabolites.R) describes the full FI and physical-phenotype GWAS and sample/model settings; those comments do not by themselves prove the same people or model were used for the modified FI.

The [data request draft](../DATA_REQUEST_DRAFT.md) now names the relevant files and requests a matched full-FI comparator. **No request was sent.** This audit does not supply a modified-FI GWAS, genetic correlation, or evidence of an independent sleep–frailty association.

The [repository tree](repository_tree.json), [verification code](verify.py), and [receipt](verification.json) record the pinned commit, source code hashes, file-name checks and scope. The script verifies local copies of the two public R scripts under the task's `work/` directory; the public code can be re-fetched at the pinned URLs above. A [current-release check](latest_release_check.json) on 28 October 2026 found only version 1 and no filename containing `nomental`. The repository and one release inventory were checked; the results are not an exhaustive proof of absence.
