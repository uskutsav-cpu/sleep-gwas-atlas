# Manuscript claim downgrade and resource recheck — 2026-09-25 08:22 UTC

## Manuscript wording

The abstract now describes a protocol-led review workflow rather than an assembled evidence map because the 56,117-record PubMed queue remains unscreened. The Discussion no longer characterizes the FI results, secondary latent-factor family and read-only aging-context extract as broad overlap across multiple frailty definitions. It now identifies them as descriptive patterns of global genetic correlation and states that they do not establish robustness across frailty definitions or independent confirmation. Numeric source-table values were not changed.

The revised manuscript SHA-256 is `2ef5dbc597df0b7ee0f4f0762b46f3bc7df463de6705495b4341eafe44331ced`. The quantitative-claim audit still passes all 28 selected claims (output SHA-256 `140fa8d0ee694a2ea0f5e5187681fe39c1a23ebb406529e47616fe4ab9c8f274`). The bibliography audit passes 29 citation uses across 27 entries with no missing, duplicate or uncited keys (output SHA-256 `b923a659e3cd144579f52b32a429f7ae07c70fec10e1f2eee1b3432703cf5961`; `references.bib` SHA-256 `399016abed676fcb0732aeedb8bd8979678cbd53ebaa45a80cee3cbeb95c6093`). The reporting-locator audit passes 92 checklist rows and 93 line locators with zero errors (output SHA-256 `fba861dbc97a9e7b4a80c6651dd577051e1bfb999326df9ee90637094ef35f7e`). These checks do not validate citation truth or substantive guideline compliance.

## FI×sleep LAVA and host state

The external controller still records pause `pause_20260925T060720Z_b646`, desired concurrency four, `PAUSED_AFTER_WORKERS_EXITED`, 8,176/29,940 receipts, 21,764 pause holds, zero workers and no `runner.lock`. The independent receipt checkpoint at 06:23:41 UTC verified all 8,176 receipt identities/checksums and exact pending-slot holds with zero duplicate claims/events. Its frozen analysis-lock SHA-256 is `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256 is `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Those large receipt sets were not rescanned at 08:22.

A read-only 08:22 UTC host snapshot found four unrelated Brain6 R workers and an unrelated SMFVI Python analysis active, with no LAVA worker. The eight-core host showed 0.45% CPU idle, load average 8.91, 133 MiB physically unused RAM, 3,457 MiB compressed memory, 30% system-wide free memory, and 7,635.44/11,264 MiB swap used. Internal disk sampling ranged approximately 408–934 MiB/s; the external SSD had 1.6 TiB free. The LAVA pause and all receipts/claims remain unchanged. No LAVA workers were launched because CPU and internal I/O were saturated and physical headroom was very low.

The full test suite and integrated preflight were not rerun at this load. Overall project readiness remains NO-GO; review screening, independent paired replication, alternate frailty source eligibility, exact cohort overlap and LAVA family completeness/QC remain unresolved.
