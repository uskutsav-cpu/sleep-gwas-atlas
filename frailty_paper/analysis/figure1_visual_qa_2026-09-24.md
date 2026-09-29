# Figure 1 visual and provenance QA

Review date: 2026-09-24.

## Rebuild

Rebuilt `figure1_evidence_framework` with the pinned Python 3.11.11 environment after detecting stale provenance for three committed inputs: the Fried FFS source-access audit, the physical-component source audit, and the targeted-resource audit. The Fried update confirms that the UK Biobank polygenic score is not full GWAS data or sleep-pair replication. The physical-component update resolves operational field definitions but leaves release-specific GWAS provenance and QC unresolved. The HFRS audit documents publisher hit tables but no complete genome-wide statistics or verified custom endpoint file.

The HFRS panel now says that publisher hit tables are available while full summary statistics and the exact custom endpoint file remain unverified. No analysis eligibility, estimates, or review counts changed.

## Evidence and counts

- All 10 recorded input hashes match current files.
- The script and source-data TSV hashes match the provenance record.
- All 3 output hashes (SVG, PDF, PNG) match the provenance record.
- Review snapshot: 60,989 PubMed records; 4,897 linked duplicates; 56,092 unique records queued; zero manual-database records; screening remains not started.
- Analysis snapshot: 12 FI pairs (9 inherited q≤0.05 under the all-396 family); 84 latent-factor pairs (48 q≤0.05 within the separate 84-pair family); 72 aging-context pairs (26 inherited q≤0.05 under the all-396 family).

## PDF visual inspection

Rendered the vector PDF at 144 dpi for inspection. `pdfinfo` reports a single page at 1500 × 1125 pt. The four panels, labels, counts, caveats, and updated HFRS status are visible; no clipping, overlap, or rasterization defects were observed. The figure remains a provisional PubMed-only framework, not a completed PRISMA flow or final result figure.
