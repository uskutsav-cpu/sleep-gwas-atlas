# Five-track PLACO exploratory dependency classification

Analysis label: **EXPLORATORY**. This branch starts from the receipt-bound
five-track PLACO candidate set (25 pair-specific intervals, 2,686 candidate
variants) and the separate 20-group geographic reconciliation. It does not
promote a LAVA result or modify the frozen shared-locus evidence tiers.

The existing dependency chain in `extensions/brain6/README.md`,
`extensions/brain6/brain6/pipeline.py`, the follow-up and molecular templates,
and `brain6/config/shared_locus_evidence_tiers_v1.json` is: dense GWAS and
source metadata -> PLACO and locus/LD evidence -> signed effect/LD alignment ->
trait-specific fine-mapping -> trait-trait coloc -> QTL coloc -> functional
context -> integrated evidence. The frozen confirmatory route additionally
requires a complete QC-passing LAVA family before local-rg promotion and
final region tiers. The canonical family failed that gate.

| Step | Class | Scientific boundary and present prerequisite |
|---|:---:|---|
| Verify five-track PLACO receipts, candidate identities, coordinates, and 25-to-20 membership | A | Reuse unchanged pair-level family and geographic outputs; grouping is not an independent-signal count. |
| Reconcile candidates with independent genotype LD and report unmatched alleles/discordant clumps | A | Diagnostic LD may be computed without LAVA; it cannot rewrite the frozen UKB clumping rule. |
| Attach already available pair-level replication and source-sensitivity evidence | A | Label global replication separately from locus-specific replication; no global result validates an individual region. |
| Compare candidate variants/regions with dated published GWAS loci and genes | A | Match build, coordinates, trait, and cohort overlap; catalogue occurrence is context, not independent replication. |
| Align two GWAS effect signs to verified allele orientation and signed LD | B | Exploratory direction consistency only; require exact allele and source/hash checks, and retain ambiguous or missing cases. |
| Trait-specific SuSiE fine-mapping | B | Exploratory only if both dense GWAS, per-variant sample-size semantics, ancestry-matched signed LD, effect alignment, LD coverage, and a reviewed native runtime satisfy the method contract. A PLACO clump alone is insufficient. |
| Trait-trait coloc | B | Exploratory only after both traits have valid regional inputs; multi-signal coloc additionally depends on sound trait-specific fine-mapping. Record prior sensitivity and all failed regions. |
| Full-cis eQTL/sQTL indexing and GWAS-QTL coloc | B | Exploratory only with source-verified, build/allele-compatible, tissue-relevant QTL and reviewed sample-size/prior settings. Nearest-gene context does not substitute for QTL evidence. |
| Gene and regulatory annotation | B | Descriptive coordinate-based context may be attached to every region; gene prioritization requires fine-mapping/coloc or another independent functional link. Annotation cannot be used to assign the frozen region tiers. |
| Tissue, cell-type, and pathway analyses | B | Exploratory only when source panels, a defensible variant/gene weighting scheme, background universe, and multiplicity family are fixed; candidate overlap or nearest-gene lists alone do not support enrichment claims. |
| Cross-disorder candidate comparison | B | Geographic overlap and LD can be described; shared causal signals or recurrent mechanisms require stronger evidence. |
| Further source-verified multi-trait LAVA rescue pilots | B | Separate prospective sensitivity branch; keep sample-size semantics and source selection outcome-blind, with distinct configs and receipts. |
| Canonical local-rg, 12,475-slot bivariate LAVA family, and local sharing claims | C | Blocked by immutable v3 `FAILED_QC_NOT_PROMOTED` (3,720/17,465 NOT_RUN; frozen maximum 873). |
| Final shared-locus tiers, confirmatory cross-layer categories, and LAVA-confirmed region promotion | C | Frozen tier rule requires the complete eligible LAVA family; downstream biological context cannot change this gate. |
| Confirmatory mechanisms, causal genes/variants, therapeutic targets, or publication-ready validated locus claims | C | Require all relevant independent statistical and biological gates; none follows from PLACO association alone. |

Classes: **A** can be run on PLACO candidates without LAVA confirmation;
**B** is exploratory/diagnostic only and still requires its own method inputs;
**C** remains `BLOCKED_LAVA`. A and B outputs must carry their own provenance
and cannot be inserted into the frozen tier decision as confirming evidence.
