# Track B fresh literature audit

Audit date: 2026-09-01

This audit searched each frozen pair using the exact phenotypes, documented synonyms, and all method families requested in the Track B specification. Search-engine discovery was followed by inspection of PubMed, PMC, or publisher primary-source records. A negative search row means only **no exact report was identified by this protocol on this date**; it is not proof that no paper exists.

## Pair A — snoring / parental lifespan

The global relationship is not novel: Campos et al. directly reported a negative correlation between snoring and “parents age at death” (`rg=-0.1279`, `SE=0.0633`, `P=0.0432`). It uses the same or overlapping snoring GWAS and an unresolved/different parental-lifespan source, so it is prior global evidence rather than independent replication. Independent or partly independent parental-lifespan cohorts exist, but this audit did not identify public, dense, signed, closely matched summary statistics sufficient for the requested standardized LDSC replication. No exact local-to-cellular mechanistic chain was identified.

## Pair B — insomnia / ADHD

The global relationship and generic shared-genetic claim are already well established. Carpena et al. reported direct global correlation, shared gene-level signals and MR. Demontis et al. included insomnia in the ADHD GWAS correlation analysis. Xue et al. reported conjunction-FDR/ASSET loci and cortical GABAergic/glutamatergic cell enrichments across insomnia and psychiatric disorders. Zu et al. subsequently reported `rg=0.32`, 12 PLACO loci, four coloc signals with `PP4>0.90`, and brain/cell enrichment. These studies materially narrow novelty: a defensible Track B advance would require independent replication and stronger signal-aware, cell-specific QTL/chromatin convergence, not another generic pleiotropy claim.

## Positive control — insomnia / frailty

Song et al. already reported global LDSC, a rho-HESS locus at 3p21.31, MTAG/CPASSOC shared variants, broad brain-tissue enrichment, bulk GTEx SMR and bidirectional MR. The control is therefore appropriate for pipeline concordance but not novelty. This audit did not identify an exact signal-aware, independently cell-resolved mechanism.

## Evidence classification counts

- `DIRECT_PRIOR_GLOBAL`: 6
- `DIRECT_PRIOR_LOCAL`: 1
- `PRIOR_CELL_MECHANISM`: 2
- `PRIOR_SHARED_LOCUS`: 5
- `RELATED_ONLY`: 8
- `UNCERTAIN`: 1
- `UNDERREPORTED`: 3

## Guardrails

- Literature findings do not change any downstream threshold.
- Dataset overlap is carried forward as a limitation, not silently relabeled as replication.
- Shared-locus association is not colocalization; bulk eQTL/SMR is not a cell-specific mediator.
- Cell enrichment is not evidence that a cell type is causal.
- The 2026 literature means Pair B novelty must be judged against already published PLACO, coloc, and cell-enrichment results.
