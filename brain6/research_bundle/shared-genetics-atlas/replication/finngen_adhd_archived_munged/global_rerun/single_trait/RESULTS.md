# Single-trait heritability recheck

With the archived ADHD-derived file and the same reference and LDSC installation as the cross-trait rerun, standalone `--h2` prints **h²=0.1892 (SE 0.0238)** under the archived file’s N convention on 1,156,825 reference-merged SNPs. That exactly matches the historical aggregate's displayed h² and SE. The cross-trait run instead reported 0.1904 (SE 0.0250) on 1,108,865 allele-compatible markers. A numerical match supports a plausible explanation of the historical summary, but cannot establish which input bytes or command the historical run actually used.

The standalone log reports a two-step estimator with cutoff 30; the default cross-trait log does not. A subsequent [explicit two-step cross-trait sensitivity](../two_step/RESULTS.md) shows that the estimation setting and marker support both differ. Therefore, SNP-count differences alone do not prove the cause of the h² difference. The ADHD source identity, source processing, N coding, Finnish reference fit and participant overlap remain unverified; this is not liability-scale heritability or an independent replication.

The [protocol](PROTOCOL.md), [run script](run.sh), [LDSC output](adhd_h2.log), [console](console.log), and [receipt](receipt.json) preserve the run and input hashes.

A [sample-size audit](../n_scale/RESULTS.md) finds that N matches the historical case–control effective sample size; the printed h² changes sharply under a synthetic constant N rescaling. It is not a validated disease-scale heritability estimate.
