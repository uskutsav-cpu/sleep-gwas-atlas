# Pre-outcome SUPERGNOVA reference correction

The first `insomnia__adhd` process launched with the 22-chromosome combined
PLINK panel. At the local covariance stage, code review found that upstream
`ld.getBlockLefts(coords, 1)` requires sorted genetic positions. The combined
panel's interpolated cM resets on each new chromosome, so that input violates
the implementation's LD-score window assumption. This is an engineering input
error independent of any candidate association outcome.

The process was terminated with signal 15 before its `to_csv` call. Exit status
was 143; only its log was present, at
`/Volumes/Extreme SSD/brain6-work/brain6-alternative-local-validation-v1/runs_v1/insomnia__adhd.aborted_combined_reference.log`
(SHA-256 `80c7af05869a2fab68196e17f59a392232c97eebc47777f056ddf044cb3ff600`).
No local output or success receipt existed. The log records 883,054 shared
SNPs, global heritability and entry into the local covariance stage but no
local statistic. The log is preserved as an aborted attempt under a new name.

The corrected reference uses the upstream-supported `@` chromosome pattern,
splitting the same prepared genotype records and genetic map by chromosome.
The GWAS sources, total N values, 1,693 non-MHC LD blocks, five-pair family,
four workers, statistical method, and family-wide threshold do not change.
The amended protocol and all split reference hashes are committed before any
real local covariance output is read or generated.
