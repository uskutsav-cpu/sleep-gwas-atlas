# Chromosome 5 regional original-release orientation audit

The original GRCh37 source was read by verified HTTP ranges for the prespecified chromosome 5 window and compared with the physical-function values in both existing matched colocalization files after pinned-chain positive-strand liftover and exact non-palindromic allele-pair matching.

| Matched regional file | Rows checked | Concordant | Opposite effect allele | Discrepancies or missing unique original keys |
|---|---:|---:|---:|---:|
| Insomnia | 2,451 | 2,451 | 2,451 | 0 |
| ADHD | 2,360 | 2,360 | 2,360 | 0 |

For every matched row, flipping the original effect allele produced the harmonized beta and complementary effect-allele frequency, and SE and p agreed within the scripted numeric tolerances. These checks strengthen confidence in the **physical-function source orientation for the retained rows**. They do not validate the frozen insomnia/ADHD effects, independently establish participant overlap, rule out other signals, or establish colocalization. The source was queried only over the selected region; its prior full-file MD5 receipt does not establish whole-file identity during this later query or verify the full separate harmonized release.

The [audit script](audit.py), [protocol](PROTOCOL.md), and [receipt](receipt.json) record source URL, byte range, chain hash, matched-file hashes, row counts and any issues. The [regional results](../RESULTS.md) remain prior-sensitive and inconclusive.
