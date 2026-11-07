# Removing either marker of the out-of-range LD pair

We repeated the selected FI–insomnia locus 464 calculation at the original 99% LD-retention threshold, removing one marker at a time from the pair whose saved reference entry exceeded one. All other GWAS inputs, reference, overlap file, processing settings and seed were held fixed. Both replays retained 2,740 markers and 216 components.

| Support | FI genetic-variance diagonal | FI status |
|---|---:|---|
| Original 2,741 markers | −4.191×10⁻⁹ | Failed |
| Remove rs62262671 | −5.903×10⁻¹⁰ | Failed |
| Remove rs55754265 | −1.696×10⁻⁹ | Failed |

Removing either one individually did **not** resolve the negative FI variance in this calculation. Both estimates moved closer to zero, but the starting and resulting estimates are close to the failure boundary. The comparison does not establish the role of the out-of-range entry: removing a marker also changes the LD decomposition, other numerical issues may remain, and this is one post-selected locus. It does not validate the reference, justify excluding either marker, or show that the pair has no influence. No bivariate genetic correlation was estimated, and the original family quality gates remain failed.

The [pre-run protocol](PROTOCOL.md), [R replay](run.R), [results](results.json), and [arithmetic verification](verify.py) retain the exact tests. The tests compare against the [previously retained failure](../../RESULTS.md), and their reference and provenance limitations apply.
