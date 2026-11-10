# Why the Liu atlas subtype and union peaks can disagree

For rs9536347, the three posted subtype intervals overlap **nearby union intervals**, although the variant itself lies outside those union intervals. This makes the two observations geometrically compatible; it does not establish the true regulatory effect of the allele.

| Subtype file | Subtype interval | Intersecting union interval | Overlap in bases, assuming half-open intervals |
|---|---|---|---:|
| Inh_ENOX2_SPHKAP | chr13:53069441–53069941 | chr13:53069133–53069633 | 192 |
| Oli_DPP10 | chr13:53069787–53070287 | chr13:53070074–53070574 | 213 |
| Oli_OPALIN | chr13:53069738–53070238 | chr13:53070074–53070574 | 164 |

The candidate position used was chr13:53069876 (GRCh38, conditional on the source-file assembly). The result is robust to a one-base start-coordinate convention at this interior point, although reported intersection lengths can shift by a base. The downloaded union file contains **1,010,153 intervals**, matching the count reported in [Liu et al.](https://doi.org/10.1016/j.cell.2025.06.031). This is a consistency check, not a provider checksum or proof that the posted file is exactly the one used for every figure.

The paper says it used ArchR `addReproduciblePeakSet` to derive subtype and union sets. [ArchR's documentation](https://www.archrproject.com/archive/1.0.1/bookdown/the-iterative-overlap-peak-merging-procedure.html) describes merging fixed-width peaks by iterative overlap removal, retaining a selected peak and removing directly overlapping competitors. This is a **plausible explanation** for a subtype interval containing a variant while an overlapping selected union interval does not contain it. We have not recovered the actual selection log or demonstrated which input peak won; the particular cause of this example is therefore not proven. The atlas files' individual assembly and exact export coordinate conventions still require confirmation.

The [interval audit](stage6_liu_union_subtype_audit.csv) is reproducible with `python3 outputs/stage6_liu_union_subtype_audit.py`; the [receipt](stage6_liu_union_subtype_audit.json) records hashes for the locally retained union and previously extracted subtype matches. No test of cell-type enrichment, allele-specific activity, insomnia, frailty, participant-level prediction or genetic causality was performed. The broader [candidate overlap review](stage6_liu_brain_ccre_candidate_review.md) documents the 22-marker screen and limitations.
