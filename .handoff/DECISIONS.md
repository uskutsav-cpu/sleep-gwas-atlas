# DECISIONS (scientific and engineering)

**2026-08-01 — h² is computed on the OBSERVED scale.** Every `pop_prev` in the
registry is `UNKNOWN` because none had a citation (rule 1: an uncited
prevalence is worse than a blank). Z = h²/SE is invariant to the liability
conversion, so the QC gate is unaffected. Confirmed empirically: shortsleep Z
went 23.7 → 23.6 and longsleep 13.6 → 13.8 across the change. Liability
conversion is deferred until the mentor supplies sourced prevalences.

**2026-08-01 — chronotype uses the UKB-only stratum.** GCST007576 (EUR
449,734) chosen over GCST007565, which carries the access-restricted 23andMe
component. Consistent with the project's own known-traps list.

**2026-08-01 — raw inputs are evicted after munging.** Munged sumstats are
~8 MB; raw + harmonized are ~700 MB per trait. The registry keeps url + sha256
+ bytes, so every input is exactly reproducible. Validated in practice:
re-downloaded traits returned byte-identical hashes.

**2026-08-01 — a documented project prior was found to be WRONG.** `CLAUDE.md`
states shortsleep and "especially longsleep" are expected to fail the power
gate. Measured, both PASS (Z = 23.6 and 13.8). longsleep fails the separate
MiXeR threshold (N_eff×h² = 3,567 < 12,000), which is a different question.
The prior should be corrected, not repeated.

**2026-08-01 — secondary sleep and aging traits are kept in a separate panel.**
frailty, healthspan and parental lifespan are aging traits, not sleep. They are
not substitutes for missing core anchors and must not inflate the primary
sleep × disease family.

**OPEN — healthspan / parental-lifespan sign inconsistency.** healthspan
correlates positively with BMI and depression while parental lifespan
correlates negatively with the same two. If both were coded "longer is
healthier" the signs would agree. At least one is coded as risk rather than
duration. **Must be resolved from source documentation before any rg sign in
the secondary panel is interpreted.**
