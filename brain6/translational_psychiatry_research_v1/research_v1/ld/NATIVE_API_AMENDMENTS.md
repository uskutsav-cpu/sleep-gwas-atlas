# Native interface corrections; scientific gates unchanged

The initial diagnostic attempt used the current documentation's `$table`
field. The pinned susieR 0.14.2 implementation returns `$conditional_dist`.
The missing field triggered a hard failure, never QC passage. The full
initial failure log is preserved as
`native_RSS_consistency_initial_interface_failure.log`.

The pinned large-N implementation requires omitting N. Explicit n=Inf
would enter its finite-N (n-1)/(z²+n-2) adjustment and produce Inf/Inf.
The corrected omission implements the prespecified N-infinity sensitivity;
it does not change the primary finite-N analysis or any gate. The 37 initial warnings comprise 36 optimizer messages from the invalid
Inf call plus the ggplot2 deprecation warning. A deliberate excluded
reproduction retains all 36 optimizer messages in
`initial_invalid_N_infinity_reproduced_warnings.tsv`; that invalid call is
never used as a scientific QC/pass result. The corrected run
retains the expected large-N approximation warning and a ggplot2
`aes_string()` deprecation warning. Neither was interpreted as calibrated
scientific evidence.

A following output-writing attempt hit R's newline/else parsing rule after
printing both valid diagnostics. Its full serializer/control-flow failure
log is preserved as
`native_RSS_consistency_initial_serializer_failure.log`. Braces repaired
that code path; the final native diagnostic execution exits successfully.
A strict adapter now requires a finite 5-column, p-row `conditional_dist`
table. Three synthetic native regressions verify the corrected API,
large-N omission and hard failure for a missing pinned diagnostic field.
Ten source-free regressions also cover signed allele parity, palindromic
exclusion, BGZF CRC/truncation, exact HTTP ranges and independent posterior
arithmetic. None is counted as clinical or biological validation.

Large diagnostic RDS files initially retained ggplot environments and full
LD eigenvectors. A serialization-only cleanup retains every kriging table
entry, s estimate, variant key and eigenvalue extrema in compact local RDS;
full matrix/factor hashes and exact regeneration code preserve the reference.
All source-derived diagnostic tables stay local under source restrictions.
Original native fine-mapping fit RDS objects and coloc outputs are unchanged.
