# Interrupted source stream

On 28 October 2026 the first questionnaire-insomnia comparison was interrupted: LDSC exited with macOS `Bus error: 10` (exit 138). Immediately afterward `/Volumes/Extreme SSD` and `/Volumes/Brain6Work` were absent, and `diskutil list external` did not show the external physical SSD. The runtime, reference and an input were on those volumes. The source stream stopped before completion, no checksum receipt was written, and the zero-byte log contains no genetic-correlation result. This observation is consistent with a drive disconnection but does not independently establish the cause of the bus error.

Do not use the partial stream or report a result. Once the SSD is visible again, verify input and reference paths and free capacity, then rerun the prespecified full stream and both analyses. The publisher checksum must match on each completed stream before results can be interpreted. Do not delete files to make space without explicit authorization.

Follow-up: the SSD mounted and the workspace image was reattached. The questionnaire and registry inputs matched their previously recorded SHA-256 values. The separate Brain6NHANES volume still returned an I/O error, so the retry uses a temporary directory on Brain6Work and distinct output names. The SSD had about 105 MiB free and Brain6Work about 32 MiB; no large archive is being saved.
