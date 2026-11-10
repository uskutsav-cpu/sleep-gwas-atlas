# BBJ insomnia extraction status: complete

On 28 October 2026 exec session 50944 finished with shell pipeline exit code **0**. The checked-range autosomal gzip stream passed its gzip CRC and length check and scanned 13,435,997 rows. The complete ZIP was not retained and no publisher full-file checksum was verified. The [candidate review](stage6_bbj_insomnia_candidate_review.md), [receipt](stage6_bbj_insomnia_extract_receipt.json), and aligned table document the result. The earlier interrupted partial files remain labelled partial and should not be used for absence or completeness claims.

The preparation script was exercised on synthetic subsets to check sign orientation, the explicit integrity confirmation flag, missing pairs, the low-count flag and frequency alignment. Its first full run revealed an extra paper-control record; the script was corrected to exclude that recognized control from the score candidates and then completed. Testing does not independently validate source conventions, reference population fit, sample independence or prediction.

The [rare-variant warning](stage6_bbj_insomnia_partial_rare_variant_qc.md) arose in the partial fragment; the exact rs12524473 row was subsequently recovered in the complete stream and the warning is confirmed at the source-row level. Individual-level score validation and a sleep-excluded FI comparison remain outstanding.
