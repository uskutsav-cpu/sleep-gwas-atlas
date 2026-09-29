# Fried Frailty Score Figshare browser recheck

Recheck date: 2026-09-23. This records a second ordinary-route attempt to obtain the Fried Frailty Score (FFS) summary-statistics file described in `fried_ffs_source_access_audit_2026-09-23.md`.

## Attempt

Opened the article-linked Figshare share page (`https://figshare.com/s/6683396c68807fe4e729`) in the Codex in-app browser. The visible page described `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`, displayed a size of 609.63 MB, showed the item as privately shared, and exposed a “Download (609.63 MB)” link. The page's cookie dialog initially covered the page; selected “Reject all” and then clicked the visible download link. The browser page remained on the share page without a visible completion or error state. A local search of the expected Downloads and project data locations found no file modified during the attempt.

## Result and handling

No completed download was verified. Do not treat the local zero-byte placeholder, if present, as data. The prior audit's direct request had returned an HTTP 202 WAF challenge; this browser interaction did not resolve that access problem. No CAPTCHA was solved, no warning bypassed, and no credentials or forms were used. Keep the FFS source unavailable pending a normal, verified download route. If the data provider's anti-bot/access restriction remains, request access or assistance through the provider's ordinary contact route rather than bypassing it.

## Alternate-host check

Checked the Knowledge Portal Network listing for dataset `Ye2023_Frailty_EU`, which identifies the Fried Frailty Score GWAS (386,565 European-ancestry participants) and exposes a “download” link. Following that link resolves to the same Figshare share item and the same 609.63 MB file; it is not a separate mirror or alternative access route. No additional download was attempted through a guessed endpoint. [Knowledge Portal listing](https://www.kp4cd.org/index.php/node/1238).
