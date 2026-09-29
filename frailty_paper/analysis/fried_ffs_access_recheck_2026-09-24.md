# Fried Frailty Score Figshare access recheck

**Checked:** 2026-09-24 01:39 UTC (2026-09-23 20:39 CDT)
**Source:** the article-linked Figshare share item for `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv` (Figshare file 38980967).

## Ordinary route result

Sent a header-only request to the same official file URL already cited by the share item:

```text
curl -sSIL --max-time 30 'https://figshare.com/ndownloader/files/38980967?private_link=6683396c68807fe4e729'
```

The response was HTTP/2 202 with `x-amzn-waf-action: challenge`, `content-length: 0`, and `content-type: text/html; charset=UTF-8`. This ordinary route did not return file bytes. No alternate host or endpoint was queried; no challenge, CAPTCHA, login, or access control was bypassed.

## Local artifact state

At 2026-09-23 20:38 CDT, the Downloads item `/Users/swethasunilkumar/Downloads/Unconfirmed 591651.crdownload` was still 6,296,188 bytes with its last modification at 2026-09-23 15:06:15 CDT. The registered project target on the external SSD remained zero bytes (last modification 2026-09-22 17:44:43 CDT). Neither file is a verified complete dataset. The partial was left untouched.

## Conclusion

The linked file remains unavailable for harmonization or analysis through the verified ordinary route. Preserve the existing access blocker. Resume only if the provider's normal share route yields a complete file or an authorized provider supplies a copy; then verify the exact byte count, checksum, schema, genome build, and effect semantics before registration as acquired data.

Related source and prior ordinary-browser attempts are documented in `fried_ffs_source_access_audit_2026-09-23.md`, `fried_ffs_figshare_recheck_2026-09-23.md`, and `fried_ffs_browser_recheck_2026-09-23.md`.
