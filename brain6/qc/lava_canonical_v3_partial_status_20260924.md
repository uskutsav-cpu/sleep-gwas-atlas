# Canonical LAVA v3 partial status snapshot

Snapshot time: 2026-09-24 03:01 UTC
Run ID: `d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b`
Analysis ID: `brain6-lava-canonical-v3`

This is a partial-run diagnostic, not a completed family audit or promotable LAVA result. At snapshot creation, the checksum-bound receipt verifier accepted all 786 observed receipts, representing 5,502 of the 17,465 intended trait-by-locus cells; no observed receipt was invalid. The production process remained active. There was no `latest_audit.json`, `family_decision.json`, or complete canonical family result.

## Status by trait among observed receipts

| Trait | Verified cells | TESTED | NOT_RUN | FAILED | NOT_RUN among observed |
|---|---:|---:|---:|---:|---:|
| ADHD | 786 | 738 | 48 | 0 | 6.11% |
| Bipolar disorder | 786 | 759 | 27 | 0 | 3.44% |
| Insomnia | 786 | 578 | 208 | 0 | 26.46% |
| Long sleep | 786 | 390 | 396 | 0 | 50.38% |
| Major depressive disorder | 786 | 607 | 179 | 0 | 22.77% |
| Parkinson's disease | 786 | 608 | 178 | 0 | 22.65% |
| Schizophrenia | 786 | 769 | 17 | 0 | 2.16% |
| **All traits** | **5,502** | **4,449** | **1,053** | **0** | **19.14%** |

Every observed `NOT_RUN` cell has reason `LOW_LOCAL_H2_UNDERPOWERED`; the observed cells contain no execution failures. This is a scientific local-h² eligibility result, not a runner failure. The partial snapshot already has 1,053 untested cells, whereas the frozen family permits at most `floor(0.05 × 17,465) = 873`. Because completed receipts are immutable and the remaining family can only add untested cells, the canonical family cannot pass its frozen untested-cell QC gate, even if every remaining cell is tested. Keep the threshold unchanged and do not promote downstream local-rg results. Continue the active run to complete auditable status accounting.

These per-trait rates describe only the first 786 receipt-verified loci. They are not full-family rates, pairwise local-rg results, or conclusions about biological effects.

## Integrity bindings

- Frozen family config SHA256: `f30c46f1eaaaf50d4697cf8af97f3acd44bef2ad8e067c016f3b2568b424b417`
- Frozen execution config SHA256: `f37a212ff44472004549095a3b17d28da0d130e2ab2a9ea507feb87312b65d27`
- Canonical Python coordinator SHA256: `79b66315996e855be0e25f11bc2cb54c59532a1d833c0228262094779f764d1b`
- Canonical R worker script SHA256 (receipt-bound): `fb97cecf41faa7e4361af94873623f7904f103cd8993cba5610894ed7310a492`
- LAVA reference provenance SHA256: `35ab371935b9c260dab21be248f4907b693565787ef9f62334e7d8030204270e`
- Worker count in the frozen execution configuration: **4**
- Receipt identity and per-cell integrity: independently verified for all 786 receipts in this snapshot.

The failed v2 baseline and roundoff run remain separate and untouched.
