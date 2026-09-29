#!/usr/bin/env python3
"""Audit frozen discovery-extension outputs without rerunning analyses."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def row_count(path: Path) -> int:
    with path.open(encoding='utf-8', errors='replace') as f:
        return max(0, sum(1 for _ in f) - 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='.')
    args = ap.parse_args()
    root = Path(args.repo).resolve()
    ext = root / 'discovery_extension'
    checkpoint = json.loads((ext / 'core_checkpoint.json').read_text())
    report = json.loads((ext / 'provenance/final_report.json').read_text())

    core_present: list[str] = []
    core_missing: list[str] = []
    core_mismatch: list[str] = []
    for rel, expected in checkpoint['artifact_hashes_sha256'].items():
        path = root / rel
        if rel == 'environment/tool_versions.tsv':
            pinned = subprocess.check_output(
                ['git', 'show', f"{checkpoint['git_commit']}:{rel}"],
                cwd=root, text=True,
            ).splitlines()
            current = path.read_text(encoding='utf-8').splitlines()
            if current[:len(pinned)] == pinned:
                core_present.append(f'{rel} (append-only prefix policy)')
            else:
                core_mismatch.append(f'{rel} (pinned rows changed or reordered)')
        elif not path.is_file():
            core_missing.append(rel)
        elif sha256(path) != expected:
            core_mismatch.append(rel)
        else:
            core_present.append(rel)

    ancestor = subprocess.run(
        ['git', 'merge-base', '--is-ancestor', checkpoint['git_commit'], 'HEAD'],
        cwd=root, check=False,
    ).returncode == 0

    output_checks: list[tuple[str, str, str, bool]] = []
    for name, record in report['outputs'].items():
        path = root / record['path']
        actual = sha256(path) if path.is_file() else 'MISSING'
        output_checks.append((name, record['path'], record['sha256'], actual == record['sha256']))

    extension_sources = [
        ('100-trait rerun h2', 'discovery_extension/results/ldsc/extension_trait_readiness.tsv',
         json.loads((ext / 'provenance/prioritization_post_replication.json').read_text())['extension_h2_sha256']),
        ('1200-pair rg family', 'discovery_extension/results/ldsc/extension_rg_matrix.tsv',
         json.loads((ext / 'provenance/ldsc_collation.json').read_text())['output_sha256']),
        ('217-pair replication family', 'discovery_extension/results/replication/replication_results.tsv',
         json.loads((ext / 'provenance/replication_results.json').read_text())['output_sha256']),
    ]
    for name, rel, expected in extension_sources:
        path = root / rel
        actual = sha256(path) if path.is_file() else 'MISSING'
        output_checks.append((name, rel, expected, actual == expected))

    n_h2 = row_count(root / extension_sources[0][1])
    n_rg = row_count(root / extension_sources[1][1])
    n_rep = row_count(root / extension_sources[2][1])
    n_top = row_count(root / report['outputs']['top_novel_discoveries']['path'])
    verified_outputs = sum(check[3] for check in output_checks)
    report_lines = [
        '# Frozen discovery-extension audit',
        '',
        'Read-only, checksum-based review. No extension analysis was rerun.',
        '',
        '## Core checkpoint',
        '',
        f"Checkpoint commit: `{checkpoint['git_commit']}`; ancestor of current `HEAD`: **{ancestor}**.",
        f"Pinned artifacts verified under their exact or documented append-only rules: **{len(core_present)}**.",
        f"Pinned artifacts missing: **{len(core_missing)}**; mismatches: **{len(core_mismatch)}**.",
        'The core checkpoint is **not fully re-verifiable from this checkout** because required files are absent. '
        'This qualifies the older extension status report, which recorded the checkpoint as passing at generation time.',
        '',
        'Missing pinned core artifacts:',
    ]
    report_lines += [f'- `{x}`' for x in core_missing] or ['- None']
    report_lines += ['', 'Pinned core mismatches:', *([f'- `{x}`' for x in core_mismatch] or ['- None'])]
    report_lines += ['', '## Frozen extension outputs', '',
        f"Reported outputs matching their provenance hashes: **{verified_outputs}/{len(output_checks)}**.",
        f'Result row counts: {n_h2} h2 traits; {n_rg} rg pairs; {n_rep} replication candidates; {n_top} top replicated discoveries.',
        f"Frozen report records {report['extension_fdr_hit_count']} extension-FDR-significant pairs among {report['executed_test_count']} planned tests.",
        f"Extension readiness: `{report['publication_readiness']}`.",
        '', 'Checksums:', '']
    report_lines += [f'- `{name}` — `{rel}`: **{"PASS" if passed else "FAIL"}**' for name, rel, _, passed in output_checks]
    report_lines += ['', '## Reuse boundary', '',
        'The existing extension reports qualified global discovery and independent replication. '
        'Local correlation, pleiotropy, fine-mapping, colocalization, mechanism and causal interpretation remain blocked upstream; '
        'the extension report encodes these as `NA_BLOCKED_UPSTREAM`, not as zero findings. '
        'This exploratory 100-trait family remains separate from the 45-trait core atlas and from the current frailty analysis.', '']
    text = '\n'.join(report_lines)
    out = ext / 'analysis/frozen_extension_audit.md'
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.read_text(encoding='utf-8') != text:
        raise SystemExit(f'Refusing to replace non-identical audit: {out}')
    out.write_text(text, encoding='utf-8')
    print(f'Wrote {out.relative_to(root)}; core missing={len(core_missing)}, core mismatches={len(core_mismatch)}, extension output checks={verified_outputs}/{len(output_checks)}')


if __name__ == '__main__':
    main()
