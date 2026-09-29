#!/usr/bin/env python3
"""Plot frozen sleep–frailty global-rg estimates with their distinct families."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


SLEEP_ORDER = ["insomnia", "sleepdur", "shortsleep", "longsleep", "chronotype", "sleepiness", "napping", "snoring", "sleep_apnea", "sleep_efficiency", "accel_sleep_duration", "sleep_timing"]
SLEEP_LABELS = ["Insomnia", "Sleep duration", "Short sleep", "Long sleep", "Chronotype", "Sleepiness", "Napping", "Snoring", "Sleep apnea", "Actigraphy efficiency", "Actigraphy duration", "Actigraphy timing"]
FACTORS = ["frailty", "frailty_general", "frailty_factor_1", "frailty_factor_2", "frailty_factor_3", "frailty_factor_4", "frailty_factor_5", "frailty_factor_6"]
FACTOR_LABELS = ["FI\nPrimary", "General factor", "F1: Limited\nSocial support", "F2: Unhealthy\nLifestyle", "F3\nMultimorbidity", "F4\nMetabolic problems", "F5\nPoorer cognition", "F6\nDisability"]
INPUTS = [
    "frailty_paper/analysis/frozen_atlas_frailty_global_rg.tsv",
    "frailty_paper/analysis/frozen_fi_correction_sensitivity.tsv",
    "frailty_paper/results/frailty_v1/rg_sleep_latent_factors.tsv",
    "frailty_paper/analysis/replication_resource_audit.tsv",
    "frailty_paper/config/secondary_frailty_factors_v1.tsv",
    "frailty_paper/config/secondary_sleep_factor_panel_v1.tsv",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_immutable(path: Path, data: bytes, replace: bool) -> None:
    if path.exists() and path.read_bytes() == data:
        return
    if path.exists() and not replace:
        raise SystemExit(f"Refusing to overwrite non-identical output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def blend(value: float) -> tuple[str, str]:
    """Return accessible diverging-cell fill and foreground colors."""
    v = max(-0.7, min(0.7, value)) / 0.7
    neutral = (241, 244, 247)
    target = (184, 101, 57) if v < 0 else (20, 92, 114)
    t = abs(v)
    rgb = tuple(round(n + (c - n) * t) for n, c in zip(neutral, target))
    fill = "#" + "".join(f"{c:02X}" for c in rgb)
    return fill, ("#FFFFFF" if abs(v) > 0.72 else "#182230")


def svg_text(x: float, y: float, text: str, *, size: int = 16, weight: str = "400", fill: str = "#182230", anchor: str = "start") -> str:
    return f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}px" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{esc(text)}</text>'


def build_svg(matrix: dict[tuple[str, str], dict[str, str]], counts: dict[str, int]) -> bytes:
    w, h = 1900, 1290
    a = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<defs><linearGradient id="rg" x1="0%" x2="100%"><stop offset="0%" stop-color="#B86539"/><stop offset="50%" stop-color="#F1F4F7"/><stop offset="100%" stop-color="#145C72"/></linearGradient></defs>',
        '<rect width="100%" height="100%" fill="#F5F7FA"/>',
        svg_text(44, 54, "Genome-wide genetic correlation across sleep and frailty", size=31, weight="700"),
        svg_text(44, 83, "Primary Frailty Index plus seven latent dimensions; estimates are grouped by their frozen correction family", size=17, fill="#475467"),
    ]
    # Matrix panel
    a += ['<rect x="35" y="120" width="1370" height="930" rx="16" fill="#FFFFFF" stroke="#D9E0E8" stroke-width="1.5"/>',
          svg_text(64, 160, "A", size=24, weight="700", fill="#145C72"),
          svg_text(102, 160, "Sleep × frailty rg estimates", size=21, weight="700"),
          svg_text(64, 188, "Cell value = rg; * indicates q≤0.05 within that column’s own family", size=14, fill="#667085"),
          svg_text(64, 210, "FI q-values: frozen all-396 atlas family   |   Latent-factor q-values: fixed 84-pair sensitivity family", size=14, fill="#667085")]
    x0, y0, cw, ch = 365, 270, 124, 51
    for j, label in enumerate(FACTOR_LABELS):
        x = x0 + j * cw + cw / 2
        lines = label.split("\n")
        for k, line in enumerate(lines):
            a.append(svg_text(x, 244 + 19 * k, line, size=13, weight="700", fill="#344054", anchor="middle"))
    for i, trait in enumerate(SLEEP_ORDER):
        y = y0 + i * ch
        a.append(svg_text(335, y + 30, SLEEP_LABELS[i], size=15, fill="#344054", anchor="end"))
        for j, factor in enumerate(FACTORS):
            row = matrix[trait, factor]
            rg = float(row["rg"])
            sig = float(row["q"]) <= 0.05
            fill, text_fill = blend(rg)
            x = x0 + j * cw
            a.append(f'<rect x="{x}" y="{y}" width="{cw-3}" height="{ch-3}" rx="5" fill="{fill}" stroke="#FFFFFF" stroke-width="1"/>')
            a.append(svg_text(x + (cw - 3) / 2, y + 30, f"{rg:.2f}{'*' if sig else ''}", size=14, weight="700" if sig else "400", fill=text_fill, anchor="middle"))
    legend_y = 915
    a += [f'<rect x="{x0}" y="{legend_y}" width="330" height="16" rx="8" fill="url(#rg)"/>',
          svg_text(x0, legend_y + 38, "−0.7", size=13, fill="#475467"),
          svg_text(x0 + 165, legend_y + 38, "0", size=13, fill="#475467", anchor="middle"),
          svg_text(x0 + 330, legend_y + 38, "+0.7", size=13, fill="#475467", anchor="end"),
          svg_text(x0 + 365, legend_y + 13, "Genetic correlation (rg)", size=14, fill="#475467"),
          svg_text(64, 1004, "No direct tests compare estimates between FI and latent factors; differences in color/significance are descriptive.", size=13, fill="#667085")]

    # Integrity and replication panel
    a += ['<rect x="1430" y="120" width="435" height="930" rx="16" fill="#FFFFFF" stroke="#D9E0E8" stroke-width="1.5"/>',
          svg_text(1458, 160, "B", size=24, weight="700", fill="#145C72"),
          svg_text(1496, 160, "Robustness and replication", size=20, weight="700")]
    cards = [
        ("PRIMARY FI", f"{counts['fi_bh']}/12", f"pass inherited all-396 BH q≤0.05; {counts['fi_bonf']}/12 also pass same-family Bonferroni sensitivity", "#EAF3F7", "#145C72"),
        ("LATENT FACTORS", f"{counts['latent_bh']}/84", "pass BH q≤0.05 in the fixed secondary family; sensitivity-only, exact overlap unknown", "#EAF4F1", "#23745F"),
        ("INDEPENDENT PAIRWISE REPLICATION", "None established", "Audited candidates do not provide an eligible independent sleep–frailty rg pair", "#FFF5E8", "#9A5A12"),
        ("COHORT OVERLAP", "Unresolved", "UK Biobank overlap is expected in 11/12 primary pairs and in latent sources; exact intersections unknown", "#F0F2F5", "#475467"),
    ]
    y = 205
    for title, value, detail, fill, color in cards:
        a.append(f'<rect x="1458" y="{y}" width="379" height="163" rx="12" fill="{fill}" stroke="#D0D7E0" stroke-width="1.2"/>')
        a.append(svg_text(1480, y + 31, title, size=13, weight="700", fill=color))
        a.append(svg_text(1480, y + 69, value, size=24, weight="700", fill=color))
        lines = detail.split("\n")
        words = detail.split()
        wrapped: list[str] = []
        line = ""
        for word in words:
            if line and len(line) + len(word) + 1 > 39:
                wrapped.append(line)
                line = word
            else:
                line = f"{line} {word}".strip()
        if line:
            wrapped.append(line)
        for k, line in enumerate(wrapped):
            a.append(svg_text(1480, y + 101 + 18 * k, line, size=13, fill="#475467"))
        y += 177
    a += [svg_text(1458, 960, "LDSC rg is genome-wide covariance,", size=13, fill="#475467"),
          svg_text(1458, 981, "not causation or a causal mechanism.", size=13, fill="#475467")]
    a += ['<rect x="35" y="1080" width="1830" height="145" rx="14" fill="#EEF1F4" stroke="#CCD4DD"/>',
          svg_text(62, 1118, "INTERPRETATION", size=15, weight="700", fill="#344054"),
          svg_text(220, 1118, "The FI values reuse frozen atlas results; the seven-factor family is secondary and sensitivity-only.", size=15, fill="#475467"),
          svg_text(220, 1147, "The q-value denominators differ by column. No plotted pattern establishes independent replication, local sharing, or which frailty dimension drives overlap.", size=14, fill="#475467"),
          svg_text(220, 1180, "See the source-data TSV for estimates, SEs, P values, family labels, overlap status and interpretation limits.", size=14, fill="#475467"),
          "</svg>"]
    return ("\n".join(a) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--replace-outputs", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    outdir = (args.outdir or repo / "frailty_paper/analysis").resolve()
    inputs = {rel: repo / rel for rel in INPUTS}
    missing = [rel for rel, p in inputs.items() if not p.is_file()]
    if missing:
        raise SystemExit(f"Missing locked input(s): {missing}")

    fi_rows = read_tsv(inputs[INPUTS[0]])
    bonf_rows = read_tsv(inputs[INPUTS[1]])
    latent_rows = read_tsv(inputs[INPUTS[2]])
    replication_rows = read_tsv(inputs[INPUTS[3]])
    factor_panel = read_tsv(inputs[INPUTS[4]])
    sleep_panel = read_tsv(inputs[INPUTS[5]])
    if len(fi_rows) != 12 or len(bonf_rows) != 12 or len(latent_rows) != 84:
        raise SystemExit("The frozen global-rg family sizes changed; update and re-audit before plotting")
    if len({r["sleep_trait"] for r in fi_rows}) != 12:
        raise SystemExit("FI rows do not contain 12 unique sleep traits")
    if any(r["family_denominator"] != "84" or r["interpretation_status"] != "SENSITIVITY_ONLY" for r in latent_rows):
        raise SystemExit("Latent family lock/status mismatch")
    if [r["trait_id"] for r in factor_panel] != FACTORS[1:]:
        raise SystemExit("Latent factor configuration/order mismatch")
    if [r["trait_id"] for r in sleep_panel if r["domain"] == "sleep"] != SLEEP_ORDER:
        raise SystemExit("Latent sleep panel configuration/order mismatch")
    if any(r["bonferroni_family_size"] != "396" for r in bonf_rows):
        raise SystemExit("Bonferroni sensitivity denominator mismatch")
    if any(r["pairwise_replication_eligibility"].startswith("ELIGIBLE") for r in replication_rows):
        raise SystemExit("Replication eligibility changed; revise the figure before use")

    matrix: dict[tuple[str, str], dict[str, str]] = {}
    data_rows: list[dict[str, str]] = []
    for row in fi_rows:
        matrix[row["sleep_trait"], "frailty"] = {
            "rg": row["global_rg"], "se": row["global_rg_se"], "p": row["global_rg_p"],
            "q": row["global_rg_fdr_all_396"], "family": "frozen atlas all-396 BH",
            "denominator": "396", "status": "PRIMARY_FROZEN_ATLAS_REUSE",
            "overlap": "UK Biobank expected for 11/12 pairs; exact participant intersections unknown",
        }
    for row in latent_rows:
        matrix[row["sleep_trait"], row["disease_trait"]] = {
            "rg": row["rg"], "se": row["se"], "p": row["p_value"], "q": row["q_value"],
            "family": "secondary sleep x latent factors BH", "denominator": "84",
            "status": row["interpretation_status"], "overlap": row["cohort_overlap_status"],
        }
    expected = {(sleep, factor) for sleep in SLEEP_ORDER for factor in FACTORS}
    if set(matrix) != expected:
        raise SystemExit(f"Matrix mismatch: missing={sorted(expected-set(matrix))}; unexpected={sorted(set(matrix)-expected)}")
    for sleep in SLEEP_ORDER:
        for factor in FACTORS:
            row = matrix[sleep, factor]
            data_rows.append({"row_type": "estimate", "sleep_trait": sleep, "frailty_definition": factor,
                              "rg": row["rg"], "se": row["se"], "p_value": row["p"], "q_value": row["q"],
                              "q_family": row["family"], "family_denominator": row["denominator"],
                              "status": row["status"], "cohort_overlap": row["overlap"]})
    fi_bh = sum(float(r["global_rg_fdr_all_396"]) <= 0.05 for r in fi_rows)
    fi_bonf = sum(r["bonferroni_all_396_significant_at_0.05"] == "TRUE" for r in bonf_rows)
    latent_bh = sum(float(r["q_value"]) <= 0.05 for r in latent_rows)
    counts = {"fi_bh": fi_bh, "fi_bonf": fi_bonf, "latent_bh": latent_bh}
    if counts != {"fi_bh": 9, "fi_bonf": 9, "latent_bh": 48}:
        raise SystemExit(f"Correction results changed: {counts}")
    summary_items = [
        ("FI BH significant", str(fi_bh), "frozen all-396 BH", "396"),
        ("FI same-family Bonferroni significant", str(fi_bonf), "P×396≤0.05 sensitivity", "396"),
        ("Latent BH significant", str(latent_bh), "secondary family BH; sensitivity-only", "84"),
        ("Independent pairwise replication", "None established", "resource audit; no eligible pairwise sleep–frailty rg", "NA"),
    ]
    for label, value, family, denominator in summary_items:
        data_rows.append({"row_type": "summary", "sleep_trait": "NA", "frailty_definition": label,
                          "rg": "NA", "se": "NA", "p_value": "NA", "q_value": value,
                          "q_family": family, "family_denominator": denominator,
                          "status": "DESCRIPTIVE_STATUS", "cohort_overlap": "exact participant intersections unknown"})

    columns = ["row_type", "sleep_trait", "frailty_definition", "rg", "se", "p_value", "q_value", "q_family", "family_denominator", "status", "cohort_overlap"]
    tsv_data = "\t".join(columns) + "\n" + "".join("\t".join(r[c].replace("\t", " ") for c in columns).rstrip("\t") + "\n" for r in data_rows)
    out_tsv = outdir / "figure2_global_rg_matrix_source_data.tsv"
    write_immutable(out_tsv, tsv_data.encode(), args.replace_outputs)

    svg = build_svg(matrix, counts)
    out_svg = outdir / "figure2_global_rg_matrix.svg"
    out_pdf = outdir / "figure2_global_rg_matrix.pdf"
    out_png = outdir / "figure2_global_rg_matrix.png"
    write_immutable(out_svg, svg, args.replace_outputs)
    pdf_cmd = ["rsvg-convert", "-f", "pdf", "-o", str(out_pdf), str(out_svg)]
    png_cmd = ["rsvg-convert", "-w", "4000", "-f", "png", "-b", "#F5F7FA", "-o", str(out_png), str(out_svg)]
    if args.replace_outputs or not out_pdf.exists():
        subprocess.run(pdf_cmd, check=True)
    if args.replace_outputs or not out_png.exists():
        subprocess.run(png_cmd, check=True)
    script = Path(__file__).resolve()
    prov = {
        "schema_version": "frailty_paper_figure2_global_rg.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable, script.relative_to(repo).as_posix(), *sys.argv[1:]],
        "working_directory": str(repo),
        "script": script.relative_to(repo).as_posix(), "script_sha256": sha256(script),
        "python": platform.python_version(),
        "rsvg_convert_executable": subprocess.check_output(["which", "rsvg-convert"], text=True).strip(),
        "rsvg_convert_version": subprocess.check_output(["rsvg-convert", "--version"], text=True).strip(),
        "inputs": {rel: sha256(path) for rel, path in inputs.items()},
        "source_data_tsv": out_tsv.relative_to(repo).as_posix(),
        "source_data_tsv_sha256": sha256(out_tsv),
        "families": {"FI": {"pairs": 12, "q_family": "frozen all-396 BH", "bh_q_le_0_05": fi_bh, "same_family_bonferroni_p_times_396_le_0_05": fi_bonf}, "latent": {"pairs": 84, "q_family": "fixed 84-pair BH", "bh_q_le_0_05": latent_bh, "interpretation": "SENSITIVITY_ONLY"}},
        "independent_pairwise_replication": "NONE_ESTABLISHED_IN_AUDITED_CANDIDATES",
        "outputs": {path.relative_to(repo).as_posix(): sha256(path) for path in (out_svg, out_pdf, out_png)},
        "render_commands": [pdf_cmd, png_cmd],
        "rendering": {"pdf": "vector SVG rendering", "png_width_px": 4000},
        "limitations": ["FI and latent estimates use distinct correction families and are not formally compared.", "Exact participant intersections remain unknown.", "Latent-factor estimates are sensitivity-only, not independent replication.", "Global genetic correlation is not causation and does not identify loci or mechanisms."],
    }
    out_prov = outdir / "figure2_global_rg_matrix.provenance.json"
    write_immutable(out_prov, (json.dumps(prov, indent=2, sort_keys=True)+"\n").encode(), args.replace_outputs)
    print(json.dumps({"pairs": len(matrix), "counts": counts, "outputs": prov["outputs"], "provenance": str(out_prov)}, indent=2))


if __name__ == "__main__":
    main()
