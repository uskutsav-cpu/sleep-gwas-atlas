"""Manifest-driven orchestration for the locked atlas-v1.0 panel."""
import csv


configfile: "config/workflow.yaml"

PANEL = config["panel"]
PANEL_LOCK = config["panel_lock"]
SOURCES = config["sources"]
SCHEMAS = config["schemas"]
PYTHON = config["python"]
LDSC_PYTHON = config["ldsc_python"]
LDSC_DIR = config["ldsc_dir"]
SELECTED = config.get("phase0_traits", [])
H2_SCALE = config.get("h2_scale", "liability")

with open(PANEL, newline="", encoding="utf-8") as handle:
    PANEL_ROWS = list(csv.DictReader(handle, delimiter="\t"))
RAW_BY_TRAIT = {row["trait_id"]: f"data/raw/{row['raw_file']}" for row in PANEL_ROWS}
PANEL_IDS = set(RAW_BY_TRAIT)
unknown = set(SELECTED).difference(PANEL_IDS)
if unknown:
    raise ValueError(f"phase0_traits contains IDs outside the locked panel: {sorted(unknown)}")


rule all:
    input:
        "results/tables/analysis_panel_provenance.tsv",
        "results/tables/source_readiness.tsv",


rule panel_contract:
    input:
        manifest=PANEL,
        lock=PANEL_LOCK,
    output:
        "results/tables/analysis_panel_provenance.tsv",
    shell:
        "{PYTHON} scripts/00_validate_panel.py --manifest {input.manifest} "
        "--lock {input.lock} --write-provenance {output}"


rule source_readiness:
    input:
        provenance="results/tables/analysis_panel_provenance.tsv",
        manifest=PANEL,
        lock=PANEL_LOCK,
        sources=SOURCES,
        schemas=SCHEMAS,
    output:
        "results/tables/source_readiness.tsv",
    shell:
        "{PYTHON} scripts/10_phase0_audit.py --config {input.manifest} "
        "--lock {input.lock} --sources {input.sources} --out {output} --strict"
        " --schemas {input.schemas}"


rule smoke_test:
    input:
        manifest=PANEL,
        lock=PANEL_LOCK,
    output:
        "results/_smoketest/SMOKE_TEST_OK",
    shell:
        "PYTHON_BIN={PYTHON} bash scripts/run_smoke_test.sh"


rule phase0:
    input:
        provenance="results/tables/analysis_panel_provenance.tsv",
        raw=[RAW_BY_TRAIT[trait] for trait in SELECTED],
    output:
        touch("results/logs/workflow_phase0.complete"),
    params:
        traits=" ".join(SELECTED),
    shell:
        "test -n '{params.traits}' || "
        "(echo 'ERROR: set phase0_traits in config/workflow.yaml' >&2; exit 1); "
        "PYTHON_BIN={PYTHON} LDSC_PYTHON={LDSC_PYTHON} LDSC_DIR={LDSC_DIR} "
        "bash scripts/02_munge.sh {params.traits}"


rule h2:
    input:
        phase0="results/logs/workflow_phase0.complete",
        munged=expand("data/munged/{trait}.sumstats.gz", trait=SELECTED),
    output:
        "results/tables/h2_summary.tsv",
    params:
        traits=" ".join(SELECTED),
        scale_flag="" if H2_SCALE == "liability" else "--observed-scale",
    shell:
        "test -n '{params.traits}' || "
        "(echo 'ERROR: set phase0_traits in config/workflow.yaml' >&2; exit 1); "
        "PYTHON_BIN={PYTHON} LDSC_PYTHON={LDSC_PYTHON} LDSC_DIR={LDSC_DIR} "
        "H2_OUT={output} bash scripts/03_h2_qc.sh {params.scale_flag} {params.traits}"


rule phase1_rg:
    input:
        h2="results/tables/h2_summary.tsv",
        munged=expand("data/munged/{trait}.sumstats.gz", trait=SELECTED),
    output:
        rg="results/tables/rg_matrix.tsv",
        inclusion="results/tables/phase1_inclusion.tsv",
    shell:
        "PYTHON_BIN={PYTHON} LDSC_PYTHON={LDSC_PYTHON} LDSC_DIR={LDSC_DIR} "
        "RG_OUT={output.rg} INCLUSION_OUT={output.inclusion} "
        "bash scripts/04_rg.sh --h2 {input.h2}"
