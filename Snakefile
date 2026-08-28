"""Manifest-driven orchestration for the locked atlas-v1.0 panel."""
import csv


configfile: "config/workflow.yaml"

PANEL = config["panel"]
PANEL_LOCK = config["panel_lock"]
SOURCES = config["sources"]
SCHEMAS = config["schemas"]
VARIANT_MAPPINGS = config["variant_mappings"]
LIFTOVER_PLANS = config["liftover_plans"]
PYTHON = config["python"]
LDSC_PYTHON = config["ldsc_python"]
LDSC_DIR = config["ldsc_dir"]
RSCRIPT = config["rscript"]
EUR_LD_DIR = config["eur_ld_dir"]
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
        variant_mappings=VARIANT_MAPPINGS,
        liftover_plans=LIFTOVER_PLANS,
    output:
        "results/tables/source_readiness.tsv",
    shell:
        "{PYTHON} scripts/10_phase0_audit.py --config {input.manifest} "
        "--lock {input.lock} --sources {input.sources} --out {output} --strict"
        " --schemas {input.schemas}"
        " --variant-mappings {input.variant_mappings}"
        " --liftover-plans {input.liftover_plans}"


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


rule full_covariance:
    input:
        panel=PANEL,
        lock=PANEL_LOCK,
        munged=expand("data/munged/{trait}.sumstats.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
        ld_scores=expand(EUR_LD_DIR + "/{chromosome}.l2.ldscore.gz", chromosome=range(1, 23)),
        ld_m=expand(EUR_LD_DIR + "/{chromosome}.l2.M_5_50", chromosome=range(1, 23)),
    output:
        covariance="results/tables/ldsc_covariance_45x45.tsv",
        correlations="results/tables/ldsc_genetic_correlation_45x45.tsv",
        intercepts="results/tables/ldsc_intercept_45x45.tsv",
        sampling="results/tables/ldsc_sampling_covariance_1035x1035.tsv.gz",
        pairs="results/tables/ldsc_covariance_pairs.tsv",
        scales="results/tables/ldsc_covariance_trait_scales.tsv",
        metadata="results/tables/ldsc_covariance_metadata.tsv",
        structure="results/tables/ldsc_covariance_structure.rds",
    shell:
        "R_BIN={RSCRIPT} bash scripts/25_genomicsem_covariance.sh "
        "--panel {input.panel} --munged-dir data/munged --ld-dir {EUR_LD_DIR} "
        "--out-dir results/tables --log-prefix results/logs/genomicsem/ldsc_45_trait"
