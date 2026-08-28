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
LAVA_REFERENCE_PREFIX = config["lava_reference_prefix"]
LAVA_LOCUS_FILE = config["lava_locus_file"]
MIXER_POLICY = config["mixer_policy"]
PLEIOTROPY_POLICY = config["pleiotropy_policy"]
DOWNSTREAM_POLICY = config["downstream_policy"]
ATLAS_SCHEMA = config["atlas_schema"]
SELECTED = config.get("phase0_traits", [])
H2_SCALE = config.get("h2_scale", "liability")

with open(PANEL, newline="", encoding="utf-8") as handle:
    PANEL_ROWS = list(csv.DictReader(handle, delimiter="\t"))
RAW_BY_TRAIT = {row["trait_id"]: f"data/raw/{row['raw_file']}" for row in PANEL_ROWS}
PANEL_IDS = set(RAW_BY_TRAIT)
SLEEP_TRAITS = [row["trait_id"] for row in PANEL_ROWS if row["domain"] == "sleep"]
NON_SLEEP_TRAITS = [row["trait_id"] for row in PANEL_ROWS if row["domain"] != "sleep"]
PLEIOTROPY_PAIRS = [f"{sleep}__{other}" for sleep in SLEEP_TRAITS for other in NON_SLEEP_TRAITS]
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
        diagnostics="results/tables/ldsc_covariance_diagnostics.tsv",
        structure="results/tables/ldsc_covariance_structure.rds",
    shell:
        "R_BIN={RSCRIPT} bash scripts/25_genomicsem_covariance.sh "
        "--panel {input.panel} --munged-dir data/munged --ld-dir {EUR_LD_DIR} "
        "--out-dir results/tables --log-prefix results/logs/genomicsem/ldsc_45_trait"


rule genomic_sem_trait_qc:
    input:
        panel=PANEL,
        h2="results/tables/h2_summary.tsv",
        pairs="results/tables/ldsc_covariance_pairs.tsv",
    output:
        "results/tables/genomicsem_trait_inclusion.tsv",
    shell:
        "{PYTHON} scripts/27_genomicsem_trait_qc.py --panel {input.panel} "
        "--phase1-h2 {input.h2} --covariance-pairs {input.pairs} --out {output}"


rule chromosome_split_inputs:
    input:
        panel=PANEL,
        inclusion="results/tables/genomicsem_trait_inclusion.tsv",
        munged=expand("data/munged/{trait}.sumstats.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
    output:
        provenance="results/tables/chromosome_split_provenance.json",
        odd=expand("data/munged_chromosome_split/odd/{trait}.sumstats.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
        even=expand("data/munged_chromosome_split/even/{trait}.sumstats.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
    shell:
        "{PYTHON} scripts/28_prepare_chromosome_split.py"


rule genomic_sem_discovery_odd:
    input:
        inclusion="results/tables/genomicsem_trait_inclusion.tsv",
        split="results/tables/chromosome_split_provenance.json",
    output:
        structure="results/tables/genomicsem_discovery_odd.rds",
        metadata="results/tables/genomicsem_discovery_odd_metadata.tsv",
    shell:
        "{RSCRIPT} scripts/28_chromosome_split_covariance.R --split odd"


rule genomic_sem_validation_even:
    input:
        inclusion="results/tables/genomicsem_trait_inclusion.tsv",
        split="results/tables/chromosome_split_provenance.json",
    output:
        structure="results/tables/genomicsem_validation_even.rds",
        metadata="results/tables/genomicsem_validation_even_metadata.tsv",
    shell:
        "{RSCRIPT} scripts/28_chromosome_split_covariance.R --split even "
        "--munged-dir data/munged_chromosome_split/even --ld-dir ref/eur_w_ld_chr_even "
        "--out {output.structure} --metadata {output.metadata} "
        "--log-prefix results/logs/genomicsem/validation_even"


rule genomic_sem_model:
    input:
        discovery="results/tables/genomicsem_discovery_odd.rds",
        validation="results/tables/genomicsem_validation_even.rds",
        inclusion="results/tables/genomicsem_trait_inclusion.tsv",
    output:
        efa="results/tables/genomic_sem_efa_models.tsv",
        fits="results/tables/genomic_sem_model_fit.tsv",
        loadings="results/tables/genomic_sem_factor_loadings.tsv",
        syntax="results/tables/genomic_sem_model_syntax.tsv",
        diagnostics="results/tables/genomic_sem_split_diagnostics.tsv",
    shell:
        "{RSCRIPT} scripts/29_genomicsem_model.R"


rule lava_inputs:
    input:
        panel=PANEL,
        policy="config/lava_analysis_policy.json",
        intercepts="results/tables/ldsc_intercept_45x45.tsv",
        rg="results/tables/rg_matrix.tsv",
        munged=expand("data/munged/{trait}.sumstats.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
    output:
        info="results/tables/lava_input_info.tsv",
        overlap="results/tables/lava_sample_overlap.txt",
        pairs="results/tables/lava_pair_manifest.tsv",
        provenance="results/tables/lava_input_provenance.tsv",
        runtime="results/tables/lava_runtime_policy.tsv",
        diagnostics="results/tables/lava_input_diagnostics.json",
    shell:
        "{PYTHON} scripts/31_prepare_lava.py"


rule lava:
    input:
        info="results/tables/lava_input_info.tsv",
        overlap="results/tables/lava_sample_overlap.txt",
        pairs="results/tables/lava_pair_manifest.tsv",
        provenance="results/tables/lava_input_provenance.tsv",
        runtime="results/tables/lava_runtime_policy.tsv",
        locus=LAVA_LOCUS_FILE,
        reference=expand(
            LAVA_REFERENCE_PREFIX + "_chr{chromosome}.{suffix}",
            chromosome=range(1, 23), suffix=["info", "bcor"],
        ),
    output:
        status="results/tables/lava_locus_status.tsv",
        univariate="results/tables/lava_univariate.tsv",
        bivariate="results/tables/lava_bivariate.tsv",
    shell:
        "{RSCRIPT} scripts/33_run_lava.R && "
        "{PYTHON} scripts/34_validate_lava.py --quiet"


rule mixer_preflight:
    input:
        panel=PANEL,
        policy=MIXER_POLICY,
        prefilters="config/hm3_prefilter_plans.tsv",
        harmonized=expand("data/harmonized/{trait}.harmonized.tsv.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
        qc=expand("data/harmonized/{trait}.qc.txt", trait=[row["trait_id"] for row in PANEL_ROWS]),
    output:
        report="results/tables/mixer_preflight.json",
        traits="results/tables/mixer_input_readiness.tsv",
    shell:
        "{PYTHON} scripts/35_mixer_preflight.py --report-only"


rule pleiotropy_preflight:
    input:
        panel=PANEL,
        policy=PLEIOTROPY_POLICY,
        sources="config/pleiotropy_reference_sources.tsv",
        prefilters="config/hm3_prefilter_plans.tsv",
        harmonized=expand("data/harmonized/{trait}.harmonized.tsv.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
        qc=expand("data/harmonized/{trait}.qc.txt", trait=[row["trait_id"] for row in PANEL_ROWS]),
    output:
        report="results/tables/pleiotropy_preflight.json",
        traits="results/tables/pleiotropy_input_readiness.tsv",
    shell:
        "{PYTHON} scripts/42_pleiotropy_preflight.py --report-only"


rule pleiotropy_pairs:
    input:
        panel=PANEL,
        policy=PLEIOTROPY_POLICY,
        prefilters="config/hm3_prefilter_plans.tsv",
        rg="results/tables/rg_matrix.tsv",
        harmonized=expand("data/harmonized/{trait}.harmonized.tsv.gz", trait=[row["trait_id"] for row in PANEL_ROWS]),
        qc=expand("data/harmonized/{trait}.qc.txt", trait=[row["trait_id"] for row in PANEL_ROWS]),
    output:
        manifest="results/tables/pleiotropy_pair_manifest.tsv",
        lock="results/tables/pleiotropy_pair_manifest.lock.json",
    shell:
        "{PYTHON} scripts/43_prepare_pleiotropy_pairs.py --report-only"


rule pleiotropy_pair_input:
    input:
        policy=PLEIOTROPY_POLICY,
    output:
        pair="results/pleiotropy/inputs/{pair_id}.tsv.gz",
        provenance="results/pleiotropy/inputs/{pair_id}.provenance.json",
    shell:
        "{PYTHON} scripts/44_materialize_pleiotropy_pair.py {wildcards.pair_id} --materialize"


rule placo_task:
    input:
        pair="results/pleiotropy/inputs/{pair_id}.tsv.gz",
        provenance="results/pleiotropy/inputs/{pair_id}.provenance.json",
    output:
        task="results/pleiotropy/tasks/{pair_id}.tsv",
        lock="results/pleiotropy/tasks/{pair_id}.lock.tsv",
    shell:
        "{PYTHON} scripts/45_prepare_placo_task.py {wildcards.pair_id}"


rule placo_pair:
    input:
        task="results/pleiotropy/tasks/{pair_id}.tsv",
        lock="results/pleiotropy/tasks/{pair_id}.lock.tsv",
    output:
        hits="results/pleiotropy/placo/{pair_id}.hits.tsv",
        summary="results/pleiotropy/placo/{pair_id}.summary.tsv",
    shell:
        "{RSCRIPT} scripts/46_run_placo_pair.R {input.task} {input.lock} --execute"


rule pleiofdr_trait:
    input:
        source="data/harmonized/{trait}.harmonized.tsv.gz",
        qc="data/harmonized/{trait}.qc.txt",
        template="ref/pleiofdr/9545380.ref",
        policy=PLEIOTROPY_POLICY,
    output:
        mat="data/pleiofdr/{trait}.mat",
        provenance="data/pleiofdr/{trait}.provenance.json",
    shell:
        "{PYTHON} scripts/47_prepare_pleiofdr_trait.py {wildcards.trait} --materialize"


rule conjfdr_task:
    input:
        sleep=lambda wildcards: "data/pleiofdr/" + wildcards.pair_id.split("__", 1)[0] + ".mat",
        non_sleep=lambda wildcards: "data/pleiofdr/" + wildcards.pair_id.split("__", 1)[1] + ".mat",
        template="ref/pleiofdr/9545380.ref",
        reference="ref/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat",
        patch="patches/pleiofdr-enable-overlap.patch",
    output:
        config="results/pleiotropy/conjfdr_tasks/{pair_id}.config.txt",
        task="results/pleiotropy/conjfdr_tasks/{pair_id}.tsv",
        lock="results/pleiotropy/conjfdr_tasks/{pair_id}.lock.tsv",
    shell:
        "{PYTHON} scripts/48_prepare_conjfdr_task.py {wildcards.pair_id}"


rule conjfdr_pair:
    input:
        task="results/pleiotropy/conjfdr_tasks/{pair_id}.tsv",
        lock="results/pleiotropy/conjfdr_tasks/{pair_id}.lock.tsv",
    output:
        completion="results/pleiotropy/conjfdr/{pair_id}/atlas_completion.tsv",
    shell:
        "{PYTHON} scripts/49_run_conjfdr_pair.py {input.task} {input.lock} --execute"


rule pleiotropy:
    input:
        placo_hits=expand("results/pleiotropy/placo/{pair_id}.hits.tsv", pair_id=PLEIOTROPY_PAIRS),
        placo_summaries=expand("results/pleiotropy/placo/{pair_id}.summary.tsv", pair_id=PLEIOTROPY_PAIRS),
        conjfdr=expand("results/pleiotropy/conjfdr/{pair_id}/atlas_completion.tsv", pair_id=PLEIOTROPY_PAIRS),
        locus=LAVA_LOCUS_FILE,
    output:
        placo="results/tables/placo_loci.tsv",
        conjfdr="results/tables/conjfdr_loci.tsv",
        shared="results/atlas/shared_loci.tsv",
        provenance="results/atlas/shared_loci.provenance.json",
    shell:
        "{PYTHON} scripts/50_collate_pleiotropy.py"


rule atlas_core:
    input:
        panel=PANEL,
        lock=PANEL_LOCK,
        h2="results/tables/h2_summary.tsv",
        rg="results/tables/rg_matrix.tsv",
        policy=DOWNSTREAM_POLICY,
    output:
        traits="results/atlas/traits.tsv",
        pairs="results/atlas/trait_pairs.tsv",
        provenance="results/atlas/core.provenance.json",
    shell:
        "{PYTHON} scripts/51_build_atlas_core.py"


rule validate_integrated_atlas:
    input:
        policy=DOWNSTREAM_POLICY,
        schema=ATLAS_SCHEMA,
        tables=expand(
            "results/atlas/{name}",
            name=[
                "traits.tsv", "trait_pairs.tsv", "loci.tsv", "variants.tsv",
                "genes.tsv", "regulatory_elements.tsv", "cell_types.tsv",
                "pathways.tsv", "causal_tests.tsv", "edges.tsv",
            ],
        ),
    output:
        touch("results/atlas/ATLAS_SCHEMA_OK"),
    shell:
        "{PYTHON} scripts/52_validate_integrated_atlas.py --quiet"


rule validate_robustness:
    input:
        policy=DOWNSTREAM_POLICY,
        summary="results/tables/robustness_summary.tsv",
    output:
        touch("results/tables/ROBUSTNESS_OK"),
    shell:
        "{PYTHON} scripts/53_validate_robustness.py --quiet"


rule atlas_release:
    input:
        atlas="results/atlas/ATLAS_SCHEMA_OK",
        robustness="results/tables/ROBUSTNESS_OK",
    output:
        touch("releases/ATLAS_V1_RELEASE_OK"),
    shell:
        "{PYTHON} scripts/54_build_release.py --execute && "
        "{PYTHON} scripts/55_validate_release.py --quiet"
