"""Manifest-driven orchestration for the locked atlas-v1.0 panel."""
import csv
import json


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
FINE_MAPPING_POLICY = config["fine_mapping_policy"]
MOLECULAR_POLICY = config.get("molecular_policy", "config/molecular_analysis_policy.json")
INTERPRETATION_POLICY = config.get("interpretation_policy", "config/interpretation_analysis_policy.json")
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

with open(INTERPRETATION_POLICY, encoding="utf-8") as handle:
    INTERPRETATION_SPEC = json.load(handle)
with open(INTERPRETATION_SPEC["screen_registry_v4"]["component_manifest"], encoding="utf-8") as handle:
    SCREEN_SOURCE_BUNDLE = json.load(handle)
SCREEN_SOURCE_PATHS = [component["path"] for component in SCREEN_SOURCE_BUNDLE["components"]]
with open(INTERPRETATION_SPEC["hocomoco_v14"]["component_manifest"], encoding="utf-8") as handle:
    HOCOMOCO_SOURCE_BUNDLE = json.load(handle)
HOCOMOCO_SOURCE_PATHS = [component["path"] for component in HOCOMOCO_SOURCE_BUNDLE["components"]]
with open(INTERPRETATION_SPEC["abc_2021"]["component_manifest"], encoding="utf-8") as handle:
    ABC_SOURCE_BUNDLE = json.load(handle)
ABC_SOURCE_PATHS = [component["path"] for component in ABC_SOURCE_BUNDLE["components"]]
with open(INTERPRETATION_SPEC["pchic_2016"]["component_manifest"], encoding="utf-8") as handle:
    PCHIC_SOURCE_BUNDLE = json.load(handle)
PCHIC_SOURCE_PATHS = [component["path"] for component in PCHIC_SOURCE_BUNDLE["components"]]
with open(INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"], encoding="utf-8") as handle:
    FUMA_SOURCE_BUNDLE = json.load(handle)
FUMA_SOURCE_PATHS = [component["path"] for component in FUMA_SOURCE_BUNDLE["components"]]
FUMA_MATRIX_PATHS = [
    f"{INTERPRETATION_SPEC['fuma_scrna']['matrix_cache_dir']}/{dataset['dataset_id']}.txt"
    for dataset in FUMA_SOURCE_BUNDLE["datasets"]
]
MAGMA_REFERENCE_PATHS = [
    f"{INTERPRETATION_SPEC['fuma_scrna']['reference_dir']}/{member['name']}"
    for member in FUMA_SOURCE_BUNDLE["reference_members"]
]


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


rule fine_mapping_preflight:
    input:
        panel=PANEL,
        policy=FINE_MAPPING_POLICY,
        downstream=DOWNSTREAM_POLICY,
        sources="config/fine_mapping_sources.tsv",
        references="config/fine_mapping_method_references.tsv",
    output:
        report="results/tables/fine_mapping_preflight.json",
        traits="results/tables/fine_mapping_input_readiness.tsv",
    shell:
        "{PYTHON} scripts/56_finemapping_preflight.py --report-only"


checkpoint fine_mapping_loci:
    input:
        preflight="results/tables/fine_mapping_preflight.json",
        shared="results/atlas/shared_loci.tsv",
        policy=FINE_MAPPING_POLICY,
    output:
        manifest="results/tables/fine_mapping_locus_manifest.tsv",
        lock="results/tables/fine_mapping_locus_manifest.lock.json",
    shell:
        "{PYTHON} scripts/57_prepare_finemapping_loci.py"


rule fine_mapping_input:
    input:
        manifest="results/tables/fine_mapping_locus_manifest.tsv",
        lock="results/tables/fine_mapping_locus_manifest.lock.json",
        preflight="results/tables/fine_mapping_preflight.json",
        references=expand(
            LAVA_REFERENCE_PREFIX + "_chr{chromosome}.{suffix}",
            chromosome=range(1, 23), suffix=["info", "bcor"],
        ),
        extracted="ref/lava/ukb_v1.1/extracted_manifest.tsv",
    output:
        summary1="data/fine_mapping/{shared_locus_id}/sleep.tsv.gz",
        summary2="data/fine_mapping/{shared_locus_id}/non_sleep.tsv.gz",
        order="data/fine_mapping/{shared_locus_id}/variants.tsv",
        ld="data/fine_mapping/{shared_locus_id}/ld.tsv.gz",
        task="results/fine_mapping/tasks/{shared_locus_id}.tsv",
        task_lock="results/fine_mapping/tasks/{shared_locus_id}.lock.json",
    shell:
        "{PYTHON} scripts/58_materialize_finemapping_locus.py {wildcards.shared_locus_id} --materialize"


rule fine_mapping_locus:
    input:
        task="results/fine_mapping/tasks/{shared_locus_id}.tsv",
        task_lock="results/fine_mapping/tasks/{shared_locus_id}.lock.json",
    output:
        variants="results/fine_mapping/runs/{shared_locus_id}/variants.tsv",
        credible="results/fine_mapping/runs/{shared_locus_id}/credible_sets.tsv",
        coloc="results/fine_mapping/runs/{shared_locus_id}/colocalization.tsv",
        shared="results/fine_mapping/runs/{shared_locus_id}/shared_variant_posteriors.tsv",
        diagnostics="results/fine_mapping/runs/{shared_locus_id}/diagnostics.tsv",
        provenance="results/fine_mapping/runs/{shared_locus_id}/provenance.json",
    shell:
        "{PYTHON} scripts/59_run_finemapping_locus.py {wildcards.shared_locus_id} --execute"


def fine_mapping_run_provenance(wildcards):
    checkpoint_output = checkpoints.fine_mapping_loci.get(**wildcards).output.manifest
    with open(checkpoint_output, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return expand(
        "results/fine_mapping/runs/{shared_locus_id}/provenance.json",
        shared_locus_id=[row["shared_locus_id"] for row in rows],
    )


rule fine_mapping:
    input:
        runs=fine_mapping_run_provenance,
        manifest="results/tables/fine_mapping_locus_manifest.tsv",
        lock="results/tables/fine_mapping_locus_manifest.lock.json",
    output:
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        credible="results/tables/fine_mapping_credible_sets.tsv",
        diagnostics="results/tables/fine_mapping_diagnostics.tsv",
        coloc="results/tables/trait_trait_colocalization.tsv",
        provenance="results/atlas/fine_mapping.provenance.json",
    shell:
        "{PYTHON} scripts/60_collate_finemapping.py"


rule molecular_metadata:
    input:
        policy=MOLECULAR_POLICY,
    output:
        eqtl_metadata="ref/molecular/eqtl_catalogue_r7/dataset_metadata_r7.tsv",
        eqtl_paths="ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths.tsv",
        gtex_paths="ref/molecular/eqtl_catalogue_r7/tabix_ftp_paths_imported.tsv",
        genes="ref/molecular/gtex_v8/gencode.v26.GRCh38.genes.gtf",
    shell:
        "bash scripts/61_fetch_molecular_metadata.sh"


rule molecular_preflight:
    input:
        policy=MOLECULAR_POLICY,
        sources="config/molecular_source_registry.tsv",
        metadata=rules.molecular_metadata.output,
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        trait_coloc="results/tables/trait_trait_colocalization.tsv",
        fine_mapping_provenance="results/atlas/fine_mapping.provenance.json",
    output:
        readiness="results/tables/molecular_input_readiness.tsv",
        report="results/tables/molecular_preflight.json",
    shell:
        "{PYTHON} scripts/61_molecular_preflight.py"


checkpoint molecular_search_plan:
    input:
        preflight="results/tables/molecular_preflight.json",
        policy=MOLECULAR_POLICY,
        sources="config/molecular_source_registry.tsv",
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        fine_mapping_provenance="results/atlas/fine_mapping.provenance.json",
    output:
        plan="results/tables/molecular_search_plan.tsv",
        lock="results/tables/molecular_search_plan.lock.json",
    shell:
        "{PYTHON} scripts/62_prepare_molecular_search_plan.py"


rule molecular_source_search:
    input:
        plan="results/tables/molecular_search_plan.tsv",
        lock="results/tables/molecular_search_plan.lock.json",
        variants="results/atlas/variants.tsv",
    output:
        normalized="results/molecular/search/{search_task_id}/normalized_qtl.tsv.gz",
        provenance="results/molecular/search/{search_task_id}/provenance.json",
    shell:
        "{PYTHON} scripts/63_run_molecular_search.py {wildcards.search_task_id} --execute"


def molecular_search_provenance(wildcards):
    plan = checkpoints.molecular_search_plan.get(**wildcards).output.plan
    with open(plan, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return expand(
        "results/molecular/search/{search_task_id}/provenance.json",
        search_task_id=[row["search_task_id"] for row in rows],
    )


checkpoint molecular_features:
    input:
        searches=molecular_search_provenance,
        plan="results/tables/molecular_search_plan.tsv",
        lock="results/tables/molecular_search_plan.lock.json",
    output:
        coverage="results/tables/molecular_search_coverage.tsv",
        manifest="results/tables/molecular_feature_manifest.tsv",
        exclusions="results/tables/molecular_feature_exclusions.tsv",
        lock="results/tables/molecular_feature_manifest.lock.json",
    shell:
        "{PYTHON} scripts/64_lock_molecular_features.py"


rule molecular_coloc_input:
    input:
        manifest="results/tables/molecular_feature_manifest.tsv",
        lock="results/tables/molecular_feature_manifest.lock.json",
    output:
        task="results/molecular/tasks/{comparison_id}.tsv",
        task_lock="results/molecular/tasks/{comparison_id}.lock.json",
        trait="data/molecular_coloc/{comparison_id}/trait.tsv.gz",
        molecular="data/molecular_coloc/{comparison_id}/molecular.tsv.gz",
        variants="data/molecular_coloc/{comparison_id}/variants.tsv",
        ld="data/molecular_coloc/{comparison_id}/ld.tsv.gz",
    shell:
        "{PYTHON} scripts/65_materialize_molecular_coloc.py {wildcards.comparison_id} --materialize"


rule molecular_coloc_run:
    input:
        task="results/molecular/tasks/{comparison_id}.tsv",
        task_lock="results/molecular/tasks/{comparison_id}.lock.json",
        trait="data/molecular_coloc/{comparison_id}/trait.tsv.gz",
        molecular="data/molecular_coloc/{comparison_id}/molecular.tsv.gz",
        variants="data/molecular_coloc/{comparison_id}/variants.tsv",
        ld="data/molecular_coloc/{comparison_id}/ld.tsv.gz",
    output:
        variants="results/molecular/coloc_runs/{comparison_id}/variants.tsv",
        credible="results/molecular/coloc_runs/{comparison_id}/credible_sets.tsv",
        coloc="results/molecular/coloc_runs/{comparison_id}/colocalization.tsv",
        shared="results/molecular/coloc_runs/{comparison_id}/shared_variant_posteriors.tsv",
        diagnostics="results/molecular/coloc_runs/{comparison_id}/diagnostics.tsv",
        provenance="results/molecular/coloc_runs/{comparison_id}/provenance.json",
    shell:
        "{PYTHON} scripts/66_run_molecular_coloc.py {wildcards.comparison_id} --execute"


def molecular_coloc_provenance(wildcards):
    manifest = checkpoints.molecular_features.get(**wildcards).output.manifest
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return expand(
        "results/molecular/coloc_runs/{comparison_id}/provenance.json",
        comparison_id=[row["comparison_id"] for row in rows],
    )


checkpoint twas_phi_inventory:
    input:
        policy=MOLECULAR_POLICY,
    output:
        inventory="results/tables/twas_phi_model_inventory.tsv",
        lock="results/tables/twas_phi_model_inventory.lock.json",
        snapshots=expand("results/sources/twas_phi_model_inventory/page_{page:02d}.html", page=range(1, 6)),
    shell:
        "{PYTHON} scripts/67_lock_twas_model_inventory.py --execute"


rule twas_model_download:
    input:
        inventory="results/tables/twas_phi_model_inventory.tsv",
        lock="results/tables/twas_phi_model_inventory.lock.json",
    output:
        "results/tables/twas_phi_model_downloads.lock.json",
    params:
        acknowledgement="--acknowledge-large-download" if config.get("acknowledge_large_downloads", False) else "",
    shell:
        "bash scripts/67_fetch_twas_resources.sh --models {params.acknowledgement}"


rule metaxcan_runtime:
    input:
        policy=MOLECULAR_POLICY,
    output:
        archive=".molecular-env/source_archives/MetaXcan_v0.8.1.tar.gz",
        entrypoint=".molecular-env/MetaXcan/software/SPrediXcan.py",
        python=".molecular-env/python/bin/python",
    shell:
        "bash scripts/67_fetch_twas_resources.sh --runtime"


rule twas_model_index:
    input:
        inventory="results/tables/twas_phi_model_inventory.tsv",
        inventory_lock="results/tables/twas_phi_model_inventory.lock.json",
        downloads="results/tables/twas_phi_model_downloads.lock.json",
    output:
        registry="results/tables/twas_model_registry.tsv",
        variants="results/tables/twas_model_variants/GTEx_v8_ELASTIC_NET_PHI_eQTL.tsv.gz",
        phi_exclusions="results/tables/twas_model_phi_exclusions.tsv",
        lock="results/tables/twas_model_registry.lock.json",
    shell:
        "{PYTHON} scripts/68_index_twas_models.py --materialize"


rule twas_preflight:
    input:
        models=rules.twas_model_index.output,
        runtime=rules.metaxcan_runtime.output,
        fine_mapping="results/atlas/fine_mapping.provenance.json",
    output:
        readiness="results/tables/twas_input_readiness.tsv",
        report="results/tables/twas_preflight.json",
    shell:
        "{PYTHON} scripts/61_molecular_preflight.py --readiness-out {output.readiness} --out {output.report}"


checkpoint twas_manifest:
    input:
        preflight="results/tables/twas_preflight.json",
        models="results/tables/twas_model_registry.tsv",
        models_lock="results/tables/twas_model_registry.lock.json",
        traits="results/atlas/traits.tsv",
    output:
        mapping="results/tables/twas_gwas_mapping_manifest.tsv",
        eligibility="results/tables/twas_trait_eligibility.tsv",
        runs="results/tables/twas_run_manifest.tsv",
        lock="results/tables/twas_run_manifest.lock.json",
    shell:
        "{PYTHON} scripts/69_prepare_twas_manifest.py --preflight {input.preflight}"


rule twas_gwas_mapping:
    input:
        manifest="results/tables/twas_gwas_mapping_manifest.tsv",
        lock="results/tables/twas_run_manifest.lock.json",
    output:
        data="data/twas/{model_family}/{trait_id}.tsv.gz",
        lock="data/twas/{model_family}/{trait_id}.tsv.lock.json",
    shell:
        "{PYTHON} scripts/70_materialize_twas_gwas.py {wildcards.trait_id} {wildcards.model_family} --materialize"


def twas_mapping_locks(wildcards):
    manifest = checkpoints.twas_manifest.get(**wildcards).output.mapping
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row["mapped_gwas_lock_path"] for row in rows]


def twas_run_dependencies(wildcards):
    manifest = checkpoints.twas_manifest.get().output.runs
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["run_id"] == wildcards.run_id]
    if len(selected) != 1:
        raise ValueError(f"unknown TWAS run ID: {wildcards.run_id}")
    row = selected[0]
    return [row["mapped_gwas_path"], row["mapped_gwas_lock_path"]]


rule twas_run:
    input:
        dependencies=twas_run_dependencies,
        manifest="results/tables/twas_run_manifest.tsv",
        lock="results/tables/twas_run_manifest.lock.json",
        runtime=".molecular-env/python/bin/python",
    output:
        result="results/twas/runs/{run_id}/spredixcan.csv",
        provenance="results/twas/runs/{run_id}/provenance.json",
    shell:
        "{PYTHON} scripts/71_run_twas.py {wildcards.run_id} --execute"


def twas_run_provenance(wildcards):
    manifest = checkpoints.twas_manifest.get(**wildcards).output.runs
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return expand(
        "results/twas/runs/{run_id}/provenance.json",
        run_id=[row["run_id"] for row in rows],
    )


rule twas:
    input:
        mappings=twas_mapping_locks,
        runs=twas_run_provenance,
        manifest="results/tables/twas_run_manifest.tsv",
        lock="results/tables/twas_run_manifest.lock.json",
        eligibility="results/tables/twas_trait_eligibility.tsv",
    output:
        results="results/tables/twas.tsv",
        coverage="results/tables/twas_coverage.tsv",
        provenance="results/tables/twas.provenance.json",
    shell:
        "{PYTHON} scripts/72_collate_twas.py"


rule molecular_integration:
    input:
        qtl_runs=molecular_coloc_provenance,
        twas=rules.twas.output,
        feature_manifest="results/tables/molecular_feature_manifest.tsv",
        feature_lock="results/tables/molecular_feature_manifest.lock.json",
        search_coverage="results/tables/molecular_search_coverage.tsv",
        search_plan="results/tables/molecular_search_plan.tsv",
        search_plan_lock="results/tables/molecular_search_plan.lock.json",
        models=rules.twas_model_index.output,
        trait_coloc="results/tables/trait_trait_colocalization.tsv",
        loci="results/atlas/loci.tsv",
        traits="results/atlas/traits.tsv",
        genes_annotation="ref/molecular/gtex_v8/gencode.v26.GRCh38.genes.gtf",
    output:
        molecular_coloc="results/tables/molecular_qtl_colocalization.tsv",
        all_coloc="results/tables/colocalization.tsv",
        coverage="results/tables/molecular_locus_coverage.tsv",
        evidence="results/tables/molecular_evidence.tsv",
        genes="results/atlas/genes.tsv",
        provenance="results/atlas/molecular.provenance.json",
    shell:
        "{PYTHON} scripts/73_collate_molecular.py"


rule interpretation_preflight:
    input:
        policy=INTERPRETATION_POLICY,
        downstream=DOWNSTREAM_POLICY,
        sources="config/interpretation_source_registry.tsv",
        references="config/interpretation_method_references.tsv",
        source_manifests=[
            INTERPRETATION_SPEC["screen_registry_v4"]["component_manifest"],
            INTERPRETATION_SPEC["hocomoco_v14"]["component_manifest"],
            INTERPRETATION_SPEC["abc_2021"]["component_manifest"],
            INTERPRETATION_SPEC["pchic_2016"]["component_manifest"],
            INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"],
        ],
        source_files=(
            SCREEN_SOURCE_PATHS + HOCOMOCO_SOURCE_PATHS + ABC_SOURCE_PATHS
            + PCHIC_SOURCE_PATHS + FUMA_SOURCE_PATHS
        ),
        molecular=rules.molecular_integration.output,
        traits="results/atlas/traits.tsv",
        pairs="results/atlas/trait_pairs.tsv",
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        genes="results/atlas/genes.tsv",
    output:
        readiness="results/tables/interpretation_source_readiness.tsv",
        report="results/tables/interpretation_preflight.json",
    shell:
        "{PYTHON} scripts/74_interpretation_preflight.py --report-only"


rule abc_overlap_cache:
    input:
        policy=INTERPRETATION_POLICY,
        sources="config/interpretation_source_registry.tsv",
        manifest=INTERPRETATION_SPEC["abc_2021"]["component_manifest"],
        source=ABC_SOURCE_PATHS,
        variants="results/atlas/variants.tsv",
        genes="results/atlas/genes.tsv",
    output:
        cache=INTERPRETATION_SPEC["abc_2021"]["cache_path"],
        provenance=INTERPRETATION_SPEC["abc_2021"]["cache_provenance_path"],
    shell:
        "{PYTHON} scripts/84_prepare_abc_overlap_cache.py"


rule pchic_overlap_cache:
    input:
        policy=INTERPRETATION_POLICY,
        sources="config/interpretation_source_registry.tsv",
        manifest=INTERPRETATION_SPEC["pchic_2016"]["component_manifest"],
        source=PCHIC_SOURCE_PATHS,
        variants="results/atlas/variants.tsv",
        genes="results/atlas/genes.tsv",
    output:
        cache=INTERPRETATION_SPEC["pchic_2016"]["cache_path"],
        provenance=INTERPRETATION_SPEC["pchic_2016"]["cache_provenance_path"],
    shell:
        "{PYTHON} scripts/86_prepare_pchic_overlap_cache.py"


rule fuma_scrna_matrices:
    input:
        policy=INTERPRETATION_POLICY,
        manifest=INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"],
        source=FUMA_SOURCE_PATHS,
    output:
        matrices=FUMA_MATRIX_PATHS,
        provenance=INTERPRETATION_SPEC["fuma_scrna"]["matrix_cache_provenance_path"],
    shell:
        "{PYTHON} scripts/88_materialize_fuma_resources.py --matrices"


rule magma_eur_reference:
    input:
        policy=INTERPRETATION_POLICY,
        manifest=INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"],
        archive=INTERPRETATION_SPEC["fuma_scrna"]["reference_archive_path"],
    output:
        members=MAGMA_REFERENCE_PATHS,
        provenance=INTERPRETATION_SPEC["fuma_scrna"]["reference_extract_provenance_path"],
    params:
        acknowledgement=(
            "--acknowledge-large-extract"
            if config.get("acknowledge_large_extracts", False) else ""
        ),
    shell:
        "{PYTHON} scripts/88_materialize_fuma_resources.py --reference {params.acknowledgement}"


rule magma_gene_annotation:
    input:
        policy=INTERPRETATION_POLICY,
        manifest=INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"],
        reference=rules.magma_eur_reference.output,
        binary=INTERPRETATION_SPEC["fuma_scrna"]["magma_binary_path"],
        genes=INTERPRETATION_SPEC["fuma_scrna"]["gene_location_path"],
    output:
        annotation=INTERPRETATION_SPEC["fuma_scrna"]["gene_annotation_prefix"] + ".genes.annot",
        log=INTERPRETATION_SPEC["fuma_scrna"]["gene_annotation_prefix"] + ".log",
        provenance=INTERPRETATION_SPEC["fuma_scrna"]["gene_annotation_provenance_path"],
    shell:
        "{PYTHON} scripts/89_prepare_magma_annotation.py"


rule magma_gene_results:
    input:
        policy=INTERPRETATION_POLICY,
        manifest=INTERPRETATION_SPEC["fuma_scrna"]["component_manifest"],
        annotation=rules.magma_gene_annotation.output,
        reference=rules.magma_eur_reference.output,
        binary=INTERPRETATION_SPEC["fuma_scrna"]["magma_binary_path"],
        gwas=INTERPRETATION_SPEC["fuma_scrna"]["gwas_path_template"],
    output:
        raw=INTERPRETATION_SPEC["fuma_scrna"]["gene_results_dir"] + "/{trait_id}/{trait_id}.genes.raw",
        result=INTERPRETATION_SPEC["fuma_scrna"]["gene_results_dir"] + "/{trait_id}/{trait_id}.genes.out",
        log=INTERPRETATION_SPEC["fuma_scrna"]["gene_results_dir"] + "/{trait_id}/{trait_id}.log",
        provenance=INTERPRETATION_SPEC["fuma_scrna"]["gene_results_dir"] + "/{trait_id}/provenance.json",
    shell:
        "{PYTHON} scripts/90_prepare_magma_gene_results.py {wildcards.trait_id}"


checkpoint interpretation_tasks:
    input:
        preflight="results/tables/interpretation_preflight.json",
        policy=INTERPRETATION_POLICY,
        sources="config/interpretation_source_registry.tsv",
        references="config/interpretation_method_references.tsv",
        traits="results/atlas/traits.tsv",
        pairs="results/atlas/trait_pairs.tsv",
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        genes="results/atlas/genes.tsv",
    output:
        manifest="results/tables/interpretation_task_manifest.tsv",
        lock="results/tables/interpretation_task_manifest.lock.json",
    shell:
        "{PYTHON} scripts/75_prepare_interpretation_tasks.py"


def interpretation_task_dependencies(wildcards):
    manifest = checkpoints.interpretation_tasks.get(**wildcards).output.manifest
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [
        row for row in rows
        if row["task_id"] == wildcards.task_id and row["analysis_family"] == wildcards.family
    ]
    if len(selected) != 1:
        raise ValueError(f"unknown interpretation task ID: {wildcards.task_id}")
    task = selected[0]
    if task["source_id"] == INTERPRETATION_SPEC["fuma_scrna"]["source_id"]:
        gene_prefix = (
            f"{INTERPRETATION_SPEC['fuma_scrna']['gene_results_dir']}/"
            f"{task['trait_id']}/{task['trait_id']}"
        )
        return [
            *FUMA_MATRIX_PATHS,
            INTERPRETATION_SPEC["fuma_scrna"]["matrix_cache_provenance_path"],
            gene_prefix + ".genes.raw", gene_prefix + ".genes.out", gene_prefix + ".log",
            f"{INTERPRETATION_SPEC['fuma_scrna']['gene_results_dir']}/{task['trait_id']}/provenance.json",
        ]
    if task["source_id"] == INTERPRETATION_SPEC["pchic_2016"]["source_id"]:
        if task["domain"] in set(INTERPRETATION_SPEC["pchic_2016"]["available_domains"]):
            return [
                rules.pchic_overlap_cache.output.cache,
                rules.pchic_overlap_cache.output.provenance,
            ]
        return [
            INTERPRETATION_SPEC["pchic_2016"]["component_manifest"],
            *PCHIC_SOURCE_PATHS,
        ]
    if task["source_id"] == INTERPRETATION_SPEC["abc_2021"]["source_id"]:
        return [
            rules.abc_overlap_cache.output.cache,
            rules.abc_overlap_cache.output.provenance,
        ]
    if task["source_id"] == INTERPRETATION_SPEC["hocomoco_v14"]["source_id"]:
        screen_ids = set(INTERPRETATION_SPEC["screen_registry_v4"]["source_ids"])
        screen_tasks = [
            row for row in rows
            if row["analysis_family"] == "regulatory" and row["locus_id"] == task["locus_id"]
            and row["domain"] == task["domain"] and row["source_id"] in screen_ids
        ]
        if len(screen_tasks) != len(screen_ids):
            raise ValueError(f"motif task lacks exact SCREEN dependencies: {task['task_id']}")
        return [
            INTERPRETATION_SPEC["hocomoco_v14"]["component_manifest"],
            INTERPRETATION_SPEC["regulatory_build_harmonization"]["chain_path"],
            *HOCOMOCO_SOURCE_PATHS,
            *[row["normalized_result_path"] for row in screen_tasks],
            *[row["provenance_path"] for row in screen_tasks],
        ]
    if task["source_id"] in set(INTERPRETATION_SPEC["screen_registry_v4"]["source_ids"]):
        return [
            INTERPRETATION_SPEC["screen_registry_v4"]["component_manifest"],
            INTERPRETATION_SPEC["regulatory_build_harmonization"]["chain_path"],
            *SCREEN_SOURCE_PATHS,
        ]
    if task["source_id"] == INTERPRETATION_SPEC["promoter_mapping"]["source_id"]:
        return [
            INTERPRETATION_SPEC["regulatory_build_harmonization"]["chain_path"],
            "ref/molecular/gtex_v8/gencode.v26.GRCh38.genes.gtf",
        ]
    if task["method"] == "regulatory_overlap":
        return [row["provenance_path"] for row in rows if row["analysis_family"] == "regulatory"]
    if task["method"] == "QTL_colocalization":
        return ["results/tables/molecular_evidence.tsv", "results/atlas/molecular.provenance.json"]
    return []


rule interpretation_task:
    input:
        dependencies=interpretation_task_dependencies,
        manifest="results/tables/interpretation_task_manifest.tsv",
        lock="results/tables/interpretation_task_manifest.lock.json",
    output:
        result="results/interpretation/runs/{family}/{task_id}/normalized.tsv",
        provenance="results/interpretation/runs/{family}/{task_id}/provenance.json",
    shell:
        "{PYTHON} scripts/76_run_interpretation_task.py {wildcards.task_id} --execute"


def interpretation_task_provenance(wildcards):
    manifest = checkpoints.interpretation_tasks.get(**wildcards).output.manifest
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row["provenance_path"] for row in rows]


rule interpretation:
    input:
        task_provenance=interpretation_task_provenance,
        manifest="results/tables/interpretation_task_manifest.tsv",
        lock="results/tables/interpretation_task_manifest.lock.json",
    output:
        regulatory="results/atlas/regulatory_elements.tsv",
        cells="results/atlas/cell_types.tsv",
        pathways="results/atlas/pathways.tsv",
        causal="results/atlas/causal_tests.tsv",
        coverage="results/tables/interpretation_coverage.tsv",
        provenance="results/atlas/interpretation.provenance.json",
    shell:
        "{PYTHON} scripts/77_collate_interpretation.py"


rule atlas_edges:
    input:
        interpretation=rules.interpretation.output,
        traits="results/atlas/traits.tsv",
        pairs="results/atlas/trait_pairs.tsv",
        loci="results/atlas/loci.tsv",
        variants="results/atlas/variants.tsv",
        genes="results/atlas/genes.tsv",
    output:
        edges="results/atlas/edges.tsv",
        conclusions="results/tables/major_conclusions.tsv",
        provenance="results/atlas/edges.provenance.json",
    shell:
        "{PYTHON} scripts/78_build_atlas_edges.py"


checkpoint robustness_tasks:
    input:
        policy=INTERPRETATION_POLICY,
        downstream=DOWNSTREAM_POLICY,
        conclusions="results/tables/major_conclusions.tsv",
        edges="results/atlas/edges.tsv",
    output:
        manifest="results/tables/robustness_task_manifest.tsv",
        lock="results/tables/robustness_task_manifest.lock.json",
    shell:
        "{PYTHON} scripts/79_prepare_robustness_tasks.py"


rule robustness_task:
    input:
        manifest="results/tables/robustness_task_manifest.tsv",
        lock="results/tables/robustness_task_manifest.lock.json",
    output:
        result="results/robustness/runs/{task_id}/normalized.tsv",
        provenance="results/robustness/runs/{task_id}/provenance.json",
    shell:
        "{PYTHON} scripts/80_run_robustness_task.py {wildcards.task_id} --execute"


def robustness_task_provenance(wildcards):
    manifest = checkpoints.robustness_tasks.get(**wildcards).output.manifest
    with open(manifest, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row["provenance_path"] for row in rows]


rule robustness:
    input:
        task_provenance=robustness_task_provenance,
        manifest="results/tables/robustness_task_manifest.tsv",
        lock="results/tables/robustness_task_manifest.lock.json",
    output:
        summary="results/tables/robustness_summary.tsv",
        provenance="results/tables/robustness.provenance.json",
    shell:
        "{PYTHON} scripts/81_collate_robustness.py"


rule validate_integrated_atlas:
    input:
        policy=DOWNSTREAM_POLICY,
        schema=ATLAS_SCHEMA,
        molecular_coverage="results/tables/molecular_locus_coverage.tsv",
        interpretation_policy=INTERPRETATION_POLICY,
        interpretation_coverage="results/tables/interpretation_coverage.tsv",
        interpretation_provenance="results/atlas/interpretation.provenance.json",
        edges_provenance="results/atlas/edges.provenance.json",
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
        summary=rules.robustness.output.summary,
        provenance=rules.robustness.output.provenance,
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
