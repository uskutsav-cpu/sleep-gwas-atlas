"""Stable, component-scoped policy fingerprints for CATlas cache provenance."""
from __future__ import annotations

import hashlib
import json


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def reference_policy_sha256(policy: dict[str, object]) -> str:
    spec = policy["catlas_adult_v4"]
    keys = (
        "source_id", "component_manifest", "component_manifest_sha256",
        "expected_selected_cells", "selected_cells_by_domain", "cell_type_id_prefix",
        "analysis_reference_prefix", "variant_cache_path", "variant_cache_provenance_path",
        "variant_cache_fields", "minimum_maf", "ld_pruning_r2", "ld_pruning_window_kb",
    )
    causal = policy["causal_inference"]
    return canonical_sha256({
        "catlas_reference": {key: spec[key] for key in keys},
        "build_harmonization": policy["regulatory_build_harmonization"],
        "runtime_manifest": {
            "path": causal["component_manifest"],
            "sha256": causal["component_manifest_sha256"],
        },
    })


def trait_policy_sha256(policy: dict[str, object]) -> str:
    spec = policy["catlas_adult_v4"]
    keys = (
        "gwas_path_template", "variant_cache_path", "variant_cache_provenance_path",
        "trait_cache_path_template", "trait_cache_provenance_path_template",
        "variant_cache_fields", "trait_cache_fields", "minimum_trait_variant_coverage",
    )
    return canonical_sha256({
        "reference_policy_sha256": reference_policy_sha256(policy),
        "catlas_trait": {key: spec[key] for key in keys},
    })
