from validate_current_outputs import profile_loo_tables_equivalent


def test_profile_loo_comparison_allows_only_serialization_scale_float_drift():
    frozen = [{
        "claim_id": "profile_a_b",
        "omitted_sleep_trait": "sleepiness",
        "n_sleep_traits": "11",
        "leave_one_out_r": "0.967280655333",
    }]
    rerendered = [{**frozen[0], "leave_one_out_r": "0.96728065533297"}]
    assert profile_loo_tables_equivalent(frozen, rerendered)


def test_profile_loo_comparison_rejects_material_or_structural_changes():
    frozen = [{"claim_id": "profile_a_b", "leave_one_out_r": "0.967280655333"}]
    changed_value = [{"claim_id": "profile_a_b", "leave_one_out_r": "0.9672806553"}]
    changed_label = [{"claim_id": "profile_a_c", "leave_one_out_r": "0.967280655333"}]
    assert not profile_loo_tables_equivalent(frozen, changed_value)
    assert not profile_loo_tables_equivalent(frozen, changed_label)
