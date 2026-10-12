"""Candidate identity exclusion is independent of training-label availability."""

from dataclasses import replace

import pytest

from selection_landscape.model import enumerate_candidates
from test_model_prediction import synthetic_model


def test_enumeration_excludes_all_observed_identities_not_just_training_pairs():
    model = replace(synthetic_model(), training_pairs=("A1V__C2S",))
    observed = [
        {"mutation_1": "A1V", "mutation_2": "C2S"},
        {"mutation_1": "C2S", "mutation_2": "E4Q"},
    ]
    rows = enumerate_candidates(model, observed)
    assert {row["candidate_id"] for row in rows} == {
        "A1V__D3N", "A1V__E4Q", "C2S__D3N", "D3N__E4Q"
    }
    assert next(row["sequence"] for row in rows if row["candidate_id"] == "D3N__E4Q") == "ACNQ"


def test_same_position_substitutions_are_never_combined():
    model = replace(synthetic_model(), mutation_effects={"A1V": .4, "A1T": .2, "C2S": .2})
    assert {row["candidate_id"] for row in enumerate_candidates(model, [])} == {
        "A1T__C2S", "A1V__C2S"
    }


def test_training_identities_cannot_be_omitted_from_exclusion_list():
    model = replace(synthetic_model(), training_pairs=("A1V__C2S",))
    with pytest.raises(ValueError, match="training"):
        enumerate_candidates(model, [])
