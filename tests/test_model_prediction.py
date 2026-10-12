"""Public prediction behavior for explicitly synthetic, single-parent variants."""

import pytest

from selection_landscape.model import CalibratedAdditiveModel, predict


def synthetic_model():
    return CalibratedAdditiveModel(
        wild_type_sequence="ACDE",
        wild_type_label=0.5,
        epsilon=0.01,
        mutation_effects={"A1V": 0.4, "C2S": 0.2, "D3N": -0.1, "E4Q": -0.3},
        intercept=0.1,
        slope=0.8,
        label_name="synthetic_score",
    )


def test_frozen_synthetic_parameters_produce_known_scores_and_order():
    candidates = [
        {"mutation_1": "C2S", "mutation_2": "E4Q", "sequence": "ASDQ"},
        {"mutation_1": "A1V", "mutation_2": "D3N", "sequence": "VCNE"},
        {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE"},
    ]
    rows = predict(synthetic_model(), candidates)
    assert [row["candidate_id"] for row in rows] == [
        "A1V__C2S", "A1V__D3N", "C2S__E4Q"
    ]
    assert [row["predicted_label"] for row in rows] == pytest.approx(
        [0.6410674063348171, 0.5841905229354073, 0.5049998333399998], abs=1e-15
    )
    assert [row["rank"] for row in rows] == [1, 2, 3]


@pytest.mark.parametrize("record", [
    {"mutation_1": "C2S", "mutation_2": "A1V", "sequence": "VSDE"},
    {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSNE"},
    {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE", "label": 0.9},
    {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE", "candidate_id": "incorrect"},
    {"mutation_1": "A1V", "mutation_2": "E4V", "sequence": "VCDV"},
])
def test_prediction_rejects_wrong_identity_extra_mutations_and_labels(record):
    with pytest.raises(ValueError):
        predict(synthetic_model(), [record])


def test_duplicate_candidates_are_rejected():
    row = {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE"}
    with pytest.raises(ValueError, match="duplicate"):
        predict(synthetic_model(), [row, row])


def test_numeric_model_parameters_have_one_float_contract():
    model = CalibratedAdditiveModel("ACDE", "0.5", "0.01", {"A1V": ".4", "C2S": ".2"},
                                   "0.1", "0.8")
    assert all(isinstance(v, float) for v in
               [model.wild_type_label, model.epsilon, model.intercept, model.slope])
    assert predict(model, [{"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE"}])[0][
        "predicted_label"] == pytest.approx(.6410674063348171)
