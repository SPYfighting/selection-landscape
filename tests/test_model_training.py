"""Training labels are caller-supplied measurements, never implicit NGS scores."""

import pytest

from selection_landscape.model import TrainingData, fit_model, predict


def synthetic_records():
    wild_type = {"sequence": "ACDE", "label": 0.5}
    singles = [
        {"mutation": "A1V", "sequence": "VCDE", "label": 0.5},
        {"mutation": "C2S", "sequence": "ASDE", "label": 0.5},
        {"mutation": "D3N", "sequence": "ACNE", "label": 0.5},
        {"mutation": "E4Q", "sequence": "ACDQ", "label": 0.5},
    ]
    doubles = [
        {"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE", "label": 0.2, "weight": 1},
        {"mutation_1": "D3N", "mutation_2": "E4Q", "sequence": "ACNQ", "label": 0.8, "weight": 3},
    ]
    return wild_type, singles, doubles


def test_generic_labels_and_explicit_weights_fit_the_known_weighted_mean():
    data = TrainingData.from_records(*synthetic_records())
    model = fit_model(data, epsilon=0.01, label_name="synthetic_activity")
    result = predict(model, [{"mutation_1": "A1V", "mutation_2": "D3N", "sequence": "VCNE"}])
    # All input single-mutant effects are zero. The independently known least-
    # squares optimum of the two weighted observations is (0.2 + 3*0.8)/4 = 0.65.
    assert result[0]["predicted_label"] == pytest.approx(0.65, abs=2e-5)
    assert model.label_name == "synthetic_activity"
    assert model.training_pairs == ("A1V__C2S", "D3N__E4Q")


@pytest.mark.parametrize("invalid", [-0.01, 1.01, float("nan"), float("inf")])
def test_invalid_labels_are_rejected_without_normalizing(invalid):
    wt, singles, doubles = synthetic_records()
    singles[0]["label"] = invalid
    with pytest.raises(ValueError, match="label"):
        TrainingData.from_records(wt, singles, doubles)


def test_undeclared_weak_tier_cannot_silently_enter_training():
    wt, singles, doubles = synthetic_records()
    doubles[0]["tier"] = "weak"
    with pytest.raises(ValueError, match="unexpected"):
        TrainingData.from_records(wt, singles, doubles)


@pytest.mark.parametrize("weight", [0, -1, float("nan"), float("inf")])
def test_double_weights_must_be_positive_finite_values(weight):
    wt, singles, doubles = synthetic_records()
    doubles[0]["weight"] = weight
    with pytest.raises(ValueError, match="weight"):
        TrainingData.from_records(wt, singles, doubles)


@pytest.mark.parametrize("epsilon", [0, 0.5, -0.1, 1e-20, float("nan")])
def test_clipping_precision_must_be_explicit_and_valid(epsilon):
    with pytest.raises(ValueError, match="epsilon"):
        fit_model(TrainingData.from_records(*synthetic_records()), epsilon=epsilon)


def test_only_a_canonical_single_change_matches_a_single_mutant_record():
    wt, singles, doubles = synthetic_records()
    singles[0]["sequence"] = "VSDE"
    with pytest.raises(ValueError, match="single|sequence"):
        TrainingData.from_records(wt, singles, doubles)
