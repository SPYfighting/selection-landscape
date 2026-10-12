"""Evaluation operates on matched identities and resamples position-pair groups."""

import pytest

from selection_landscape.model import evaluate_predictions


def synthetic_evaluation():
    ids = ["A1V__C2S", "A1T__C2S", "A1V__D3N", "A1T__D3N"]
    predictions = [dict(candidate_id=i, predicted_label=p)
                   for i, p in zip(ids, [.9, .7, .4, .2])]
    labels = [dict(candidate_id=i, label=y, group_id=g)
              for i, y, g in zip(ids, [.8, .6, .2, .4], ["1_2", "1_2", "1_3", "1_3"])]
    return predictions, labels


def test_evaluation_known_metrics_and_repeatable_group_bootstrap():
    predictions, labels = synthetic_evaluation()
    result = evaluate_predictions(predictions, labels, bootstrap_replicates=20, seed=7)
    metrics = {r["metric"]: r["value"] for r in result["metrics"]}
    assert metrics["spearman"] == pytest.approx(.8)
    assert metrics["mae"] == pytest.approx(.15)
    assert metrics["top10_percent_ndcg"] == 1.0
    assert result["top_k"] == 1
    assert result["n_groups"] == 2
    assert len(result["bootstrap_metrics"]) == 80
    assert result == evaluate_predictions(predictions, list(reversed(labels)),
                                          bootstrap_replicates=20, seed=7)


def test_evaluation_never_silently_intersects_missing_labels():
    predictions, labels = synthetic_evaluation()
    with pytest.raises(ValueError, match="identities"):
        evaluate_predictions(predictions, labels[:-1])


def test_evaluation_rejects_groups_that_split_a_position_pair():
    predictions, labels = synthetic_evaluation()
    labels[1]["group_id"] = "different"
    with pytest.raises(ValueError, match="position"):
        evaluate_predictions(predictions, labels)


def test_constant_correlations_and_zero_relevance_have_defined_zero_values():
    predictions, labels = synthetic_evaluation()
    for row in labels:
        row["label"] = 0.0
    metrics = {r["metric"]: r["value"] for r in evaluate_predictions(predictions, labels)["metrics"]}
    assert metrics["spearman"] == metrics["pearson"] == metrics["top10_percent_ndcg"] == 0.0
