"""Label-scale evaluation and bootstrap resampling of mutation-position pairs."""

import math
import re

import numpy as np
from scipy.stats import pearsonr, spearmanr

from .validation import ModelError, bounded_label, columns


def _position_pair(identifier):
    match = re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]([1-9][0-9]*)[ACDEFGHIKLMNPQRSTVWY]__"
                         r"[ACDEFGHIKLMNPQRSTVWY]([1-9][0-9]*)[ACDEFGHIKLMNPQRSTVWY]", identifier)
    if match is None or int(match[1]) >= int(match[2]):
        raise ModelError("candidate_id must encode a canonical double mutation")
    return int(match[1]), int(match[2])


def _metrics(target, prediction, identifiers):
    k = max(1, math.ceil(.10 * len(target)))
    discount = 1.0 / np.log2(np.arange(2, k + 2))
    predicted_order = np.lexsort((identifiers, -prediction))[:k]
    ideal_order = np.lexsort((identifiers, -target))[:k]
    ideal = float(np.sum(target[ideal_order] * discount))
    ndcg = float(np.sum(target[predicted_order] * discount)) / ideal if ideal > 0 else 0.0
    constant = np.ptp(target) == 0 or np.ptp(prediction) == 0
    return {
        "spearman": 0.0 if constant else float(spearmanr(target, prediction)[0]),
        "top10_percent_ndcg": ndcg,
        "mae": float(np.mean(np.abs(target - prediction))),
        "pearson": 0.0 if constant else float(pearsonr(target, prediction)[0]),
    }, k


def evaluate_predictions(predictions, labels, *, bootstrap_replicates=0, seed=20260825):
    """Evaluate matched predictions and labels; optionally resample position pairs.

    Predictions require candidate_id and predicted_label. Labels require
    candidate_id, label and group_id. Each group must correspond to exactly one
    mutation-position pair. No unmatched rows are dropped.
    """
    if isinstance(bootstrap_replicates, bool) or not isinstance(bootstrap_replicates, int) or bootstrap_replicates < 0:
        raise ModelError("bootstrap_replicates must be a nonnegative integer")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ModelError("seed must be a nonnegative integer")
    predicted = {}
    for row in predictions:
        columns(row, {"candidate_id", "predicted_label"}, {"mutation_1", "mutation_2", "sequence", "rank"})
        identifier = str(row["candidate_id"])
        _position_pair(identifier)
        if identifier in predicted:
            raise ModelError("duplicate prediction identity")
        predicted[identifier] = bounded_label(row["predicted_label"])
    observed, pair_to_group, group_to_pair = {}, {}, {}
    for row in labels:
        columns(row, {"candidate_id", "label", "group_id"})
        identifier = str(row["candidate_id"])
        position_pair = _position_pair(identifier)
        group = row["group_id"]
        if not isinstance(group, str) or not group.strip():
            raise ModelError("group_id must be a nonempty string")
        if identifier in observed:
            raise ModelError("duplicate label identity")
        if (position_pair in pair_to_group and pair_to_group[position_pair] != group) or (
                group in group_to_pair and group_to_pair[group] != position_pair):
            raise ModelError("group_id must uniquely identify a mutation-position pair")
        pair_to_group[position_pair], group_to_pair[group] = group, position_pair
        observed[identifier] = (bounded_label(row["label"]), group)
    if set(predicted) != set(observed):
        raise ModelError("prediction and label identities must match exactly")
    if len(predicted) < 2:
        raise ModelError("evaluation requires at least two predictions")
    # Preserve prediction row order, while matching labels by identity.
    ids = np.asarray(list(predicted), dtype=str)
    target = np.asarray([observed[i][0] for i in ids], dtype=float)
    scores = np.asarray([predicted[i] for i in ids], dtype=float)
    metrics, k = _metrics(target, scores, ids)
    groups = sorted(group_to_pair)
    by_group = {g: np.asarray([n for n, i in enumerate(ids) if observed[i][1] == g], dtype=int)
                for g in groups}
    rng = np.random.default_rng(seed)
    bootstrap = []
    for replicate in range(bootstrap_replicates):
        draws = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_group[g] for g in draws])
        values, _ = _metrics(target[indices], scores[indices], ids[indices])
        bootstrap.extend({"replicate_id": replicate, "metric": name, "value": value}
                         for name, value in values.items())
    summary = []
    for name, value in metrics.items():
        ci = [None, None]
        if bootstrap_replicates:
            ci = np.quantile([r["value"] for r in bootstrap if r["metric"] == name],
                             [.025, .975], method="linear").tolist()
        summary.append({"metric": name, "value": value, "ci_low": ci[0], "ci_high": ci[1]})
    return {"n_samples": len(ids), "top_k": k, "n_groups": len(groups),
            "bootstrap_replicates": bootstrap_replicates, "seed": seed,
            "metrics": summary, "bootstrap_metrics": bootstrap}
