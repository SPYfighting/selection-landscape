"""Calibrated additive prediction from observed single-mutant labels."""

from dataclasses import dataclass
from typing import Iterable, Mapping

import numpy as np
from scipy.optimize import minimize

from .validation import (ModelError, bounded_label, finite_number, mutation_identity, pair_sequence,
                         valid_epsilon, valid_sequence, validated_pairs, validated_training)


@dataclass(frozen=True)
class TrainingData:
    """WT, measured singles and explicitly selected training doubles."""

    wild_type: Mapping[str, object]
    single_mutants: tuple[Mapping[str, object], ...]
    observed_doubles: tuple[Mapping[str, object], ...]

    def __post_init__(self):
        wt, singles, doubles = validated_training(self.wild_type, self.single_mutants,
                                                 self.observed_doubles)
        object.__setattr__(self, "wild_type", wt)
        object.__setattr__(self, "single_mutants", singles)
        object.__setattr__(self, "observed_doubles", doubles)

    @classmethod
    def from_records(cls, wild_type, single_mutants, observed_doubles):
        return cls(dict(wild_type), tuple(dict(row) for row in single_mutants),
                   tuple(dict(row) for row in observed_doubles))


@dataclass(frozen=True)
class CalibratedAdditiveModel:
    """A fixed single-parent model; labels have a common 0–1 scale."""

    wild_type_sequence: str
    wild_type_label: float
    epsilon: float
    mutation_effects: Mapping[str, float]
    intercept: float
    slope: float
    label_name: str = "label"
    training_pairs: tuple[str, ...] = ()

    def __post_init__(self):
        valid_sequence(self.wild_type_sequence)
        object.__setattr__(self, "wild_type_label", bounded_label(self.wild_type_label))
        object.__setattr__(self, "epsilon", valid_epsilon(self.epsilon))
        object.__setattr__(self, "intercept", finite_number(self.intercept, "intercept"))
        object.__setattr__(self, "slope", finite_number(self.slope, "slope"))
        if self.slope < 0:
            raise ModelError("slope must be nonnegative")
        if not isinstance(self.label_name, str) or not self.label_name.strip():
            raise ModelError("label_name must be a nonempty description")
        if not self.mutation_effects:
            raise ModelError("model must contain measured single-mutant effects")
        effects = {}
        for mutation, value in self.mutation_effects.items():
            mutation_identity(mutation, self.wild_type_sequence)
            effects[mutation] = finite_number(value, "mutation effect")
        object.__setattr__(self, "mutation_effects", effects)
        pairs = tuple(self.training_pairs)
        if len(set(pairs)) != len(pairs):
            raise ModelError("duplicate training pair identity")
        for identifier in pairs:
            parts = str(identifier).split("__")
            if len(parts) != 2:
                raise ModelError("training pair identity must contain two mutation IDs")
            pair_sequence(parts[0], parts[1], self.wild_type_sequence, effects)
        object.__setattr__(self, "training_pairs", pairs)


def _sigmoid(values: np.ndarray) -> np.ndarray:
    positive = values >= 0
    output = np.empty_like(values, dtype=float)
    output[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exponential = np.exp(values[~positive])
    output[~positive] = exponential / (1.0 + exponential)
    return output


def fit_model(
    data: TrainingData, *, epsilon: float, label_name: str = "label"
) -> CalibratedAdditiveModel:
    """Fit two sigmoid calibration parameters using weighted squared error."""
    epsilon = valid_epsilon(epsilon)
    # Recheck caller-owned records at the operation boundary.
    data = TrainingData.from_records(data.wild_type, data.single_mutants, data.observed_doubles)
    wild = float(np.clip(float(data.wild_type["label"]), epsilon, 1.0 - epsilon))
    labels = np.clip(np.asarray([float(r["label"]) for r in data.single_mutants]),
                     epsilon, 1.0 - epsilon)
    wild_array = np.array([wild])
    effects = np.log(labels / (1.0 - labels)) - np.log(wild_array / (1.0 - wild_array))[0]
    effect_map = dict(sorted((str(row["mutation"]), float(value))
                            for row, value in zip(data.single_mutants, effects)))
    x = np.asarray([effect_map[str(r["mutation_1"])] + effect_map[str(r["mutation_2"])]
                    for r in data.observed_doubles], dtype=float)
    y = np.asarray([float(r["label"]) for r in data.observed_doubles], dtype=float)
    weights = np.asarray([float(r.get("weight", 1.0)) for r in data.observed_doubles], dtype=float)

    def objective(parameters):
        prediction = _sigmoid(parameters[0] + parameters[1] * x)
        return float(np.average((prediction - y) ** 2, weights=weights))

    result = minimize(objective, np.array([0.0, 1.0]), method="L-BFGS-B",
                      bounds=[(None, None), (0.0, None)])
    if not result.success:
        raise ModelError(f"calibration failed: {result.message}")
    return CalibratedAdditiveModel(
        wild_type_sequence=str(data.wild_type["sequence"]),
        wild_type_label=float(data.wild_type["label"]), epsilon=float(epsilon),
        mutation_effects=effect_map, intercept=float(result.x[0]), slope=float(result.x[1]),
        label_name=label_name,
        training_pairs=tuple(sorted(f"{r['mutation_1']}__{r['mutation_2']}"
                                    for r in data.observed_doubles)),
    )


def predict(
    model: CalibratedAdditiveModel, candidates: Iterable[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Score double mutants; break exact score ties by canonical candidate ID."""
    rows = []
    latent = []
    for record in validated_pairs(candidates, model.wild_type_sequence, model.mutation_effects):
        first, second = str(record["mutation_1"]), str(record["mutation_2"])
        rows.append({
            "candidate_id": f"{first}__{second}",
            "mutation_1": first,
            "mutation_2": second,
            "sequence": str(record["sequence"]),
        })
        latent.append(model.mutation_effects[first] + model.mutation_effects[second])
    scores = _sigmoid(model.intercept + model.slope * np.asarray(latent, dtype=float))
    for row, score in zip(rows, scores):
        row["predicted_label"] = float(score)
    rows.sort(key=lambda row: (-row["predicted_label"], row["candidate_id"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    return rows
