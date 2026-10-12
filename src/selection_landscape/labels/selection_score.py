"""Calculate a three-round, dual-chain NGS-derived selection score.

The published score is a selection signal. It is neither an affinity
measurement nor a probability of experimental success.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import math
from pathlib import Path
from typing import Any, Iterator, Mapping


ROUNDS = (1, 2, 3)
CHAINS = ("Crick", "Watson")


def _rows(path: Path) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    source_name = path.name[:-3] if path.suffix == ".gz" else path.name
    delimiter = "\t" if source_name.endswith(".tsv") else ","
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def _yes(value: Any) -> bool:
    word = str(value).strip().lower()
    if word in {"true", "1", "yes"}:
        return True
    if word in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean: {value!r}")


def _count(value: Any) -> int:
    number = float(value)
    if not math.isfinite(number) or number < 0 or not number.is_integer():
        raise ValueError(f"molecule_count must be a nonnegative integer: {value!r}")
    return int(number)


def _percentiles(values: list[float]) -> list[float]:
    """Ascending, one-based average-rank percentiles on the primary universe."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    result = [0.0] * n
    start = 0
    while start < n:
        end = start + 1
        while end < n and values[order[end]] == values[order[start]]:
            end += 1
        percentile = ((start + 1 + end) / 2) / n
        for position in range(start, end):
            result[order[position]] = percentile
        start = end
    return result


def _empirical_midrank(value: float, values: list[float]) -> float:
    left = bisect.bisect_left(values, value)
    right = bisect.bisect_right(values, value)
    return (left + right + 1.0) / (2.0 * len(values))


def _validate_profile(profile: Mapping[str, Any]) -> None:
    expected = {
        "primary_minimum_pooled_molecules": 10,
        "low_weight_minimum_pooled_molecules": 5,
        "low_weight_maximum_pooled_molecules": 9,
        "minimum_molecules_per_chain_cell": 1,
        "transition_weights": {"r1_to_r2": 0.4, "r2_to_r3": 0.6},
        "weak_transition_weight": 0.5,
        "round3_abundance_weight": 0.2,
        "decline_tolerance": 0.7,
        "decline_penalty_strength": 0.2,
    }
    for key, value in expected.items():
        if profile.get(key) != value:
            raise ValueError(f"unsupported selection-score profile value: {key}")


def _read_eligibility(path: Path) -> set[str]:
    allowed: set[str] = set()
    seen: set[str] = set()
    for row in _rows(path):
        aa = str(row.get("aa_sequence") or row.get("aa_sequence_1_to_120") or "")
        if not aa or aa in seen:
            raise ValueError("eligibility must contain unique, nonempty amino-acid sequences")
        seen.add(aa)
        if _yes(row["model_eligible"]) and _yes(row["sequence_qc_eligible"]):
            allowed.add(aa)
    return allowed


def _aggregate(path: Path, eligible: set[str]) -> tuple[dict[str, dict[tuple[int, str], int]], dict[tuple[int, str], int]]:
    totals = {(round_id, chain): 0 for round_id in ROUNDS for chain in CHAINS}
    sample_ids: dict[tuple[int, str], str] = {}
    counts: dict[str, dict[tuple[int, str], int]] = {}
    parameter_names: set[str] = set()
    for row in _rows(path):
        parameter_names.add(row.get("param_name", ""))
        if len(parameter_names) > 1:
            raise ValueError("counts contain multiple reconstruction parameter sets; supply one parameter set matching eligibility")
        round_id = _count(row["round"])
        chain = row["chain"]
        if round_id not in ROUNDS or chain not in CHAINS:
            raise ValueError("expected exactly three rounds and Crick/Watson chains")
        key = (round_id, chain)
        sample_id = row.get("sample_id") or f"R{round_id}_{chain}"
        if key in sample_ids and sample_ids[key] != sample_id:
            raise ValueError(f"multiple samples for round/chain {key}")
        sample_ids[key] = sample_id
        if not _yes(row["productive"]):
            continue
        molecule_count = _count(row["molecule_count"])
        totals[key] += molecule_count
        aa = row["aa_sequence"]
        if not aa:
            raise ValueError("productive count has no amino-acid sequence")
        if aa in eligible:
            counts.setdefault(aa, {(r, c): 0 for r in ROUNDS for c in CHAINS})[key] += molecule_count
    if set(sample_ids) != set(totals) or any(value <= 0 for value in totals.values()):
        raise ValueError("all six round/chain samples must have positive productive totals")
    return counts, totals


def _score(row: dict[str, Any], p12: float, p23: float, p3: float, profile: Mapping[str, Any]) -> None:
    weights = profile["transition_weights"]
    transition_mean = float(weights["r1_to_r2"]) * p12 + float(weights["r2_to_r3"]) * p23
    weak_weight = float(profile["weak_transition_weight"])
    trajectory = (1.0 - weak_weight) * transition_mean + weak_weight * min(p12, p23)
    round3_weight = float(profile["round3_abundance_weight"])
    base = (1.0 - round3_weight) * trajectory + round3_weight * p3
    tolerance = float(profile["decline_tolerance"])
    d12 = 0.0 if row["enrich_1_2"] >= tolerance else 1.0 - row["enrich_1_2"] / tolerance
    d23 = 0.0 if row["enrich_2_3"] >= tolerance else 1.0 - row["enrich_2_3"] / tolerance
    penalty = float(weights["r1_to_r2"]) * d12 + float(weights["r2_to_r3"]) * d23
    row.update({
        "percentile_enrich_1_2": p12, "percentile_enrich_2_3": p23,
        "percentile_round3_abundance": p3,
        "minimum_transition_percentile": min(p12, p23),
        "selection_score_base": base, "decline_penalty": penalty,
        "label": max(0.0, min(1.0, base - float(profile["decline_penalty_strength"]) * penalty)),
    })


def calculate_selection_labels(
    counts_path: Path,
    eligibility_path: Path,
    profile: Mapping[str, Any],
    *,
    target_id: str,
    lead_aa: str,
) -> list[dict[str, Any]]:
    """Return ranked primary and low-weight labels from molecule counts.

    Sample normalization includes every productive molecule, including
    sequences that fail ranking eligibility. Only eligible amino-acid sequences
    enter the ranked populations. The model is not imported or called.
    """
    if not target_id or not lead_aa:
        raise ValueError("target_id and lead_aa are required")
    _validate_profile(profile)
    eligible = _read_eligibility(eligibility_path)
    counts, totals = _aggregate(counts_path, eligible)
    primary: list[dict[str, Any]] = []
    low_weight: list[dict[str, Any]] = []
    for aa in sorted(eligible):
        aa_counts = counts.get(aa, {})
        pooled: dict[int, int] = {}
        balanced: dict[int, float] = {}
        row: dict[str, Any] = {"target_id": target_id, "aa_sequence": aa}
        for round_id in ROUNDS:
            crick = aa_counts.get((round_id, "Crick"), 0)
            watson = aa_counts.get((round_id, "Watson"), 0)
            pooled[round_id] = crick + watson
            crick_frequency = crick / totals[(round_id, "Crick")]
            watson_frequency = watson / totals[(round_id, "Watson")]
            balanced[round_id] = math.sqrt(crick_frequency * watson_frequency)
            row[f"molecule_r{round_id}_crick"] = crick
            row[f"molecule_r{round_id}_watson"] = watson
            row[f"molecule_r{round_id}_pooled"] = pooled[round_id]
            row[f"balanced_abundance_r{round_id}"] = balanced[round_id]
        minimum = min(pooled.values())
        if any(aa_counts.get((r, c), 0) < 1 for r in ROUNDS for c in CHAINS):
            continue
        if minimum < int(profile["low_weight_minimum_pooled_molecules"]):
            continue
        e12 = balanced[2] / balanced[1]
        e23 = balanced[3] / balanced[2]
        row.update({
            "enrich_1_2": e12, "enrich_2_3": e23,
            "log2_enrich_1_2": math.log2(e12),
            "log2_enrich_2_3": math.log2(e23),
            "minimum_pooled_molecules": minimum,
            "total_pooled_molecules": sum(pooled.values()),
        })
        if minimum >= int(profile["primary_minimum_pooled_molecules"]):
            row.update({"tier": "primary", "weight": 1.0})
            primary.append(row)
        else:
            row.update({"tier": "low_weight", "weight": minimum / 10})
            low_weight.append(row)
    if not primary:
        raise ValueError("no primary labels satisfy the frozen molecule support rule")
    for field, destination in (
        ("log2_enrich_1_2", "p12"),
        ("log2_enrich_2_3", "p23"),
        ("balanced_abundance_r3", "p3"),
    ):
        values = [float(row[field]) for row in primary]
        ranks = _percentiles(values)
        for row, rank in zip(primary, ranks):
            row[destination] = rank
        ordered = sorted(values)
        for row in low_weight:
            row[destination] = _empirical_midrank(float(row[field]), ordered)
    for row in primary + low_weight:
        _score(row, row.pop("p12"), row.pop("p23"), row.pop("p3"), profile)
    def sort_key(row: Mapping[str, Any]) -> tuple[float, float, int, int, str]:
        return (-float(row["label"]), -float(row["minimum_transition_percentile"]),
                -int(row["minimum_pooled_molecules"]), -int(row["total_pooled_molecules"]),
                str(row["aa_sequence"]))
    primary.sort(key=sort_key)
    low_weight.sort(key=sort_key)
    for population in (primary, low_weight):
        for index, row in enumerate(population, start=1):
            row["rank"] = index
    return primary + low_weight
