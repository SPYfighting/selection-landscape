"""Observable behavior of the published NGS-derived selection label."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from selection_landscape.labels.selection_score import calculate_selection_labels


def _write(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_selection_score_cli_uses_six_samples_and_frozen_primary_ranks(tmp_path: Path) -> None:
    counts = tmp_path / "counts.csv"
    rows = []
    for round_id, aa_count, ab_count in ((1, 10, 10), (2, 20, 10), (3, 30, 10)):
        for chain in ("Crick", "Watson"):
            for aa, molecule_count in (("AAA", aa_count), ("ABA", ab_count)):
                rows.append(
                    {"sample_id": f"r{round_id}_{chain}", "round": round_id,
                     "chain": chain, "aa_sequence": aa,
                     "molecule_count": molecule_count, "productive": True}
                )
    _write(counts, list(rows[0]), rows)
    eligibility = tmp_path / "eligibility.csv"
    _write(eligibility, ["aa_sequence", "model_eligible", "sequence_qc_eligible"], [
        {"aa_sequence": aa, "model_eligible": True, "sequence_qc_eligible": True}
        for aa in ("AAA", "ABA")
    ])
    profile = tmp_path / "profile.json"
    profile.write_text(json.dumps({
        "primary_minimum_pooled_molecules": 10,
        "low_weight_minimum_pooled_molecules": 5,
        "low_weight_maximum_pooled_molecules": 9,
        "minimum_molecules_per_chain_cell": 1,
        "transition_weights": {"r1_to_r2": 0.4, "r2_to_r3": 0.6},
        "weak_transition_weight": 0.5,
        "round3_abundance_weight": 0.2,
        "decline_tolerance": 0.7,
        "decline_penalty_strength": 0.2,
    }), encoding="utf-8")
    target = tmp_path / "target.json"
    target.write_text(json.dumps({"target_id": "synthetic", "lead_aa": "AAA"}), encoding="utf-8")
    output = tmp_path / "labels.csv"
    command = [sys.executable, "-m", "selection_landscape.labels", "--counts", str(counts),
               "--eligibility", str(eligibility), "--profile", str(profile),
               "--target", str(target), "--output", str(output)]
    completed = subprocess.run(command, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr
    with output.open(encoding="utf-8", newline="") as handle:
        labels = list(csv.DictReader(handle))
    assert [(row["aa_sequence"], row["rank"]) for row in labels] == [("AAA", "1"), ("ABA", "2")]
    assert float(labels[0]["label"]) == 1.0
    assert abs(float(labels[1]["label"]) - 0.4961904761904762) < 1e-12
    assert all(row["tier"] == "primary" for row in labels)


def test_label_calculation_rejects_mixed_reconstruction_parameters(tmp_path: Path) -> None:
    examples = Path(__file__).parents[1] / "examples" / "ngs"
    with (examples / "synthetic_counts.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    other_rows = [dict(row, param_name="alternative") for row in rows]
    counts = tmp_path / "mixed.csv"
    _write(counts, list(rows[0]), rows + other_rows)
    eligibility = tmp_path / "eligibility.csv"
    _write(eligibility, ["aa_sequence", "model_eligible", "sequence_qc_eligible"], [
        {"aa_sequence": aa, "model_eligible": True, "sequence_qc_eligible": True}
        for aa in sorted({row["aa_sequence"] for row in rows})
    ])
    profile = json.loads((examples / "selection_profile.json").read_text())
    target = json.loads((examples / "target.json").read_text())
    with pytest.raises(ValueError, match="parameter"):
        calculate_selection_labels(counts, eligibility, profile, **target)
