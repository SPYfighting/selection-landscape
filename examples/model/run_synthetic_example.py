"""Run a synthetic file workflow against an explicitly defined generating model."""

import argparse
import json
import math
from pathlib import Path

from selection_landscape.model.cli import main
from selection_landscape.model.io import read_csv


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True, help="A new output directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    def src(name):
        return str(source / name)
    def out(name):
        return str(args.output / name)
    commands = [
        ["train", "--wild-type", src("synthetic_wild_type.csv"), "--singles", src("synthetic_single_mutants.csv"),
         "--doubles", src("synthetic_training_doubles.csv"), "--epsilon", "0.01", "--label-name", "synthetic_score",
         "--output", out("synthetic_model.json")],
        ["enumerate", "--model", out("synthetic_model.json"), "--observed-pairs", src("synthetic_all_observed_pairs.csv"),
         "--output", out("synthetic_candidates.csv")],
        ["predict", "--model", out("synthetic_model.json"), "--candidates", out("synthetic_candidates.csv"),
         "--output", out("synthetic_predictions.csv")],
        ["evaluate", "--predictions", out("synthetic_predictions.csv"), "--labels", src("synthetic_evaluation_labels.csv"),
         "--bootstrap-replicates", "2000", "--seed", "20260825", "--output", out("synthetic_metrics.json"),
         "--bootstrap-output", out("synthetic_bootstrap_metrics.csv")],
    ]
    for command in commands:
        if main(command):
            raise SystemExit(2)
    predictions = read_csv(args.output / "synthetic_predictions.csv")
    definition = json.loads((source / "synthetic_definition.json").read_text())
    assert len(predictions) == 4
    assert [r["candidate_id"] for r in predictions] == definition["expected_candidate_order"]
    effects = definition["single_mutant_effects"]
    for row in predictions:
        latent = effects[row["mutation_1"]] + effects[row["mutation_2"]]
        known_score = 1.0 / (1.0 + math.exp(-(definition["intercept"] + definition["slope"] * latent)))
        assert abs(float(row["predicted_label"]) - known_score) < definition["fitted_prediction_absolute_tolerance"]


if __name__ == "__main__":
    run()
