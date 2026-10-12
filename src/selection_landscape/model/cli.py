"""Command-line training, candidate generation, prediction and evaluation."""

import argparse
import json
from pathlib import Path
import sys

from . import (TrainingData, enumerate_candidates, evaluate_predictions, fit_model,
               load_model, predict, save_model)
from .io import read_csv, write_csv
from .validation import ModelError

_CANDIDATE_FIELDS = ["candidate_id", "mutation_1", "mutation_2", "sequence"]


def _parser():
    parser = argparse.ArgumentParser(description="Calibrated additive prediction from bounded labels.")
    commands = parser.add_subparsers(dest="command", required=True)
    train = commands.add_parser("train", help="Fit sigmoid calibration from measured singles and doubles")
    train.add_argument("--wild-type", required=True, help="One-row CSV: sequence,label")
    train.add_argument("--singles", required=True, help="CSV: mutation,sequence,label")
    train.add_argument("--doubles", required=True, help="CSV: mutation_1,mutation_2,sequence,label[,weight]")
    train.add_argument("--epsilon", required=True, type=float, help="Fixed logit-clipping epsilon, strictly between 0 and 0.5")
    train.add_argument("--label-name", default="label")
    predict_parser = commands.add_parser("predict", help="Predict supplied double-mutant candidates")
    predict_parser.add_argument("--model", required=True)
    predict_parser.add_argument("--candidates", required=True)
    enumerate_parser = commands.add_parser("enumerate", help="Enumerate unseen pairs of measured single mutations")
    enumerate_parser.add_argument("--model", required=True)
    enumerate_parser.add_argument("--observed-pairs", required=True,
                                  help="Complete observed identity CSV, including held-out pairs; no labels")
    evaluate_parser = commands.add_parser("evaluate", help="Evaluate predictions with position-pair bootstrap intervals")
    evaluate_parser.add_argument("--predictions", required=True)
    evaluate_parser.add_argument("--labels", required=True, help="CSV: candidate_id,label,group_id")
    evaluate_parser.add_argument("--bootstrap-replicates", type=int, default=0)
    evaluate_parser.add_argument("--seed", type=int, default=20260825)
    evaluate_parser.add_argument("--bootstrap-output", help="Optional CSV containing each resampled metric")
    for command in [train, predict_parser, enumerate_parser, evaluate_parser]:
        command.add_argument("--output", required=True, help="New output path; existing files are never replaced")
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        outputs = [args.output]
        if args.command == "evaluate" and args.bootstrap_output:
            outputs.append(args.bootstrap_output)
        if len({Path(path).resolve() for path in outputs}) != len(outputs):
            raise ModelError("output paths must be different")
        for path in outputs:
            if Path(path).exists():
                raise ModelError(f"output already exists: {path}")
            if not Path(path).parent.is_dir():
                raise ModelError(f"output parent directory does not exist: {Path(path).parent}")
        if args.command == "train":
            wild_type = read_csv(args.wild_type)
            if len(wild_type) != 1:
                raise ModelError("wild-type CSV must contain exactly one record")
            data = TrainingData.from_records(wild_type[0], read_csv(args.singles), read_csv(args.doubles))
            model = fit_model(data, epsilon=args.epsilon, label_name=args.label_name)
            save_model(model, args.output)
            summary = {"single_mutants": len(model.mutation_effects), "training_doubles": len(model.training_pairs)}
        elif args.command == "predict":
            rows = predict(load_model(args.model), read_csv(args.candidates))
            write_csv(args.output, rows, _CANDIDATE_FIELDS + ["predicted_label", "rank"])
            summary = {"predictions": len(rows)}
        elif args.command == "enumerate":
            rows = enumerate_candidates(load_model(args.model), read_csv(args.observed_pairs))
            write_csv(args.output, rows, _CANDIDATE_FIELDS)
            summary = {"candidates": len(rows)}
        else:
            result = evaluate_predictions(read_csv(args.predictions), read_csv(args.labels),
                                          bootstrap_replicates=args.bootstrap_replicates, seed=args.seed)
            bootstrap = result.pop("bootstrap_metrics")
            with Path(args.output).open("x", encoding="utf-8") as handle:
                json.dump(result, handle, indent=2, allow_nan=False)
                handle.write("\n")
            if args.bootstrap_output:
                write_csv(args.bootstrap_output, bootstrap, ["replicate_id", "metric", "value"])
            summary = {"n_samples": result["n_samples"], "bootstrap_replicates": result["bootstrap_replicates"]}
        print(json.dumps({"output": str(args.output), **summary}, allow_nan=False))
        return 0
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
