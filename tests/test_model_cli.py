"""The file workflow preserves explicit labels and never overwrites outputs."""

import csv
import json

from selection_landscape.model.cli import main
from test_model_training import synthetic_records


def write_table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_train_enumerate_predict_evaluate_file_workflow(tmp_path, capsys):
    wt, singles, doubles = synthetic_records()
    for name, rows in [("wt", [wt]), ("singles", singles), ("doubles", doubles),
                       ("observed", [{k: r[k] for k in ["mutation_1", "mutation_2"]} for r in doubles])]:
        write_table(tmp_path / f"{name}.csv", rows)
    def path(name):
        return str(tmp_path / name)
    assert main(["train", "--wild-type", path("wt.csv"), "--singles", path("singles.csv"),
                 "--doubles", path("doubles.csv"), "--epsilon", "0.01",
                 "--output", path("model.json")]) == 0
    assert main(["enumerate", "--model", path("model.json"), "--observed-pairs", path("observed.csv"),
                 "--output", path("candidates.csv")]) == 0
    assert main(["predict", "--model", path("model.json"), "--candidates", path("candidates.csv"),
                 "--output", path("predictions.csv")]) == 0
    with (tmp_path / "predictions.csv").open(encoding="utf-8-sig") as handle:
        predictions = list(csv.DictReader(handle))
    assert len(predictions) == 4
    assert all(abs(float(r["predicted_label"]) - .65) < 2e-5 for r in predictions)
    write_table(tmp_path / "labels.csv", [dict(candidate_id=r["candidate_id"], label=.5,
                                              group_id=r["candidate_id"]) for r in predictions])
    assert main(["evaluate", "--predictions", path("predictions.csv"), "--labels", path("labels.csv"),
                 "--bootstrap-replicates", "20", "--seed", "7", "--output", path("metrics.json"),
                 "--bootstrap-output", path("bootstrap.csv")]) == 0
    assert json.loads((tmp_path / "metrics.json").read_text())["n_samples"] == 4
    assert main(["predict", "--model", path("model.json"), "--candidates", path("candidates.csv"),
                 "--output", path("predictions.csv")]) == 2
    assert "already exists" in capsys.readouterr().err


def test_duplicate_csv_headers_are_rejected(tmp_path, capsys):
    (tmp_path / "bad.csv").write_text("sequence,label,label\nACDE,.5,.7\n")
    assert main(["train", "--wild-type", str(tmp_path / "bad.csv"), "--singles", "unused",
                 "--doubles", "unused", "--epsilon", ".01", "--output", str(tmp_path / "out")]) == 2
    assert "header" in capsys.readouterr().err
