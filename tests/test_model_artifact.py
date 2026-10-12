"""Portable JSON artifacts preserve scores and reject silent replacement."""

import json

import pytest

from selection_landscape.model import load_model, predict, save_model
from selection_landscape.model import CalibratedAdditiveModel
from test_model_prediction import synthetic_model


def test_json_roundtrip_preserves_prediction_and_rejects_overwrite(tmp_path):
    path = tmp_path / "synthetic-model.json"
    model = synthetic_model()
    candidates = [{"mutation_1": "A1V", "mutation_2": "C2S", "sequence": "VSDE"}]
    save_model(model, path)
    assert predict(load_model(path), candidates) == predict(model, candidates)
    with pytest.raises(FileExistsError):
        save_model(model, path)


def test_changed_artifact_parameters_fail_integrity_check(tmp_path):
    path = tmp_path / "synthetic-model.json"
    save_model(synthetic_model(), path)
    content = json.loads(path.read_text())
    content["model"]["calibration"]["slope"] = 2.0
    path.write_text(json.dumps(content))
    with pytest.raises(ValueError, match="checksum"):
        load_model(path)


@pytest.mark.parametrize("pairs", [("C2S__A1V",), ("A1V__A1V",), ("A1V__C2S", "A1V__C2S")])
def test_training_pair_provenance_requires_canonical_unique_identities(pairs):
    with pytest.raises(ValueError):
        CalibratedAdditiveModel("ACDE", .5, .01, {"A1V": .4, "C2S": .2}, .1, .8, training_pairs=pairs)
