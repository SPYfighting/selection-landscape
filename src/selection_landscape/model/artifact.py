"""Versioned JSON storage with a checksum over the complete model payload."""

import hashlib
import json
from pathlib import Path

from .core import CalibratedAdditiveModel
from .validation import ModelError, columns


def _digest(payload):
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True, allow_nan=False)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def save_model(model: CalibratedAdditiveModel, path: str | Path) -> None:
    """Write a new artifact. An existing destination is never overwritten."""
    payload = {
        "kind": "calibrated_additive",
        "wild_type": {"sequence": model.wild_type_sequence, "label": model.wild_type_label},
        "epsilon": model.epsilon,
        "mutation_effects": dict(sorted(model.mutation_effects.items())),
        "calibration": {"method": "sigmoid", "intercept": model.intercept, "slope": model.slope},
        "label": {"name": model.label_name, "domain": "unit_interval", "direction": "higher_is_better"},
        "training_pairs": list(model.training_pairs),
    }
    envelope = {"schema_version": 1, "model": payload, "model_sha256": _digest(payload)}
    serialized = json.dumps(envelope, sort_keys=True, indent=2, allow_nan=False) + "\n"
    with Path(path).open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(serialized)


def load_model(path: str | Path) -> CalibratedAdditiveModel:
    """Load only the supported additive schema after checking its checksum."""
    try:
        envelope = json.loads(Path(path).read_text(encoding="utf-8"))
        columns(envelope, {"schema_version", "model", "model_sha256"})
        if envelope["schema_version"] != 1:
            raise ModelError("unsupported model schema_version")
        payload = envelope["model"]
        if _digest(payload) != envelope["model_sha256"]:
            raise ModelError("model checksum mismatch")
        columns(payload, {"kind", "wild_type", "epsilon", "mutation_effects", "calibration", "label", "training_pairs"})
        if payload["kind"] != "calibrated_additive":
            raise ModelError("only calibrated additive models are supported")
        columns(payload["wild_type"], {"sequence", "label"})
        calibration, label = payload["calibration"], payload["label"]
        columns(calibration, {"method", "intercept", "slope"})
        columns(label, {"name", "domain", "direction"})
        if calibration["method"] != "sigmoid":
            raise ModelError("only sigmoid calibration is supported")
        if label["domain"] != "unit_interval" or label["direction"] != "higher_is_better":
            raise ModelError("model labels must use a common higher-is-better unit interval")
        if not isinstance(payload["training_pairs"], list):
            raise ModelError("training_pairs must be a list")
        return CalibratedAdditiveModel(
            wild_type_sequence=payload["wild_type"]["sequence"],
            wild_type_label=payload["wild_type"]["label"], epsilon=payload["epsilon"],
            mutation_effects=payload["mutation_effects"],
            intercept=calibration["intercept"], slope=calibration["slope"],
            label_name=label["name"], training_pairs=tuple(payload["training_pairs"]),
        )
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ModelError(f"invalid model artifact: {exc}") from exc
