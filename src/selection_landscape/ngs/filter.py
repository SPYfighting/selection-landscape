"""Derive clone eligibility from UMI molecule counts, independently of model training."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path


ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    name = path.name[:-3] if path.suffix == ".gz" else path.name
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t" if name.endswith(".tsv") else ",")


def _integer(raw: str) -> int:
    value = float(raw)
    if not math.isfinite(value) or value < 0 or not value.is_integer():
        raise ValueError(f"invalid count: {raw!r}")
    return int(value)


def _yes(raw: str) -> bool:
    value = str(raw).lower()
    if value in {"true", "1", "yes"}:
        return True
    if value in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean: {raw!r}")


def filter_counts(counts_path: Path, output_path: Path, *, parameter_name: str, length: int, start_motif: str, end_motif: str) -> dict[str, int]:
    """Write one eligibility row per amino-acid sequence for a three-round assay."""
    if output_path.exists():
        raise FileExistsError(f"output already exists: {output_path}")
    if length <= 0 or not start_motif or not end_motif:
        raise ValueError("positive length and both sequence motifs are required")
    clones: dict[str, dict] = {}
    n_rows = 0
    for row in _rows(counts_path):
        if row.get("param_name") != parameter_name:
            continue
        aa = row.get("aa_sequence", "")
        if not aa:
            continue
        round_id = _integer(row["round"])
        if round_id not in (1, 2, 3):
            raise ValueError("round must be 1, 2 or 3")
        clone = clones.setdefault(aa, {"rounds": defaultdict(int), "productive": False, "reason": row.get("productive_reason", "")})
        count = _integer(row["molecule_count"])
        clone["rounds"][round_id] += count
        productive = _yes(row["productive"])
        if productive and not clone["productive"]:
            clone["productive"] = True
            clone["reason"] = row.get("productive_reason", "")
        n_rows += 1
    if n_rows == 0:
        raise ValueError(f"no count rows for parameter set {parameter_name!r}")
    fields = ["aa_sequence", "aa_length", "productive", "sequence_qc_eligible", "sequence_qc_reason", "label_eligible", "label_reason", "forward_eligible_r2", "forward_reason", "model_eligible", "model_reason", "mol_r1", "mol_r2", "mol_r3"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_sequence = n_label = n_forward = 0
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for aa in sorted(clones):
            clone = clones[aa]
            if not clone["productive"]:
                reason = clone["reason"] or "not_productive"
            elif len(aa) != length:
                reason = "aa_length"
            elif "*" in aa:
                reason = "stop_codon"
            elif any(residue not in ALPHABET for residue in aa):
                reason = "unresolved_amino_acid"
            elif not aa.startswith(start_motif):
                reason = "bad_start"
            elif not aa.endswith(end_motif):
                reason = "bad_end"
            else:
                reason = "ok"
            seq_ok = reason == "ok"
            rounds = clone["rounds"]
            label_ok = seq_ok and all(rounds[r] >= 1 for r in (1, 2, 3))
            forward_ok = seq_ok and all(rounds[r] >= 1 for r in (1, 2))
            label_reason = "parent_false" if not seq_ok else ("ok" if label_ok else "missing_round_coverage")
            forward_reason = "parent_false" if not seq_ok else ("ok" if forward_ok else "missing_r1_r2_coverage")
            writer.writerow({"aa_sequence": aa, "aa_length": len(aa), "productive": clone["productive"], "sequence_qc_eligible": seq_ok, "sequence_qc_reason": reason, "label_eligible": label_ok, "label_reason": label_reason, "forward_eligible_r2": forward_ok, "forward_reason": forward_reason, "model_eligible": label_ok, "model_reason": label_reason, "mol_r1": rounds[1], "mol_r2": rounds[2], "mol_r3": rounds[3]})
            n_sequence += seq_ok
            n_label += label_ok
            n_forward += forward_ok
    return {"input_rows": n_rows, "amino_acid_sequences": len(clones), "sequence_qc_eligible": n_sequence, "label_eligible": n_label, "forward_eligible_r2": n_forward, "model_eligible": n_label}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--parameter-name", required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--start-motif", required=True)
    parser.add_argument("--end-motif", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = filter_counts(args.counts, args.output, parameter_name=args.parameter_name, length=args.length, start_motif=args.start_motif, end_motif=args.end_motif)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
