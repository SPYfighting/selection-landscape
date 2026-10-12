"""Command-line entry for NGS-derived selection labels."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .selection_score import calculate_selection_labels


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--eligibility", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    profile = json.loads(args.profile.read_text(encoding="utf-8"))
    target = json.loads(args.target.read_text(encoding="utf-8"))
    rows = calculate_selection_labels(
        args.counts, args.eligibility, profile,
        target_id=str(target["target_id"]), lead_aa=str(target["lead_aa"]),
    )
    if args.output.exists():
        raise ValueError(f"output already exists: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({
        "primary": sum(row["tier"] == "primary" for row in rows),
        "low_weight": sum(row["tier"] == "low_weight" for row in rows),
        "output": str(args.output),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
