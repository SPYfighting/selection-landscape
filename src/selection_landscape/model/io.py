"""Lossless floating-point text and strict CSV shape checks."""

import csv
from pathlib import Path

from .validation import ModelError


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames
        if not header or any(not field for field in header) or len(set(header)) != len(header):
            raise ModelError("CSV must have a nonempty, unique header")
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ModelError(f"CSV row {line} does not match the header")
            rows.append(row)
    return rows


def write_csv(path, rows, fields):
    # UTF-8 BOM also allows direct opening in Excel without changing cell text.
    with Path(path).open("x", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
