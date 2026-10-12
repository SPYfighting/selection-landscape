from pathlib import Path

import pytest

from selection_landscape.ngs.rebuild import paired_fastq_records


def test_paired_fastq_rejects_truncated_mate(tmp_path: Path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.write_text("@pair/1\nACGT\n+\nIIII\n", encoding="ascii")
    second.write_text("@pair/2\nACGT\n+\n", encoding="ascii")
    with pytest.raises(ValueError, match="truncated or unpaired"):
        list(paired_fastq_records(first, second))


def test_paired_fastq_rejects_mismatched_identifiers(tmp_path: Path):
    first = tmp_path / "first.fastq"
    second = tmp_path / "second.fastq"
    first.write_text("@a/1\nACGT\n+\nIIII\n", encoding="ascii")
    second.write_text("@b/2\nACGT\n+\nIIII\n", encoding="ascii")
    with pytest.raises(ValueError, match="identifiers differ"):
        list(paired_fastq_records(first, second))
