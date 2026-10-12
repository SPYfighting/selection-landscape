import csv
import subprocess
import sys
from pathlib import Path


def test_filter_counts_cli(tmp_path: Path):
    counts = tmp_path / "counts.csv"
    with counts.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["param_name", "aa_sequence", "round", "productive", "productive_reason", "molecule_count", "read_count"])
        writer.writeheader()
        for aa, rounds in [("EVQACDVTVSS", (2, 3, 4)), ("EVQACDVTVST", (2, 3, 4)), ("EVQACDVTVSS*", (2, 3, 4)), ("EVQCCDVTVSS", (2, 3, 0))]:
            for r, count in enumerate(rounds, 1):
                if count:
                    writer.writerow(dict(param_name="strict", aa_sequence=aa, round=r, productive="true", productive_reason="ok", molecule_count=count, read_count=count))
    output = tmp_path / "eligible.csv"
    result = subprocess.run([sys.executable, "-m", "selection_landscape.ngs.filter", "--counts", str(counts), "--parameter-name", "strict", "--length", "11", "--start-motif", "EVQ", "--end-motif", "VTVSS", "--output", str(output)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    with output.open(newline="") as handle:
        rows = {r["aa_sequence"]: r for r in csv.DictReader(handle)}
    assert rows["EVQACDVTVSS"]["model_eligible"] == "True"
    assert rows["EVQACDVTVST"]["sequence_qc_reason"] == "bad_end"
    assert rows["EVQACDVTVSS*"]["sequence_qc_reason"] == "aa_length"
    assert rows["EVQCCDVTVSS"]["label_reason"] == "missing_round_coverage"
