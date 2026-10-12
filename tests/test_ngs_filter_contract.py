import subprocess
import sys
from pathlib import Path


def run_filter(counts: Path, output: Path):
    return subprocess.run([
        sys.executable, "-m", "selection_landscape.ngs.filter", "--counts", str(counts),
        "--parameter-name", "strict", "--length", "11", "--start-motif", "EVQ",
        "--end-motif", "VTVSS", "--output", str(output),
    ], capture_output=True, text=True)


def write_counts(path: Path, parameter="strict", encoding="utf-8"):
    path.write_text(
        "param_name,aa_sequence,round,productive,productive_reason,molecule_count\n"
        + "".join(f"{parameter},EVQACDVTVSS,{r},true,ok,12\n" for r in (1, 2, 3)),
        encoding=encoding,
    )


def test_filter_refuses_to_replace_an_existing_result(tmp_path):
    counts, output = tmp_path / "counts.csv", tmp_path / "result.csv"
    write_counts(counts)
    output.write_text("existing scientific result\n")
    completed = run_filter(counts, output)
    assert completed.returncode != 0
    assert output.read_text() == "existing scientific result\n"


def test_filter_rejects_an_absent_parameter_without_creating_empty_output(tmp_path):
    counts, output = tmp_path / "counts.csv", tmp_path / "result.csv"
    write_counts(counts, parameter="different")
    assert run_filter(counts, output).returncode != 0
    assert not output.exists()


def test_filter_accepts_utf8_bom_and_emits_machine_readable_result(tmp_path):
    import json
    counts, output = tmp_path / "counts.csv", tmp_path / "result.csv"
    write_counts(counts, encoding="utf-8-sig")
    completed = run_filter(counts, output)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["amino_acid_sequences"] == 1
