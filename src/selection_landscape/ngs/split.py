"""Competitive paired-read assignment to a primary or competing reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterator


def _single_fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    records = [line for line in lines if line.startswith(">")]
    if len(records) != 1:
        raise ValueError("each reference must contain exactly one sequence")
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    if not sequence:
        raise ValueError("reference sequence is empty")
    return sequence.upper()


def _read_pair(first: Any, second: Any) -> tuple[list[str], list[str]] | None:
    r1 = [first.readline() for _ in range(4)]
    r2 = [second.readline() for _ in range(4)]
    if not r1[0] and not r2[0]:
        return None
    if any(not line for line in r1 + r2):
        raise ValueError("paired FASTQ is truncated")
    if not r1[0].startswith("@") or not r2[0].startswith("@") or not r1[2].startswith("+") or not r2[2].startswith("+"):
        raise ValueError("invalid FASTQ record")
    if len(r1[1].rstrip()) != len(r1[3].rstrip()) or len(r2[1].rstrip()) != len(r2[3].rstrip()):
        raise ValueError("sequence and quality lengths differ")
    return r1, r2


def _qname(header: str) -> str:
    name = header[1:].split()[0]
    return name.removesuffix("/1").removesuffix("/2")


def _sam_groups(path: Path) -> Iterator[tuple[str, list[str]]]:
    with path.open(encoding="utf-8") as handle:
        current_name = ""
        group: list[str] = []
        for line in handle:
            if line.startswith("@"):
                continue
            name = line.split("\t", 1)[0]
            if group and name != current_name:
                yield current_name, group
                group = []
            current_name = name
            group.append(line)
        if group:
            yield current_name, group


def classify_pair(scores: dict[str, int | None], *, min_pair_score: int, min_score_margin: int) -> str:
    """Classify paired alignment scores using the frozen score and margin rule."""
    usable = {name: score for name, score in scores.items() if score is not None}
    if not usable:
        return "unassigned"
    ordered = sorted(usable.items(), key=lambda item: (-int(item[1]), item[0]))
    best_name, best_score = ordered[0]
    if int(best_score) < min_pair_score:
        return "unassigned"
    if len(ordered) > 1 and int(best_score) - int(ordered[1][1]) < min_score_margin:
        return "ambiguous"
    return f"{best_name}-only"


def _paired_scores(records: list[str]) -> dict[str, int | None]:
    per_mate: dict[str, dict[str, int | None]] = {
        "R1": {"primary": None, "competitor": None},
        "R2": {"primary": None, "competitor": None},
    }
    for line in records:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            raise ValueError("malformed SAM alignment")
        flag = int(fields[1])
        reference = fields[2]
        if flag & 4 or reference not in ("primary", "competitor"):
            continue
        mate = "R1" if flag & 64 else ("R2" if flag & 128 else "")
        if not mate:
            continue
        score = next((int(field[5:]) for field in fields[11:] if field.startswith("AS:i:")), None)
        if score is None:
            continue
        old = per_mate[mate][reference]
        per_mate[mate][reference] = score if old is None else max(old, score)
    result: dict[str, int | None] = {}
    for reference in ("primary", "competitor"):
        values = [per_mate[mate][reference] for mate in ("R1", "R2") if per_mate[mate][reference] is not None]
        result[reference] = sum(values) if values else None
    return result


def split_mixed_fastq(
    sample_manifest: Path,
    primary_reference: Path,
    competitor_reference: Path,
    output_dir: Path,
    *,
    bowtie2: str = "bowtie2",
    bowtie2_build: str = "bowtie2-build",
    min_pair_score: int = 60,
    min_score_margin: int = 12,
    max_pairs: int | None = None,
    threads: int = 1,
) -> list[dict[str, Any]]:
    """Write target-assigned FASTQ while conserving every input pair."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("output directory must be empty")
    output_dir.mkdir(parents=True, exist_ok=True)
    combined = output_dir / "competitive_reference.fa"
    combined.write_text(
        f">primary\n{_single_fasta(primary_reference)}\n"
        f">competitor\n{_single_fasta(competitor_reference)}\n",
        encoding="utf-8",
    )
    index = output_dir / "competitive_reference"
    subprocess.run([bowtie2_build, str(combined), str(index)], check=True,
                   stdout=(output_dir / "index.stdout.log").open("w"),
                   stderr=(output_dir / "index.stderr.log").open("w"))
    with sample_manifest.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle))
    if not samples:
        raise ValueError("sample manifest is empty")
    identifiers = [sample["sample_id"] for sample in samples]
    if len(set(identifiers)) != len(identifiers) or any(not item or item in {".", ".."} or "/" in item or "\\" in item for item in identifiers):
        raise ValueError("sample identifiers must be unique single path components")
    results: list[dict[str, Any]] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        r1 = Path(sample["r1_path"])
        r2 = Path(sample["r2_path"])
        if not r1.is_absolute():
            r1 = sample_manifest.parent / r1
        if not r2.is_absolute():
            r2 = sample_manifest.parent / r2
        sample_dir = output_dir / sample_id
        sample_dir.mkdir()
        sam = sample_dir / "alignments.sam"
        command = [bowtie2, "--very-sensitive-local", "--reorder", "--no-mixed", "--no-discordant",
                   "-p", str(threads), "-k", "2", "-x", str(index), "-1", str(r1), "-2", str(r2), "-S", str(sam)]
        if max_pairs is not None:
            command.extend(["-u", str(max_pairs)])
        with (sample_dir / "bowtie2.stdout.log").open("w") as stdout, (sample_dir / "bowtie2.stderr.log").open("w") as stderr:
            subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        labels = ("primary-only", "competitor-only", "ambiguous", "unassigned")
        handles = {
            label: (
                gzip.open(sample_dir / f"{label}_R1.fq.gz", "wt", encoding="utf-8"),
                gzip.open(sample_dir / f"{label}_R2.fq.gz", "wt", encoding="utf-8"),
            )
            for label in labels
        }
        counts: Counter[str] = Counter()
        try:
            with gzip.open(r1, "rt", encoding="utf-8") as first, gzip.open(r2, "rt", encoding="utf-8") as second:
                groups = iter(_sam_groups(sam))
                while max_pairs is None or sum(counts.values()) < max_pairs:
                    pair = _read_pair(first, second)
                    if pair is None:
                        break
                    first_record, second_record = pair
                    name = _qname(first_record[0])
                    if _qname(second_record[0]) != name:
                        raise ValueError("FASTQ mate names differ")
                    group_name, records = next(groups)
                    if group_name != name:
                        raise ValueError(f"SAM and FASTQ order differs: {group_name} != {name}")
                    label = classify_pair(_paired_scores(records), min_pair_score=min_pair_score,
                                          min_score_margin=min_score_margin)
                    handles[label][0].writelines(first_record)
                    handles[label][1].writelines(second_record)
                    counts[label] += 1
        finally:
            for pair_handles in handles.values():
                pair_handles[0].close()
                pair_handles[1].close()
        results.append({"sample_id": sample_id, "round": sample["round"], "chain": sample["chain"],
                        "pairs": sum(counts.values()), **{key: counts[key] for key in labels}})
    report = output_dir / "split_report.csv"
    with report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    (output_dir / "split_config.json").write_text(json.dumps({
        "min_pair_score": min_pair_score, "min_score_margin": min_score_margin,
        "max_pairs": max_pairs, "threads": threads,
    }, indent=2) + "\n", encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--primary-reference", type=Path, required=True)
    parser.add_argument("--competitor-reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bowtie2", default="bowtie2")
    parser.add_argument("--bowtie2-build", default="bowtie2-build")
    parser.add_argument("--min-pair-score", type=int, default=60)
    parser.add_argument("--min-score-margin", type=int, default=12)
    parser.add_argument("--max-pairs", type=int)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.max_pairs is not None and args.max_pairs <= 0 or args.threads <= 0:
        parser.error("max-pairs and threads must be positive")
    rows = split_mixed_fastq(args.samples, args.primary_reference, args.competitor_reference,
                             args.output_dir, bowtie2=args.bowtie2, bowtie2_build=args.bowtie2_build,
                             min_pair_score=args.min_pair_score, min_score_margin=args.min_score_margin,
                             max_pairs=args.max_pairs, threads=args.threads)
    print(json.dumps({"samples": len(rows), "report": str(args.output_dir / "split_report.csv")}))


if __name__ == "__main__":
    main()
