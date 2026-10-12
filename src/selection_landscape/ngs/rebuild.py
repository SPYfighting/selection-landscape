"""Indel-aware UMI consensus and VHH molecule counting from paired FASTQ.

The reconstruction function preserves the historical CIGAR vote, per-UID
consensus, minimum ORF coverage, translation, and molecule counting rules.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Optional

from .vhh import aa_productive, alignment_orf_votes, find_anchor_uid, translate_orf
from .read_flow import reconcile_alignment_flow

def resolve_tool(project_root: Path, configured: str) -> Path:
    path = Path(configured)
    if path.is_absolute() and path.is_file():
        return path
    if not path.is_absolute() and (project_root / path).is_file():
        return project_root / path
    discovered = shutil.which(configured)
    if discovered:
        return Path(discovered)
    raise FileNotFoundError(f"required command not found: {configured}")


def _read_name(header: str) -> str:
    fields = header[1:].split()
    if not fields:
        raise ValueError("FASTQ read identifier is empty")
    name = fields[0]
    # Bowtie2 removes one terminal mate suffix when writing the SAM QNAME.
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    if not name:
        raise ValueError("FASTQ read identifier is empty after removing its mate suffix")
    return name


def paired_fastq_records(r1: Path, r2: Path):
    """Yield validated paired FASTQ records; reject truncation and mate mismatch."""
    opener1 = gzip.open if r1.suffix == ".gz" else open
    opener2 = gzip.open if r2.suffix == ".gz" else open
    with opener1(r1, "rt", encoding="ascii") as first, opener2(r2, "rt", encoding="ascii") as second:
        index = 0
        seen_names: set[str] = set()
        while True:
            rec1 = [first.readline() for _ in range(4)]
            rec2 = [second.readline() for _ in range(4)]
            if not any(rec1) and not any(rec2):
                break
            index += 1
            for mate, rec in (("R1", rec1), ("R2", rec2)):
                if not all(rec):
                    raise ValueError(f"truncated or unpaired FASTQ at pair {index}, {mate}")
                header, sequence, separator, quality = (line.rstrip("\r\n") for line in rec)
                if not header.startswith("@") or not separator.startswith("+") or not sequence or len(sequence) != len(quality):
                    raise ValueError(f"invalid FASTQ record at pair {index}, {mate}")
            id1 = _read_name(rec1[0])
            id2 = _read_name(rec2[0])
            if id1 != id2:
                raise ValueError(f"paired FASTQ identifiers differ at pair {index}")
            if id1 in seen_names:
                raise ValueError(f"duplicate paired FASTQ identifier at pair {index}: {id1}")
            seen_names.add(id1)
            yield tuple(line.rstrip("\r\n") for line in rec1), tuple(line.rstrip("\r\n") for line in rec2)


def rebuild_sample(
    *,
    sample: dict[str, str],
    r1: Path,
    r2: Path,
    param: dict[str, Any],
    config: dict[str, Any],
    project_root: Path,
    work_dir: Path,
    max_reads: Optional[int],
) -> dict[str, Any]:
    """PE-align trimmed reads; UID from terminal anchor; ORF by per-UID position consensus."""

    anchor = config["uid"]["anchor"]
    uid_len = int(config["uid"]["uid_length"])
    max_mismatch = int(param["anchor_mismatch_max"])
    orf = config["orf"]
    orf_start = int(orf["nt_start_zero_based"])
    orf_end = int(orf["nt_end_exclusive"])
    min_cov = float(param.get("min_orf_coverage", config.get("min_orf_coverage", 0.9)))
    bowtie2 = resolve_tool(project_root, config["tools"]["bowtie2"])
    samtools = resolve_tool(project_root, config["tools"]["samtools"])
    index = str(project_root / config["bowtie2_index"])

    work_dir.mkdir(parents=True, exist_ok=True)
    uid_map = work_dir / "read_uid.tsv"
    sorted_bam = work_dir / "aligned.sorted.bam"
    resume_bam = sorted_bam.exists() and (work_dir / "aligned.sorted.bam.bai").exists()
    resume_uid = uid_map.exists() and uid_map.stat().st_size > 0

    # When max_reads is set, write a named subset FASTQ for reproducible pilot sizing.
    # Full runs align the gzipped trimmed FASTQs directly (no uncompressed rewrite).
    use_subset = max_reads is not None
    pe_r1 = work_dir / "subset_R1.fastq" if use_subset else r1
    pe_r2 = work_dir / "subset_R2.fastq" if use_subset else r2

    uid_found = 0
    reads_in = 0
    if resume_uid and resume_bam and not use_subset:
        # Resume after OOM/kill: keep existing UID map + BAM.
        with uid_map.open(encoding="utf-8") as handle:
            next(handle, None)
            for line in handle:
                if line.strip():
                    uid_found += 1
        reads_in = uid_found  # lower bound; exact reads_in recorded from prior bowtie stderr if needed
        print(f"[progress] resume UID+BAM {sample['sample_id']} uid_rows={uid_found}", file=sys.stderr, flush=True)
    else:
        with uid_map.open("w", encoding="utf-8", newline="") as um:
            writer = csv.writer(um, delimiter="\t")
            writer.writerow(["read_id", "uid", "orientation", "anchor_mismatch", "mate"])
            if use_subset:
                f1 = pe_r1.open("w", encoding="utf-8")
                f2 = pe_r2.open("w", encoding="utf-8")
            else:
                f1 = f2 = None
            try:
                for rec1, rec2 in paired_fastq_records(r1, r2):
                    h1_id, s1, _, q1 = rec1
                    _, s2, _, q2 = rec2
                    reads_in += 1
                    if use_subset:
                        read_id = f"r{reads_in - 1}"
                    else:
                        read_id = _read_name(h1_id)
                    hit1 = find_anchor_uid(s1, anchor, uid_len, max_mismatch)
                    hit2 = None if hit1 else find_anchor_uid(s2, anchor, uid_len, max_mismatch)
                    hit = hit1 or hit2
                    if hit:
                        uid_found += 1
                        writer.writerow(
                            [
                                read_id,
                                hit["uid"],
                                hit["orientation"],
                                hit["anchor_mismatch"],
                                "R1" if hit1 else "R2",
                            ]
                        )
                    if use_subset:
                        assert f1 is not None and f2 is not None
                        f1.write(f"@{read_id}\n{s1}\n+\n{q1 or 'I' * len(s1)}\n")
                        f2.write(f"@{read_id}\n{s2}\n+\n{q2 or 'I' * len(s2)}\n")
                    if max_reads is not None and reads_in >= max_reads:
                        break
            finally:
                if f1 is not None:
                    f1.close()
                if f2 is not None:
                    f2.close()

    if uid_found == 0:
        return {
            "sample_id": sample["sample_id"],
            "status": "FAILED",
            "reason_code": "no_uid",
            "reads_in": reads_in,
            "uid_found": 0,
            "aligned": 0,
            "productive": 0,
            "unique_aa": 0,
            "lead_recovered": False,
            "clones": [],
            "param_name": param["name"],
        }

    mode = param["bowtie2_mode"].split()
    cmd = [
        str(bowtie2),
        "-x",
        index,
        "-1",
        str(pe_r1),
        "-2",
        str(pe_r2),
        "-p",
        str(config.get("threads", 8)),
        *mode,
    ]
    if not resume_bam:
        # Stream SAM→BAM→sort to avoid multi-GB SAM on full libraries.
        bt = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
        view = subprocess.Popen(
            [str(samtools), "view", "-bS", "-"],
            stdin=bt.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if bt.stdout is not None:
            bt.stdout.close()
        sort_proc = subprocess.run(
            [str(samtools), "sort", "-o", str(sorted_bam), "-"],
            stdin=view.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if view.stdout is not None:
            view.stdout.close()
        bt_stderr = bt.stderr.read().decode("utf-8", errors="replace") if bt.stderr else ""
        view_err = view.stderr.read().decode("utf-8", errors="replace") if view.stderr else ""
        bt_rc = bt.wait()
        view_rc = view.wait()
        (work_dir / "bowtie2.stderr").write_text(
            bt_stderr + ("\n" + view_err if view_err else ""), encoding="utf-8"
        )
        if bt_rc != 0 or view_rc != 0 or sort_proc.returncode != 0:
            return {
                "sample_id": sample["sample_id"],
                "status": "FAILED",
                "reason_code": "bowtie2_or_samtools_failed",
                "reads_in": reads_in,
                "uid_found": uid_found,
                "aligned": 0,
                "productive": 0,
                "unique_aa": 0,
                "lead_recovered": False,
                "clones": [],
                "param_name": param["name"],
            }
        subprocess.run([str(samtools), "index", str(sorted_bam)], check=True)
    else:
        print(f"[progress] resume BAM {sample['sample_id']}", file=sys.stderr, flush=True)

    # Compact numpy allele matrices per UID (uint16[orf_len,5]) — far smaller than nested Counters.
    import numpy as np

    read_uid: dict[str, str] = {}
    with uid_map.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            read_uid[row["read_id"]] = row["uid"]

    base_index = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
    base_chars = "ACGT-"
    orf_len = orf_end - orf_start
    uid_mats: dict[str, Any] = {}
    uid_hits: Counter = Counter()
    aligned = 0
    aligned_without_orf_vote = 0
    aligned_with_orf_vote = 0
    indel_events = 0
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("pysam is required for UMI reconstruction") from exc
    with pysam.AlignmentFile(str(sorted_bam), "rb") as af:
        for aln in af.fetch(until_eof=True):
            if aln.is_unmapped or aln.mapping_quality < int(param.get("min_mapq", 0)):
                continue
            uid = read_uid.get(aln.query_name, "")
            if not uid:
                continue
            aligned += 1
            votes, events = alignment_orf_votes(
                aln.query_sequence or "",
                list(aln.cigartuples or []),
                int(aln.reference_start),
                orf_start,
                orf_end,
            )
            indel_events += len(events)
            if not votes:
                aligned_without_orf_vote += 1
                continue
            aligned_with_orf_vote += 1
            mat = uid_mats.get(uid)
            if mat is None:
                mat = np.zeros((orf_len, 5), dtype=np.uint16)
                uid_mats[uid] = mat
            uid_hits[uid] += 1
            for pos, base in votes.items():
                bi = base_index.get(base)
                if bi is None:
                    continue
                i = pos - orf_start
                if 0 <= i < orf_len and mat[i, bi] < 65535:
                    mat[i, bi] += 1

    nt_uids: dict[str, set[str]] = defaultdict(set)
    nt_reads: Counter = Counter()
    nt_indel: dict[str, bool] = {}
    min_cov_n = int(math.ceil(min_cov * orf_len))
    for uid, mat in uid_mats.items():
        covered = int(np.count_nonzero(mat.sum(axis=1)))
        if covered < min_cov_n:
            continue
        chars: list[str] = []
        has_indel = False
        for i in range(orf_len):
            row = mat[i]
            if row.sum() == 0:
                chars.append("N")
                continue
            bi = int(row.argmax())
            ch = base_chars[bi]
            chars.append(ch)
            if ch == "-":
                has_indel = True
        nt = "".join(chars)
        if nt.count("N") > orf_len - min_cov_n:
            continue
        nt_uids[nt].add(uid)
        nt_reads[nt] += int(uid_hits[uid])
        nt_indel[nt] = nt_indel.get(nt, False) or has_indel
    uid_mats.clear()

    clones = []
    productive = 0
    lead = orf["lead_aa"]
    lead_recovered = False
    for nt, uids in nt_uids.items():
        aa = translate_orf(nt, frame=int(orf["frame_zero_based"]))
        ok, reason = aa_productive(aa, orf["start_motif"], orf["end_motif"], int(orf["aa_length"]))
        if ok:
            productive += 1
        if aa == lead:
            lead_recovered = True
        clones.append(
            {
                "sample_id": sample["sample_id"],
                "round": int(sample["round"]),
                "chain": sample["chain"],
                "nt_sequence": nt,
                "aa_sequence": aa,
                "aa_length": len(aa),
                "read_count": int(nt_reads[nt]),
                "molecule_count": len(uids),
                "productive": ok,
                "productive_reason": reason,
                "has_indel_gaps": bool(nt_indel.get(nt, False)),
                "uid_ok": True,
                "strand_conflict": False,
                "start_motif": orf["start_motif"],
                "end_motif": orf["end_motif"],
            }
        )

    released_clone_read_count = sum(int(clone["read_count"]) for clone in clones)
    read_flow = reconcile_alignment_flow(
        aligned_records=aligned,
        records_without_orf_vote=aligned_without_orf_vote,
        clone_table_read_count=released_clone_read_count,
    )
    if read_flow["records_with_orf_vote"] != aligned_with_orf_vote:
        raise ValueError("internal alignment-flow vote count mismatch")

    return {
        "sample_id": sample["sample_id"],
        "round": int(sample["round"]),
        "chain": sample["chain"],
        "status": "SUCCESS",
        "reason_code": "ok",
        "reads_in": reads_in,
        "uid_found": uid_found,
        "aligned": aligned,
        "aligned_without_orf_vote": aligned_without_orf_vote,
        "aligned_with_orf_vote": aligned_with_orf_vote,
        "aligned_not_released_after_orf_vote": read_flow[
            "records_not_released_after_vote"
        ],
        "released_clone_read_count": released_clone_read_count,
        "aligned_to_clone_read_conservation_ok": read_flow[
            "conservation_ok"
        ],
        "indel_events": indel_events,
        "productive": productive,
        "molecule_count": sum(len(uids) for uids in nt_uids.values()),
        "unique_nt": len(nt_uids),
        "unique_aa": len({c["aa_sequence"] for c in clones if c["productive"]}),
        "lead_recovered": lead_recovered,
        "clones": clones,
        "param_name": param["name"],
    }


def count_fastq(
    sample_manifest: Path,
    reference_fasta: Path,
    configuration: dict[str, Any],
    output_dir: Path,
    *,
    parameter_name: str,
    bowtie2: str = "bowtie2",
    bowtie2_build: str = "bowtie2-build",
    samtools: str = "samtools",
    max_reads: int | None = None,
    threads: int = 1,
) -> list[dict[str, Any]]:
    """Count UMI consensus molecules for each sample in a CSV manifest."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("output directory must be empty")
    output_dir.mkdir(parents=True, exist_ok=True)
    params = [item for item in configuration["parameter_sets"] if item["name"] == parameter_name]
    if len(params) != 1:
        raise ValueError(f"exactly one named reconstruction parameter set is required: {parameter_name}")
    index_dir = output_dir / "reference_index"
    index_dir.mkdir()
    reference_copy = index_dir / "reference.fa"
    shutil.copyfile(reference_fasta, reference_copy)
    index = index_dir / "ref"
    with (index_dir / "bowtie2_build.stdout.log").open("w") as stdout, (index_dir / "bowtie2_build.stderr.log").open("w") as stderr:
        subprocess.run([bowtie2_build, str(reference_copy), str(index)], check=True,
                       stdout=stdout, stderr=stderr)
    config = dict(configuration)
    config["bowtie2_index"] = str(index.relative_to(output_dir))
    config["tools"] = {"bowtie2": bowtie2, "samtools": samtools}
    config["threads"] = threads
    delimiter = "\t" if sample_manifest.suffix == ".tsv" else ","
    with sample_manifest.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter=delimiter))
    if not samples:
        raise ValueError("sample manifest is empty")
    identifiers = [sample["sample_id"] for sample in samples]
    if len(set(identifiers)) != len(identifiers) or any(not item or item in {".", ".."} or "/" in item or "\\" in item for item in identifiers):
        raise ValueError("sample identifiers must be unique single path components")
    summaries: list[dict[str, Any]] = []
    clones: list[dict[str, Any]] = []
    for sample in samples:
        first = Path(sample["r1_path"])
        second = Path(sample["r2_path"])
        if not first.is_absolute():
            first = sample_manifest.parent / first
        if not second.is_absolute():
            second = sample_manifest.parent / second
        result = rebuild_sample(
            sample=sample, r1=first, r2=second, param=params[0], config=config,
            project_root=output_dir, work_dir=output_dir / "work" / sample["sample_id"],
            max_reads=max_reads,
        )
        if result["status"] != "SUCCESS":
            raise RuntimeError(f"reconstruction failed for {sample['sample_id']}: {result['reason_code']}")
        for clone in result["clones"]:
            clones.append({**clone, "param_name": parameter_name})
        summaries.append({key: value for key, value in result.items() if key != "clones"})
    fields = ["sample_id", "round", "chain", "param_name", "nt_sequence", "aa_sequence",
              "aa_length", "read_count", "molecule_count", "productive", "productive_reason",
              "has_indel_gaps"]
    with (output_dir / "molecule_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(clones)
    (output_dir / "sample_summary.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
    return summaries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--parameter-name", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bowtie2", default="bowtie2")
    parser.add_argument("--bowtie2-build", default="bowtie2-build")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--max-reads", type=int)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.max_reads is not None and args.max_reads <= 0 or args.threads <= 0:
        parser.error("max-reads and threads must be positive")
    config = json.loads(args.config.read_text(encoding="utf-8"))
    summaries = count_fastq(
        args.samples, args.reference, config, args.output_dir,
        parameter_name=args.parameter_name, bowtie2=args.bowtie2,
        bowtie2_build=args.bowtie2_build, samtools=args.samtools,
        max_reads=args.max_reads, threads=args.threads,
    )
    print(json.dumps({"samples": len(summaries), "molecule_counts": str(args.output_dir / "molecule_counts.tsv")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
