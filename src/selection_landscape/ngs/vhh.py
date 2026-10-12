"""UID extraction and indel-aware VHH NT/AA reconstruction."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def find_anchor_uid(
    sequence: str,
    anchor: str,
    uid_length: int,
    max_mismatch: int,
    terminal_window: int = 50,
) -> Optional[dict[str, Any]]:
    """Return a UID adjacent to the configured anchor.

    Observed PE layout:
    - R2 often forward-sense 5' of amplicon.
    - R1 often RC of 3' end: RC(UID) then RC(anchor) then RC(body).
    Search only terminal windows to avoid coding-sequence false anchors.
    """

    seq = sequence.upper()
    anchor = anchor.upper()
    rc_anchor = revcomp(anchor)

    # Forward: anchor then UID near 3' terminus.
    fwd_start = max(0, len(seq) - terminal_window - uid_length)
    fwd_end = max(0, len(seq) - len(anchor) - uid_length + 1)
    best: Optional[dict[str, Any]] = None
    for i in range(fwd_start, fwd_end):
        window = seq[i : i + len(anchor)]
        mm = hamming(window, anchor)
        if mm <= max_mismatch:
            uid = seq[i + len(anchor) : i + len(anchor) + uid_length]
            if len(uid) == uid_length and set(uid) <= set("ACGT"):
                cand = {
                    "orientation": "forward",
                    "anchor_start": i,
                    "uid": uid,
                    "body": seq[:i],
                    "anchor_mismatch": mm,
                }
                if best is None or mm < best["anchor_mismatch"]:
                    best = cand

    # RC layout near 5': strip leading N junk, then RC(UID)+RC(anchor).
    lead_n = len(seq) - len(seq.lstrip("N"))
    core = seq[lead_n:]
    rc_end = min(terminal_window, max(0, len(core) - len(rc_anchor) + 1))
    for i in range(uid_length, rc_end):
        window = core[i : i + len(rc_anchor)]
        mm = hamming(window, rc_anchor)
        if mm <= max_mismatch:
            uid_rc = core[i - uid_length : i]
            if len(uid_rc) == uid_length and set(uid_rc) <= set("ACGT"):
                uid = revcomp(uid_rc)
                cand = {
                    "orientation": "reverse_complement",
                    "anchor_start": lead_n + i,
                    "uid": uid,
                    "body": revcomp(core[i + len(rc_anchor) :]),
                    "anchor_mismatch": mm,
                }
                if best is None or mm < best["anchor_mismatch"]:
                    best = cand
    return best


def cigar_to_query_ref_sequence(
    query_seq: str,
    cigar_tuples: Sequence[tuple[int, int]],
    ref_start: int,
) -> tuple[str, list[dict[str, Any]]]:
    """Rebuild reference-oriented NT allowing insertions/deletions.

    Returns aligned_ref_span_query_chars for positions covering the reference
    and a list of indel events.
    """

    # pysam: 0=M,1=I,2=D,3=N,4=S,5=H,7==,8=X
    q_pos = 0
    r_pos = ref_start
    rebuilt: dict[int, str] = {}
    events: list[dict[str, Any]] = []
    for op, length in cigar_tuples:
        if op in (0, 7, 8):  # match/equal/diff
            for _ in range(length):
                rebuilt[r_pos] = query_seq[q_pos]
                q_pos += 1
                r_pos += 1
        elif op == 1:  # insertion in query
            insert = query_seq[q_pos : q_pos + length]
            events.append({"type": "insertion", "ref_pos": r_pos, "seq": insert})
            q_pos += length
        elif op == 2:  # deletion in query
            for _ in range(length):
                rebuilt[r_pos] = "-"
                events.append({"type": "deletion", "ref_pos": r_pos, "seq": "-"})
                r_pos += 1
        elif op == 4:  # soft clip
            q_pos += length
        elif op == 5:  # hard clip
            continue
        elif op == 3:  # ref skip
            r_pos += length
        else:
            continue
    if not rebuilt:
        return "", events
    min_r = min(rebuilt)
    max_r = max(rebuilt)
    chars = [rebuilt.get(pos, "N") for pos in range(min_r, max_r + 1)]
    return "".join(chars), events


def alignment_orf_votes(
    query_seq: str,
    cigar_tuples: Sequence[tuple[int, int]],
    ref_start: int,
    orf_start: int,
    orf_end: int,
) -> tuple[dict[int, str], list[dict[str, Any]]]:
    """Map one alignment onto ORF reference coordinates (base or '-')."""

    q_pos = 0
    r_pos = ref_start
    votes: dict[int, str] = {}
    events: list[dict[str, Any]] = []
    for op, length in cigar_tuples:
        if op in (0, 7, 8):
            for _ in range(length):
                if orf_start <= r_pos < orf_end:
                    votes[r_pos] = query_seq[q_pos]
                q_pos += 1
                r_pos += 1
        elif op == 1:
            insert = query_seq[q_pos : q_pos + length]
            if orf_start <= r_pos < orf_end:
                events.append({"type": "insertion", "ref_pos": r_pos, "seq": insert})
            q_pos += length
        elif op == 2:
            for _ in range(length):
                if orf_start <= r_pos < orf_end:
                    votes[r_pos] = "-"
                    events.append({"type": "deletion", "ref_pos": r_pos, "seq": "-"})
                r_pos += 1
        elif op == 4:
            q_pos += length
        elif op == 3:
            r_pos += length
        elif op == 5:
            continue
        else:
            continue
    return votes, events


def consensus_orf_from_position_votes(
    position_counters: Mapping[int, Counter],
    orf_start: int,
    orf_end: int,
    min_coverage: float = 0.9,
) -> Optional[dict[str, Any]]:
    """Majority-vote ORF NT string from per-position base counters."""

    chars: list[str] = []
    covered = 0
    indel_gaps = False
    for pos in range(orf_start, orf_end):
        counter = position_counters.get(pos)
        if not counter:
            chars.append("N")
            continue
        base, _ = counter.most_common(1)[0]
        chars.append(base)
        covered += 1
        if base == "-":
            indel_gaps = True
    length = orf_end - orf_start
    frac = covered / length if length else 0.0
    if frac < min_coverage:
        return None
    nt = "".join(chars)
    if "N" in nt.replace("-", ""):
        # Allow sparse interior N only if coverage gate already passed on observed positions;
        # reject unresolved IUPAC for translation fidelity.
        if nt.count("N") > int((1.0 - min_coverage) * length):
            return None
    return {"nt_sequence": nt, "coverage_fraction": frac, "has_indel_gaps": indel_gaps}


def translate_orf(nt: str, frame: int = 0) -> str:
    from Bio.Seq import Seq

    clean = nt.replace("-", "")
    seq = clean[frame:]
    seq = seq[: len(seq) - (len(seq) % 3)]
    if not seq:
        return ""
    return str(Seq(seq).translate(to_stop=False))


def aa_productive(aa: str, start_motif: str, end_motif: str, expected_len: int) -> tuple[bool, str]:
    if not aa:
        return False, "empty_aa"
    if "*" in aa:
        return False, "stop_codon"
    if len(aa) != expected_len:
        return False, f"aa_length_{len(aa)}"
    if not aa.startswith(start_motif):
        return False, "bad_start_motif"
    if not aa.endswith(end_motif):
        return False, "bad_end_motif"
    return True, "ok"


@dataclass
class RebuildStats:
    reads_in: int = 0
    uid_found: int = 0
    aligned: int = 0
    productive: int = 0
    unique_nt: int = 0
    unique_aa: int = 0
    lead_recovered: bool = False


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(columns)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in cols})


def majority_nt_from_uid_families(
    uid_to_nts: Mapping[str, Counter],
) -> dict[str, Any]:
    """Collapse UID families: each UID votes for its majority NT; then aggregate AA clones."""

    uid_nt: dict[str, str] = {}
    for uid, counter in uid_to_nts.items():
        if not counter:
            continue
        uid_nt[uid] = counter.most_common(1)[0][0]
    nt_uids: dict[str, set[str]] = defaultdict(set)
    nt_reads: Counter = Counter()
    for uid, nt in uid_nt.items():
        nt_uids[nt].add(uid)
        nt_reads[nt] += sum(uid_to_nts[uid].values())
    rows = []
    for nt, uids in nt_uids.items():
        rows.append(
            {
                "nt_sequence": nt,
                "read_count": int(nt_reads[nt]),
                "molecule_count": len(uids),
            }
        )
    return {"uid_assigned": len(uid_nt), "clones": rows}
