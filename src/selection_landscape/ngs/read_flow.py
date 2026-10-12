"""Exact accounting helpers for alignment-record flow into VHH clone tables."""

from __future__ import annotations

from typing import Any, Sequence


def cigar_has_orf_vote(
    cigar_tuples: Sequence[tuple[int, int]],
    ref_start: int,
    orf_start: int,
    orf_end: int,
) -> bool:
    """Return whether the rebuild engine emits at least one ORF vote."""
    ref_position = ref_start
    for operation, length in cigar_tuples:
        if operation in (0, 7, 8, 2):
            ref_end = ref_position + length
            if ref_position < orf_end and ref_end > orf_start:
                return True
            ref_position = ref_end
        elif operation == 3:
            ref_position += length
        elif operation in (1, 4, 5):
            continue
    return False


def reconcile_alignment_flow(
    *,
    aligned_records: int,
    records_without_orf_vote: int,
    clone_table_read_count: int,
) -> dict[str, Any]:
    """Partition accepted alignment records into exact downstream outcomes."""
    values = (
        aligned_records,
        records_without_orf_vote,
        clone_table_read_count,
    )
    if any(not isinstance(value, int) or value < 0 for value in values):
        raise ValueError("alignment-flow counts must be non-negative integers")
    records_with_orf_vote = aligned_records - records_without_orf_vote
    records_not_released_after_vote = (
        records_with_orf_vote - clone_table_read_count
    )
    if records_with_orf_vote < 0 or records_not_released_after_vote < 0:
        raise ValueError("alignment-flow counts are not monotonic")
    return {
        "aligned_records": aligned_records,
        "records_without_orf_vote": records_without_orf_vote,
        "records_with_orf_vote": records_with_orf_vote,
        "records_not_released_after_vote": records_not_released_after_vote,
        "clone_table_read_count": clone_table_read_count,
        "conservation_ok": (
            aligned_records
            == records_without_orf_vote
            + records_not_released_after_vote
            + clone_table_read_count
        ),
        "reason_definitions": {
            "records_without_orf_vote": (
                "accepted alignment has no match/deletion position inside "
                "the frozen ORF interval"
            ),
            "records_not_released_after_vote": (
                "alignment contributed an ORF vote, but its UID consensus "
                "did not pass the frozen coverage/unresolved-base release gate"
            ),
            "clone_table_read_count": (
                "alignment is represented in a released UID consensus clone"
            ),
        },
    }
