from selection_landscape.ngs.split import _paired_scores, classify_pair


def _sam(name: str, flag: int, reference: str, score: int) -> str:
    return f"{name}\t{flag}\t{reference}\t1\t60\t4M\t=\t1\t0\tACGT\tIIII\tAS:i:{score}\n"


def test_competitive_pair_scores_and_margin():
    records = [_sam("read", 64, "primary", 48), _sam("read", 128, "primary", 50),
               _sam("read", 64, "competitor", 40), _sam("read", 128, "competitor", 39)]
    scores = _paired_scores(records)
    assert scores == {"primary": 98, "competitor": 79}
    assert classify_pair(scores, min_pair_score=60, min_score_margin=12) == "primary-only"
    assert classify_pair({"primary": 90, "competitor": 80}, min_pair_score=60, min_score_margin=12) == "ambiguous"
    assert classify_pair({"primary": 59, "competitor": None}, min_pair_score=60, min_score_margin=12) == "unassigned"
