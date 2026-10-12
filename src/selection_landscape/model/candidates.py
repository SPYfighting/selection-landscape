"""Enumerate new two-position combinations of measured single mutations."""

from itertools import combinations
from typing import Iterable, Mapping

from .core import CalibratedAdditiveModel
from .validation import ModelError, columns, mutation_identity, pair_sequence


def enumerate_candidates(
    model: CalibratedAdditiveModel, observed_pairs: Iterable[Mapping[str, object]]
) -> list[dict[str, object]]:
    """Exclude a caller-supplied, complete, label-free list of observed pairs.

    Include observed identities even if they were held out, had low support or
    otherwise supplied no training label. The training set alone is insufficient.
    """
    known = model.mutation_effects
    observed = set()
    for record in observed_pairs:
        columns(record, {"mutation_1", "mutation_2"}, {"pair_id"})
        first, second = str(record["mutation_1"]), str(record["mutation_2"])
        pair_sequence(first, second, model.wild_type_sequence, known)
        identifier = f"{first}__{second}"
        if "pair_id" in record and record["pair_id"] != identifier:
            raise ModelError("observed pair_id does not match its mutations")
        if identifier in observed:
            raise ModelError(f"duplicate observed pair: {identifier}")
        observed.add(identifier)
    if not set(model.training_pairs) <= observed:
        raise ModelError("observed identity list omits training pairs")
    ordered = sorted(known, key=lambda mutation: (mutation_identity(mutation, model.wild_type_sequence)[0], mutation))
    rows = []
    for first, second in combinations(ordered, 2):
        p1 = mutation_identity(first, model.wild_type_sequence)[0]
        p2 = mutation_identity(second, model.wild_type_sequence)[0]
        if p1 == p2 or f"{first}__{second}" in observed:
            continue
        rows.append({"candidate_id": f"{first}__{second}", "mutation_1": first,
                     "mutation_2": second,
                     "sequence": pair_sequence(first, second, model.wild_type_sequence, known)})
    return rows
