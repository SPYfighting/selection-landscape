"""Sequence identity and bounded-label checks at public model boundaries."""

import math
import re


class ModelError(ValueError):
    """Model inputs or an artifact violate the documented contract."""


_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_MUTATION = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])$")


def finite_number(value, field):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ModelError(f"{field} must be a finite number") from exc
    if not math.isfinite(number):
        raise ModelError(f"{field} must be a finite number")
    return number


def bounded_label(value):
    number = finite_number(value, "label")
    if not 0.0 <= number <= 1.0:
        raise ModelError("label must lie in [0, 1]; transform labels explicitly before training")
    return number


def valid_epsilon(value):
    number = finite_number(value, "epsilon")
    if not 0.0 < number < 0.5:
        raise ModelError("epsilon must lie strictly between 0 and 0.5")
    if 1.0 - number == 1.0:
        raise ModelError("epsilon is too small to distinguish 1-epsilon from 1 in floating-point arithmetic")
    return number


def valid_sequence(value):
    if not isinstance(value, str) or not value or not set(value) <= _AMINO_ACIDS:
        raise ModelError("sequence must contain only the 20 standard uppercase amino acids")
    return value


def columns(record, required, optional=()):
    keys = set(record)
    missing = set(required) - keys
    unexpected = keys - set(required) - set(optional)
    if missing or unexpected:
        raise ModelError(f"missing fields {sorted(missing)}; unexpected fields {sorted(unexpected)}")


def mutation_identity(mutation, wild_type):
    match = _MUTATION.fullmatch(str(mutation))
    if match is None:
        raise ModelError(f"invalid mutation: {mutation}")
    wild, position, changed = match.group(1), int(match.group(2)), match.group(3)
    if position > len(wild_type) or wild_type[position - 1] != wild or wild == changed:
        raise ModelError(f"mutation does not match the WT sequence: {mutation}")
    return position, changed


def pair_sequence(first, second, wild_type, known_mutations):
    if first not in known_mutations or second not in known_mutations:
        raise ModelError("pair refers to an unknown single mutation")
    position_1, changed_1 = mutation_identity(first, wild_type)
    position_2, changed_2 = mutation_identity(second, wild_type)
    if position_1 >= position_2:
        raise ModelError("pair must use different positions in canonical increasing position order")
    sequence = list(wild_type)
    sequence[position_1 - 1], sequence[position_2 - 1] = changed_1, changed_2
    return "".join(sequence)


def validated_pairs(records, wild_type, known_mutations, *, labelled=False):
    required = {"mutation_1", "mutation_2", "sequence"}
    optional = {"candidate_id"}
    if labelled:
        required.add("label")
        optional = {"pair_id", "weight"}
    output, seen = [], set()
    for source in records:
        row = dict(source)
        columns(row, required, optional)
        first, second = str(row["mutation_1"]), str(row["mutation_2"])
        identifier = f"{first}__{second}"
        rebuilt = pair_sequence(first, second, wild_type, known_mutations)
        if str(row["sequence"]) != rebuilt:
            raise ModelError(f"pair sequence has incorrect or additional substitutions: {identifier}")
        id_field = "pair_id" if labelled else "candidate_id"
        if id_field in row and row[id_field] != identifier:
            raise ModelError(f"{id_field} does not match its mutations")
        if identifier in seen:
            raise ModelError(f"duplicate pair: {identifier}")
        seen.add(identifier)
        normalized = {id_field: identifier, "mutation_1": first, "mutation_2": second,
                      "sequence": rebuilt}
        if labelled:
            normalized["label"] = bounded_label(row["label"])
            weight = finite_number(row.get("weight", 1.0), "weight")
            if weight <= 0:
                raise ModelError("weight must be positive")
            normalized["weight"] = weight
        output.append(normalized)
    return output


def validated_training(wild_type, single_mutants, observed_doubles):
    columns(wild_type, {"sequence", "label"})
    wt = {"sequence": valid_sequence(wild_type["sequence"]),
          "label": bounded_label(wild_type["label"])}
    singles, known = [], set()
    for source in single_mutants:
        columns(source, {"mutation", "sequence", "label"})
        mutation = str(source["mutation"])
        position, changed = mutation_identity(mutation, wt["sequence"])
        expected = list(wt["sequence"])
        expected[position - 1] = changed
        if source["sequence"] != "".join(expected):
            raise ModelError("single-mutant sequence must encode exactly the declared change")
        if mutation in known:
            raise ModelError(f"duplicate single mutation: {mutation}")
        known.add(mutation)
        singles.append({"mutation": mutation, "sequence": source["sequence"],
                        "label": bounded_label(source["label"])})
    doubles = validated_pairs(observed_doubles, wt["sequence"], known, labelled=True)
    if len(doubles) < 2:
        raise ModelError("calibration requires at least two observed double mutants")
    return wt, tuple(singles), tuple(doubles)
