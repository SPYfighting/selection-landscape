# NGS-derived selection score

This module turns a three-round, two-chain molecule count table into a continuous **NGS-derived selection score** named `label`. It is independent of model training and inference. The score summarizes measured selection trajectories; it is not a binding affinity or a probability of experimental success.

## Count filtering

The filter groups count rows by amino-acid sequence. A clone passes sequence quality control when it is productive, has the configured length and start/end motifs, contains no stop codon, and uses only the 20 canonical amino acids. It is label eligible when it has at least one molecule in each of rounds 1, 2 and 3. The output records `sequence_qc_eligible`, `label_eligible`, `forward_eligible_r2`, `model_eligible`, the reason codes and each round's pooled molecule count. `forward_eligible_r2` uses only rounds 1 and 2; it is an audit field and is not used by the published three-round score.

```bash
python -m selection_landscape.ngs.filter \
  --counts molecule_counts.tsv \
  --parameter-name strict \
  --length 126 --start-motif EVQ --end-motif VTVSS \
  --output eligibility.csv
```

Replace length and motifs with those of the studied target. The module accepts CSV or TSV counts, optionally gzip-compressed, with or without a UTF-8 byte-order mark. An existing output or an absent parameter set is rejected. It does not import the model. The label calculator requires counts from one reconstruction parameter set, matching the eligibility calculation; a table mixing parameter sets is rejected.

## Score calculation

For each of the six round/chain samples, the denominator is **all productive molecules**, including clones that fail ranking eligibility. An eligible clone's frequency in each chain equals its molecule count divided by that denominator. For round *r*, balanced abundance is the geometric mean of its Crick and Watson frequencies. The two enrichment ratios are balanced abundance in round 2 divided by round 1 and in round 3 divided by round 2.

The primary ranking set requires at least one molecule in each of the six chain/round cells and at least 10 pooled molecules in **each** round. Clones with 5–9 pooled molecules in their weakest round form a separate low-weight set. Log₂ enrichment for both transitions and round-3 balanced abundance are converted to ascending average-rank percentiles using the primary set. Low-weight percentiles are mapped onto those same frozen primary distributions.

Let `P12`, `P23` and `P3` be these percentiles. The transition average is `0.4 P12 + 0.6 P23`; trajectory is half this average plus half `min(P12, P23)`; base score is `0.8 trajectory + 0.2 P3`. Each enrichment ratio below 0.7 incurs a proportional decline penalty, combined at weights 0.4 and 0.6. The final `label` is the base minus 0.2 times that penalty, clipped to [0, 1]. The exact parameters are in `examples/ngs/selection_profile.json`.

```bash
python -m selection_landscape.labels \
  --counts molecule_counts.tsv \
  --eligibility eligibility.csv \
  --profile selection_profile.json \
  --target target.json \
  --output labels.csv
```

`target.json` contains a `target_id` and `lead_aa`; both are explicit assay metadata. The output contains `label`, rank, tier, support counts, balanced abundances, enrichment ratios and intermediate score components. The model consumes a separately prepared canonical training CSV and has no import path into these NGS modules.

## Synthetic example and frozen-data validation

The files in `examples/ngs/` are **synthetic**. They demonstrate two eligible clones, one primary and one low-weight:

```bash
python -m selection_landscape.ngs.filter \
  --counts examples/ngs/synthetic_counts.csv --parameter-name strict \
  --length 11 --start-motif EVQ --end-motif VTVSS \
  --output /tmp/example_eligibility.csv
python -m selection_landscape.labels \
  --counts examples/ngs/synthetic_counts.csv \
  --eligibility /tmp/example_eligibility.csv \
  --profile examples/ngs/selection_profile.json \
  --target examples/ngs/target.json \
  --output /tmp/example_labels.csv
```

For a private frozen full assay, the public filter reproduced all checked amino-acid identities, eligibility decisions, reason codes and per-round count fields. Passing its output to the public label calculator yielded the same complete CSV as using the historical eligibility table; the checked labels, ranks and enrichment components also matched the historical score output exactly. These checks start from the frozen molecule count table. They do not establish full-FASTQ byte-for-byte equivalence because the historical full-run reconstruction source snapshot is unavailable. Detailed sample sizes and comparisons remain in the private verification record.
