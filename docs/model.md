# Calibrated additive model

The model predicts a bounded experimental or selection label for a new double mutant, using measured single mutants of one parent sequence and observed double mutants for calibration. Label calculation is a separate operation: the model never reads sequencing counts, round names or a particular enrichment formula.

This package implements final fitting and inference for the selected calibrated additive model. It does not repeat historical architecture selection, compare protein language models, or enable a residual branch.

## Label contract

WT, single mutants and calibration doubles must use the same label definition, scale and experimental context. Labels must be finite numbers in [0, 1], with larger values representing the desired direction. Out-of-range values are rejected; they are never automatically normalized. Raw counts, enrichment ratios, dissociation constants and unbounded activities therefore require an explicit, scientifically justified conversion outside this module. A bounded label does not by itself establish biological comparability.

Sequences contain the 20 standard uppercase amino acids. Mutation identifiers use the parent residue, one-based position and substituted residue, for example `A1V`. A double mutant must contain exactly its two declared substitutions, at different positions, listed in increasing position order. Its canonical identity is `A1V__C2S`. Unmeasured single substitutions and variants with more than two substitutions are outside this implementation's prediction domain.

## Calculation

Let \(y_0\) be the WT label and \(y_i\) the observed label for single substitution \(i\). For a fixed, explicitly supplied \(\epsilon\), define

\[
\tilde y=\min(1-\epsilon,\max(\epsilon,y)),\qquad
a_i=\log\frac{\tilde y_i}{1-\tilde y_i}
    -\log\frac{\tilde y_0}{1-\tilde y_0}.
\]

The prediction for substitutions \(i,j\) is

\[
u_{ij}=a_i+a_j,\qquad
\widehat y_{ij}=\frac{1}{1+\exp[-(\beta_0+\beta_1u_{ij})]}.
\]

Only \(\beta_0\) and \(\beta_1\ge0\) are fitted to the observed double-mutant labels, by minimizing their weighted mean squared error. Single-mutant effects remain fixed. Fitting uses SciPy's L-BFGS-B implementation, initialized at \((0,1)\), with an unbounded intercept and a nonnegative slope. Omitted double-mutant weights are one; supplied weights must be finite and positive. The supplied epsilon must satisfy \(0<\epsilon<0.5\) and remain distinguishable from zero at the upper clipping boundary in float arithmetic. It is saved with the model, rather than recomputed from a later dataset.

This is a monotonic calibration of additive single-mutant effects, with two fitted parameters. It has no neural-network encoder or pair-specific interaction parameters. A predicted selection label is not a binding-affinity measurement. Different candidate combinations can have exactly equal scores; ties are ordered by canonical candidate ID.

## CSV files

Headers are exact and case-sensitive. Extra columns are rejected at model boundaries. UTF-8 files with or without a BOM are accepted; generated CSVs include a BOM for direct opening in Excel. Floating-point values are written without decimal rounding.

| Input | Required columns | Optional columns |
|---|---|---|
| WT (one row) | `sequence,label` | — |
| Singles | `mutation,sequence,label` | — |
| Calibration doubles | `mutation_1,mutation_2,sequence,label` | `pair_id,weight` |
| All observed pair identities | `mutation_1,mutation_2` | `pair_id` |
| Prediction candidates | `mutation_1,mutation_2,sequence` | `candidate_id` |
| Evaluation labels | `candidate_id,label,group_id` | — |

The calibration CSV contains only records intended for training. The complete observed-pair identity CSV is separate and should include training, validation, test and other already observed pairs, without their labels. Enumeration excludes every supplied observed identity and requires that all training identities are included. It combines measured substitutions at distinct positions; it does not generate an unrestricted sequence library.

Prediction output contains `candidate_id,mutation_1,mutation_2,sequence,predicted_label,rank`, sorted from highest score to lowest. A JSON model stores the WT, epsilon, single-mutant effects, calibration parameters, label description and training-pair identities, together with a SHA-256 integrity checksum. This checksum detects accidental modification; it is not a digital signature.

## Python interface

```python
from selection_landscape.model import (
    TrainingData, fit_model, enumerate_candidates, predict,
    save_model, load_model, evaluate_predictions,
)

data = TrainingData.from_records(wild_type, singles, training_doubles)
model = fit_model(data, epsilon=0.001, label_name="normalized_activity")
save_model(model, "model.json")
candidates = enumerate_candidates(model, all_observed_pair_identities)
predictions = predict(load_model("model.json"), candidates)
evaluation = evaluate_predictions(
    held_out_predictions, held_out_labels,
    bootstrap_replicates=2000, seed=20260825,
)
```

Inputs in this example are dictionaries or lists of dictionaries with the CSV fields above. Numeric CSV text and finite numeric values are normalized to Python floats. Labels are checked again at the fit boundary. Held-out prediction identities must match held-out labels exactly; evaluation does not silently take an intersection.

## Command-line interface

After installation, use `selection-landscape-model`, or equivalently `python -m selection_landscape.model.cli`.

```sh
selection-landscape-model train --wild-type wt.csv --singles singles.csv \
  --doubles training_doubles.csv --epsilon 0.001 --label-name normalized_activity \
  --output model.json
selection-landscape-model enumerate --model model.json \
  --observed-pairs all_observed_pairs.csv --output candidates.csv
selection-landscape-model predict --model model.json \
  --candidates candidates.csv --output predictions.csv
selection-landscape-model evaluate --predictions held_out_predictions.csv \
  --labels held_out_labels.csv --bootstrap-replicates 2000 --seed 20260825 \
  --output metrics.json --bootstrap-output bootstrap_metrics.csv
```

Output paths must not already exist and their parent directories must exist. Existing model files and results are never replaced. A runnable synthetic demonstration is provided in [`examples/model`](../examples/model/README.md).

## Evaluation and reproducibility

Evaluation reports Spearman correlation, Pearson correlation, mean absolute error and linear-gain NDCG at \(k=\lceil0.1N\rceil\). NDCG uses the supplied label directly as relevance, with discount \(1/\log_2(r+1)\); it does not apply exponential gain. Both predicted and ideal ranks break ties by canonical candidate ID. A zero ideal gain gives NDCG zero. A constant observed or predicted vector gives correlation zero by convention.

For bootstrap intervals, `group_id` must identify the unordered pair of mutated positions: substitutions at the same two positions belong to the same group. Each group must represent exactly one position pair, and each position pair must belong to exactly one group. The algorithm samples the sorted set of groups with replacement, draws as many groups as the original set contains, and retains every row of every sampled group, including repeated draws. It recomputes all four metrics and the top-10% cutoff in each replicate. The 95% interval is the linear-interpolated 2.5th–97.5th percentile. The optional bootstrap CSV preserves every replicate's values.

Bootstrap intervals describe uncertainty within the evaluated dataset. They do not substitute for a prospective test or establish generalization to another parent protein or assay. Architecture selection and its audit evidence are separate from this final-fit implementation.

For strict result preservation, reuse the frozen model artifact and retain full score precision, candidate identities and tie-breaking rules. Refitting the same objective with different numerical-library versions or floating-point parsing can change the final low-order digits and reorder nearly tied candidates. A refit should be evaluated as a separate artifact before replacing any established results. Record the input hashes, package versions, epsilon, weights, split identities and model checksum for every scientific run.
