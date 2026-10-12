# Selection Landscape

CPU-based tools for processing targeted sequencing selections, calculating selection labels, and predicting new combinations of measured protein substitutions.

The workflow has three independent stages:

1. **NGS processing:** assign paired reads to competing references and reconstruct productive sequences and UMI molecule counts.
2. **Label calculation:** filter molecule-count tables and calculate a bounded selection score from three rounds and two sequencing chains.
3. **Combination prediction:** use measured single-mutant labels and observed double mutants to fit a calibrated additive model and rank new double mutants.

The model accepts generic, finite `label` values in **[0, 1]**, with higher values representing the desired direction. It does not import the NGS or label-calculation modules. WT, single mutants and calibration doubles must share the same label definition and scale. Labels outside this range are rejected rather than automatically normalized.

This repository contains MIT-licensed source code and synthetic examples. It does not include experimental reads, study datasets, or fitted models from experiments.

## Installation

Python **3.11 or later** is required. From a source checkout:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The model requires NumPy and SciPy and runs on a CPU. No GPU, pretrained protein model or embedding download is required. For raw-read processing and tests, install the optional dependencies:

```sh
python -m pip install -e ".[ngs,test]"
```

Raw-read processing also requires the external **Bowtie2**, **bowtie2-build** and **samtools** executables on `PATH`, or explicit executable paths supplied to the relevant commands. The `ngs` extra installs Biopython and pysam. Filtering existing count tables and calculating labels do not require these alignment tools.

[`requirements-tested.txt`](requirements-tested.txt) records the Python 3.13 environment used for verification. These exact versions are distinct from the supported dependency ranges in [`pyproject.toml`](pyproject.toml).

## Run the synthetic model example

Run this command from the repository root. The output directory must be new.

```sh
python examples/model/run_synthetic_example.py --output example-model-output
```

The example uses artificial labels generated from explicitly specified positive and negative single-mutant effects. It trains on two synthetic double mutants, excludes those identities, and predicts the four remaining combinations. It checks the predicted order and agreement with independently specified generating scores.

| Output | Contents |
|---|---|
| `synthetic_model.json` | Single-mutant effects, calibration parameters and integrity checksum |
| `synthetic_candidates.csv` | Unobserved double-mutant identities and sequences |
| `synthetic_predictions.csv` | Predicted labels and ranks |
| `synthetic_metrics.json` | Correlations, MAE, top-10% NDCG and bootstrap intervals |
| `synthetic_bootstrap_metrics.csv` | All metrics from 2,000 group-bootstrap replicates |

See the [example definition and expected scores](examples/model/README.md). This constructed example demonstrates software behavior, not experimental performance.

For other inputs, the model command exposes four operations:

```sh
selection-landscape-model train --help
selection-landscape-model enumerate --help
selection-landscape-model predict --help
selection-landscape-model evaluate --help
```

The model combines fixed WT-relative single-mutant effects and fits two sigmoid calibration parameters using observed double-mutant labels. Candidate enumeration combines measured substitutions at different positions and excludes a supplied complete list of observed pair identities. Prediction input contains identities and sequences, without labels. The [model documentation](docs/model.md) defines the mathematics, CSV schemas, Python interface and evaluation procedure.

## Run the synthetic count-to-label example

The included count table is synthetic and can be processed without raw reads or an aligner. Use output names that do not already exist:

```sh
python -m selection_landscape.ngs.filter \
  --counts examples/ngs/synthetic_counts.csv \
  --parameter-name strict \
  --length 11 --start-motif EVQ --end-motif VTVSS \
  --output example-label-output/eligibility.csv

python -m selection_landscape.labels \
  --counts examples/ngs/synthetic_counts.csv \
  --eligibility example-label-output/eligibility.csv \
  --profile examples/ngs/selection_profile.json \
  --target examples/ngs/target.json \
  --output example-label-output/labels.csv
```

`eligibility.csv` records sequence-quality and round-coverage decisions. `labels.csv` contains the selection label, rank, support tier, abundances, enrichment ratios and intermediate score components. The example produces one primary and one low-weight label. Selection labels describe sequencing-based selection trajectories; they are not binding affinities or experimental-success probabilities.

The label output is an auditable analysis table. Before model fitting, prepare the explicit WT, single-mutant and selected calibration-double CSVs described in [model.md](docs/model.md); the model does not automatically select rows from the full label table.

## Process paired FASTQ files

Supply your own reads, reference FASTAs and a CSV manifest with `sample_id,round,chain,r1_path,r2_path`. Competitive assignment and molecule reconstruction are separate commands:

```sh
python -m selection_landscape.ngs.split \
  --samples samples.csv \
  --primary-reference primary.fa \
  --competitor-reference competitor.fa \
  --output-dir split-output --threads 1

python -m selection_landscape.ngs.rebuild \
  --samples primary_samples.csv \
  --reference primary.fa \
  --config reconstruction.json \
  --parameter-name strict \
  --output-dir count-output --threads 1
```

After assignment, prepare `primary_samples.csv` to point to each sample's primary-only paired FASTQ files. Define the actual UID anchor, coding-region coordinates and productivity rules in `reconstruction.json`. The provided [configuration template](examples/ngs/rebuild_config_template.json) uses synthetic values and must be adapted to the assay.

Assignment produces paired FASTQ classes and `split_report.csv`. Reconstruction produces `molecule_counts.tsv` and `sample_summary.json`, including read-flow conservation checks. Pass the count table through filtering and label calculation. See [NGS processing](docs/ngs.md) for input contracts, external-tool options, pilot modes and verification boundaries.

## Tests and documentation

```sh
python -m pytest -q
```

Tests cover label calculations, sequence identity, calibrated prediction, evaluation, artifact integrity and file interfaces. The synthetic examples can also be run independently of the test suite. Save input hashes, configurations, tool versions, model checksums and output files for scientific runs. Reuse a frozen model artifact when exact established predictions must be preserved; refitting under different numerical environments can reorder nearly tied candidates.

- [Model definition, schemas and evaluation](docs/model.md)
- [Selection-score calculation](docs/labels.md)
- [Paired FASTQ and molecule reconstruction](docs/ngs.md)
- [Synthetic model inputs and expected results](examples/model/README.md)
- [Contributing](CONTRIBUTING.md)

## License and citation

Code and included synthetic examples are distributed under the [MIT License](LICENSE). Dependencies retain their respective licenses. Researchers supply their own data and are responsible for its permitted use and distribution.

Please cite the software version and repository URL using the metadata in [CITATION.cff](CITATION.cff). The repository location is [SPYfighting/selection-landscape](https://github.com/SPYfighting/selection-landscape). No paper citation or archival DOI is assigned in this source snapshot.
