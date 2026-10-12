# Paired FASTQ to UMI molecule counts

This module processes a three-round, two-chain selection assay. It has no model dependency. Raw reads and reference sequences are supplied by the researcher and are not included in this repository.

## Inputs and tools

- Paired gzip-compressed FASTQ files with matching read identifiers and complete four-line records.
- One primary and one competing single-record FASTA for competitive assignment; a primary FASTA for molecule reconstruction.
- A CSV sample manifest with `sample_id,round,chain,r1_path,r2_path`. Paths may be absolute or relative to the manifest. Sample identifiers must be unique.
- Bowtie2 (verified with 2.5.3), samtools (verified with 1.6), Python with NumPy, pysam and Biopython. Specify the executable paths through CLI options if they are not on `PATH`.

The split stage aligns each read pair against a combined reference with Bowtie2 `--very-sensitive-local --reorder --no-mixed --no-discordant -k 2`. For each reference, it sums the best alignment score from each mapped mate. A pair is assigned to the best reference when its summed score is at least 60 and, when an alternative alignment is available, exceeds the other reference by at least 12; otherwise it is ambiguous or unassigned. The thresholds are command-line parameters. All four classes are written as paired FASTQ files, and `split_report.csv` records pair conservation.

```bash
python -m selection_landscape.ngs.split \
  --samples samples.csv \
  --primary-reference primary.fa \
  --competitor-reference competitor.fa \
  --output-dir split_run \
  --threads 1
```

For molecule reconstruction, prepare a second manifest whose `r1_path` and `r2_path` point to the `primary-only` files for each sample. The JSON template in `examples/ngs/rebuild_config_template.json` shows the required fields; **its synthetic motif and coordinates are illustrative, not assay defaults**. Record the actual reference, UID anchor, ORF boundaries, parameter set and tool versions for each study. The output directory must be new or empty; the public CLI never silently resumes old alignments. It rejects truncated FASTQ records, unequal mate counts, mismatched identifiers, and unequal sequence/quality lengths. Terminal `/1` and `/2` mate suffixes are normalized consistently with Bowtie2. Duplicate normalized identifiers are rejected; this check keeps a set of read names in memory, so its memory use increases with input pair count.

```bash
python -m selection_landscape.ngs.rebuild \
  --samples primary_samples.csv \
  --reference primary.fa \
  --config reconstruction.json \
  --parameter-name strict \
  --output-dir count_run \
  --threads 1
```

The reconstruction preserves the assay's terminal-anchor UID search, paired alignment, CIGAR-aware ORF vote, per-UID consensus, translation and productive-sequence checks. `molecule_counts.tsv` contains sample, round, chain, nucleotide and amino-acid sequences, read and molecule counts, and productivity. `sample_summary.json` records UID recovery, alignment flow, released molecules and conservation checks. Diagnostics go to standard error or log files; standard output is a compact JSON result.

The output of this stage enters `selection_landscape.ngs.filter`, followed by the independent label module described in [labels.md](labels.md). `--max-pairs` on the split CLI and `--max-reads` on the reconstruction CLI are explicitly sized pilot modes; omit both for a complete analysis. Do not combine pilot counts with full-run labels.

## Verification boundary

The public split and reconstruction code passed a limited real-read smoke test, including pair conservation and alignment-to-clone read conservation. This is an executable check, not evidence that historical full-FASTQ outputs are byte-identical: the exact source snapshots of the historical full reconstruction run are unavailable. The frozen full count table and downstream filter/label outputs were checked separately as described in [labels.md](labels.md). Detailed sample results remain in the private verification record.
