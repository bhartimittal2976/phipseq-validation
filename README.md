# PhIP-Seq FASTQ validation

This repository demonstrates a reproducible FASTQ-to-candidate-region
PhIP-Seq analysis using the shallow Pan-CoV example dataset bundled with
[PhIP-Flow](https://github.com/matsengrp/phip-flow).

This is an exploratory computational validation. It is not a biological or
clinical validation study.

## Dataset

The analysis uses the eight example FASTQ files supplied with PhIP-Flow
V1.12:

- two input-library controls;
- two bead-only controls;
- four antibody-containing experimental samples.

The peptide library contains 10,047 Pan-CoV and control peptide records.
PhIP-Flow supplies the example sample table, FASTQ files, and peptide table.

## Workflow

The analysis performs:

1. Sample and peptide-table validation.
2. FASTQ alignment to the designed oligonucleotide library with Bowtie.
3. SAM-to-peptide count conversion and alignment QC.
4. Counts-per-million and size-factor normalization.
5. Fold enrichment relative to input-library controls.
6. Z-score modelling relative to bead-only controls.
7. Adjacent overlapping-tile detection.
8. Candidate-region visualization.

The final Nextflow execution completed successfully with no failed, aborted,
pending, or running tasks.

## Requirements

- macOS or Linux
- Docker
- Java 17 or later
- Nextflow 25.10.4
- Python 3

PhIP-Flow V1.12 was run with Nextflow 25.10.4 because Nextflow 26's stricter
parser rejects legacy top-level syntax in that workflow release.

## Install the project-local Nextflow version

From the repository root:

```bash
mkdir -p tools
cd tools
curl -s https://get.nextflow.io | NXF_VER=25.10.4 bash
mv nextflow nextflow-25.10.4
chmod +x nextflow-25.10.4
cd ..
```

The downloaded executable is ignored by Git.

## Run PhIP-Flow

```bash
mkdir -p run-pan-cov
cd run-pan-cov

../tools/nextflow-25.10.4 run matsengrp/phip-flow \
  -r V1.12 \
  -profile docker \
  --run_edgeR false \
  --run_BEER false \
  --run_cpm_enr_workflow true \
  --run_zscore_fit_predict true \
  -resume
```

This processes the bundled FASTQs through alignment, peptide counting,
normalization, fold enrichment, and bead-background Z-score modelling.
edgeR and BEER are optional analyses and were not used for this validation.

## Call adjacent candidate regions

The library uses 39-amino-acid peptide tiles beginning 20 positions apart.
Adjacent tiles therefore share a 19-amino-acid sequence.

From the repository root, generate the strict exploratory calls:

```bash
python3 scripts/call_candidate_regions.py \
  --wide-dir run-pan-cov/results/wide_data \
  --out run-pan-cov/results/candidate_regions-z15.csv \
  --min-z 15 \
  --min-count 5 \
  --tile-step 20
```

Generate the SVG report:

```bash
python3 scripts/plot_candidate_regions.py \
  --wide-dir run-pan-cov/results/wide_data \
  --candidates run-pan-cov/results/candidate_regions-z15.csv \
  --out run-pan-cov/results/candidate_regions-z15.svg
```

Both scripts use only the Python standard library.

## Strict candidate calls

At `Z >= 15` and raw count `>= 5`, the analysis identified three collapsed
adjacent-tile candidate regions:

- Sample 4: NL63 nucleocapsid, shared region 61-79.
- Sample 4: SARS nucleocapsid, shared region 241-259.
- Sample 5: HKU1 ORF8-like protein, shared region 21-39.

Repeated strain annotations and shared 1a/1ab sequences are collapsed so that
redundant library records are not counted as independent discoveries.

## Outputs

The repository retains these compact final outputs:

- `run-pan-cov/results/candidate_regions.csv`
- `run-pan-cov/results/candidate_regions-z15.csv`
- `run-pan-cov/results/candidate_regions-z15.svg`

Regenerable Nextflow work files, logs, wide matrices, binary datasets,
downloaded executables, raw FASTQ files, and machine-specific resolved
configuration are excluded from Git.

## Limitations

The example FASTQs contain only a shallow subset of reads and are intended for
workflow validation. Large fold-enrichment values can result from low peptide
representation in the input-library controls. The candidate regions must not
be interpreted as validated antibody epitopes without independent experimental
evidence.
