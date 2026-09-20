# Reproducible PhIP-Seq validation

This repository demonstrates a small, reproducible computational validation of
Phage ImmunoPrecipitation Sequencing (PhIP-Seq) analysis. It uses the synthetic
example dataset distributed with
[`openvax/phipkit`](https://github.com/openvax/phipkit).

![PhIP-Seq computational validation infographic](assets/phipseq-validation-linkedin-infographic.png)

## Result

The complete workflow reproduced the supplied reference results:

- 811/811 peptide enrichment-score rows matched
- 67/67 expected overlapping peptide-pair hits were recovered
- 11/11 expected antigen-region calls were recovered
- Four nonredundant antigen regions were identified across two synthetic samples
- A three-page antigen visualization was generated for ORF7b, Spike and Envelope

## Workflow

```text
Peptide annotations + clone-count matrix + reference proteins
                            |
                            v
             all-against-all peptide BLAST
                            |
                            v
        enrichment versus bead-only negative controls
                            |
                            v
           significant overlapping-peptide pairs
                            |
                            v
              reference-protein region calls
                            |
                            v
             visualization + automated validation
```

This repository begins with an analysis-ready clone-count matrix. In a complete
experiment, phage-library construction, antibody immunoprecipitation, PCR,
sequencing, read alignment and clone counting occur before this stage.

## Requirements

- macOS or Linux
- [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- Git

The workflow was initially validated on Apple Silicon macOS and uses
platform-resolved packages from conda-forge and bioconda, allowing the same
environment recipe to be used on Linux. The plotting code in `phipkit` 0.0.2
requires `pandas<2`; the provided environment records this compatibility
constraint.

## Installation

```bash
micromamba create -f environment.yml
micromamba activate phipseq
```

## Reproduce the analysis

```bash
bash scripts/run_analysis.sh
```

The script performs six analysis stages, regenerates all files under
`results/`, and runs the automated reference comparison. A successful run ends
with:

```text
PASS: 811 enrichment-score rows match
PASS: 67/67 expected peptide-pair hits recovered
PASS: 11/11 expected antigen-region calls recovered
VALIDATION SUCCESSFUL
```

To validate existing results without rerunning the analysis:

```bash
python scripts/validate_results.py
```

## Main outputs

- `results/scores.csv`: enrichment scores relative to bead-only controls
- `results/hits.csv`: significant overlapping peptide-pair hits
- `results/antigens.csv`: consolidated antigen-region calls
- `results/antigens.pdf`: three-page antigen-region visualization

## Interpretation

The nonredundant synthetic calls are:

| Sample | Antigen | Region | Sequence | Supporting clones |
|---|---|---:|---|---:|
| `sample_a` | ORF7b | 20-36 | `VLIMLIIFWFSLELQD` | 2 |
| `sample_a` | Spike | 740-748 | `YICGDSTE` | 8 |
| `sample_b` | Envelope | 44-52 | `NIVNVSLV` | 13 |
| `sample_b` | Envelope | 8-20 | `TGTLIVNSVLLF` | 3 |

## Scope and limitations

This is a computational reproducibility exercise using synthetic data. It does
not validate phage-library construction, immunoprecipitation, sequencing,
patient antibodies, binding affinity, clinical significance or biological
specificity. PhIP-Seq primarily captures linear peptide epitopes and generally
does not reproduce conformational epitopes or post-translational modifications.

## Data and software provenance

The input and expected-output files under `data/` are from the Apache-2.0
licensed [`openvax/phipkit` example dataset](https://github.com/openvax/phipkit/tree/master/example-data).
The analysis uses revision
`267a9b0bb3b79b3ff416e07ce390e43b84a02be4` of `phipkit`.

## References

- Larman HB et al. *Autoantigen discovery with a synthetic human peptidome.*
  Nature Biotechnology (2011). <https://doi.org/10.1038/nbt.1856>
- Mohan D et al. *PhIP-Seq characterization of serum antibodies using
  oligonucleotide-encoded peptidomes.* Nature Protocols (2018).
  <https://doi.org/10.1038/s41596-018-0025-6>
- Galloway J et al. *phippery: a software suite for PhIP-Seq data analysis.*
  Bioinformatics Advances (2023). <https://doi.org/10.1093/bioadv/vbad104>
