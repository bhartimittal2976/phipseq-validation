#!/usr/bin/env bash
set -euo pipefail

# Resolve the project directory regardless of where this script is launched.
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
data_dir="${project_dir}/data"
results_dir="${project_dir}/results"

mkdir -p "${results_dir}"

echo "[1/6] Finding overlaps among displayed peptide clones"
phipkit-blast "${data_dir}/annotations.csv" \
  --out "${results_dir}/blast.all_against_all.csv"

echo "[2/6] Mapping displayed peptides to complete reference proteins"
phipkit-blast "${data_dir}/annotations.csv" \
  --reference "${data_dir}/proteins.fasta" \
  --out "${results_dir}/blast.reference.csv"

echo "[3/6] Calculating enrichment relative to bead-only controls"
phipkit-score "${data_dir}/counts.csv" \
  --out "${results_dir}/scores.csv"

echo "[4/6] Calling significant overlapping-peptide pairs"
phipkit-call-hits \
  "${results_dir}/blast.all_against_all.csv" \
  "${results_dir}/scores.csv" \
  --counts "${data_dir}/counts.csv" \
  --fdr 0.1 \
  --out "${results_dir}/hits.csv"

echo "[5/6] Consolidating peptide hits into antigen regions"
phipkit-call-antigens \
  "${results_dir}/blast.reference.csv" \
  "${results_dir}/hits.csv" \
  --out "${results_dir}/antigens.csv"

echo "[6/6] Creating the antigen-region PDF"
phipkit-plot-antigens \
  "${results_dir}/blast.reference.csv" \
  "${results_dir}/hits.csv" \
  "${results_dir}/antigens.csv" \
  --out "${results_dir}/antigens.pdf"

python "${project_dir}/scripts/validate_results.py"
