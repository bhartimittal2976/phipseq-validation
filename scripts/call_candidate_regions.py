#!/usr/bin/env python3
"""Call adjacent, overlapping PhIP-Seq candidate regions from PhIP-Flow output."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def load_matrix(path):
    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle)
        samples = next(reader)[1:]
        matrix = {
            row[0]: [float(value) for value in row[1:]]
            for row in reader
        }
    return samples, matrix


def longest_suffix_prefix(left, right):
    for length in range(min(len(left), len(right)), 0, -1):
        if left[-length:] == right[:length]:
            return left[-length:]
    return ""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find adjacent enriched peptide tiles and collapse redundant regions."
    )
    parser.add_argument("--wide-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-z", type=float, default=10.0)
    parser.add_argument("--min-count", type=float, default=3.0)
    parser.add_argument("--tile-step", type=int, default=20)
    return parser.parse_args()


def main():
    args = parse_args()
    samples, counts = load_matrix(args.wide_dir / "data_counts.csv.gz")
    enrichment_samples, enrichment = load_matrix(
        args.wide_dir / "data_enrichment.csv.gz"
    )
    zscore_samples, zscores = load_matrix(args.wide_dir / "data_zscore.csv.gz")
    if samples != enrichment_samples or samples != zscore_samples:
        raise ValueError("Matrix sample columns do not match")

    with gzip.open(
        args.wide_dir / "data_peptide_annotation_table.csv.gz", "rt"
    ) as handle:
        annotations = {row["peptide_id"]: row for row in csv.DictReader(handle)}

    with gzip.open(
        args.wide_dir / "data_sample_annotation_table.csv.gz", "rt"
    ) as handle:
        sample_annotations = {row["sample_id"]: row for row in csv.DictReader(handle)}

    output = {}
    for sample_index, sample_id in enumerate(samples):
        sample_meta = sample_annotations[sample_id]
        if sample_meta["control_status"] != "empirical":
            continue

        proteins = defaultdict(dict)
        for peptide_id, values in zscores.items():
            zscore = values[sample_index]
            count = counts[peptide_id][sample_index]
            if zscore < args.min_z or count < args.min_count:
                continue

            annotation = annotations[peptide_id]
            protein_key = (annotation["Virus"], annotation["Full_name"])
            start = int(annotation["Prot_Start"])
            previous = proteins[protein_key].get(start)
            if previous is None or zscore > previous[0]:
                proteins[protein_key][start] = (zscore, peptide_id)

        for (virus, full_name), starts in proteins.items():
            ordered = sorted(starts)
            for start1, start2 in zip(ordered, ordered[1:]):
                if start2 - start1 != args.tile_step:
                    continue

                _, peptide1 = starts[start1]
                _, peptide2 = starts[start2]
                annotation1 = annotations[peptide1]
                annotation2 = annotations[peptide2]
                overlap = longest_suffix_prefix(annotation1["Prot"], annotation2["Prot"])
                if not overlap:
                    continue

                region_start = start2
                region_end = region_start + len(overlap) - 1
                collapse_key = (sample_id, virus, overlap)
                record = output.get(collapse_key)
                if record is None:
                    record = {
                        "sample_id": sample_id,
                        "library_batch": sample_meta["library_batch"],
                        "virus": virus,
                        "candidate_region_start": region_start,
                        "candidate_region_end": region_end,
                        "shared_sequence": overlap,
                        "overlap_length": len(overlap),
                        "supporting_annotations": set(),
                        "peptide_ids": set(),
                        "tile_starts": set(),
                        "counts": [],
                        "zscores": [],
                        "library_fold_enrichment": [],
                    }
                    output[collapse_key] = record

                record["supporting_annotations"].add(full_name)
                record["peptide_ids"].update((peptide1, peptide2))
                record["tile_starts"].update((start1, start2))
                record["counts"].extend(
                    (counts[peptide1][sample_index], counts[peptide2][sample_index])
                )
                record["zscores"].extend(
                    (zscores[peptide1][sample_index], zscores[peptide2][sample_index])
                )
                record["library_fold_enrichment"].extend(
                    (
                        enrichment[peptide1][sample_index],
                        enrichment[peptide2][sample_index],
                    )
                )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample_id",
        "library_batch",
        "virus",
        "candidate_region_start",
        "candidate_region_end",
        "shared_sequence",
        "overlap_length",
        "supporting_annotations",
        "peptide_ids",
        "tile_starts",
        "max_count",
        "min_zscore",
        "max_zscore",
        "min_library_fold_enrichment",
        "max_library_fold_enrichment",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in sorted(
            output.values(),
            key=lambda item: (int(item["sample_id"]), item["virus"], item["candidate_region_start"]),
        ):
            writer.writerow(
                {
                    "sample_id": record["sample_id"],
                    "library_batch": record["library_batch"],
                    "virus": record["virus"],
                    "candidate_region_start": record["candidate_region_start"],
                    "candidate_region_end": record["candidate_region_end"],
                    "shared_sequence": record["shared_sequence"],
                    "overlap_length": record["overlap_length"],
                    "supporting_annotations": ";".join(sorted(record["supporting_annotations"])),
                    "peptide_ids": ";".join(sorted(record["peptide_ids"], key=int)),
                    "tile_starts": ";".join(map(str, sorted(record["tile_starts"]))),
                    "max_count": f"{max(record['counts']):.0f}",
                    "min_zscore": f"{min(record['zscores']):.4f}",
                    "max_zscore": f"{max(record['zscores']):.4f}",
                    "min_library_fold_enrichment": f"{min(record['library_fold_enrichment']):.4f}",
                    "max_library_fold_enrichment": f"{max(record['library_fold_enrichment']):.4f}",
                }
            )

    print(f"Wrote {len(output)} collapsed candidate regions to {args.out}")


if __name__ == "__main__":
    main()
