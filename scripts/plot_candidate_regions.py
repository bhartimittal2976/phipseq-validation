#!/usr/bin/env python3
"""Render collapsed PhIP-Seq candidate regions as a dependency-free SVG."""

import argparse
import csv
import gzip
import html
import math
from pathlib import Path


NAVY = "#12304A"
BLUE = "#197A96"
TEAL = "#2AA495"
CORAL = "#F27652"
PALE = "#EAF5F7"
GRID = "#CADCE2"
TEXT = "#163047"


def load_matrix(path):
    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle)
        samples = next(reader)[1:]
        values = {row[0]: [float(x) for x in row[1:]] for row in reader}
    return samples, values


def esc(value):
    return html.escape(str(value))


def text(x, y, value, size=18, weight=400, fill=TEXT, anchor="start"):
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}">{esc(value)}</text>'
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wide-dir", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    samples, counts = load_matrix(args.wide_dir / "data_counts.csv.gz")
    _, enrichment = load_matrix(args.wide_dir / "data_enrichment.csv.gz")
    _, zscores = load_matrix(args.wide_dir / "data_zscore.csv.gz")

    with gzip.open(
        args.wide_dir / "data_peptide_annotation_table.csv.gz", "rt"
    ) as handle:
        annotations = {row["peptide_id"]: row for row in csv.DictReader(handle)}
    with args.candidates.open() as handle:
        candidates = list(csv.DictReader(handle))

    width = 1800
    panel_height = 260
    height = 150 + panel_height * len(candidates) + 70
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#F7FBFC"/>',
        f'<rect x="35" y="28" width="1730" height="90" rx="22" fill="{NAVY}"/>',
        text(70, 77, "PhIP-Seq candidate regions from FASTQ data", 34, 700, "white"),
        text(70, 104, "Adjacent-tile sensitivity analysis: Z ≥ 15, count ≥ 5", 18, 400, "#CFEAF0"),
    ]

    for index, candidate in enumerate(candidates):
        sample_id = candidate["sample_id"]
        sample_index = samples.index(sample_id)
        peptide_ids = candidate["peptide_ids"].split(";")
        by_start = {}
        for peptide_id in peptide_ids:
            start = int(annotations[peptide_id]["Prot_Start"])
            if start not in by_start or zscores[peptide_id][sample_index] > zscores[by_start[start]][sample_index]:
                by_start[start] = peptide_id
        starts = sorted(by_start)
        representative = [by_start[starts[0]], by_start[starts[-1]]]

        top = 145 + index * panel_height
        svg.extend(
            [
                f'<rect x="35" y="{top}" width="1730" height="230" rx="18" fill="white" stroke="{GRID}" stroke-width="2"/>',
                text(65, top + 40, f"Sample {sample_id} · {candidate['virus']}", 24, 700, NAVY),
                text(65, top + 70, candidate["supporting_annotations"].replace(";", " · "), 16, 400, BLUE),
                text(65, top + 105, f"Shared region {candidate['candidate_region_start']}–{candidate['candidate_region_end']}", 19, 700, CORAL),
                f'<rect x="65" y="{top + 120}" width="470" height="43" rx="10" fill="{PALE}"/>',
                text(300, top + 149, candidate["shared_sequence"], 23, 700, NAVY, "middle"),
            ]
        )

        for tile_index, peptide_id in enumerate(representative):
            annotation = annotations[peptide_id]
            start = int(annotation["Prot_Start"])
            end = start + len(annotation["Prot"]) - 1
            x = 610 + tile_index * 560
            count = counts[peptide_id][sample_index]
            fold = enrichment[peptide_id][sample_index]
            zscore = zscores[peptide_id][sample_index]
            svg.extend(
                [
                    text(x, top + 36, f"Tile {tile_index + 1}: {start}–{end}", 19, 700, TEXT),
                    f'<rect x="{x}" y="{top + 55}" width="500" height="34" rx="8" fill="{BLUE if tile_index == 0 else TEAL}" opacity="0.92"/>',
                    text(x + 250, top + 79, annotation["Prot"], 15, 700, "white", "middle"),
                    text(x, top + 120, f"Raw count: {count:.0f}", 18, 700, TEXT),
                    text(x, top + 150, f"Library fold enrichment: {fold:.2f}", 18, 700, TEXT),
                    text(x, top + 180, f"Bead-background Z-score: {zscore:.2f}", 18, 700, TEXT),
                ]
            )

        # A compact visual encoding of the minimum supporting strength.
        min_z = float(candidate["min_zscore"])
        bar_width = min(320, 8 * min_z)
        svg.extend(
            [
                text(65, top + 205, "Minimum supporting Z", 15, 700, TEXT),
                f'<rect x="245" y="{top + 188}" width="320" height="20" rx="10" fill="{PALE}"/>',
                f'<rect x="245" y="{top + 188}" width="{bar_width:.1f}" height="20" rx="10" fill="{CORAL}"/>',
                text(575, top + 205, f"{min_z:.2f}", 15, 700, CORAL),
            ]
        )

    svg.extend(
        [
            text(40, height - 30, "Exploratory computational calls from a shallow bundled test dataset; not biological or clinical validation.", 16, 700, NAVY),
            "</svg>",
        ]
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(svg))
    print(f"Wrote {len(candidates)} candidate panels to {args.out}")


if __name__ == "__main__":
    main()
