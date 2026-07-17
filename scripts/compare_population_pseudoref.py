#!/usr/bin/env python3
"""Summarize a blind comparison of a GeneMiner2 population reference and Ref2.

The external FASTA is comparison-only.  The PAF supplied to this program must
be generated after the population reference is finalized, for example:

  minimap2 -x asm5 --secondary=no Ref2.fasta population_reference.fasta > out.paf
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name: str | None = None
    sequence_length = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = sequence_length
                name = line[1:].split()[0]
                sequence_length = 0
            else:
                sequence_length += len(line)
    if name is not None:
        lengths[name] = sequence_length
    return lengths


def normalized_locus(name: str) -> str:
    # Probe-derived GM2 loci use v1__uce-123; Ref2 uses v1_uce-123.
    return name.replace("__", "_")


def best_alignments(paf: Path) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    with paf.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, block, mapq = fields[:12]
            candidate = {
                "query": query,
                "query_length": int(qlen),
                "query_aligned_bases": int(qend) - int(qstart),
                "target": target,
                "target_length": int(tlen),
                "target_aligned_bases": int(tend) - int(tstart),
                "matches": int(matches),
                "block_length": int(block),
                "mapq": int(mapq),
                "strand": strand,
            }
            previous = best.get(query)
            if previous is None or (candidate["matches"], candidate["block_length"], candidate["mapq"]) > (
                previous["matches"], previous["block_length"], previous["mapq"],
            ):
                best[query] = candidate
    return best


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--population-reference", type=Path, required=True)
    parser.add_argument("--paper-reference", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    generated = fasta_lengths(args.population_reference)
    paper = fasta_lengths(args.paper_reference)
    paper_by_locus = {normalized_locus(name): name for name in paper}
    paf_rows = best_alignments(args.paf)

    rows: list[dict[str, object]] = []
    for query, qlength in sorted(generated.items()):
        alignment = paf_rows.get(query)
        row: dict[str, object] = {
            "population_locus": query,
            "normalized_locus": normalized_locus(query),
            "population_length": qlength,
            "paper_same_locus_present": normalized_locus(query) in paper_by_locus,
            "paper_same_locus": paper_by_locus.get(normalized_locus(query), ""),
            "best_paper_locus": "",
            "best_query_coverage": "",
            "best_identity": "",
            "best_mapq": "",
            "best_is_same_locus": "",
        }
        if alignment:
            row.update(
                {
                    "best_paper_locus": alignment["target"],
                    "best_query_coverage": f"{alignment['query_aligned_bases'] / qlength:.6f}",
                    "best_identity": f"{alignment['matches'] / alignment['block_length']:.6f}" if alignment["block_length"] else "",
                    "best_mapq": alignment["mapq"],
                    "best_is_same_locus": normalized_locus(query) == normalized_locus(str(alignment["target"])),
                }
            )
        rows.append(row)

    aligned = [row for row in rows if row["best_paper_locus"]]
    high_confidence = [
        row for row in aligned
        if float(str(row["best_query_coverage"])) >= 0.8 and float(str(row["best_identity"])) >= 0.9
    ]
    same_locus = [row for row in aligned if row["best_is_same_locus"] is True]
    metrics = [
        ("population_loci", len(generated)),
        ("population_total_bases", sum(generated.values())),
        ("paper_ref2_loci", len(paper)),
        ("paper_ref2_total_bases", sum(paper.values())),
        ("normalized_locus_intersection", len(set(map(normalized_locus, generated)) & set(paper_by_locus))),
        ("population_loci_with_any_paf_alignment", len(aligned)),
        ("population_loci_high_confidence_best_alignment", len(high_confidence)),
        ("population_loci_best_alignment_same_normalized_locus", len(same_locus)),
    ]
    with (args.output_dir / "pseudoref_vs_ref2_metrics.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in metrics:
            handle.write(f"{metric}\t{value}\n")
    write_tsv(
        args.output_dir / "pseudoref_vs_ref2_locus_alignment.tsv",
        list(rows[0]) if rows else ["population_locus"],
        rows,
    )


if __name__ == "__main__":
    main()
