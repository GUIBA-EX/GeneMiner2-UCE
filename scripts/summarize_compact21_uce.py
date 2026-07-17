#!/usr/bin/env python3
"""Summarize a compact GeneMiner2 UCE validation run by sample and population.

This script is deliberately read-only with respect to the run directory.  It
collects the per-sample UCE assembly/rescue tables written by GeneMiner2 and
creates transparent TSV summaries for the 21-sample validation panel.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def truth(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "y"}


def read_populations(path: Path) -> dict[str, str]:
    populations: dict[str, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            populations[row["sample_id"]] = row["paper_population"]
    return populations


def sample_id_from_dir(path: Path, populations: dict[str, str]) -> str | None:
    # Formal GeneMiner2 directories have a numeric prefix, e.g. 1_Dk40.
    suffix = path.name.split("_", 1)[-1].lower()
    for sample in populations:
        if sample.lower() == suffix:
            return sample
    return None


def write_tsv(path: Path, header: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--populations", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    populations = read_populations(args.populations)
    sample_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []

    for summary in sorted(args.run_dir.glob("*/uce_assembly_summary.csv")):
        sample = sample_id_from_dir(summary.parent, populations)
        if sample is None:
            continue
        rows = list(csv.DictReader(summary.open(newline="")))
        accepted = sum(truth(row.get("accepted")) for row in rows)
        status_counts = Counter((row.get("status") or "").strip() for row in rows)
        low_quality = sum(truth(row.get("low_quality")) for row in rows)
        sample_rows.append(
            {
                "sample_id": sample,
                "paper_population": populations[sample],
                "assembly_rows": len(rows),
                "accepted_loci": accepted,
                "accepted_fraction": f"{accepted / len(rows):.6f}" if rows else "",
                "low_quality_loci": low_quality,
                "rescue_summary_present": (summary.parent / "uce_rescue_summary.csv").is_file(),
            }
        )
        for status, count in sorted(status_counts.items()):
            status_rows.append(
                {
                    "sample_id": sample,
                    "paper_population": populations[sample],
                    "status": status,
                    "loci": count,
                }
            )

    sample_rows.sort(key=lambda row: (str(row["paper_population"]), str(row["sample_id"])))
    population_totals: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        population_totals[str(row["paper_population"])].append(row)
    population_rows = []
    for population, rows in sorted(population_totals.items()):
        accepted = [int(row["accepted_loci"]) for row in rows]
        population_rows.append(
            {
                "paper_population": population,
                "n_samples": len(rows),
                "mean_accepted_loci": f"{sum(accepted) / len(accepted):.2f}",
                "min_accepted_loci": min(accepted),
                "max_accepted_loci": max(accepted),
            }
        )

    write_tsv(
        args.output_dir / "assembly_by_sample.tsv",
        [
            "sample_id", "paper_population", "assembly_rows", "accepted_loci",
            "accepted_fraction", "low_quality_loci", "rescue_summary_present",
        ],
        sample_rows,
    )
    write_tsv(
        args.output_dir / "assembly_status_by_sample.tsv",
        ["sample_id", "paper_population", "status", "loci"],
        status_rows,
    )
    write_tsv(
        args.output_dir / "assembly_by_population.tsv",
        ["paper_population", "n_samples", "mean_accepted_loci", "min_accepted_loci", "max_accepted_loci"],
        population_rows,
    )


if __name__ == "__main__":
    main()
