#!/usr/bin/env python3
"""Conservative, assembly-free profiling of recruited Symbiodiniaceae ITS2 reads."""

import argparse
import csv
import math
import os
import re
import subprocess
from collections import Counter, defaultdict

UNNAMED = re.compile(r"^\d+_([A-I])$")
NAMED = re.compile(r"^([A-I])(?:\d|$)")
CIGAR_QUERY = re.compile(r"(\d+)([MIDNSHP=X])")

def parse_reference_header(header):
    """Extract the legacy clade and named SymPortal DIV aliases from a header."""
    clades, named = [], []
    for alias in header.split("/"):
        unnamed = UNNAMED.match(alias)
        if unnamed:
            clades.append(unnamed.group(1))
            continue
        match = NAMED.match(alias)
        if match:
            clades.append(match.group(1))
            named.append(alias)
    return (clades[0] if len(set(clades)) == 1 else "", "/".join(sorted(set(named))))

def query_coverage(cigar, sequence):
    if cigar == "*" or sequence == "*" or not sequence:
        return 0.0
    consumed = sum(int(length) for length, op in CIGAR_QUERY.findall(cigar) if op in "MIS=X")
    return consumed / len(sequence)

def classify_fragment(records, ref_info, score_delta, min_coverage):
    records = [record for record in records if record[2] >= min_coverage]
    if not records:
        return "unclassified", "unclassified", "", "", ""
    best_score = max(record[1] for record in records)
    refs = sorted({record[0] for record in records if record[1] >= best_score - score_delta})
    labels = {ref_info[ref][1] for ref in refs if ref_info[ref][1]}
    clades = {ref_info[ref][0] for ref in refs if ref_info[ref][0]}
    if len(labels) == 1 and len(clades) == 1:
        return "div", next(iter(labels)), next(iter(clades)), str(best_score), ";".join(refs)
    if len(clades) == 1:
        clade = next(iter(clades))
        return "clade_only", clade, clade, str(best_score), ";".join(refs)
    return "ambiguous", "ambiguous", "", str(best_score), ";".join(refs)

def clade_support(records, ref_info, score_delta, min_coverage):
    """Allocate one fragment across its plausible clades.

    References are first collapsed to clades, so a clade with more database
    entries does not gain support merely from reference redundancy.
    """
    records = [record for record in records if record[2] >= min_coverage]
    if not records:
        return {}
    best_score = max(record[1] for record in records)
    clade_scores = {}
    for ref, score, _ in records:
        clade = ref_info[ref][0]
        if clade and score >= best_score - score_delta:
            clade_scores[clade] = max(clade_scores.get(clade, float("-inf")), score)
    if not clade_scores:
        return {}
    weights = {clade: math.exp((score - best_score) / 5.0)
               for clade, score in clade_scores.items()}
    total = sum(weights.values())
    return {clade: weight / total for clade, weight in weights.items()}

def passes_specificity_filter(records, ref_info, row, min_coverage, min_score_margin):
    """Require a classified clade to outrank the nearest alternative clade."""
    if row[0] not in {"div", "clade_only"}:
        return False, ""
    winner = row[2]
    clade_scores = {}
    for ref, score, coverage in records:
        clade = ref_info[ref][0]
        if clade and coverage >= min_coverage:
            clade_scores[clade] = max(clade_scores.get(clade, float("-inf")), score)
    winner_score = clade_scores.get(winner)
    if winner_score is None:
        return False, ""
    alternatives = [score for clade, score in clade_scores.items() if clade != winner]
    if not alternatives:
        return True, "NA"
    margin = winner_score - max(alternatives)
    return margin >= min_score_margin, str(margin)

def main():
    parser = argparse.ArgumentParser(description="Classify recruited ITS2 fragments without assembly")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--reads", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--minimap2", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--score-delta", type=int, default=5)
    parser.add_argument("--min-query-coverage", type=float, default=0.85)
    parser.add_argument("--min-clade-score-margin", type=int, default=10, help="Minimum winning-clade AS margin for the specificity filter")
    args = parser.parse_args()
    if args.score_delta < 0:
        parser.error("--score-delta must be non-negative")
    if args.min_clade_score_margin < 0:
        parser.error("--min-clade-score-margin must be non-negative")
    if not 0 < args.min_query_coverage <= 1:
        parser.error("--min-query-coverage must be in (0, 1]")
    ref_info = {}
    with open(args.reference) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                ref_info[name] = parse_reference_header(name)
    command = [args.minimap2, "-ax", "sr", "-t", str(max(args.threads, 1)), "-N", "50", "--secondary=yes", args.reference, *args.reads]
    sam = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE).stdout
    by_read = defaultdict(list)
    for line in sam.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        by_read[fields[0]]
        if fields[2] == "*" or fields[2] not in ref_info:
            continue
        score = next((int(field[5:]) for field in fields[11:] if field.startswith("AS:i:")), None)
        if score is not None:
            by_read[fields[0]].append((fields[2], score, query_coverage(fields[5], fields[9])))
    os.makedirs(args.output, exist_ok=True)
    counts = Counter()
    strict_clade_counts = Counter()
    fractional_clade_counts = Counter()
    passed_counts = Counter()
    with open(os.path.join(args.output, "its2_fragment_assignments.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["fragment", "classification_level", "label", "clade", "best_alignment_score", "candidate_references", "specificity_filter", "runner_up_clade_score_margin"])
        for fragment in sorted(by_read):
            row = classify_fragment(by_read[fragment], ref_info, args.score_delta, args.min_query_coverage)
            counts[row[:3]] += 1
            passed, margin = passes_specificity_filter(by_read[fragment], ref_info, row, args.min_query_coverage, args.min_clade_score_margin)
            if passed:
                passed_counts[row[:3]] += 1
                if row[0] in {"div", "clade_only"}:
                    strict_clade_counts[row[2]] += 1
                for clade, support in clade_support(by_read[fragment], ref_info, args.score_delta, args.min_query_coverage).items():
                    fractional_clade_counts[clade] += support
            writer.writerow([fragment, *row, "pass" if passed else "fail", margin])
    with open(os.path.join(args.output, "its2_profile.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["classification_level", "label", "clade", "fragments_all", "fragments_passing_specificity_filter"])
        for level, label, clade in sorted(set(counts) | set(passed_counts)):
            writer.writerow([level, label, clade, counts[(level, label, clade)], passed_counts[(level, label, clade)]])
    strict_total = sum(strict_clade_counts.values())
    fractional_total = sum(fractional_clade_counts.values())
    with open(os.path.join(args.output, "its2_clade_composition.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["clade", "strict_fragments", "strict_proportion", "fractional_its2_support", "fractional_its2_proportion"])
        for clade in "ABCDEFGHI":
            strict = strict_clade_counts[clade]
            fractional = fractional_clade_counts[clade]
            writer.writerow([clade, strict,
                             f"{strict / strict_total:.8f}" if strict_total else "NA",
                             f"{fractional:.6f}",
                             f"{fractional / fractional_total:.8f}" if fractional_total else "NA"])

if __name__ == "__main__":
    main()
