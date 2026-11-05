#!/usr/bin/env python3
"""Set every sample's trait profile side by side: one row per sample, one column per
trait, so a batch of strains can be compared at a glance.

Runs as rule traits_summary (shared/90_report.smk), once per run, after rule
trait_profile (shared/40_annotation.smk) has written every sample's
{sample}_trait_profile.tsv. Nothing is decided again here: each cell is copied from
that sample's profile, so the reasons stay in the per-sample files.

What it reads
-------------
  --profiles     every sample's {sample}_trait_profile.tsv (one row per trait, in
                 the order of trait_rules_v1.tsv; columns sample, trait_id,
                 decision and more). The sample name is read from the sample
                 column, not from the file name.

What it writes
--------------
  --out-summary  09.report/traits_summary.tsv: a sample column, then one column per
                 trait_id in the order of trait_rules_v1.tsv; each cell is called,
                 partial or not called. One row per sample, in the order the
                 profiles were given: the rule passes them in the order of
                 SAMPLES, the sample names BacFlux read from the input file
                 names, sorted by name.
  --out-mqc      09.report/traits_summary_mqc.tsv: the same grid as numbers, 2
                 called, 1 partial, 0 not called, under '#' lines that tell
                 MultiQC how to draw it (MULTIQC_HEADER below). Rule multiqc reads
                 it as custom content (MultiQC's name for a table made by another
                 program, which it finds by the _mqc in the file name) and draws a
                 heatmap: MultiQC colours a heatmap only from numbers, while a
                 person reading the table wants the words, hence two files.

All profiles must list the same traits in the same order, which they do when one
run of the workflow wrote them; anything else stops the run, because a grid of
mismatched columns would put calls under the wrong trait.

Standard library only (Python 3.12 or later). Rule traits_summary runs it in
workflow/envs/python.yaml.
"""

import argparse
import os
import sys


# The three decisions trait_profile.py writes, and the number each becomes in the
# MultiQC heatmap: the stronger the call, the darker the cell.
DECISION_SCORES = {"called": 2, "partial": 1, "not called": 0}

# The profile columns this script reads, found by name.
PROFILE_COLUMNS_NEEDED = ["sample", "trait_id", "decision"]

# The '#' lines at the top of the MultiQC file, read by MultiQC as YAML. The id is
# one no MultiQC module uses. xcats_samples is off because the columns are traits,
# not samples; the rows are samples. The colour scale is fixed from 0 to 2 so a
# partial call always has the middle colour: in a batch with no called trait, the
# scale would otherwise stretch to 0-1 and draw partial calls in the colour of a
# call. MultiQC offers to cluster rows and columns, but shows the table order first.
# legend is off: MultiQC's default colour bar runs from 0 to 2 and suggests a
# continuous scale for what are three categories. MultiQC 1.33 writes the value
# in each cell when the grid has fewer than 400 cells (with about fifty traits, up
# to eight samples); for a larger batch, the description names the colour of each
# call, and the hover text of a cell gives its value. The colours named there are
# those of colstops, the last line.
MULTIQC_HEADER = [
    "# id: 'bacflux_traits'",
    "# section_name: 'Trait profile'",
    "# description: 'One row per strain, one column per trait: 2 = called (dark "
    "blue), 1 = partial (light blue), 0 = not called (light grey). A call means the "
    "set of genes behind the trait is present, not that the strain shows it. The "
    "genes and the reason behind each call are in the sample folders under "
    "04.annotation/traits/.'",
    "# plot_type: 'heatmap'",
    "# pconfig:",
    "#     id: 'bacflux_traits_heatmap'",
    "#     title: 'BacFlux: trait profile'",
    "#     xlab: 'Trait'",
    "#     ylab: 'Sample'",
    "#     zlab: 'Call (2 called, 1 partial, 0 not called)'",
    "#     min: 0",
    "#     max: 2",
    "#     tt_decimals: 0",
    "#     square: false",
    "#     xcats_samples: false",
    "#     ycats_samples: true",
    "#     legend: false",
    "#     colstops: [[0, '#f0f0f0'], [0.5, '#9ecae1'], [1, '#08519c']]",
]


def read_profile(path):
    """Read one {sample}_trait_profile.tsv.

    Takes in: a profile written by trait_profile.py.
    Returns:  (sample, [(trait_id, decision), ...]) in the file's order.
    Stops the run when the file has no rows, lacks a needed column, has a row
    with more or fewer fields than the header, names more than one sample, or
    holds a decision other than the three known ones: each means the file did
    not come from trait_profile.py as expected.
    """
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n").rstrip("\r") for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"{path}: empty file, not a trait profile.")
    header = lines[0].split("\t")
    missing = [name for name in PROFILE_COLUMNS_NEEDED if name not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}; not a "
                         "trait profile.")
    sample_index = header.index("sample")
    trait_index = header.index("trait_id")
    decision_index = header.index("decision")

    samples = set()
    calls = []
    for number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"{path}: line {number} has {len(fields)} fields, the "
                             f"header has {len(header)}.")
        decision = fields[decision_index]
        if decision not in DECISION_SCORES:
            raise ValueError(f"{path}: line {number}: unknown decision '{decision}'.")
        samples.add(fields[sample_index])
        calls.append((fields[trait_index], decision))
    if not calls:
        raise ValueError(f"{path}: no trait rows.")
    if len(samples) != 1:
        raise ValueError(f"{path}: names {len(samples)} samples "
                         f"({', '.join(sorted(samples))}); a profile holds one.")
    return samples.pop(), calls


def build_grid(paths):
    """Read every profile and line them up under one list of traits.

    Takes in: the profile paths, in the order the rule passes them (SAMPLES,
              sorted by name).
    Returns:  (trait_ids, [(sample, [decision per trait]), ...]).
    Stops the run when two profiles name the same sample, or list different
    traits or the same traits in a different order.
    """
    trait_ids = None
    first_path = None
    grid = []
    seen = {}
    for path in paths:
        sample, calls = read_profile(path)
        if sample in seen:
            raise ValueError(f"{path}: sample {sample} is also in {seen[sample]}.")
        seen[sample] = path
        these_traits = [trait_id for trait_id, _decision in calls]
        if trait_ids is None:
            trait_ids = these_traits
            first_path = path
        elif these_traits != trait_ids:
            raise ValueError(f"{path}: its traits differ from those of {first_path} "
                             "(different trait table?); the grid would mix up "
                             "columns.")
        grid.append((sample, [decision for _trait_id, decision in calls]))
    return trait_ids, grid


def write_lines(path, lines):
    """Write text lines to path, each ending in a newline, creating the folder if
    needed. main writes both outputs with it."""
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("".join(line + "\n" for line in lines))


def main(argv=None):
    """Command-line entry point used by rule traits_summary.

    Reads every profile into one grid (build_grid), writes the grid in words and
    the same grid as numbers under MULTIQC_HEADER, and prints one line for the
    log.
    """
    parser = argparse.ArgumentParser(
        description=("Set every sample's trait profile side by side: one row per "
                     "sample, one column per trait.")
    )
    parser.add_argument("--profiles", nargs="+", required=True,
                        help="Every sample's {sample}_trait_profile.tsv.")
    parser.add_argument("--out-summary", required=True,
                        help="Output: the grid with the words called / partial / "
                             "not called.")
    parser.add_argument("--out-mqc", required=True,
                        help="Output: the same grid as 2 / 1 / 0 with a MultiQC "
                             "header, drawn as a heatmap.")
    args = parser.parse_args(argv)

    trait_ids, grid = build_grid(args.profiles)

    # The grid in words: a header of trait_ids, then one line per sample.
    summary_lines = ["\t".join(["sample"] + trait_ids)]
    for sample, decisions in grid:
        summary_lines.append("\t".join([sample] + decisions))
    write_lines(args.out_summary, summary_lines)

    # The same grid for MultiQC: its '#' lines first, then the header, then each
    # decision as its number (DECISION_SCORES).
    mqc_lines = list(MULTIQC_HEADER)
    mqc_lines.append("\t".join(["sample"] + trait_ids))
    for sample, decisions in grid:
        scores = [str(DECISION_SCORES[decision]) for decision in decisions]
        mqc_lines.append("\t".join([sample] + scores))
    write_lines(args.out_mqc, mqc_lines)

    # One line for the log, so it shows what was combined.
    print(f"Trait summary: {len(grid)} sample(s) x {len(trait_ids)} trait(s) "
          f"written to {args.out_summary} and {args.out_mqc}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
