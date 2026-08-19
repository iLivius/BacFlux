#!/usr/bin/env python3
"""
Score the EBI mobilome-annotation-pipeline (MAP) and BacFlux's own ICE/IME
caller against ONE curated pilot set, using the same "recovered" rule the
phase-7 benchmark harness uses: an element counts as recovered if any call
overlaps its curated interval at all.

Why a second, narrower scorer exists alongside compare_ice_calls.py:
compare_ice_calls.py attributes every call to whichever curated element it
overlaps most, across the whole of ground_truth.tsv. That is the right view for
inspecting boundaries, but it mixes ICEs into an IME comparison, because these
chromosomes carry both. This script answers the narrower question the pilot was
designed for - "of the N elements we deliberately curated, how many did each
caller find?" - and nothing else.

Usage:
    python3 score_pilot_elements.py <pilot_set.tsv> <map_run_name>
    e.g.  python3 score_pilot_elements.py \
              /media/data/.../phase7_benchmark/ime_pilot_set.tsv  ime12

Inputs
  pilot_set.tsv : the curated pilot table (element, accession, start, end, ...)
  MAP calls     : runs/<map_run_name>/results/<sample>/prediction/icefinder2lite/<sample>_ices.tsv
  BacFlux calls : phase7_benchmark/work/<sample>/ice_elements.tsv
  where <sample> is the accession with '.' replaced by '_', which is how both
  the benchmark harness and our MAP samplesheets name things.

Output: a per-element table on stdout plus the recovered counts for each caller.
"""

import csv
import os
import sys

BENCH = "/media/data/antonielli_dir/BacFlux_v2_validation/phase7_benchmark"
EBI = "/media/data/antonielli_dir/BacFlux_v2_validation/ebi_map"


def map_calls(sample, run_name):
    """MAP's predicted elements for one genome, as (start, end, type) tuples.

    MAP stores coordinates as "start..end" in the ice_location column. Its
    contigs are renamed to contig_N, but renaming does not shift coordinates,
    and every pilot genome here is a single replicon, so the numbers are
    directly comparable to the curated ones.
    """
    path = f"{EBI}/runs/{run_name}/results/{sample}/prediction/icefinder2lite/{sample}_ices.tsv"
    calls = []
    if not os.path.exists(path):
        return calls
    for row in csv.DictReader(open(path), delimiter="\t"):
        location = row.get("ice_location", "")
        if ".." not in location:
            continue
        start_text, end_text = location.split("..", 1)
        calls.append((int(start_text), int(end_text), row.get("ice_type", "")))
    return calls


def bacflux_calls(sample):
    """BacFlux's predicted elements for the same genome, from the phase-7 harness."""
    path = f"{BENCH}/work/{sample}/ice_elements.tsv"
    calls = []
    if not os.path.exists(path):
        return calls
    for row in csv.DictReader(open(path), delimiter="\t"):
        calls.append((int(row["start"]), int(row["end"]), row.get("element_type", "")))
    return calls


def best_overlap(calls, true_start, true_end):
    """The single call overlapping the curated element most, and by how much."""
    best_bp = 0
    best_call = None
    for call_start, call_end, call_type in calls:
        shared = max(0, min(call_end, true_end) - max(call_start, true_start) + 1)
        if shared > best_bp:
            best_bp, best_call = shared, (call_start, call_end, call_type)
    return best_bp, best_call


def main():
    pilot_path = sys.argv[1]
    run_name = sys.argv[2]

    rows = list(csv.DictReader(open(pilot_path), delimiter="\t"))
    print(f"{'element':<26}{'accession':<14}{'len':>8}  {'EBI MAP':<24}{'BacFlux':<24}")

    map_found = 0
    bacflux_found = 0
    for row in rows:
        accession = row["accession"]
        sample = accession.replace(".", "_")
        true_start, true_end = int(row["start"]), int(row["end"])
        true_len = true_end - true_start + 1

        map_bp, map_call = best_overlap(map_calls(sample, run_name), true_start, true_end)
        bac_bp, bac_call = best_overlap(bacflux_calls(sample), true_start, true_end)

        map_text = f"{100 * map_bp / true_len:.1f}% {map_call[2][:10]}" if map_call else "-- no call"
        bac_text = f"{100 * bac_bp / true_len:.1f}% {bac_call[2][:10]}" if bac_call else "-- no call"
        map_found += map_bp > 0
        bacflux_found += bac_bp > 0

        print(f"{row['element']:<26}{accession:<14}{true_len:>8}  {map_text:<24}{bac_text:<24}")

    total = len(rows)
    print(f"\nRECOVERED (any overlap):  EBI MAP {map_found}/{total}    BacFlux {bacflux_found}/{total}")


if __name__ == "__main__":
    main()
