#!/usr/bin/env python3
"""
Put the EBI mobilome-annotation-pipeline (MAP) ICE/IME calls next to BacFlux's
own calls and next to ICEberg's curated coordinates, for the same genomes.

Why this exists
---------------
BacFlux's ICE caller (workflow/scripts/mobilome/conjscan_to_ice.py) and MAP's
ICEfinder2-lite are two independent implementations that happen to run the SAME
CONJScan/ICEscan v2.0.1 HMM profile set. So when they disagree about where an
element starts and ends, the disagreement is in the logic wrapped around the
models, not in the evidence the models saw. That makes MAP a useful second
opinion on our boundary errors: if MAP misses the same edge we miss, the
curated coordinate is the thing to question; if MAP nails it and we do not, the
error is ours.

Inputs (all read-only, nothing here writes into the repo or the benchmark)
  1. MAP     : runs/<run>/results/<sample>/prediction/icefinder2lite/<sample>_ices.tsv
               one row per predicted ICE; coordinates live in `ice_location` as
               "start..end" on MAP's renamed contigs (contig_1, contig_2, ...).
  2. MAP     : runs/<run>/results/<sample>/preprocessing/<sample>_contigID.map
               two columns, ">contig_N <original name>", used to translate MAP's
               renamed contigs back to the accessions everything else uses.
               Renaming changes names only, never sequence, so coordinates are
               already in the same frame - no offset correction is needed.
  3. BacFlux : phase7_benchmark/work/<sample>/ice_elements.tsv
  4. Truth   : phase7_benchmark/ground_truth.tsv, ICEberg's curated coordinates.

Output
  A TSV on stdout, one row per (genome, caller, element), plus overlap and
  boundary-offset columns measured against the curated element. Consumed by a
  human reading the head-to-head, not by another script.
"""

import csv
import os
import sys

BENCH = "/media/data/antonielli_dir/BacFlux_v2_validation/phase7_benchmark"
EBI = "/media/data/antonielli_dir/BacFlux_v2_validation/ebi_map"


def read_contig_name_map(path):
    """Translate MAP's renamed contigs back to the original accessions.

    MAP's RENAME step rewrites every contig header to contig_1, contig_2, ...
    and records the mapping in this file. We need the reverse direction so that
    a MAP call on "contig_1" can be reported against, say, "CP042858.1".

    Returns {map_name: original_name}; empty dict if the file is absent.
    """
    names = {}
    if not os.path.exists(path):
        return names
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                map_name = fields[0].lstrip(">")
                names[map_name] = fields[1]
    return names


def read_map_ices(sample, run_name):
    """Read MAP's ICE predictions for one genome.

    MAP writes `ice_location` as "start..end" (1-based inclusive, GenBank
    style). We split that into integers so the offsets against the curated
    coordinates are a plain subtraction.
    """
    base = os.path.join(EBI, "runs", run_name, "results", sample)
    ices_path = os.path.join(base, "prediction", "icefinder2lite", f"{sample}_ices.tsv")
    name_map = read_contig_name_map(
        os.path.join(base, "preprocessing", f"{sample}_contigID.map")
    )

    calls = []
    if not os.path.exists(ices_path):
        return calls

    with open(ices_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            location = row.get("ice_location", "")
            if ".." not in location:
                continue
            start_text, end_text = location.split("..", 1)
            contig = row.get("contig", "")
            calls.append(
                {
                    "caller": "EBI_MAP",
                    "element_type": row.get("ice_type", ""),
                    "contig": name_map.get(contig, contig),
                    "start": int(start_text),
                    "end": int(end_text),
                    # MAP reports "-" when vmatch found no flanking direct repeat.
                    "boundary_evidence": row.get("direct_repeats", "-"),
                    "extra": f"relaxase={row.get('relaxase_type', '-')};"
                             f"mpf={row.get('mating_pair_formation_systems', '-')};"
                             f"near_rna={row.get('close_to_RNA', '-')}",
                }
            )
    return calls


def read_bacflux_ices(sample):
    """Read BacFlux's own ICE/IME calls for the same genome.

    Produced by the phase-7 benchmark harness (process_one.sh -> conjscan_to_ice.py).
    """
    path = os.path.join(BENCH, "work", sample, "ice_elements.tsv")
    calls = []
    if not os.path.exists(path):
        return calls

    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            calls.append(
                {
                    "caller": "BacFlux",
                    "element_type": row.get("element_type", ""),
                    "contig": row.get("contig", ""),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    # Ours records HOW the edge was set: tRNA-anchored, de novo
                    # repeat search, or nothing found.
                    "boundary_evidence": f"{row.get('boundary_method', '-')}"
                                         f"/att{row.get('att_length_bp', '0')}bp",
                    "extra": f"conf={row.get('confidence', '-')};"
                             f"relaxase={row.get('relaxase_type', '-')};"
                             f"mpf={row.get('mpf_type', '-')}",
                }
            )
    return calls


def read_ground_truth():
    """ICEberg's curated element coordinates, keyed by accession.

    One accession can carry more than one curated element, so the value is a list.
    """
    truth = {}
    with open(os.path.join(BENCH, "ground_truth.tsv")) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            truth.setdefault(row["accession"], []).append(
                {
                    "kind": row["kind"],
                    "element": row["element"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "organism": row["organism"],
                }
            )
    return truth


def overlap_bp(a_start, a_end, b_start, b_end):
    """Base pairs shared by two intervals; 0 if they do not touch."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def main():
    run_name = sys.argv[1] if len(sys.argv) > 1 else "pilot3"
    # sample -> the accession ICEberg curates it under. The benchmark harness
    # turns '.' into '_' for directory names, so we carry both forms.
    samples = sys.argv[2:] or ["CP042858_1", "NC_000964_3", "U15027"]

    truth = read_ground_truth()
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "sample", "caller", "element_type", "contig", "start", "end", "length_bp",
            "boundary_evidence", "extra",
            "matched_curated_element", "curated_start", "curated_end", "curated_len",
            "overlap_bp", "pct_of_curated", "start_offset", "end_offset",
        ]
    )

    for sample in samples:
        # The benchmark harness names its work directories after the accession
        # with '.' replaced by '_', so "CP002888.1" becomes "CP002888_1" and
        # "NC_013316.1" becomes "NC_013316_1". Going back is ambiguous, because
        # RefSeq accessions legitimately contain an underscore ("NC_013316").
        # Only the LAST underscore is ever the version separator, so restore
        # that one and try both spellings against the curated table.
        candidates = [sample]
        if "_" in sample and sample.rsplit("_", 1)[1].isdigit():
            head, version = sample.rsplit("_", 1)
            candidates.append(f"{head}.{version}")
        curated = []
        for candidate in candidates:
            if candidate in truth:
                curated = truth[candidate]
                break

        calls = read_map_ices(sample, run_name) + read_bacflux_ices(sample)

        if not calls:
            writer.writerow([sample, "NO_CALLS", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""])
            continue

        for call in calls:
            # Attribute each call to whichever curated element it overlaps most.
            best = None
            best_overlap = 0
            for entry in curated:
                shared = overlap_bp(call["start"], call["end"], entry["start"], entry["end"])
                if shared > best_overlap:
                    best, best_overlap = entry, shared

            if best is None:
                writer.writerow(
                    [
                        sample, call["caller"], call["element_type"], call["contig"],
                        call["start"], call["end"], call["end"] - call["start"] + 1,
                        call["boundary_evidence"], call["extra"],
                        "NONE", "", "", "", 0, "", "", "",
                    ]
                )
                continue

            curated_len = best["end"] - best["start"] + 1
            writer.writerow(
                [
                    sample, call["caller"], call["element_type"], call["contig"],
                    call["start"], call["end"], call["end"] - call["start"] + 1,
                    call["boundary_evidence"], call["extra"],
                    best["element"], best["start"], best["end"], curated_len,
                    best_overlap, f"{100.0 * best_overlap / curated_len:.1f}",
                    # Positive = our/MAP edge sits inside the curated element.
                    call["start"] - best["start"], call["end"] - best["end"],
                ]
            )


if __name__ == "__main__":
    main()
