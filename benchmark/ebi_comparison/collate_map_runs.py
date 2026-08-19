#!/usr/bin/env python3
"""
Collate the EBI mobilome-annotation-pipeline (MAP) benchmark runs into one
per-genome status table.

Why this exists
---------------
MAP was run over the 40 BacFlux benchmark genomes in four separate Nextflow
runs (pilot3, ime12, negative, ice15), because they were launched at different
times as the host became free. Each run writes its own Nextflow execution
trace. This script reads those traces plus the published results tree and
answers, per genome: did MAP finish, how long did it take, and where did its
ICE/IME calls land.

The one non-obvious point about MAP's output
--------------------------------------------
MAP's ICEfinder2-lite subworkflow is a chain:

    HMMSCAN -> PRESCAN_TO_FASTA -> {VMATCH, MACSYFINDER, BLASTP_PROKKA}
            -> REFINE_BOUNDARIES -> <sample>_ices.tsv

Only REFINE_BOUNDARIES writes the *_ices.tsv table. If PRESCAN_TO_FASTA finds
no candidate region, subworkflows/local/icefinder2lite.nf filters that sample
out of the channel and the later steps never run for it, so NO FILE IS WRITTEN
AT ALL. A missing *_ices.tsv therefore means "MAP looked and called nothing",
which is a real result, not a crash. This script distinguishes the two by
checking that PRESCAN_TO_FASTA completed for the sample.

Inputs : runs/<run>/results/pipeline_info/execution_trace_*.txt  (Nextflow trace)
         runs/<run>/results/<sample>/prediction/icefinder2lite/<sample>_ices.tsv
Output : a TSV on stdout, one row per genome, consumed by the benchmark writeup.
"""

import csv
import glob
import os
import re
from datetime import datetime

EBI_MAP = "/media/data/antonielli_dir/BacFlux_v2_validation/ebi_map"
BENCH = "/media/data/antonielli_dir/BacFlux_v2_validation/phase7_benchmark"

# The runs to scan, in the order they were launched. A genome can appear in
# more than one run (NC_000964.3 and U15027 were in the 3-genome smoke test as
# well as their real set); later runs win, since they are the completed ones.
RUNS = ["pilot3", "ime12", "negative", "ice15"]

# MAP samplesheet names replace '.' with '_' (CP042858.1 -> CP042858_1),
# because Nextflow uses the sample id in file names.
def to_map_id(accession):
    return accession.replace(".", "_")


def load_sets():
    """Read the three benchmark set definitions -> {accession: [kinds]}.

    Two accessions (NC_004668.1, CP048437.1) carry both a curated ICE and a
    curated IME, which is why 18+12+12 rows collapse to 40 unique genomes.
    """
    membership = {}
    for fname, kind in [
        ("pilot_set.tsv", "ICE"),
        ("ime_pilot_set.tsv", "IME"),
        ("negative_set.tsv", "NEG"),
    ]:
        with open(os.path.join(BENCH, fname)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                acc = row["accession"]
                membership.setdefault(acc, [])
                if kind not in membership[acc]:
                    membership[acc].append(kind)
    return membership


def parse_duration(text):
    """Nextflow duration strings ('2s', '1m 4s', '1h 2m 3s', '350ms') -> seconds."""
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)\s*(ms|s|m|h|d)", text):
        value = float(value)
        total += value * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    return total


def scan_run(run):
    """Read one run's newest execution trace.

    Returns {sample: {...}} recording, per sample, which pipeline stages ran,
    whether anything failed, and the wall-clock span from the first task
    submitted for that sample to the last one finishing. That span overlaps
    between samples because Nextflow runs them concurrently - it is a per-genome
    elapsed time, NOT a share of the total, and the two do not sum.
    """
    traces = sorted(glob.glob(
        os.path.join(EBI_MAP, "runs", run, "results", "pipeline_info",
                     "execution_trace_*.txt")))
    if not traces:
        return {}
    # Several traces exist when a run was restarted; the largest is the real one.
    trace = max(traces, key=os.path.getsize)

    samples = {}
    with open(trace) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["name"]
            match = re.search(r"\(([^)]*)\)$", name)
            if not match:
                continue
            sample = match.group(1)
            # RENAME emits per-contig-size FASTAs, so some task tags are
            # '<sample>_100kb_contigs.1' rather than the bare sample id.
            sample = re.sub(r"_(1|5|100)kb_contigs\.\d+$", "", sample)
            if not sample or sample.isdigit():
                continue  # MULTIQC / DUMPSOFTWAREVERSIONS are tagged '(1)'

            rec = samples.setdefault(sample, {
                "stages": set(), "failed": [], "start": None, "end": None,
                "cached": False})
            process = name.split(" (")[0].split(":")[-1]
            if row["status"] in ("COMPLETED", "CACHED"):
                rec["stages"].add(process)
            else:
                rec["failed"].append(f"{process}:{row['status']}")

            if row["status"] == "CACHED":
                # A cached task keeps the submit time and duration of the run
                # that first produced it, which can be hours earlier. Including
                # those would make this genome look far slower than it was, so
                # they are excluded from the timing and just flagged.
                rec["cached"] = True
                continue

            try:
                submitted = datetime.strptime(row["submit"], "%Y-%m-%d %H:%M:%S.%f")
            except (ValueError, KeyError):
                continue
            finished = submitted.timestamp() + parse_duration(row["duration"])
            if rec["start"] is None or submitted.timestamp() < rec["start"]:
                rec["start"] = submitted.timestamp()
            if rec["end"] is None or finished > rec["end"]:
                rec["end"] = finished
    return samples


def main():
    membership = load_sets()

    # Walk the runs in launch order so a genome re-run later overwrites the
    # earlier, possibly partial, record.
    per_sample = {}
    for run in RUNS:
        for sample, rec in scan_run(run).items():
            prior = per_sample.get(sample)
            # Keep whichever record got further down the ICEfinder chain.
            if prior and "PRESCAN_TO_FASTA" in prior["rec"]["stages"] \
                     and "PRESCAN_TO_FASTA" not in rec["stages"]:
                continue
            per_sample[sample] = {"run": run, "rec": rec}

    rows = []
    for accession in sorted(membership):
        sid = to_map_id(accession)
        entry = per_sample.get(sid)
        kinds = "+".join(membership[accession])

        if entry is None:
            rows.append([accession, sid, kinds, "-", "NOT_RUN", "", "", "", ""])
            continue

        run = entry["run"]
        rec = entry["rec"]
        stages = rec["stages"]

        # Where MAP publishes the ICE/IME table for this sample.
        ices = os.path.join(EBI_MAP, "runs", run, "results", sid,
                            "prediction", "icefinder2lite", f"{sid}_ices.tsv")
        has_file = os.path.exists(ices)
        n_calls = ""
        if has_file:
            with open(ices) as fh:
                n_calls = str(max(0, sum(1 for _ in fh) - 1))

        # Decide completed vs failed, judged on the ICE branch specifically -
        # that is the only part of MAP under comparison here. Reaching
        # PRESCAN_TO_FASTA means the ICE branch ran to its decision point;
        # whether a table exists after that is biology, not a crash.
        ice_branch_ran = "PRESCAN_TO_FASTA" in stages
        broken = [f for f in rec["failed"] if "ICEFINDER" in f.upper()
                  or f.split(":")[0] in ("HMMSCAN", "PRESCAN_TO_FASTA",
                                         "REFINE_BOUNDARIES", "MACSYFINDER",
                                         "VMATCH", "RENAME", "PROKKA")]

        if broken:
            status = "FAILED"
        elif has_file or ice_branch_ran:
            status = "COMPLETED"
            if not has_file:
                n_calls = "0"
        else:
            status = "INCOMPLETE"

        if has_file:
            outcome = "ices.tsv written"
        elif status == "COMPLETED":
            outcome = "no ICE candidate at prescan; no table written"
        else:
            outcome = ";".join(rec["failed"]) or "did not reach ICEfinder"

        # Anything that failed outside the ICE branch is worth stating but does
        # not invalidate the ICE calls (e.g. a downstream GFF merge that was
        # killed when a run was stopped after ICEfinder had already finished).
        other = [f for f in rec["failed"] if f not in broken]
        if other:
            outcome += f" [non-ICE stage {','.join(other)}]"
        if rec["cached"]:
            outcome += " [some stages reused from an earlier run]"

        elapsed = ""
        if rec["start"] and rec["end"]:
            elapsed = f"{(rec['end'] - rec['start']) / 60:.1f}"

        rows.append([accession, sid, kinds, run, status, elapsed,
                     n_calls, outcome, ices if has_file else ""])

    writer = csv.writer(__import__("sys").stdout, delimiter="\t",
                        lineterminator="\n")
    writer.writerow(["accession", "map_sample_id", "set", "run", "status",
                     "elapsed_min", "n_ice_calls", "outcome", "ices_tsv_path"])
    writer.writerows(rows)


if __name__ == "__main__":
    main()
