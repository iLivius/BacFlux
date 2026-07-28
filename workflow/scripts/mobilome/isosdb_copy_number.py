#!/usr/bin/env python3
"""Estimate how many IS copies the assembly LOST, by counting reads instead of contigs.

THE PROBLEM THIS EXISTS TO MEASURE
    Insertion sequences are the single biggest cause of contig breaks in a
    short-read assembly, because multiple identical copies of one IS collapse
    into a single node in the assembly graph. So the thing the mobilome module is
    hunting is often precisely what destroyed the contig it is looking at.

    Everywhere else this module says "the located IS count is a FLOOR, not a
    count". That is honest, but it is also unquantified: a reader has no idea
    whether the floor is 2 short or 40 short. This script puts a number on it.

HOW, IN ONE SENTENCE
    Reads are immune to assembly collapse - every copy of an IS contributes its
    own reads whether or not the assembler kept them apart - so an IS present in
    five copies attracts about five times the read depth of the single-copy
    chromosome, and dividing the two gives the copy number.

        copy number  ~=  depth over the IS  /  depth over the genome as a whole

WHAT IT TAKES IN
    - BBMap covstats from mapping this sample's reads against ISOSDB, the openly
      licensed IS nucleotide database (rule isosdb_map).
    - BBMap covstats from mapping the SAME reads against the sample's own
      assembly (rule assembly_depth), which supplies the single-copy baseline.
    - ISOSDB's IS_fam_annot.txt, so results can be reported per IS FAMILY as well
      as per database entry.
    - Optionally the ISEScan table, so the located count can be compared with the
      read-based estimate - which is the whole point.

WHAT IT PRODUCES
    - a per-family summary: located copies, read-based estimate, and the delta;
    - an audit TSV, because every element excluded from the estimate needs a
      stated reason.

WHAT THIS IS NOT
    It does not change any AMR gene's mobility tier, and it must not: it says
    nothing about WHERE the extra copies are, only that they exist. It is a
    quality metric attached to the IS inventory, and it is the honest companion
    to the "floor, not a count" warning.

WHY DEPTH RATIOS AND NOT READ COUNTS
    Longer references collect more reads simply by being longer. Depth (reads x
    read length / reference length) already divides that out, which is why BBMap's
    Avg_fold is the column used rather than Plus_reads + Minus_reads.
"""

import argparse
import os
import statistics
import sys


# Below this fraction of a database entry covered by reads, the "depth" is being
# computed over a reference that is mostly untouched, and the number means
# nothing. A genuine IS present in the sample is covered end to end; a partial
# hit is usually a conserved domain shared with a different family.
DEFAULT_MIN_COVERED_PERCENT = 90.0

# Depth below this multiple of the genome baseline is treated as "not present".
# Set below 1.0 on purpose: a single-copy IS should sit at ~1x, and sampling
# noise plus mapping loss regularly pushes a real single copy to 0.6-0.8x.
DEFAULT_MIN_COPY_NUMBER = 0.5


def read_covstats(path):
    """Read a BBMap covstats file into a list of dicts.

    BBMap writes a header line beginning '#ID' and then one row per REFERENCE
    sequence. The columns used here are Avg_fold (depth) and Covered_percent
    (how much of the reference the reads actually touched).

    A missing file returns an empty list rather than raising: the read-mapping
    leg only runs in modes that have short reads, and its absence is a normal
    configuration rather than a failure.
    """
    rows = []
    if not path or not os.path.isfile(path):
        return rows

    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if len(lines) < 2:
        return rows

    header = lines[0].lstrip("#").split("\t")
    for line in lines[1:]:
        values = line.split("\t")
        rows.append(dict(zip(header, values)))
    return rows


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def genome_baseline_depth(assembly_covstats):
    """The read depth of a SINGLE-COPY region, which is the denominator.

    Takes in: covstats from mapping the reads against the sample's own assembly,
              one row per contig.
    Returns:  a float depth, or None when it cannot be established.

    Uses the LENGTH-WEIGHTED MEDIAN across contigs rather than the mean, and the
    reason matters. A bacterial assembly is one long chromosome plus a handful of
    short contigs, and those short contigs are exactly the ones with unstable
    depth - repeat collapses sit at several times the genome depth, and low-
    coverage junk sits near zero. A plain mean over contigs lets a 2 kb outlier
    weigh as much as a 5 Mb chromosome; a plain median over contigs does the
    same. Weighting by length makes the answer "the depth of a typical BASE",
    which is what a single-copy baseline should mean.
    """
    weighted = []
    for row in assembly_covstats:
        depth = to_float(row.get("Avg_fold"))
        length = to_int(row.get("Length"))
        if depth <= 0 or length <= 0:
            continue
        weighted.append((depth, length))

    if not weighted:
        return None

    weighted.sort(key=lambda item: item[0])
    total_length = sum(length for _depth, length in weighted)
    half = total_length / 2.0

    running = 0
    for depth, length in weighted:
        running += length
        if running >= half:
            return depth
    return weighted[-1][0]


def read_family_map(path):
    """ISOSDB entry ID -> IS family, from ISOSDB's IS_fam_annot.txt.

    The file is a two-column TSV, 'IS' and 'IS_fam'. Not every entry has a family
    (about 570 of 22,713 are unannotated), so a missing key is normal and those
    entries are reported under 'unassigned' rather than dropped.
    """
    families = {}
    if not path or not os.path.isfile(path):
        return families
    with open(path, encoding="utf-8", errors="replace") as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if index == 0 and fields and fields[0].strip().upper() == "IS":
                continue                      # header
            if len(fields) < 2:
                continue
            families[fields[0].strip()] = fields[1].strip()
    return families


def read_located_families(path):
    """IS family -> how many copies ISEScan actually LOCATED on the contigs.

    This is the number the read-based estimate is compared against, and the
    comparison is the deliverable: located is what survived assembly, estimated
    is what the reads say was there.
    """
    counts = {}
    if not path or not os.path.isfile(path):
        return counts

    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if len(lines) < 2:
        return counts

    header = lines[0].split("\t")
    try:
        family_index = header.index("is_family")
    except ValueError:
        try:
            family_index = header.index("family")
        except ValueError:
            return counts

    for line in lines[1:]:
        fields = line.split("\t")
        if family_index >= len(fields):
            continue
        family = fields[family_index].strip() or "unassigned"
        counts[family] = counts.get(family, 0) + 1
    return counts


def audit_row(sample, action, reason, detail, element="NA"):
    return {
        "sample": sample, "element": element,
        "action": action, "reason": reason, "detail": detail,
    }


AUDIT_COLUMNS = ["sample", "element", "action", "reason", "detail"]

SUMMARY_COLUMNS = [
    "sample", "is_family",
    "located_copies",          # what ISEScan found on the contigs
    "estimated_copies",        # what the read depth implies
    "collapse_delta",          # estimated - located; how much the assembly lost
    "n_db_entries_detected",   # how many ISOSDB entries contributed
    "max_entry_copy_number",   # the deepest single entry, for context
    "genome_baseline_depth",
]


def estimate_copies(sample, isosdb_covstats, assembly_covstats, family_map,
                    located_counts,
                    min_covered_percent=DEFAULT_MIN_COVERED_PERCENT,
                    min_copy_number=DEFAULT_MIN_COPY_NUMBER):
    """Turn read depths into a per-family copy-number estimate.

    Returns (summary_rows, audit_rows).

    A NOTE ON WHAT THE NUMBER IS AND IS NOT. ISOSDB is dereplicated at 95%
    identity, but IS families remain similar enough that one real element can
    attract reads across several database entries. BBMap is run with
    ambiguous=best so each read lands on one entry only, which stops the total
    from being multiplied - but it also means the split between near-identical
    entries is arbitrary. That is why the family total is the headline and the
    per-entry numbers are reported only as context: the family sum is robust to
    where an ambiguous read landed, the per-entry split is not.
    """
    audit_rows = []

    baseline = genome_baseline_depth(assembly_covstats)
    if baseline is None or baseline <= 0:
        audit_rows.append(audit_row(
            sample, "not_applicable", "no_genome_baseline_depth",
            "Read depth over the assembly could not be established, so there is "
            "no single-copy baseline to divide by and no copy number can be "
            "estimated. The located IS count stands on its own, still a floor."))
        return [], audit_rows

    per_family = {}
    for row in isosdb_covstats:
        entry = (row.get("ID") or "").split()[0] if row.get("ID") else ""
        if not entry:
            continue

        covered = to_float(row.get("Covered_percent"))
        depth = to_float(row.get("Avg_fold"))
        family = family_map.get(entry, "unassigned")

        if covered < min_covered_percent:
            # Depth over a reference the reads barely touched is not a depth for
            # that element - it is usually a conserved domain shared with another
            # family. Quietly averaging it in would inflate every estimate.
            audit_rows.append(audit_row(
                sample, "discarded", "database_entry_not_fully_covered",
                f"{entry} ({family}): reads cover only {covered:.1f}% of this "
                f"ISOSDB entry (needs {min_covered_percent:.0f}%), so its depth "
                "does not describe a copy of this element.",
                element=entry))
            continue

        copy_number = depth / baseline
        if copy_number < min_copy_number:
            audit_rows.append(audit_row(
                sample, "discarded", "depth_below_single_copy",
                f"{entry} ({family}): depth {depth:.1f}x is only "
                f"{copy_number:.2f}x the genome baseline ({baseline:.1f}x), below "
                f"the {min_copy_number} threshold for calling the element present.",
                element=entry))
            continue

        bucket = per_family.setdefault(family, {"copies": 0.0, "entries": 0, "max": 0.0})
        bucket["copies"] += copy_number
        bucket["entries"] += 1
        bucket["max"] = max(bucket["max"], copy_number)

    summary_rows = []
    for family in sorted(set(per_family) | set(located_counts)):
        bucket = per_family.get(family, {"copies": 0.0, "entries": 0, "max": 0.0})
        located = located_counts.get(family, 0)
        estimated = bucket["copies"]
        summary_rows.append({
            "sample": sample,
            "is_family": family,
            "located_copies": str(located),
            "estimated_copies": f"{estimated:.1f}",
            "collapse_delta": f"{estimated - located:.1f}",
            "n_db_entries_detected": str(bucket["entries"]),
            "max_entry_copy_number": f"{bucket['max']:.2f}",
            "genome_baseline_depth": f"{baseline:.1f}",
        })

    total_located = sum(located_counts.values())
    total_estimated = sum(b["copies"] for b in per_family.values())
    audit_rows.append(audit_row(
        sample, "summary", "copy_number_estimate_complete",
        f"Genome baseline depth {baseline:.1f}x. ISEScan located "
        f"{total_located} IS element(s) on the contigs; read depth against ISOSDB "
        f"implies about {total_estimated:.1f}. A positive difference is the "
        "assembly collapse this module warns about - identical IS copies merged "
        "into one contig node - and is why the located count is reported as a "
        "floor. The estimate cannot say WHERE the extra copies are."))

    return summary_rows, audit_rows


def write_tsv(path, columns, rows):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "NA")) for column in columns) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Estimate IS copy number from read depth against ISOSDB, and compare "
            "it with what ISEScan located on the contigs. Quantifies how much of "
            "the IS inventory the assembly collapsed."
        )
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--isosdb-covstats", required=True,
                        help="BBMap covstats: reads vs ISOSDB.")
    parser.add_argument("--assembly-covstats", required=True,
                        help="BBMap covstats: the SAME reads vs this sample's "
                             "assembly. Supplies the single-copy baseline.")
    parser.add_argument("--family-map", default="",
                        help="ISOSDB IS_fam_annot.txt (entry -> IS family).")
    parser.add_argument("--is-table", default="",
                        help="ISEScan table, for the located-vs-estimated comparison.")
    parser.add_argument("--min-covered-percent", type=float,
                        default=DEFAULT_MIN_COVERED_PERCENT)
    parser.add_argument("--min-copy-number", type=float,
                        default=DEFAULT_MIN_COPY_NUMBER)
    parser.add_argument("--out-table", required=True)
    parser.add_argument("--out-audit", required=True)
    args = parser.parse_args(argv)

    summary, audit = estimate_copies(
        args.sample,
        read_covstats(args.isosdb_covstats),
        read_covstats(args.assembly_covstats),
        read_family_map(args.family_map),
        read_located_families(args.is_table),
        min_covered_percent=args.min_covered_percent,
        min_copy_number=args.min_copy_number,
    )

    write_tsv(args.out_table, SUMMARY_COLUMNS, summary)
    write_tsv(args.out_audit, AUDIT_COLUMNS, audit)

    total_located = sum(to_int(r["located_copies"]) for r in summary)
    total_estimated = sum(to_float(r["estimated_copies"]) for r in summary)
    print(
        f"Sample {args.sample}: {len(summary)} IS family/families; "
        f"{total_located} located on contigs, ~{total_estimated:.0f} implied by "
        f"read depth (difference = assembly collapse; see {args.out_audit})."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
