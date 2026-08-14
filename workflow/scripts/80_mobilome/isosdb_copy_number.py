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

# EVERY `action` / `reason` PAIR THIS SCRIPT CAN WRITE. `action` says what
# happened, `reason` says why. Note that NOTHING here changes an AMR gene's
# mobility tier - this whole leg is a QC measurement of how badly the assembler
# collapsed the IS copies, reported alongside the located count.
#
#   action=discarded      an ISOSDB database entry was not counted as detected
#     database_entry_not_fully_covered  less of the entry was covered by reads
#                                       than --min-covered-percent requires, so
#                                       it is not solid evidence the element is
#                                       present
#     depth_below_single_copy           read depth over the entry came out below
#                                       one copy's worth, so it cannot support a
#                                       copy-number estimate
#
#   action=not_applicable no estimate could be made, for a reason that is NOT a
#                         failure - the located count simply stands on its own
#     no_genome_baseline_depth          assembly depth could not be established,
#                                       so there is no single-copy baseline to
#                                       divide by. Nothing can be estimated.
#     isosdb_does_not_cover_this_organism  ISOSDB has no entries matching this
#                                       genome's IS at all
#     family_absent_from_isosdb         this particular IS family is not in the
#                                       database, so its collapse cannot be
#                                       measured even though others can
#     isosdb_detected_fewer_than_located  the read-mapping leg found FEWER copies
#                                       than ISEScan located on the contigs. The
#                                       delta is reported as NA rather than as a
#                                       negative number, because a negative
#                                       "collapse" is not meaningful - it means
#                                       the database is the limiting factor here,
#                                       not the assembly.
#
#   action=summary        one closing row per family, recording that an estimate
#                         was made
#     copy_number_estimate_complete

SUMMARY_COLUMNS = [
    "sample", "is_family",
    "located_copies",          # what ISEScan found on the contigs
    "estimated_copies",        # what the read depth implies (NA if undetectable)
    "collapse_delta",          # estimated - located, or NA - see db_informative
    "db_informative",          # did ISOSDB contain this family at all?
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

        # THE CHECK THIS LEG CANNOT WORK WITHOUT.
        #
        # ISEScan finds IS elements by profile HMM, which recognises a FAMILY.
        # This leg maps reads to ISOSDB, which requires nucleotide identity to a
        # specific catalogued element. Those are not the same sensitivity, and
        # the spec (§2.2) quantifies the gap from the ISOSDB paper itself: 97.5%
        # of its transposases have protein homologs in ISfinder but only 37.9%
        # have nucleotide ones.
        #
        # So for an organism ISOSDB does not cover - which for a non-clinical,
        # environmental isolate is the NORMAL case - no entry attracts reads, the
        # estimate comes out at zero, and the delta goes NEGATIVE. Reported
        # naively that reads as "the assembly did not collapse anything, we
        # over-called", when the truth is "this database has nothing to say about
        # this organism". Measured on hybrid sample 006 (Aquipseudomonas): 17 IS
        # located by ISEScan, 1 of 22,713 ISOSDB entries covered end to end,
        # delta -7.3 on the IS3 family alone.
        #
        # Absence of a nucleotide match is not evidence of absence of copies, so
        # when nothing was detected the estimate is reported as NA rather than as
        # a number that looks like a measurement.
        # A NEGATIVE delta is never evidence of anything except incomplete
        # database coverage, so it is never reported as a delta.
        #
        # This leg can only ever say "there are AT LEAST this many more copies
        # than you located". It cannot say "there are fewer", because ISEScan's
        # profile HMMs are strictly more sensitive than nucleotide mapping to a
        # specific catalogued element - so estimated < located always means
        # ISOSDB missed some, never that ISEScan over-called. On sample 006 the
        # IS3 family had 9 elements located and exactly ONE ISOSDB entry
        # detected: db_informative is technically true, but a delta of -7.3
        # measures the database, not the assembly.
        informative = bucket["entries"] > 0 and estimated >= located
        summary_rows.append({
            "sample": sample,
            "is_family": family,
            "located_copies": str(located),
            # The estimate itself is still worth showing when anything was
            # detected - it is the delta that must not be over-read.
            "estimated_copies": f"{estimated:.1f}" if bucket["entries"] else "NA",
            "collapse_delta": f"{estimated - located:.1f}" if informative else "NA",
            "db_informative": "TRUE" if informative else "FALSE",
            "n_db_entries_detected": str(bucket["entries"]),
            "max_entry_copy_number": f"{bucket['max']:.2f}",
            "genome_baseline_depth": f"{baseline:.1f}",
        })

        if bucket["entries"] and not informative:
            audit_rows.append(audit_row(
                sample, "not_applicable", "isosdb_detected_fewer_than_located",
                f"{family}: ISEScan located {located} element(s) but ISOSDB "
                f"read depth accounts for only {estimated:.1f}, from "
                f"{bucket['entries']} database entry/entries. No collapse "
                "estimate is reported. This leg can only ever say 'at least this "
                "many MORE copies than you located' - ISEScan's profile HMMs are "
                "strictly more sensitive than nucleotide mapping to a specific "
                "catalogued element, so a shortfall always means the database "
                "missed some, never that ISEScan over-called."))

        if not bucket["entries"]:
            audit_rows.append(audit_row(
                sample, "not_applicable", "family_absent_from_isosdb",
                f"{family}: ISEScan located {located} element(s) of this family, "
                "but no ISOSDB entry for it was covered by reads, so no copy "
                "number can be estimated. This says the database does not hold "
                "this organism's version of the family at nucleotide identity - "
                "NOT that the assembly collapsed nothing. Reported as NA rather "
                "than as a negative delta, which would read like a measurement."))

    total_located = sum(located_counts.values())
    total_estimated = sum(b["copies"] for b in per_family.values())
    families_informative = sum(1 for b in per_family.values() if b["entries"] > 0)
    families_total = len(set(per_family) | set(located_counts))

    audit_rows.append(audit_row(
        sample, "summary", "copy_number_estimate_complete",
        f"Genome baseline depth {baseline:.1f}x. ISEScan located "
        f"{total_located} IS element(s) on the contigs. ISOSDB could speak to "
        f"{families_informative} of {families_total} IS family/families; where it "
        f"could, read depth implies about {total_estimated:.1f} copies. A POSITIVE "
        "difference is the assembly collapse this module warns about - identical "
        "IS copies merged into one contig node - and is why the located count is "
        "reported as a floor; the estimate cannot say WHERE the extra copies are. "
        + ("A family ISOSDB does not cover is reported NA, never as a negative "
           "delta: no nucleotide match means the database is silent about this "
           "organism, not that nothing collapsed."
           if families_informative < families_total else "")))

    if families_informative == 0:
        audit_rows.append(audit_row(
            sample, "not_applicable", "isosdb_does_not_cover_this_organism",
            f"No ISOSDB entry was covered by reads for any of the {families_total} "
            "IS family/families ISEScan found, so this leg produced no estimate at "
            "all. Expected for isolates outside ISOSDB's sampling - it is built "
            "largely from human-associated metagenomes - and it means the located "
            "IS count remains an UNQUANTIFIED floor for this sample."))

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
    n_informative = sum(1 for r in summary if r["db_informative"] == "TRUE")
    if n_informative:
        print(
            f"Sample {args.sample}: {len(summary)} IS family/families; "
            f"{total_located} located on contigs. ISOSDB could quantify collapse "
            f"for {n_informative} of them (see {args.out_table})."
        )
    else:
        print(
            f"Sample {args.sample}: {len(summary)} IS family/families, "
            f"{total_located} element(s) located on contigs, but ISOSDB could not "
            f"quantify collapse for ANY of them - the located count stays an "
            f"unquantified floor. Reasons in {args.out_audit}."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
