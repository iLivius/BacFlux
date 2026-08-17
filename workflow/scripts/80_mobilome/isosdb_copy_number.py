#!/usr/bin/env python3
"""Estimate how many IS copies the assembly LOST, by counting reads instead of
contigs.

Why the located IS count needs a number attached
------------------------------------------------
Insertion sequences are the single biggest cause of contig breaks in a
short-read assembly, because several identical copies of one IS collapse into a
single node in the assembly graph. The thing the mobilome module is hunting is
often precisely what destroyed the contig it is looking at.

Everywhere else this module says "the located IS count is a FLOOR, not a count".
That is honest, but unquantified: a reader has no idea whether the floor is 2
short or 40 short. This puts a number on it.

Reads are immune to assembly collapse — every copy of an IS contributes its own
reads whether or not the assembler kept them apart — so an IS present in five
copies attracts about five times the read depth of the single-copy chromosome,
and dividing the two gives the copy number:

    copy number  ~=  depth over the IS  /  depth over the genome as a whole

Depth and not raw read counts, because a longer reference collects more reads
simply by being longer. Depth (reads x read length / reference length) already
divides that out, which is why BBMap's Avg_fold is the column read here rather
than Plus_reads + Minus_reads.

What it reads, and what reads it
--------------------------------
  --isosdb-covstats    BBMap coverage of this sample's reads against ISOSDB, the
                       openly licensed (MIT) IS nucleotide database
                       source: rule isosdb_map
  --assembly-covstats  BBMap coverage of the SAME reads against this sample's own
                       assembly, which supplies the single-copy baseline
                       source: rule assembly_depth
  --family-map         ISOSDB's IS_fam_annot.txt, so results are reported per IS
                       FAMILY as well as per database entry
                       source: rule isosdb_db
  --is-table           the IS elements ISEScan actually located on the contigs —
                       the number the read-based estimate is set against
                       source: rule isescan_table

Out come a per-family summary (located copies, read-based estimate, delta) and an
audit TSV giving a reason for every element left out of the estimate. Nothing
downstream consumes either — they are workflow targets read by a person.

Reads are needed, so the leg runs in illumina and hybrid mode only, and only when
mobilome.isosdb.fasta_url or mobilome.isosdb.dir is set — see rule is_copy_number
in shared/80_mobilome.smk.

What this number is not
-----------------------
It changes no AMR gene's mobility tier, and must not: it says nothing about WHERE
the extra copies are, only that they exist. It is a quality metric attached to
the IS inventory, and the honest companion to the "floor, not a count" warning.

The trap is a delta that comes out NEGATIVE. When ISOSDB holds nothing matching
this organism the estimate comes out low, and a negative delta reads like "the
assembly collapsed nothing" when the truth is "this database is silent about this
genome". So it is never reported as a delta — the reasoning is in estimate_copies,
under db_informative.
"""

import argparse
import os
import statistics
import sys


# ── The two thresholds, and why they sit where they do ───────────────────────
# Both are defaults only. rule is_copy_number always passes the config values
# (mobilome.isosdb.min_covered_percent and mobilome.isosdb.min_copy_number), so
# these numbers apply when the script is run by hand.

# Below this fraction of a database entry covered by reads, the "depth" is being
# computed over a reference that is mostly untouched, and the number means
# nothing. A genuine IS present in the sample is covered end to end; a partial
# hit is usually a conserved domain shared with a different family.
DEFAULT_MIN_COVERED_PERCENT = 90.0

# Depth below this multiple of the genome baseline is treated as "not present".
# Set below 1.0 on purpose: a single-copy IS should sit at ~1x, and sampling
# noise plus mapping loss regularly pushes a real single copy to 0.6-0.8x.
DEFAULT_MIN_COPY_NUMBER = 0.5


# ── Reading BBMap's coverage table ───────────────────────────────────────────

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


# ── Small shared helpers ─────────────────────────────────────────────────────
# A cell that will not parse as a number becomes zero rather than raising, so one
# malformed row in a coverage table cannot take the whole sample down. Zero depth
# and zero length fall below every threshold here, so such a row drops out rather
# than corrupting the estimate.

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


# ── The denominator: what does single-copy sequence look like? ───────────────

def genome_baseline_depth(assembly_covstats):
    """The read depth of a SINGLE-COPY region, which is the denominator.

    Takes in: covstats from mapping the reads against the sample's own assembly,
              one row per contig.
    Returns:  a float depth, or None when it cannot be established.

    Uses the LENGTH-WEIGHTED MEDIAN across contigs rather than the mean, and the
    reason matters. A bacterial assembly is one long chromosome plus a handful of
    short contigs, and those short contigs are exactly the ones with unstable
    depth — repeat collapses sit at several times the genome depth, and low-
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


# ── ISOSDB's family annotation, and what ISEScan already located ─────────────
# The two tables that let a per-entry depth be rolled up to an IS FAMILY, which
# is the level the summary reports at. Per-entry numbers ride along as context
# only — see the ambiguous-read note in estimate_copies.

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

    This is the number the read-based estimate is set against, and the comparison
    is the deliverable: located is what survived assembly, estimated is what the
    reads say was there.

    The family column is called is_family in the tidy table isescan_to_table.py
    writes and family in ISEScan's own raw output, so both names are accepted and
    either file works. A table carrying neither column returns nothing, and every
    family is then reported with located_copies = 0: the estimate still runs, but
    the collapse delta it produces is the estimate itself and means nothing.
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


# ── What we write ────────────────────────────────────────────────────────────

def audit_row(sample, action, reason, detail, element="NA"):
    return {
        "sample": sample, "element": element,
        "action": action, "reason": reason, "detail": detail,
    }


AUDIT_COLUMNS = ["sample", "element", "action", "reason", "detail"]

# EVERY `action` / `reason` PAIR THIS SCRIPT CAN WRITE. `action` says what
# happened, `reason` says why. None of it changes an AMR gene's mobility tier —
# the whole leg is a QC measurement of how badly the assembler collapsed the IS
# copies, reported alongside the located count.
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
#                         failure — the located count simply stands on its own
#     no_genome_baseline_depth          assembly depth could not be established,
#                                       so there is no single-copy baseline to
#                                       divide by. Nothing can be estimated.
#     isosdb_does_not_cover_this_organism  no ISOSDB entry was covered by reads
#                                       for ANY family this sample has, so the
#                                       leg produced nothing at all
#     family_absent_from_isosdb         no ISOSDB entry for this one family was
#                                       covered by reads, so its collapse cannot
#                                       be measured even though other families'
#                                       can. Says the database lacks this
#                                       organism's version of the family, NOT
#                                       that nothing collapsed.
#     isosdb_detected_fewer_than_located  the read-mapping leg found FEWER copies
#                                       than ISEScan located on the contigs. The
#                                       delta is reported as NA rather than as a
#                                       negative number, because a negative
#                                       "collapse" is not meaningful — it means
#                                       the database is the limiting factor here,
#                                       not the assembly.
#
#   action=summary        ONE closing row per SAMPLE (not per family), recording
#                         what the leg managed overall
#     copy_number_estimate_complete

# The deliverable: one row per IS family, written to {sample}_is_copy_number.tsv
# and read by a person. No rule consumes it.
SUMMARY_COLUMNS = [
    "sample", "is_family",
    "located_copies",          # copies ISEScan found on the contigs
    "estimated_copies",        # copies the read depth implies; NA if none detected
    "collapse_delta",          # estimated - located; NA when db_informative is FALSE
    "db_informative",          # TRUE only when ISOSDB detected this family AND
                               # the estimate reached the located count
    "n_db_entries_detected",   # how many ISOSDB entries contributed to the estimate
    "max_entry_copy_number",   # the deepest single entry, for context
    "genome_baseline_depth",   # the single-copy denominator, identical on every row
]


# ── Depths in, per-family copy numbers out ───────────────────────────────────

def estimate_copies(sample, isosdb_covstats, assembly_covstats, family_map,
                    located_counts,
                    min_covered_percent=DEFAULT_MIN_COVERED_PERCENT,
                    min_copy_number=DEFAULT_MIN_COPY_NUMBER):
    """Turn read depths into a per-family copy-number estimate.

    Returns (summary_rows, audit_rows).

    Why the FAMILY total is the headline and the per-entry split is not. ISOSDB
    is dereplicated at 95% identity, but IS families remain similar enough that
    one real element can attract reads across several database entries. BBMap is
    run with ambiguous=best (rule isosdb_map) so each read lands on one entry
    only, which stops the total from being multiplied — but it also means the
    split between near-identical entries is arbitrary. The family sum is robust
    to where an ambiguous read landed; the per-entry numbers are context.
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
            # that element — it is usually a conserved domain shared with another
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

        # The check this leg cannot work without.
        #
        # ISEScan finds IS elements by profile HMM, which recognises a FAMILY.
        # This leg maps reads to ISOSDB, which needs nucleotide identity to one
        # specific catalogued element. Those are not the same sensitivity, and the
        # spec (§2.2) quantifies the gap from the ISOSDB paper itself: 97.5% of
        # its transposases have protein homologs in ISfinder but only 37.9% have
        # nucleotide ones.
        #
        # So for an organism ISOSDB does not cover — which for a non-clinical,
        # environmental isolate is the NORMAL case — no entry attracts reads, the
        # estimate comes out at zero, and the delta goes NEGATIVE. Read naively
        # that says "the assembly collapsed nothing, we over-called", when the
        # truth is "this database has nothing to say about this organism".
        # Measured on a hybrid validation isolate: 17 IS located by
        # ISEScan, 1 of 22,713 ISOSDB entries covered end to end, delta -7.3 on
        # the IS3 family alone.
        #
        # Hence both halves of the test below. `entries > 0` alone is not enough:
        # on that same sample the IS3 family had 9 elements located and exactly
        # ONE ISOSDB entry detected, which passes `entries > 0` and would have
        # published a delta of -7.3 that measures the database, not the assembly.
        # Requiring `estimated >= located` as well is what stops it.
        #
        # Absence of a nucleotide match is not evidence of absence of copies, so
        # a family with nothing detected is reported NA rather than as a number
        # that looks like a measurement, and a negative delta is never reported as
        # a delta at all. The leg can only ever say "there are AT LEAST this many
        # more copies than you located" — never "there are fewer", because
        # estimated < located always means ISOSDB missed some, never that ISEScan
        # over-called.
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

    # One closing audit row per sample, so the audit file alone reconstructs the
    # run. `families_informative` counts a family as one ISOSDB could speak to
    # whenever ANY entry was detected — the looser of the two tests, dropping the
    # `estimated >= located` half the db_informative column also requires. So it
    # can come out higher than the db_informative=TRUE count, which is the number
    # main() prints when it finishes. The two disagreeing is normal, not a bug.
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


# ── Writing the two outputs ──────────────────────────────────────────────────

def write_tsv(path, columns, rows):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "NA")) for column in columns) + "\n")


# ── Command line ─────────────────────────────────────────────────────────────
# rule is_copy_number passes every flag; the defaults exist for a hand run. The
# closing message says whether ISOSDB could quantify the collapse at all, because
# "no estimate" is the common outcome outside ISOSDB's sampling and a reader who
# only sees the table would not know why it is empty.

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
