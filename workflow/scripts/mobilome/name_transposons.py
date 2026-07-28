#!/usr/bin/env python3
"""Name curated transposons and integrons on the contigs, so AMR genes inside
them can reach mobility tier 4.

WHY THIS EXISTS
    The mobility ladder has six rungs, and until now rung 4 was unreachable:

        4  the gene sits inside a NAMED unit transposon or integron cassette

    colocalise.py has always known how to award tier 4 - it looks for an element
    of type `unit_transposon` or `integron` that contains the gene - but nothing
    in the workflow ever produced such an element, so the branch was dead code.
    Genes that are biologically inside a named transposon surfaced at tier 3
    (pattern-matched composite) or tier 5 (just "on a plasmid") instead: not
    wrong, but less precise than the evidence allows.

    This script closes that gap. It reads a BLAST of the sample's contigs against
    TnCentral, keeps the hits good enough to call a name, and writes them in the
    element-table shape colocalise.py already consumes.

WHY A CURATED NAME IS WORTH A WHOLE TIER
    Tier 3 is INFERRED: two IS copies of the same family, the right distance
    apart, with a gene in between, so we call it a composite transposon. That
    inference has known failure modes - IS26 in particular forms translocatable
    units with its copies in DIRECT orientation, breaking the same-orientation
    rule the pattern relies on.

    A TnCentral hit is not an inference. It is a match to an element somebody
    characterised, named and deposited, whose architecture is already known - so
    it gets IS26 right for free. That is why the spec (§7) says a curated hit
    should OVERRIDE the pattern-based call rather than merely agree with it.

WHAT IT TAKES IN
    A BLAST tabular file (outfmt 6 with the columns listed in BLAST_COLUMNS),
    produced by rule mobilome_name_transposons: the sample's contigs queried
    against the TnCentral nucleotide database.

WHAT IT PRODUCES
    - an element TSV in colocalise.py's `--is-table` shape, so it is passed as
      one more element table alongside the ISEScan and CONJscan ones;
    - an audit TSV giving a reason for every hit NOT kept, per the project rule
      that filtering decisions are never silent.

WHAT IT DELIBERATELY DOES NOT DO
    TnCentral also contains plain insertion sequences (deflines starting IS...).
    Those are skipped, with an audit line: ISEScan already inventories IS
    elements on these contigs, and an IS by definition carries no passenger genes
    (spec §2.4), so it cannot put an AMR gene inside a named transposon. Emitting
    them here would double-count the same element in two tables.
"""

import argparse
import os
import sys


# The columns this script needs from `blastn -outfmt 6 ...`. The rule asks for
# exactly these, in this order, so the two must be changed together.
BLAST_COLUMNS = [
    "qseqid",    # our contig
    "sseqid",    # the TnCentral entry, e.g. Tn4401b-JX560992
    "pident",    # percent identity over the aligned part
    "length",    # alignment length
    "qstart", "qend",
    "sstart", "send",
    "evalue", "bitscore",
    "slen",      # length of the reference element - the denominator that matters
    "qlen",
]


# The naming cascade's thresholds (spec §5.4): a hit is only allowed to confer a
# NAME when it is both highly similar and nearly complete.
#
# Why coverage is measured against the SUBJECT and not the query: the question is
# "is the whole of this known element present here?", not "how much of our contig
# is covered?". A 5 kb transposon sitting in a 300 kb contig covers 1.7% of the
# query and 100% of the subject, and it is unambiguously present.
DEFAULT_MIN_IDENTITY = 90.0
DEFAULT_MIN_SUBJECT_COVERAGE = 0.80


def parse_element_name(subject_id):
    """Pull the element name out of a TnCentral defline.

    TnCentral names its entries `<NAME>-<ACCESSION>`, for example:
        Tn4401b-JX560992          -> Tn4401b
        In0-U49101                -> In0
        IS1133_Tn10_IS903B-CP000602.1 -> IS1133_Tn10_IS903B

    Split on the LAST hyphen, because element names themselves contain hyphens
    and underscores while the accession suffix does not. If there is no hyphen at
    all the whole string is taken as the name rather than dropping the hit - an
    unexpected defline shape should cost us the accession, not the element.
    """
    text = subject_id.strip()
    if "-" not in text:
        return text, ""
    name, accession = text.rsplit("-", 1)
    return name, accession


def classify_element(name):
    """Decide what KIND of thing a TnCentral name refers to.

    TnCentral mixes three sorts of entry and they do not all mean the same thing
    for the mobility ladder:

      Tn...  a transposon. Carries passenger genes - which is the entire point,
             and why it can put an AMR gene inside a named element -> tier 4.
      In...  an integron. A gene-capture platform carrying cassettes, likewise
             a named architecture around a resistance gene -> tier 4.
      IS...  a plain insertion sequence. By definition carries ONLY what it needs
             to move (spec §2.4), so it never contains a passenger AMR gene.
             ISEScan already inventories these, so naming them here would put the
             same element in two tables.

    Returns one of "unit_transposon", "integron", "insertion_sequence" or
    "unknown" - the caller decides what to do with each.
    """
    upper = name.upper()
    if upper.startswith("TN"):
        return "unit_transposon"
    if upper.startswith("IN"):
        return "integron"
    if upper.startswith("IS"):
        return "insertion_sequence"
    return "unknown"


def read_blast_hits(path):
    """Read the BLAST tabular file into a list of dicts.

    A missing or empty file is NOT an error: a genome with no curated transposon
    on it is a perfectly ordinary result, and the module degrades gracefully
    rather than hard-failing (a standing BacFlux convention).
    """
    hits = []
    if not path or not os.path.isfile(path):
        return hits

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < len(BLAST_COLUMNS):
                continue          # truncated line; skip rather than crash
            record = dict(zip(BLAST_COLUMNS, fields))
            hits.append(record)
    return hits


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


def subject_coverage(hit):
    """What fraction of the REFERENCE element this alignment covers.

    BLAST reports one HSP per line, and a real transposon hit is often split into
    several HSPs by internal indels. Using a single HSP's length therefore
    UNDERSTATES coverage - which is why merge_hits_per_element below sums the
    covered subject span across HSPs before this threshold is applied.
    """
    slen = to_int(hit.get("slen"))
    if slen <= 0:
        return 0.0
    return min(1.0, to_int(hit.get("length")) / slen)


def merge_span(intervals):
    """Total length covered by a set of (start, end) intervals, overlaps counted once.

    Needed because several BLAST HSPs against the same reference element usually
    overlap slightly; naively adding their lengths would push coverage past 100%
    and let a poor hit through the threshold.
    """
    if not intervals:
        return 0
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    covered = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start + 1
            current_start, current_end = start, end
    covered += current_end - current_start + 1
    return covered


# How far apart two HSPs against the SAME reference element may sit on the contig
# and still be treated as one copy of it, as a multiple of that element's length.
#
# This was 1.0 - a whole element length - and that was far too generous. Two
# genuine copies of a transposon sitting less than one length apart were merged
# into one element spanning both of them PLUS the chromosome in between, and
# every gene in that gap became tier-4 cargo.
#
# Measured on the KPNIH1 control: every genuine copy has a largest internal HSP
# gap of <= 140 bp, because what separates the pieces of ONE copy is an indel or
# a small internal insertion. The gaps that marked a wrongly merged pair were
# thousands of bp. So the allowance is now sized to internal indels, which is
# what it was always meant to model, with a floor for very short references.
SAME_COPY_GAP_MULTIPLE = 0.10
SAME_COPY_GAP_FLOOR_BP = 500

# The fraction of the reported interval that must actually be ALIGNED to the
# reference. This is the guard the span check below could never be.
#
# Real case from the KPNIH1 control. Tn3000 (3,235 bp) matched NZ_CP006662.2 with
# two tiny terminal inverted-repeat HSPs (84 bp and 146 bp) at ~27,200 belonging
# to a NEIGHBOURING element, plus the real copy at 29,785-32,882. The IRs pulled
# the reported interval out to 27,185-32,882, and 2,454 bp of that 5,698 bp span
# - 43% - has no alignment to Tn3000 at all; Bakta annotates it as a complete,
# unrelated IS66-family element. Terminal inverted repeats are SHARED between
# related transposons, so this is a systematic trap rather than bad luck.
#
# Measured over the whole control: five of the six reported elements have an
# aligned fraction of 1.00 and only the spurious Tn3000 falls to 0.57, so this
# threshold removes exactly the wrong one and keeps every right one.
MIN_ALIGNED_FRACTION = 0.70

# A final sanity bound on the reported element, as a multiple of the reference.
#
# Was 3.0, which was mathematically incapable of catching the case it was written
# for: with the old gap allowance two copies could sit one length apart, giving a
# span of at most 2L + L = 3L, so the guard never fired on a two-copy merge. It
# only triggered at three or more copies, where it then discarded all of them.
# With the tightened gap above, one genuine copy spans about L even with internal
# indels, so a bound just above that can actually do some work.
MAX_SPAN_MULTIPLE_OF_REFERENCE = 2.0


def cluster_hsps_by_position(group, reference_length):
    """Split HSPs against one reference element into separate COPIES.

    THE BUG THIS EXISTS TO PREVENT, because it is not hypothetical - it was found
    on the KPNIH1 positive control before this script was ever wired in:

        Tn7246 (7,325 bp reference) hit the chromosome at ~3,544,000 and again at
        ~4,480,000. Taking min(qstart) and max(qend) over all HSPs merged those
        two copies into ONE element spanning 942,502 bp - 129 times the length of
        the transposon itself. colocalise.py would then have called every AMR
        gene in that 942 kb "inside a named transposon", tier 4.

    Multiple copies of one transposon in a genome are the normal case, not an
    edge case, so this has to be handled rather than guarded against.

    Takes in: the HSPs for one (contig, reference element) pair, and the
              reference's length.
    Does:     sorts them along the contig and starts a new copy whenever the gap
              to the running end exceeds the element's own length.
    Returns:  a list of HSP lists, one per copy.
    """
    if not group:
        return []

    max_gap = max(SAME_COPY_GAP_FLOOR_BP, int(reference_length * SAME_COPY_GAP_MULTIPLE))
    ordered = sorted(group, key=lambda h: min(to_int(h["qstart"]), to_int(h["qend"])))

    clusters = [[ordered[0]]]
    running_end = max(to_int(ordered[0]["qstart"]), to_int(ordered[0]["qend"]))
    for hit in ordered[1:]:
        start = min(to_int(hit["qstart"]), to_int(hit["qend"]))
        end = max(to_int(hit["qstart"]), to_int(hit["qend"]))
        if start - running_end > max_gap:
            clusters.append([hit])          # far away: a different copy
        else:
            clusters[-1].append(hit)
        running_end = max(running_end, end)
    return clusters


def merge_hits_per_element(hits):
    """Collapse many HSPs into one candidate per COPY of a reference element.

    Takes in: raw BLAST rows.
    Does:     groups by contig and subject, splits each group into separate
              copies by position (see cluster_hsps_by_position), then for each
              copy sums the covered subject span - so a hit broken into pieces by
              indels is judged on the whole of what it covers - and takes the
              copy's query span as its coordinates.
    Returns:  a list of merged candidate dicts, one per copy.

    Why merge HSPs at all: without it, a 6 kb transposon split into three HSPs
    looks like three separate 33%-coverage hits and is thrown out by the coverage
    threshold, so the very elements most worth naming - the big, mosaic,
    clinically interesting ones - would be the ones systematically missed.

    Why merge only WITHIN a copy: see cluster_hsps_by_position.
    """
    groups = {}
    for hit in hits:
        key = (hit["qseqid"], hit["sseqid"])
        groups.setdefault(key, []).append(hit)

    merged = []
    for (contig, subject), group in groups.items():
        slen = max(to_int(h["slen"]) for h in group)
        for copy_hsps in cluster_hsps_by_position(group, slen):
            subject_intervals = [(to_int(h["sstart"]), to_int(h["send"])) for h in copy_hsps]
            query_positions = [to_int(h["qstart"]) for h in copy_hsps] + \
                              [to_int(h["qend"]) for h in copy_hsps]
            covered_bp = merge_span(subject_intervals)
            best = max(copy_hsps, key=lambda h: to_float(h["bitscore"]))

            # How much of the interval we are about to report is actually
            # ALIGNED. Terminal inverted repeats are shared between related
            # transposons, so a couple of short IR hits belonging to a
            # NEIGHBOURING element can stretch the interval far past the real
            # copy - see MIN_ALIGNED_FRACTION. Measuring it here is what makes
            # that detectable at all.
            query_intervals = [(to_int(h["qstart"]), to_int(h["qend"])) for h in copy_hsps]
            aligned_bp = merge_span(query_intervals)
            span_bp = max(query_positions) - min(query_positions) + 1

            merged.append({
                "contig": contig,
                "subject": subject,
                "start": min(query_positions),
                "end": max(query_positions),
                "aligned_bp": aligned_bp,
                "aligned_fraction": (aligned_bp / span_bp) if span_bp > 0 else 0.0,
                "identity": to_float(best["pident"]),
                "bitscore": to_float(best["bitscore"]),
                "evalue": best.get("evalue", "NA"),
                "subject_length": slen,
                "covered_bp": covered_bp,
                "coverage": (covered_bp / slen) if slen > 0 else 0.0,
                "n_hsps": len(copy_hsps),
            })
    return merged


def overlaps(a_start, a_end, b_start, b_end):
    """Do two closed intervals share at least one base?"""
    return a_start <= b_end and b_start <= a_end


def drop_redundant_overlaps(candidates):
    """Keep one named element per stretch of contig.

    TnCentral entries are nested and mosaic on purpose - Tn4401b contains Tn2,
    which contains a transposase that also appears inside other elements - so one
    genuine transposon in the assembly routinely matches half a dozen reference
    entries at once. Reporting all of them would multiply one biological element
    into several "named elements" and make the AMR table look far busier than the
    genome is.

    So: sort by bitscore, keep the best, and drop any later candidate that
    overlaps an already-kept one on the same contig. The dropped ones are
    audited, never silently discarded, so a reader can see the nesting.
    """
    kept = []
    dropped = []
    for candidate in sorted(candidates, key=lambda c: c["bitscore"], reverse=True):
        clash = None
        for existing in kept:
            if existing["contig"] != candidate["contig"]:
                continue
            if overlaps(existing["start"], existing["end"],
                        candidate["start"], candidate["end"]):
                clash = existing
                break
        if clash is None:
            kept.append(candidate)
        else:
            candidate["superseded_by"] = clash["name"]
            dropped.append(candidate)
    return kept, dropped


def audit_row(sample, action, reason, detail, contig="NA", start="NA", end="NA"):
    """One line of the decision trail.

    Same shape as the other mobilome audit files so they can be concatenated and
    read together: every filtering decision gets a row with a reason column.
    """
    return {
        "sample": sample,
        "contig": contig,
        "start": str(start),
        "end": str(end),
        "action": action,
        "reason": reason,
        "detail": detail,
    }


ELEMENT_COLUMNS = [
    "sample", "contig", "start", "end", "strand",
    "element_type", "mge_id", "mge_name",
    "identity", "subject_coverage", "aligned_fraction",
    "reference_accession", "reference_length_bp",
    "n_hsps", "confidence",
]

AUDIT_COLUMNS = ["sample", "contig", "start", "end", "action", "reason", "detail"]


def build_elements(sample, hits, min_identity, min_coverage):
    """Turn BLAST hits into named element rows, auditing everything dropped.

    Takes in: raw BLAST rows of contigs vs TnCentral.
    Does:     merge HSPs per reference element, apply the identity and coverage
              thresholds, discard the entries that cannot carry a passenger gene
              (plain IS), de-duplicate nested matches.
    Returns:  (element_rows, audit_rows).
    """
    audit_rows = []
    if not hits:
        audit_rows.append(audit_row(
            sample, "not_applicable", "no_tncentral_hits",
            "BLAST returned no hit against TnCentral for this sample, so no AMR "
            "gene can reach tier 4 by a curated name. This is a normal result: "
            "most genomes carry no characterised transposon."))
        return [], audit_rows

    candidates = merge_hits_per_element(hits)
    passing = []

    for candidate in candidates:
        name, accession = parse_element_name(candidate["subject"])
        candidate["name"] = name
        candidate["accession"] = accession
        candidate["element_type"] = classify_element(name)

        if candidate["element_type"] == "insertion_sequence":
            audit_rows.append(audit_row(
                sample, "not_applicable", "tncentral_hit_is_a_plain_is",
                f"{name}: TnCentral entry is an insertion sequence, not a "
                "transposon or integron. An IS carries only what it needs to move "
                "and so cannot hold a passenger AMR gene, and ISEScan already "
                "inventories IS elements on these contigs - naming it here would "
                "list the same element twice.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        if candidate["element_type"] == "unknown":
            audit_rows.append(audit_row(
                sample, "discarded", "tncentral_name_not_recognised",
                f"{name}: the TnCentral defline does not begin Tn, In or IS, so "
                "this script cannot tell what kind of element it is and will not "
                "guess. Check the entry by hand if the hit looks interesting.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        if candidate["identity"] < min_identity:
            audit_rows.append(audit_row(
                sample, "discarded", "identity_below_naming_threshold",
                f"{name}: {candidate['identity']:.1f}% identity is below the "
                f"{min_identity:.0f}% a curated NAME requires. A weaker match may "
                "still be a relative of this element, but calling it by this name "
                "would overstate what the sequence shows.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        # THE INTERVAL MUST ACTUALLY BE THE ELEMENT. Terminal inverted repeats
        # are shared between related transposons, so short IR hits belonging to a
        # NEIGHBOURING element cluster with the real copy and drag the reported
        # interval across DNA that has nothing to do with this transposon. On the
        # KPNIH1 control that put 2,454 bp of an unrelated IS66 element inside a
        # "Tn3000" interval - and colocalise.py would call any AMR gene in there
        # cargo of Tn3000, at high confidence, with no audit line to explain it.
        # Neither the identity, coverage nor span checks can see this: the
        # aligned parts match perfectly, there just are not enough of them.
        if candidate["aligned_fraction"] < MIN_ALIGNED_FRACTION:
            audit_rows.append(audit_row(
                sample, "discarded", "interval_mostly_unaligned",
                f"{name}: only {100 * candidate['aligned_fraction']:.0f}% of the "
                f"{candidate['end'] - candidate['start'] + 1} bp interval aligns to "
                f"the reference (needs {100 * MIN_ALIGNED_FRACTION:.0f}%). Short "
                "terminal inverted repeats shared with a neighbouring element "
                "stretch the interval over DNA that is not part of this "
                "transposon, so naming it would put unrelated genes inside a "
                "named element.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        # A copy of a 7 kb transposon should occupy about 7 kb of contig. If the
        # span is wildly larger, the HSPs were clustered into one copy when they
        # belong to several, and passing that interval on would make every AMR
        # gene inside a huge stretch of chromosome "cargo of a named transposon".
        # cluster_hsps_by_position should already have prevented this; the check
        # stays because the failure is silent and expensive.
        span_bp = candidate["end"] - candidate["start"] + 1
        span_limit = candidate["subject_length"] * MAX_SPAN_MULTIPLE_OF_REFERENCE
        if candidate["subject_length"] > 0 and span_bp > span_limit:
            audit_rows.append(audit_row(
                sample, "discarded", "element_span_implausible_for_reference",
                f"{name}: the matched region spans {span_bp} bp but the reference "
                f"element is only {candidate['subject_length']} bp "
                f"(limit {MAX_SPAN_MULTIPLE_OF_REFERENCE:.0f}x). That means "
                "separate copies of the element were joined into one interval, "
                "which would wrongly make everything between them look like cargo "
                "of a single named transposon. Dropped rather than reported.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        if candidate["coverage"] < min_coverage:
            audit_rows.append(audit_row(
                sample, "discarded", "reference_coverage_below_threshold",
                f"{name}: only {100 * candidate['coverage']:.0f}% of the "
                f"{candidate['subject_length']} bp reference element is present "
                f"(needs {100 * min_coverage:.0f}%). A fragment of a transposon is "
                "not that transposon - on a short-read assembly this usually means "
                "the element is split across contigs.",
                contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))
            continue

        passing.append(candidate)

    kept, redundant = drop_redundant_overlaps(passing)
    for candidate in redundant:
        audit_rows.append(audit_row(
            sample, "discarded", "nested_or_overlapping_tncentral_hit",
            f"{candidate['name']} overlaps the better-scoring "
            f"{candidate['superseded_by']} on the same contig and is not listed "
            "separately. TnCentral entries are deliberately nested (a large "
            "transposon contains smaller ones), so one real element matches "
            "several references; only the best-scoring one is reported.",
            contig=candidate["contig"], start=candidate["start"], end=candidate["end"]))

    element_rows = []
    for candidate in sorted(kept, key=lambda c: (c["contig"], c["start"])):
        element_rows.append({
            "sample": sample,
            "contig": candidate["contig"],
            "start": str(candidate["start"]),
            "end": str(candidate["end"]),
            # BLAST gives an orientation, but the element's own strand is not a
            # meaningful concept for a transposon in this table and colocalise
            # only uses strand for the IS-orientation tests, so leave it unset
            # rather than inventing one.
            "strand": ".",
            "element_type": candidate["element_type"],
            "mge_id": f"{candidate['contig']}|{candidate['element_type']}-"
                      f"{candidate['start']}:{candidate['end']}",
            "mge_name": candidate["name"],
            "identity": f"{candidate['identity']:.1f}",
            "subject_coverage": f"{candidate['coverage']:.3f}",
            # How much of the reported interval is really aligned. 1.000 is the
            # normal, healthy value; anything lower means the interval is padded.
            "aligned_fraction": f"{candidate['aligned_fraction']:.3f}",
            "reference_accession": candidate["accession"] or "NA",
            "reference_length_bp": str(candidate["subject_length"]),
            "n_hsps": str(candidate["n_hsps"]),
            # A curated hit at these thresholds is strong evidence of identity.
            # It says nothing about whether the element is COMPLETE in this
            # assembly, which the coverage column reports separately.
            "confidence": "high",
        })

    audit_rows.append(audit_row(
        sample, "summary", "tncentral_naming_complete",
        f"{len(element_rows)} named element(s) kept from {len(candidates)} "
        f"reference match(es) at >={min_identity:.0f}% identity and "
        f">={100 * min_coverage:.0f}% reference coverage."))

    return element_rows, audit_rows


def write_tsv(path, columns, rows):
    """Write a TSV with a header, creating the parent directory if needed.

    Always writes the header even with zero rows: an empty table with a header is
    readable by everything downstream, while a missing file breaks the join.
    """
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
            "Name curated transposons and integrons from a BLAST of the sample's "
            "contigs against TnCentral, so AMR genes inside them can reach "
            "mobility tier 4. Writes an element table in colocalise.py's shape."
        )
    )
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--blast", required=True,
                        help="BLAST tabular output (outfmt 6) of contigs vs "
                             "TnCentral, with the columns in BLAST_COLUMNS.")
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY,
                        help="Percent identity a hit needs before it may confer a "
                             f"curated name. Default {DEFAULT_MIN_IDENTITY}.")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_SUBJECT_COVERAGE,
                        help="Fraction of the REFERENCE element that must be "
                             "present. Measured against the subject, not the "
                             "query: the question is whether the whole known "
                             f"element is here. Default {DEFAULT_MIN_SUBJECT_COVERAGE}.")
    parser.add_argument("--out-table", required=True,
                        help="Named-element TSV (passed to colocalise.py as one "
                             "more --is-table).")
    parser.add_argument("--out-audit", required=True,
                        help="Decision trail: a reason for every hit not kept.")
    args = parser.parse_args(argv)

    hits = read_blast_hits(args.blast)
    elements, audit = build_elements(
        args.sample, hits, args.min_identity, args.min_coverage)

    write_tsv(args.out_table, ELEMENT_COLUMNS, elements)
    write_tsv(args.out_audit, AUDIT_COLUMNS, audit)

    n_transposon = sum(1 for row in elements if row["element_type"] == "unit_transposon")
    n_integron = sum(1 for row in elements if row["element_type"] == "integron")
    print(
        f"Sample {args.sample}: {len(elements)} named element(s) - "
        f"{n_transposon} transposon, {n_integron} integron "
        f"(from {len(hits)} BLAST HSP(s); reasons in {args.out_audit})."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
