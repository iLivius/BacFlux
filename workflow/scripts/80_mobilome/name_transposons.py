#!/usr/bin/env python3
"""Name curated transposons and integrons on the contigs, so AMR genes inside
them can reach mobility tier 4.

Why rung 4 of the ladder needs its own script
---------------------------------------------
The mobility ladder has six rungs, and rung 4 is:

    4  the gene sits inside a NAMED unit transposon or integron cassette

colocalise.py awards it whenever an element of type `unit_transposon` or
`integron` contains the gene, and this script is the only thing in the workflow
that produces such an element. Without it that branch never fires, and genes that
are biologically inside a named transposon are reported at tier 3 (pattern-matched
composite) or tier 5 (just "on a plasmid") instead: not wrong, but less precise
than the evidence allows.

The layer is opt-in, because naming needs the TnCentral database and BacFlux
never ships it — every rule here lives behind MOBILOME_NAME_ELEMENTS, which needs
the mobilome module switched on AND mobilome.tncentral.url or .dir set. Worth
knowing before reading a report: tier 4 has never actually been assigned in a
retained run — see the README's "What the benchmark does not show".

Why a curated name is worth a whole tier
----------------------------------------
Tier 3 is INFERRED: two IS copies of the same family, the right distance apart,
with a gene in between, so we call it a composite transposon. That inference has
known failure modes — IS26 in particular forms translocatable units with its
copies in DIRECT orientation, breaking the same-orientation rule the pattern
relies on.

A TnCentral hit is not an inference. It is a match to an element somebody
characterised, named and deposited, whose architecture is already known — so it
gets IS26 right for free. That is why the spec (§7) says a curated hit should
OVERRIDE the pattern-based call rather than merely agree with it.

What it takes in and what it produces
-------------------------------------
In:  a BLAST tabular file (outfmt 6, with the columns listed in BLAST_COLUMNS)
     written by rule tncentral_blast in shared/80_mobilome.smk — this sample's
     contigs queried against the TnCentral nucleotide database. That search
     applies only a loose e-value, on purpose: the naming thresholds live here,
     so every near miss arrives and can be judged and audited in one place.
Out: an element TSV in colocalise.py's --is-table shape, handed to
     amr_mge_colocalisation as one more element table alongside the ISEScan and
     CONJscan ones; plus an audit TSV giving a reason for every hit NOT kept, per
     the project rule that filtering decisions are never silent.

Four names for one step, so nobody has to guess which is which: the rule is
name_elements, the script it runs is this one, the table is NAMED_ELEMENTS_TABLE
and the log is mobilome_name_elements_*.log. Renaming any of them would change
job names and log filenames on disk for no gain, so they stay as they are.

What it deliberately does not do
--------------------------------
TnCentral also contains plain insertion sequences (deflines starting IS...).
Those are skipped, with an audit line: ISEScan already inventories IS elements on
these contigs, and an IS by definition carries no passenger genes (spec §2.4), so
it cannot put an AMR gene inside a named transposon. Emitting them here would
double-count the same element in two tables.
"""

import argparse
import os
import sys


# ── What the BLAST table gives us ────────────────────────────────────────────
# The columns this script needs from `blastn -outfmt 6 ...`. Rule tncentral_blast
# asks for exactly these, in this order, so the two must be changed together — a
# drift shows up as short lines and is audited as blast_line_unparsable.
BLAST_COLUMNS = [
    "qseqid",    # our contig
    "sseqid",    # the TnCentral entry, e.g. Tn4401b-JX560992
    "pident",    # percent identity over the aligned part
    "length",    # alignment length
    "qstart", "qend",
    "sstart", "send",
    "evalue", "bitscore",
    "slen",      # length of the reference element — the denominator that matters
    "qlen",
]


# ── How good a hit must be before it may confer a name ───────────────────────
# The naming cascade's thresholds (spec §5.4): a hit is only allowed to confer a
# NAME when it is both highly similar and nearly complete. Both numbers are the
# cascade's convention, not a biological boundary. Rule name_elements always
# passes them explicitly, from mobilome.tncentral.min_identity and
# .min_reference_coverage, which fall back to these same values in
# shared/00_common.smk; the defaults here are what a hand-run of the script gets.
#
# Why coverage is measured against the SUBJECT and not the query: the question is
# "is the whole of this known element present here?", not "how much of our contig
# is covered?". A 5 kb transposon sitting in a 300 kb contig covers 1.7% of the
# query and 100% of the subject, and it is unambiguously present.
DEFAULT_MIN_IDENTITY = 90.0
DEFAULT_MIN_SUBJECT_COVERAGE = 0.80


# ── Reading a TnCentral defline ──────────────────────────────────────────────

def parse_element_name(subject_id):
    """Pull the element name out of a TnCentral defline.

    TnCentral names its entries `<NAME>-<ACCESSION>`, for example:
        Tn4401b-JX560992          -> Tn4401b
        In0-U49101                -> In0
        IS1133_Tn10_IS903B-CP000602.1 -> IS1133_Tn10_IS903B

    Split on the LAST hyphen, because element names themselves contain hyphens
    and underscores while the accession suffix does not. The accession may also
    be empty (TnCentral ships deflines like ">Tn7246-"), so an unexpected shape
    should cost us the accession, not the element: with no hyphen at all the
    whole string is taken as the name rather than dropping the hit.
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

      Tn...  a transposon. Carries passenger genes — which is the entire point,
             and why it can put an AMR gene inside a named element -> tier 4.
      In...  an integron. A gene-capture platform carrying cassettes, likewise
             a named architecture around a resistance gene -> tier 4.
      IS...  a plain insertion sequence. By definition carries ONLY what it needs
             to move (spec §2.4), so it never contains a passenger AMR gene.
             ISEScan already inventories these, so naming them here would put the
             same element in two tables.

    The test is the first two letters of the name and nothing else, so a defline
    that starts some other way comes back "unknown" and is audited rather than
    guessed at. Returns one of "unit_transposon", "integron",
    "insertion_sequence" or "unknown" — build_elements decides what to do with
    each.
    """
    upper = name.upper()
    if upper.startswith("TN"):
        return "unit_transposon"
    if upper.startswith("IN"):
        return "integron"
    if upper.startswith("IS"):
        return "insertion_sequence"
    return "unknown"


# ── Reading the BLAST tabular file ───────────────────────────────────────────

def read_blast_hits(path):
    """Read the BLAST tabular file into a list of dicts.

    A missing or empty file is NOT an error: a genome with no curated transposon
    on it is a perfectly ordinary result, and the module degrades gracefully
    rather than hard-failing (a standing BacFlux convention).

    Returns (hits, skipped_line_numbers) so the caller can audit anything it
    could not parse instead of losing it silently.
    """
    hits = []
    skipped = []
    if not path or not os.path.isfile(path):
        return hits, skipped

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < len(BLAST_COLUMNS):
                # Skipping is right — one bad line must not lose the sample — but
                # it has to be VISIBLE. blastn does not write short lines, so if
                # this fires the file was truncated or the column list here and
                # the -outfmt in the rule have drifted apart, and the run would
                # otherwise report "no curated transposon" with total confidence.
                skipped.append(line_number)
                continue
            hits.append(dict(zip(BLAST_COLUMNS, fields)))
    return hits, skipped


# ── Measuring identity and coverage across several HSPs ──────────────────────
# to_float / to_int keep the parse forgiving: every value arrives as text from a
# BLAST line, and a field that cannot be read falls back to the default instead
# of raising, so one odd row cannot take down the sample.

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


def weighted_identity(hsps):
    """Percent identity over a whole copy, weighted by how much each HSP covers.

    BLAST reports identity per HSP. Taking the best one flatters the copy: the
    number then describes its most conserved fragment rather than the element as
    a whole. Weighting by alignment length gives "how similar is this copy to the
    reference, base for base", which is what the naming threshold is meant to ask.
    """
    total_length = sum(to_int(h["length"]) for h in hsps)
    if total_length <= 0:
        return 0.0
    weighted = sum(to_int(h["length"]) * to_float(h["pident"]) for h in hsps)
    return weighted / total_length


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


# ── Guards against joining separate copies into one element ──────────────────
# Three thresholds, all doing the same job from different angles: stop a set of
# BLAST hits being reported as one element that is really several, or as one
# element padded out with DNA that is not part of it. Each answers a failure seen
# on the positive control, and each records what set its value — a measurement,
# or the arithmetic showing the old value could never fire — because the failures
# they prevent are silent.

# How far apart two HSPs against the SAME reference element may sit on the contig
# and still be treated as one copy of it, as a multiple of that element's length.
#
# This was 1.0 — a whole element length — and that was far too generous. Two
# genuine copies of a transposon sitting less than one length apart were merged
# into one element spanning both of them PLUS the chromosome in between, and
# every gene in that gap became tier-4 cargo.
#
# Measured on the K. pneumoniae positive control: every genuine copy has a
# largest internal HSP gap of <= 140 bp, because what separates the pieces of ONE
# copy is an indel or a small internal insertion. The gaps that marked a wrongly
# merged pair were thousands of bp. So the allowance is now sized to internal
# indels, which is what it was always meant to model, with an absolute floor so
# that a very short reference still gets room for one.
SAME_COPY_GAP_MULTIPLE = 0.10
SAME_COPY_GAP_FLOOR_BP = 500

# The fraction of the reported interval that must actually be ALIGNED to the
# reference. This is the guard the span check below could never be.
#
# Real case from the K. pneumoniae positive control (the accession below pins it
# to ATCC BAA-2146). Tn3000 (3,235 bp) matched NZ_CP006662.2 with two tiny
# terminal inverted-repeat HSPs (84 bp and 146 bp) at ~27,200 belonging to a
# NEIGHBOURING element, plus the real copy at 29,785-32,882. The IRs pulled the
# reported interval out to 27,185-32,882, and 2,454 bp of that 5,698 bp span
# — 43% — has no alignment to Tn3000 at all; Bakta annotates it as a complete,
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


# ── Splitting HSPs into element copies, then merging each copy ───────────────

def cluster_hsps_by_position(group, reference_length):
    """Split HSPs against one reference element into separate COPIES.

    The bug this exists to prevent, because it is not hypothetical — it was found
    on the K. pneumoniae positive control before this script was ever wired in:

        Tn7246 (7,325 bp reference) hit the chromosome at ~3,544,000 and again at
        ~4,480,000. Taking min(qstart) and max(qend) over all HSPs merged those
        two copies into ONE element spanning 942,502 bp — 129 times the length of
        the transposon itself. colocalise.py would then have called every AMR
        gene in that 942 kb "inside a named transposon", tier 4.

    Multiple copies of one transposon in a genome are the normal case, not an
    edge case, so this has to be handled rather than guarded against.

    Takes in: the HSPs for one (contig, reference element) pair, and the
              reference's length.
    Does:     sorts them along the contig and starts a new copy whenever the gap
              to the running end is bigger than SAME_COPY_GAP_MULTIPLE of the
              reference length, or SAME_COPY_GAP_FLOOR_BP, whichever is larger.
              That allowance models an internal indel inside one copy, not the
              distance between two copies — see the constants above.
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
              copy sums the covered subject span — so a hit broken into pieces by
              indels is judged on the whole of what it covers — and takes the
              copy's query span as its coordinates.
    Returns:  a list of merged candidate dicts, one per copy.

    Why merge HSPs at all: BLAST reports one HSP per line, and a real transposon
    hit is usually broken into several by internal indels. Judging coverage on a
    single HSP's length therefore UNDERSTATES it badly — a 6 kb transposon split
    into three HSPs looks like three separate 33%-coverage hits and is thrown out
    by the coverage threshold, so the very elements most worth naming (the big,
    mosaic, clinically interesting ones) would be the ones systematically missed.
    The `subject_coverage` column written further down is computed from the
    MERGED span for exactly this reason.

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
            # merge_span counts overlapping subject intervals ONCE, so two HSPs
            # hitting the same part of the reference cannot add up to look like
            # two different parts of it. That is what keeps the summed coverage
            # honest: it measures how much of the reference is present, not how
            # many times we matched it.
            covered_bp = merge_span(subject_intervals)
            best = max(copy_hsps, key=lambda h: to_float(h["bitscore"]))

            # How much of the interval we are about to report is actually
            # ALIGNED. Terminal inverted repeats are shared between related
            # transposons, so a couple of short IR hits belonging to a
            # NEIGHBOURING element can stretch the interval far past the real
            # copy — see MIN_ALIGNED_FRACTION. Measuring it here is what makes
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
                # Length-weighted identity across the pieces of THIS copy, not
                # the identity of its best fragment. A copy made of a 3 kb HSP at
                # 99.9% and a 2 kb HSP at 85% is not a 99.9% match to the
                # reference, and reporting it as one lets sequence that would
                # never have passed the naming threshold on its own ride along
                # inside an element labelled almost perfect.
                "identity": weighted_identity(copy_hsps),
                "bitscore": to_float(best["bitscore"]),
                "evalue": best.get("evalue", "NA"),
                "subject_length": slen,
                "covered_bp": covered_bp,
                "coverage": (covered_bp / slen) if slen > 0 else 0.0,
                "n_hsps": len(copy_hsps),
            })
    return merged


# ── Keeping one named element per stretch of contig ──────────────────────────

def overlaps(a_start, a_end, b_start, b_end):
    """Do two closed intervals share at least one base?"""
    return a_start <= b_end and b_start <= a_end


def drop_redundant_overlaps(candidates):
    """Keep one named element per stretch of contig.

    TnCentral entries are nested and mosaic on purpose — Tn4401b contains Tn2,
    which contains a transposase that also appears inside other elements — so one
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


# ── What we write: the element table, the audit, and its vocabulary ──────────

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


# The element table, one row per kept copy. colocalise.py reads it as one more
# --is-table and matches columns BY NAME, not by position, so the order here is
# for the human reader; the ones it actually picks up are contig, start, end,
# strand, element_type, mge_id and mge_name. Everything else is the naming
# evidence someone needs to judge the call. All coordinates are 1-based and
# inclusive, as BLAST reports them.
#
#   sample               the --sample argument, so tables can be concatenated
#   contig               BLAST qseqid — the contig the copy sits on
#   start, end           the copy's span on the contig
#                        source: min/max query position over its HSPs
#   strand               always "." — see the note where the row is built
#   element_type         unit_transposon | integron
#                        source: classify_element() on the TnCentral name
#   mge_id               contig|element_type-start:end, the spec §9 ID format
#   mge_name             the curated name, e.g. Tn4401b
#                        source: parse_element_name() on the BLAST sseqid
#   identity             % identity of the whole copy, length-weighted over HSPs
#   subject_coverage     how much of the REFERENCE element is present, 0-1
#                        source: merged subject intervals / slen
#   aligned_fraction     how much of the reported interval really aligns, 0-1;
#                        1.000 is normal, lower means the interval is padded
#   reference_accession  the accession half of the defline, or NA
#   reference_length_bp  BLAST slen — the curated element's full length
#   n_hsps               how many HSPs were merged into this copy
#   confidence           always "high" — see the note where the row is built
ELEMENT_COLUMNS = [
    "sample", "contig", "start", "end", "strand",
    "element_type", "mge_id", "mge_name",
    "identity", "subject_coverage", "aligned_fraction",
    "reference_accession", "reference_length_bp",
    "n_hsps", "confidence",
]

# The audit table: one row per decision, kept or refused, with the reason in its
# own column. contig/start/end are "NA" on the rows that are about the file as a
# whole rather than about one hit.
AUDIT_COLUMNS = ["sample", "contig", "start", "end", "action", "reason", "detail"]

# EVERY `action` / `reason` PAIR THIS SCRIPT CAN WRITE. This is the layer that
# unlocks mobility tier 4, so a missing name is the difference between "inside
# TnX" and "just on a plasmid" — which is exactly why every refusal is recorded.
#
#   action=discarded       the BLAST hit did NOT become a named element
#     blast_line_unparsable            a line had fewer columns than expected;
#                                      means the file was truncated or the -outfmt
#                                      in the rule drifted from BLAST_COLUMNS
#     tncentral_name_not_recognised    the defline did not parse to a known kind
#                                      of element (Tn.../In.../IS...)
#     identity_below_naming_threshold  too diverged to carry the name
#     reference_coverage_below_threshold  too little of the curated element is
#                                      present. THE COMMON ONE — a fragment of a
#                                      transposon is not that transposon, and this
#                                      is why tier 4 is rarely reached in practice
#     interval_mostly_unaligned        the reported interval is padded with DNA
#                                      that does not align to the reference,
#                                      usually via terminal inverted repeats
#                                      shared with a neighbouring element
#     element_span_implausible_for_reference  the interval is far longer than the
#                                      reference element could account for
#     nested_or_overlapping_tncentral_hit  a better hit already covers this span
#
#   action=not_applicable  nothing was attempted, for a stated reason
#     no_tncentral_hits                BLAST found nothing. A normal result: most
#                                      genomes carry no characterised transposon
#     tncentral_hit_is_a_plain_is      the hit is an insertion sequence, which by
#                                      definition carries no passenger gene, so it
#                                      cannot put an AMR gene inside a named
#                                      element. ISEScan already inventories these
#
#   action=summary         one closing row recording the thresholds actually used
#     tncentral_naming_complete


# ── Applying the naming thresholds and recording every refusal ───────────────

def build_elements(sample, hits, min_identity, min_coverage, skipped_lines=None):
    """Turn BLAST hits into named element rows, auditing everything dropped.

    Takes in: raw BLAST rows of contigs vs TnCentral.
    Does:     merge HSPs per reference element, apply the identity and coverage
              thresholds, discard the entries that cannot carry a passenger gene
              (plain IS), de-duplicate nested matches.
    Returns:  (element_rows, audit_rows).

    A candidate must survive all six tests below, and it leaves at the first one
    it fails, so the audit carries one reason per hit rather than a list:

      1. plain IS               -> not_applicable, ISEScan's territory
      2. unrecognised defline   -> discarded, we will not guess the element kind
      3. identity               -> too diverged to carry the name
      4. aligned fraction       -> the interval is padded with unrelated DNA
      5. span vs reference      -> separate copies were joined into one interval
      6. reference coverage     -> only a fragment of the element is here

    The order changes only which reason is recorded, never which hits survive.
    Tests 4 and 5 are worth reading first anyway: a padded or merged interval can
    show excellent identity and coverage, because the alignments in it are real —
    there are just not enough of them to fill the span they are reported across.
    """
    audit_rows = []

    # Unparsable lines are skipped rather than fatal, but never silently: a
    # truncated BLAST file would otherwise produce "no curated transposon found"
    # with complete confidence.
    if skipped_lines:
        audit_rows.append(audit_row(
            sample, "discarded", "blast_line_unparsable",
            f"{len(skipped_lines)} BLAST line(s) had fewer than "
            f"{len(BLAST_COLUMNS)} columns and were skipped (first at line "
            f"{skipped_lines[0]}). blastn does not write short lines, so this "
            "means the file was truncated or the -outfmt in the rule no longer "
            "matches BLAST_COLUMNS here - the naming result is incomplete."))

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

        # The interval must actually be the element. Terminal inverted repeats
        # are shared between related transposons, so short IR hits belonging to a
        # NEIGHBOURING element cluster with the real copy and drag the reported
        # interval across DNA that has nothing to do with this transposon. On the
        # control run that put 2,454 bp of an unrelated IS66 element inside a
        # "Tn3000" interval — and colocalise.py would call any AMR gene in there
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


# ── Writing the two TSVs and the command line ────────────────────────────────

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
    """Command-line entry point used by rule name_elements.

    Both TSVs are always written, even when BLAST found nothing, so the rule's
    declared outputs exist and amr_mge_colocalisation can join against an empty
    table instead of special-casing a missing file.
    """
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

    hits, skipped = read_blast_hits(args.blast)
    elements, audit = build_elements(
        args.sample, hits, args.min_identity, args.min_coverage,
        skipped_lines=skipped)

    write_tsv(args.out_table, ELEMENT_COLUMNS, elements)
    write_tsv(args.out_audit, AUDIT_COLUMNS, audit)

    # Closing summary. Rule name_elements sends stdout to
    # logs/mobilome_name_elements_{sample}.log, so the counts and the path to the
    # reasons are in the log even when the element table comes out empty.
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
