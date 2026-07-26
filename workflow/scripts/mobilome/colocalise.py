#!/usr/bin/env python3
"""Place every AMR gene in its mobile-element context and give it a mobility tier.

This is the analytical heart of the mobilome module (spec `docs/mobilome_module_SPEC.md`
§7 = work package D, and §2.5 = the mobility ladder). It answers ONE question per
AMR gene:

    is this resistance gene sitting inside (or next to) a mobile genetic element,
    and if so how transferable is that element?

WHY THIS MATTERS (the biology)
------------------------------
A resistance gene that is part of the species' normal chromosomal repertoire
(an intrinsic efflux pump, a chromosomal AmpC) is a very different risk from the
same activity carried on a conjugative plasmid, because only the second one can
move into another bacterium. Regulators (EFSA) ask exactly this: is the
determinant INTRINSIC or ACQUIRED, and is it TRANSFERABLE. This script produces
the supporting evidence for that judgement. It does not produce the judgement
itself, and the output must never be labelled "EFSA-compliant".

THE MOBILITY LADDER implemented here (spec §2.5), lowest to highest:

    1  chromosomal, no mobile-element context      -> intrinsic candidate
    2  IS sitting upstream and pointing at the gene -> the IS can supply an
       outward-reading ("hybrid") promoter that raises expression. This is
       EXPRESSION MODULATION, NOT mobilisation - the gene still cannot move.
    3  gene sits between two copies of the same IS  -> a composite transposon;
       the whole block can hop within the cell.
    4  gene sits inside a named unit transposon or an integron cassette
       -> mobilisable, and with a curated architecture we can name.
    5  gene sits on a mobilisable plasmid, or inside an IME
       -> it can move to another cell, but only with a helper element.
    6  gene sits on a conjugative plasmid, or inside an ICE
       -> PREDICTED self-transmissible. Always write "predicted": the
          confirmatory experiment is a filter/broth mating assay, not software.

SPECIAL CASE, reported separately and never counted as mobilisation: an IS that
has landed INSIDE the AMR coding sequence. That usually INACTIVATES the gene
(the strain is likely susceptible despite the hit), so it gets its own
`is_inside_amr_cds` flag rather than being folded into the mobilisation call.

THE HONEST-SIGNALS RULE (spec §2.3 - read this before trusting any output)
--------------------------------------------------------------------------
On short-read (SPAdes) assemblies, IS elements are the main reason contigs
break: multiple identical IS copies collapse into one, and the assembler cannot
tell which copy went where. So the AMR gene and the IS that flanks it very often
land on DIFFERENT contigs - the exact structure we are trying to detect is what
destroyed the assembly. Consequences hard-wired here:

  * a located IS count is a FLOOR, never a count;
  * ABSENCE of a flanking IS may simply mean the contig ended, so every row
    carries `dist_to_contig_end` and `is_at_contig_boundary`, and anything
    resting on a contig boundary or on a partial hit is capped at LOW confidence;
  * published IS-detection false-discovery rates are 8-24% even on curated data,
    so confidence is always TIERED (high/medium/low) and never a bare call.

Full composite-transposon and ICE-boundary calling belongs to the long-read
workflows (BacFluxL); this script does the short-read-honest version and simply
consumes long-read element calls when they are supplied.

DATA FLOW
---------
Inputs (all plain TSV; parsed with the standard library only - no pandas, no
bedtools, the interval arithmetic is written out below so it can be read):

  --amrfinder       AMRFinderPlus 4.x TSV (rule amrfinderplus, WP-A). One row per
                    AMR/STRESS/VIRULENCE hit with contig + 1-based coordinates.
                    THIS IS THE ONLY REQUIRED INPUT - it defines the output rows.
  --is-table        mobile elements with coordinates, normally derived from
                    ISEScan (rule isescan + loaders). One row per element. An
                    optional `element_type` column lets the SAME file also carry
                    named transposons/integrons (from the TnCentral naming
                    cascade) and ICE/IME calls (BacFluxL), so tiers 3, 4 and 6
                    all read from one tidy table.
  --replicons       one row per contig: chromosome or plasmid (from Platon /
                    the plasmid-concordance step) plus, for plasmids, whether
                    CONJscan found conjugation machinery (conjugative) or only a
                    relaxase (mobilisable). Drives tiers 5 and 6.
  --contig-lengths  contig -> length, or a samtools .fai. Needed for the
                    distance-to-contig-end signal described above.

Every input EXCEPT --amrfinder may be missing, empty or header-only: ISEScan
writes no files at all when it finds no IS, and plenty of isolates have no
plasmid. Those cases degrade gracefully and are written into the audit file, so
"no evidence" is never silently confused with "no data".

Outputs:
  --out-table  one row per AMR gene, the deliverable (spec §9).
  --out-audit  the decision trail: an explicit reason for every gene that got no
               mobile-element context, every candidate structure that was
               rejected, and every confidence cap that was applied. This follows
               BacFlux's existing `contig_taxonomy_decisions.tsv` convention -
               every filtering decision must be auditable.

Consumed next by the mobilome report rule (spec §10).

Coordinates: everything here is 1-BASED AND INCLUSIVE on both ends, which is what
AMRFinderPlus and ISEScan both emit, so no conversion happens anywhere. Distances
are the number of bases strictly BETWEEN two features (adjacent features = 0).

Licensing note: this is an independent implementation. No code or expression was
taken from EBI's mobilome-annotation-pipeline (parts of which are CC BY-NC-SA and
incompatible with BacFlux's MIT licence). Only conventions are borrowed - the
Sequence-Ontology vocabulary, the "assign a CDS to an element at >=0.9 coverage"
threshold and the discard-with-reason file pattern - and conventions are facts,
not expression (spec §11).

Standalone CLI, stdlib only, exercised by test_colocalise.py with synthetic
tables and no tools or databases.
"""

import argparse
import csv
import os
import sys


# ---------------------------------------------------------------------------
# CONVENTIONS, NOT BIOLOGY
#
# Every number below is a threshold we CHOSE. None of them is a property of a
# cell. That is why the output table always reports the MEASURED distance next
# to the tier: a reader who disagrees with a threshold can re-judge a row
# without re-running anything.
# ---------------------------------------------------------------------------

# Longest block of DNA we are willing to call a composite transposon. Real
# composites run from ~2.5 kb (IS26 translocatable units) to ~25 kb; 20 kb is the
# spec's default and is overridable on the command line.
DEFAULT_MAX_COMPOSITE_SPAN_BP = 20000

# How close an upstream IS must be before we are willing to say it could be
# driving the gene from an outward-reading promoter (ladder tier 2). The classic
# examples sit very close: ISEcp1 is 42-266 bp upstream of blaCTX-M. 500 bp is a
# generous round number; beyond it, an intervening gene is likely and the
# mechanism becomes speculative.
HYBRID_PROMOTER_MAX_BP = 500

# How close any IS has to be before we bother reporting it as context at all.
# Beyond this the nearest IS is still reported in the distance column, but it no
# longer sets `mge_context`.
IS_ADJACENT_MAX_BP = 5000

# A feature this close to the end of its contig cannot be trusted to have its
# real neighbourhood assembled next to it - the contig ended, so absence of
# evidence is not evidence of absence. Triggers the boundary flag and the
# low-confidence cap.
CONTIG_END_WINDOW_BP = 1000

# Fraction of the AMR gene that must lie inside a named element before we say the
# gene is carried BY that element. Borrowed from the EBI Mobilome Annotation
# Pipeline's published thresholds (a convention, see the licensing note above).
MIN_AMR_FRACTION_INSIDE_ELEMENT = 0.9

# Fraction of the AMR gene that must be overlapped by an IS before we call the
# coding sequence disrupted (likely inactivated) rather than merely overlapping.
# Same 0.5 the spec's bedtools recipe uses (`bedtools intersect -f 0.5`).
MIN_AMR_FRACTION_OVERLAPPED_FOR_DISRUPTION = 0.5

# How close an IS must sit to the edge of a PARTIAL AMR hit before we read the
# pair as "the IS landed in this gene and split it". A transposition leaves the
# insertion flush against the interrupted sequence, usually with a short
# target-site duplication of a few bp, so this window is deliberately tight - it
# is testing for "butted up against", not "nearby". See
# find_truncating_abutting_is for why this second test is needed at all.
ABUTTING_IS_MAX_GAP_BP = 25

# Below this % coverage of the reference, an AMRFinderPlus hit is treated as a
# FRAGMENT of the gene rather than the whole thing. Only used together with the
# abutting-IS test above; on its own a low coverage means nothing about mobility.
MIN_COVERAGE_FOR_INTACT_GENE = 90.0

# IS families whose copies flank a cargo gene in DIRECT orientation. IS26 (family
# IS6) is the single most clinically important AMR architecture there is: it
# builds "translocatable units" with its copies in direct, not inverted,
# orientation. The generic composite rule below asks for the two flanking copies
# to be in the same orientation; this family is exempted from that test entirely
# so that IS26 structures are called whichever way round the two copies sit.
# Without this exception the most important architecture in the table is the one
# we would miss (spec §7, "Pitfalls").
ORIENTATION_EXEMPT_IS_FAMILIES = {"IS6", "IS26", "IS6/IS26", "IS26/IS6"}

# ISEScan family values that carry no information about WHICH element this is.
# Two elements both labelled "new" are not necessarily the same family, so a pair
# like that must not be called a composite transposon on family grounds alone.
UNINFORMATIVE_FAMILY_VALUES = {"", "NA", "N/A", "-", ".", "NEW", "UNKNOWN", "ISNCY"}

# Plain-language meaning of each rung of the ladder, written into the table so a
# biologist reading the TSV never has to look the number up. Note tier 6's
# wording: "predicted", always (spec §2.6 language discipline).
MOBILITY_TIER_LABELS = {
    "NA": "not_assessable",
    1: "intrinsic_candidate",
    2: "expression_modulation_not_mobilisation",
    3: "composite_mobilisable_within_cell",
    4: "named_element_mobilisable",
    5: "mobilisable_needs_helper",
    6: "predicted_self_transmissible",
}

# The deliverable's columns, in the order they are written. Kept as one list so
# the header and every row always agree (spec §9).
OUTPUT_COLUMNS = [
    "sample",
    "contig",
    "replicon",                  # chromosome | plasmid | unknown
    "replicon_id",               # plasmid identifier when there is one, else NA
    "amr_gene",                  # AMRFinderPlus "Element symbol"
    "amr_name",                  # AMRFinderPlus "Element name"
    "amr_element_type",          # AMR | STRESS | VIRULENCE  (AMRFinderPlus "Type")
    "amr_class",
    "amr_subclass",
    "amr_start",
    "amr_end",
    "amr_strand",
    "amrfinder_method",          # EXACTX / BLASTP / HMM / PARTIAL_CONTIG_ENDX ...
    "amr_pct_identity",
    "amr_pct_coverage",
    "amr_partial_at_contig_end", # yes | no  (Method starts with PARTIAL_CONTIG_END)
    "mge_context",               # none|is_adjacent|composite|unit_transposon|integron|ice|ime|plasmid
    "mge_id",
    "mge_name",                  # curated name when a naming cascade supplied one
    "distance_bp",               # measured bp between gene and the context element
    "orientation",               # same | opposite | unknown | NA  (IS vs gene strand)
    "is_family",
    "is_cluster",
    "n_flanking_is",             # IS on either side within the composite window
    "same_orientation",          # yes | no | NA  (the two flanking copies)
    "flanking_span_bp",          # measured span of the called composite
    "is_inside_amr_cds",         # yes | no  - likely INACTIVATION, not mobilisation
    "is_amr_overlap_bp",
    "plasmid_mobility",          # conjugative | mobilisable | non_mobilisable | unknown | NA
    "mobility_evidence",         # free text passed through from CONJscan, if given
    "contig_length",
    "dist_to_contig_end",
    "is_at_contig_boundary",     # yes | no | NA
    "spans_contigs",             # yes | no
    "mobility_tier",             # 1-6, or NA when the gene has no coordinates
    "mobility_tier_label",
    "confidence",                # high | medium | low
]

# The audit trail's columns. One row per decision, never per gene: a gene can
# generate several (no context found, plus two rejected structures, plus a
# confidence cap).
AUDIT_COLUMNS = [
    "sample",
    "contig",
    "amr_gene",
    "amr_start",
    "amr_end",
    "mobility_tier",
    "mge_context",
    "decision",   # input_missing | row_skipped | no_mge_context | tier_not_raised |
                  # evidence_rejected | confidence_capped
    "reason",     # short machine-readable slug
    "detail",     # the numbers behind the reason, in plain language
]


# ===========================================================================
# Small shared helpers: reading tables, finding columns, interval arithmetic.
# ===========================================================================


def read_tsv(path):
    """Read a tab-separated file into (header_list, list_of_row_lists).

    Returns ([], []) for a path that is None, missing, empty, or contains only
    comments. That is the graceful-degradation path the module depends on:
    ISEScan writes NO output files at all when it finds no IS elements, and an
    isolate with no plasmid may have an empty replicon table. "No file" and "no
    rows" must both mean "no data", not "crash".

    Lines starting with '#' are skipped, because several MacSyFinder/CONJscan
    derived tables carry '#' comment banners above their header.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []

    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row and not row[0].startswith("#")]

    if not rows:
        return [], []
    header = [cell.strip() for cell in rows[0]]
    data_rows = rows[1:]
    return header, data_rows


def find_column(header, accepted_names):
    """Return the index of the first column whose name matches, else None.

    Columns are matched by NAME, not by position, so an upstream tool adding or
    reordering a column does not silently shift the values we read. Several
    alternative names are accepted per field because the upstream loader may
    hand us either its own tidy names (`contig`, `start`) or ISEScan's raw ones
    (`seqID`, `isBegin`); both describe the same thing.

    Matching ignores case and surrounding spaces.
    """
    normalised_header = [name.strip().lower() for name in header]
    for candidate in accepted_names:
        wanted = candidate.strip().lower()
        if wanted in normalised_header:
            return normalised_header.index(wanted)
    return None


def cell(row, index, default="NA"):
    """Read one cell from a row, tolerating a short (ragged) row.

    Real tool output occasionally has rows with fewer fields than the header
    (trailing empty columns get dropped). Returning the default instead of
    raising means one ragged line cannot take the whole sample down; the value
    simply reads as missing, which is visible in the output.
    """
    if index is None or index >= len(row):
        return default
    value = row[index].strip()
    return value if value else default


def to_int(value):
    """Parse a coordinate into an int, or return None if it is not a number.

    Used for start/stop/length fields. None means "unusable", and the caller
    writes an audit row rather than guessing a coordinate.
    """
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def to_float(value):
    """Parse a percentage into a float, or return None if it is not a number.

    Used for AMRFinderPlus's "% Coverage of reference", which is "NA" whenever the
    hit came from an HMM rather than an alignment. None means "not stated", and
    every caller treats that as "cannot judge", never as zero.
    """
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def overlap_bp(a_start, a_end, b_start, b_end):
    """How many bases two 1-based inclusive intervals share (0 if they do not).

    Worked example: 100-200 and 150-300 share 150..200 = 51 bases.
    """
    shared = min(a_end, b_end) - max(a_start, b_start) + 1
    return shared if shared > 0 else 0


def gap_bp(a_start, a_end, b_start, b_end):
    """How many bases lie strictly BETWEEN two 1-based inclusive intervals.

    Touching or overlapping intervals give 0. Worked example: a feature ending
    at 100 and the next starting at 151 have 50 bases between them (101..150).

    This is the number reported in `distance_bp`. Reporting the measured gap
    matters more than the tier it falls into: thresholds are convention.
    """
    space = max(a_start, b_start) - min(a_end, b_end) - 1
    return space if space > 0 else 0


# ===========================================================================
# Input parsers. One per file, each returning plain Python structures.
# ===========================================================================

# AMRFinderPlus 4.x renamed most of its columns (3.x "Gene symbol" became
# "Element symbol", and so on). A parser written against the old names finds
# nothing and silently reports zero AMR genes, so we detect the old header
# explicitly and stop with a message instead.
AMRFINDER_LEGACY_NAMES = {
    "Gene symbol": "Element symbol",
    "Sequence name": "Element name",
    "Element type": "Type",
    "Protein identifier": "Protein id",
}


def parse_amrfinder(path):
    """Read the AMRFinderPlus TSV into one dict per AMR hit.

    Input: the output of the `amrfinderplus` rule (WP-A), run with `-n` (and
    normally `-p`/`-g`), so every row carries the contig it was found on plus
    1-based inclusive Start/Stop on that contig. Column names are AMRFinderPlus
    4.2.7's, verified against a real run (see docs/mobilome_wpA_ground_truth.md).

    What it produces: a list of dicts, one per hit, which becomes one row of the
    deliverable. STRESS and VIRULENCE rows are kept alongside AMR rows rather
    than dropped - a mercury-resistance (STRESS) gene inside a transposon is
    exactly the classic Tn21 cargo, and the `amr_element_type` column lets a
    reader filter to AMR only if they want to.

    Missing values are AMRFinderPlus's literal string "NA", never empty.
    """
    header, data_rows = read_tsv(path)
    if not header:
        raise ValueError(
            f"AMRFinderPlus table '{path}' is missing or empty. This input is "
            "required: it defines the rows of the report."
        )

    # Guard against being handed an AMRFinderPlus 3.x file, whose different
    # column names would otherwise parse as "no AMR genes found".
    for old_name, new_name in AMRFINDER_LEGACY_NAMES.items():
        if old_name in header and new_name not in header:
            raise ValueError(
                f"'{path}' looks like AMRFinderPlus 3.x output (found column "
                f"'{old_name}'). This module expects 4.x, which renamed it to "
                f"'{new_name}'. Re-run the amrfinderplus rule."
            )

    contig_index = find_column(header, ["Contig id"])
    start_index = find_column(header, ["Start"])
    stop_index = find_column(header, ["Stop"])
    symbol_index = find_column(header, ["Element symbol"])
    method_index = find_column(header, ["Method"])
    missing = [
        name
        for name, index in [
            ("Contig id", contig_index),
            ("Start", start_index),
            ("Stop", stop_index),
            ("Element symbol", symbol_index),
            ("Method", method_index),
        ]
        if index is None
    ]
    if missing:
        raise ValueError(
            f"AMRFinderPlus table '{path}' is missing required column(s): "
            f"{', '.join(missing)}. Header was: {header}"
        )

    # Everything else is optional: it enriches the report but the co-localisation
    # works without it, so a hand-trimmed table still runs.
    strand_index = find_column(header, ["Strand"])
    name_index = find_column(header, ["Element name"])
    type_index = find_column(header, ["Type"])
    class_index = find_column(header, ["Class"])
    subclass_index = find_column(header, ["Subclass"])
    coverage_index = find_column(header, ["% Coverage of reference"])
    identity_index = find_column(header, ["% Identity to reference"])

    amr_hits = []
    for line_number, row in enumerate(data_rows, start=2):
        amr_hits.append({
            "line_number": line_number,
            "contig": cell(row, contig_index),
            "start_raw": cell(row, start_index),
            "end_raw": cell(row, stop_index),
            "strand": cell(row, strand_index),
            "gene": cell(row, symbol_index),
            "name": cell(row, name_index),
            "element_type": cell(row, type_index),
            "amr_class": cell(row, class_index),
            "subclass": cell(row, subclass_index),
            "method": cell(row, method_index),
            "pct_coverage": cell(row, coverage_index),
            "pct_identity": cell(row, identity_index),
        })
    return amr_hits


# How the mobile-element table's `element_type` values map onto the vocabulary
# used in `mge_context`. The keys are lowercase and cover both the tidy names a
# loader would write and the Sequence-Ontology style names from the spec's §9
# output schema.
ELEMENT_TYPE_SYNONYMS = {
    "": "insertion_sequence",
    "na": "insertion_sequence",
    "is": "insertion_sequence",
    "insertion_sequence": "insertion_sequence",
    "transposon": "unit_transposon",
    "unit_transposon": "unit_transposon",
    "tn": "unit_transposon",
    "integron": "integron",
    "conjugative_integron": "ice",   # the SO term for an ICE is conjugative_integron
    "ice": "ice",
    "ime": "ime",
    # The two classes conjscan_to_ice.py reports but deliberately does NOT let
    # raise the tier (spec §8 phase 4: a conjugative region with no integrase is
    # unbounded - "report, don't call ICE"; an island has no conjugation machinery
    # at all). They must still be RECOGNISED here. While they were missing from
    # this table they parsed as element_type None, which meant they were dropped
    # from every test AND logged as an "unrecognised element_type" - so an AMR gene
    # sitting inside a predicted chromosomal conjugative region came out as
    # "intrinsic candidate" at high confidence, with a contradictory audit line
    # beside it. Recognised-but-not-tier-raising is the correct handling: they set
    # the context and cap the confidence (see the tier-1 branch), without ever
    # claiming self-transmissibility the evidence does not support.
    "conjugative_region": "conjugative_region",
    "genomic_island": "genomic_island",
    "island": "genomic_island",
    "cime": "genomic_island",        # a CIME is an integrated, non-mobile island
}

# Element types that describe a real mobile-element neighbourhood but are NOT
# allowed to raise the mobility tier on their own. Used by the tier-1 branch to
# turn "nothing found" into an honest "something is here, but it does not prove
# mobility".
CONTEXT_ONLY_ELEMENT_TYPES = {"conjugative_region", "genomic_island"}


def parse_mobile_elements(path):
    """Read the mobile-element table into one dict per element.

    WHERE IT COMES FROM: normally the ISEScan run on this sample's contigs
    (rule isescan), reshaped by the module's loader into tidy columns. The same
    file may also carry named transposons/integrons from the TnCentral naming
    cascade and, on the long-read workflows, ICE/IME calls - the optional
    `element_type` column says which is which, and rows without it are treated
    as insertion sequences (the ISEScan-only case).

    Columns read (first matching name wins, case-insensitive):
      contig        contig | seqID | sequence | seq_id      (required)
      start         start  | isBegin                         (required)
      end           end    | isEnd                           (required)
      strand        strand                        - may be EMPTY: ISEScan leaves
                                                    it blank for an IS with no
                                                    predicted transposase.
      family        family                        - IS family, e.g. IS6, IS3
      cluster       cluster                       - ISEScan's finer grouping
      element_type  element_type | type_of_element - see ELEMENT_TYPE_SYNONYMS
      complete      is_complete | complete | completeness | type
                                                  - ISEScan's `type` column is
                                                    'c' (complete) or 'p'
                                                    (partial); partial elements
                                                    are kept (we tier them
                                                    ourselves) but they lower
                                                    confidence.
      id            is_id | id | element_id | mge_id        - synthesised if absent
      name          name | element_name | mge_name          - curated name, if any

    NOTE on ISEScan's raw `.tsv`: its `type` column holds c/p, which is why
    `type` is accepted as a completeness alias and NOT as an element-type alias.
    Element type has to be spelled `element_type` to avoid that collision.

    Produces: a list of element dicts. Feeds the disruption test, the adjacency
    test, the composite-pair search and the "named element contains the gene"
    test below.
    """
    header, data_rows = read_tsv(path)
    if not header:
        return []

    contig_index = find_column(header, ["contig", "seqID", "sequence", "seq_id", "seqid"])
    start_index = find_column(header, ["start", "isBegin", "isbegin"])
    end_index = find_column(header, ["end", "isEnd", "isend"])
    if contig_index is None or start_index is None or end_index is None:
        raise ValueError(
            f"Mobile-element table '{path}' needs contig/start/end columns "
            f"(any of contig|seqID, start|isBegin, end|isEnd). Header was: {header}"
        )

    strand_index = find_column(header, ["strand"])
    # "is_family" FIRST because that is what isescan_to_table.py actually writes.
    # This used to accept only "family", which matches nothing in our own IS table,
    # so is_family came out NA on every row and the IS26/IS6 same-family tests
    # below silently never fired. Matching is exact (see find_column), so a missing
    # alias is invisible rather than an error - hence listing the real name first.
    family_index = find_column(header, ["is_family", "family"])
    cluster_index = find_column(header, ["cluster"])
    type_index = find_column(header, ["element_type", "type_of_element", "mge_type"])
    complete_index = find_column(header, ["is_complete", "complete", "completeness", "type"])
    id_index = find_column(header, ["is_id", "id", "element_id", "mge_id"])
    name_index = find_column(header, ["name", "element_name", "mge_name"])

    # Quality columns that only the ICE/IME table carries (conjscan_to_ice.py).
    # That script already decides whether an element's machinery is intact, whether
    # it straddles a contig break, and what confidence it deserves - and those
    # judgements must survive into the deliverable, or the AMR row ends up more
    # confident than the element it rests on. They are absent from the IS table,
    # where find_column simply returns None and the fields stay empty.
    machinery_index = find_column(header, ["machinery_intact"])
    spans_contigs_index = find_column(header, ["spans_contigs"])
    element_confidence_index = find_column(header, ["confidence"])
    mobility_label_index = find_column(header, ["mobility"])
    degraded_reason_index = find_column(header, ["degraded_reason"])

    elements = []
    for line_number, row in enumerate(data_rows, start=2):
        contig = cell(row, contig_index, default="")
        start = to_int(cell(row, start_index, default=""))
        end = to_int(cell(row, end_index, default=""))
        if not contig or start is None or end is None:
            sys.stderr.write(
                f"[colocalise] WARNING: skipping mobile-element row {line_number} "
                f"in '{path}' - unusable contig/start/end: {row!r}\n"
            )
            continue
        if start > end:
            # Coordinates are always forward-contig space; a reversed pair means
            # the upstream file put them the wrong way round. Swap and carry on
            # rather than lose a real element.
            start, end = end, start

        raw_type = cell(row, type_index, default="").lower()
        element_type = ELEMENT_TYPE_SYNONYMS.get(raw_type)

        elements.append({
            "line_number": line_number,
            "contig": contig,
            "start": start,
            "end": end,
            "strand": cell(row, strand_index, default=""),
            "family": cell(row, family_index, default=""),
            "cluster": cell(row, cluster_index, default=""),
            "raw_type": raw_type,
            "element_type": element_type,          # None when unrecognised
            "complete": normalise_completeness(cell(row, complete_index, default="")),
            "id": cell(row, id_index, default="") or f"{contig}_MGE_{line_number - 1}",
            "name": cell(row, name_index, default=""),
            # ICE/IME-only quality judgements; empty string for IS rows.
            "machinery_intact": cell(row, machinery_index, default=""),
            "spans_contigs": cell(row, spans_contigs_index, default=""),
            "own_confidence": cell(row, element_confidence_index, default=""),
            "mobility_label": cell(row, mobility_label_index, default=""),
            "degraded_reason": cell(row, degraded_reason_index, default=""),
        })
    return elements


def element_says_false(value):
    """Is this element-table cell an explicit FALSE?

    conjscan_to_ice.py writes TRUE/FALSE; other loaders might write no/0. Anything
    else - including an empty cell from the IS table, which has no such column - is
    treated as "not stated", so a missing column can never be read as a failure.
    """
    return str(value).strip().lower() in {"false", "no", "0"}


def element_quality_caps(element):
    """Turn an ICE/IME element's OWN quality flags into confidence caps.

    Biology and why this matters: conjscan_to_ice.py already asks the hard
    questions about each candidate - is the conjugation machinery complete, or is
    this a decayed element that can no longer move? does it straddle a contig
    break, so its extent is a guess? Those answers were being computed and then
    thrown away here, so an AMR gene sitting in a half-broken, contig-spanning
    element was still reported as tier 6 at high confidence. The element's own
    verdict must travel with it.

    Returns a list of (cap_level, reason, detail) in the same shape assess_gene's
    local `caps` list uses.
    """
    caps = []

    # Spec section 8 phase 6: anything spanning contigs is capped at low, whatever
    # else is true, because its coordinates are not established by this assembly.
    if str(element.get("spans_contigs", "")).strip().lower() in {"true", "yes", "1"}:
        caps.append((
            "low", "mge_spans_contigs",
            f"the element setting this context ({element['id']}) is reported across "
            "more than one contig, so its extent is not established by this assembly",
        ))

    # A degraded element is the classic overcall: decayed ICEs are common, and a
    # relaxase or VirB4 covering too little of its profile means the machinery
    # probably no longer works.
    if element_says_false(element.get("machinery_intact", "")):
        reason_detail = element.get("degraded_reason", "") or "machinery incomplete"
        caps.append((
            "medium", "mge_machinery_not_intact",
            f"the element setting this context ({element['id']}) has incomplete "
            f"conjugation machinery ({reason_detail}), so it may no longer be able "
            "to move itself",
        ))

    # Belt and braces: never report the AMR gene more confidently than the element
    # the call rests on, whatever the reason that element was downgraded.
    own_confidence = str(element.get("own_confidence", "")).strip().lower()
    if own_confidence in {"low", "medium"}:
        caps.append((
            own_confidence, "mge_own_confidence",
            f"the element setting this context ({element['id']}) was itself called "
            f"at {own_confidence} confidence",
        ))

    return caps


def normalise_completeness(value):
    """Turn the many spellings of complete/partial into 'complete'|'partial'|'unknown'.

    ISEScan writes 'c' or 'p'; a loader might write True/False or the words. A
    partial (fragmentary) IS call is kept - the spec says explicitly not to run
    ISEScan with --removeShortIS, because we tier partials ourselves - but it
    lowers the confidence of anything built on it.
    """
    text = str(value).strip().lower()
    if text in {"c", "complete", "true", "yes", "1"}:
        return "complete"
    if text in {"p", "partial", "false", "no", "0"}:
        return "partial"
    return "unknown"


# How the replicon table's plasmid-mobility values map onto our vocabulary. The
# distinction is the whole difference between ladder tier 5 and tier 6:
# a mobilisable plasmid carries a relaxase and can be moved by SOMEONE ELSE's
# conjugation machinery; a conjugative plasmid carries that machinery itself.
PLASMID_MOBILITY_SYNONYMS = {
    "conjugative": "conjugative",
    "conj": "conjugative",
    "self_transmissible": "conjugative",
    "self-transmissible": "conjugative",
    "mobilisable": "mobilisable",
    "mobilizable": "mobilisable",
    "mob": "mobilisable",
    "non_mobilisable": "non_mobilisable",
    "non-mobilisable": "non_mobilisable",
    "non_mobilizable": "non_mobilisable",
    "non-mobilizable": "non_mobilisable",
    "nonmobilizable": "non_mobilisable",
    "none": "non_mobilisable",
}


def parse_replicons(path):
    """Read the per-contig replicon call into {contig: {...}}.

    WHERE IT COMES FROM: Platon (via the plasmid-concordance step) says whether
    each contig is chromosome or plasmid; CONJscan/MacSyFinder says whether a
    plasmid carries conjugation machinery. Both are folded into one small table
    upstream so this script joins on one file.

    Columns read (first matching name wins):
      contig            contig | contig_id | seq_name | ID       (required)
      replicon          replicon | replicon_type | call | platon_call
                        -> chromosome | plasmid | anything else = unknown
      replicon_id       replicon_id | plasmid_id | mge_id        (defaults to contig)
      plasmid_mobility  plasmid_mobility | mobility | mob_class
      mobility_evidence mobility_evidence | evidence | mob_evidence
                        -> free text passed straight through, e.g. "MOBF;T4SS_typeF"

    A contig that is absent from this table is `unknown`, NOT chromosome. That
    matters: calling something "chromosomal, intrinsic candidate" when we simply
    never looked would be the most misleading error this module could make, so
    an unknown replicon caps confidence at medium and is written to the audit.
    """
    header, data_rows = read_tsv(path)
    if not header:
        return {}

    contig_index = find_column(header, ["contig", "contig_id", "seq_name", "ID", "id"])
    if contig_index is None:
        raise ValueError(
            f"Replicon table '{path}' needs a contig column "
            f"(contig|contig_id|seq_name|ID). Header was: {header}"
        )
    call_index = find_column(header, ["replicon", "replicon_type", "call", "platon_call"])
    replicon_id_index = find_column(header, ["replicon_id", "plasmid_id", "mge_id"])
    mobility_index = find_column(header, ["plasmid_mobility", "mobility", "mob_class"])
    evidence_index = find_column(header, ["mobility_evidence", "evidence", "mob_evidence"])

    replicons = {}
    for line_number, row in enumerate(data_rows, start=2):
        contig = cell(row, contig_index, default="")
        if not contig:
            continue

        raw_call = cell(row, call_index, default="").lower()
        if raw_call.startswith("plasmid"):
            replicon = "plasmid"
        elif raw_call.startswith("chromosome") or raw_call.startswith("chr"):
            replicon = "chromosome"
        else:
            replicon = "unknown"

        raw_mobility = cell(row, mobility_index, default="").lower()
        mobility = PLASMID_MOBILITY_SYNONYMS.get(raw_mobility, "unknown")

        if contig in replicons:
            sys.stderr.write(
                f"[colocalise] WARNING: contig '{contig}' appears more than once "
                f"in '{path}'; the last row (line {line_number}) wins.\n"
            )

        replicons[contig] = {
            "replicon": replicon,
            "replicon_id": cell(row, replicon_id_index, default="") or contig,
            "plasmid_mobility": mobility,
            "mobility_evidence": cell(row, evidence_index, default="NA"),
        }
    return replicons


def parse_contig_lengths(path):
    """Read contig lengths into {contig: length}.

    WHERE IT COMES FROM: the assembly, `02.assembly/{sample}/contigs_final.fasta`.
    Either a two-column TSV with a header (contig, length) or, because it is the
    obvious thing to have lying around, a headerless samtools `.fai` index whose
    first two fields are name and length. Both are accepted.

    WHAT IT IS FOR: the distance-to-contig-end signal. On a short-read assembly
    that number is the difference between "there is no IS next to this gene" and
    "the contig ran out before we could see one".
    """
    header, data_rows = read_tsv(path)
    if not header:
        return {}

    contig_index = find_column(header, ["contig", "contig_id", "seq_name", "name", "ID", "id"])
    length_index = find_column(header, ["length", "len", "bp", "contig_length", "size"])

    rows_to_read = data_rows
    if contig_index is None or length_index is None:
        # No recognisable header: assume the samtools .fai layout (name, length,
        # ...) and treat the "header" line we already consumed as data.
        contig_index, length_index = 0, 1
        rows_to_read = [header] + data_rows

    lengths = {}
    for row in rows_to_read:
        contig = cell(row, contig_index, default="")
        length = to_int(cell(row, length_index, default=""))
        if contig and length is not None:
            lengths[contig] = length
    return lengths


# ===========================================================================
# Biology rules. Each of these answers one small question about one gene.
# ===========================================================================


def normalise_family(value):
    """Upper-case an IS family name, or return '' when it says nothing useful.

    ISEScan writes 'new' for an element it could not assign and 'ISNCY' for
    "IS not classified yet". Two such elements are not necessarily the same
    family, so we deliberately return '' and let the caller fall back to
    ISEScan's finer `cluster` value instead of pairing them up on a non-name.
    """
    text = str(value).strip().upper()
    if text in UNINFORMATIVE_FAMILY_VALUES:
        return ""
    return text


def is_orientation_exempt(element):
    """True when this IS family is allowed to flank cargo in ANY orientation.

    See ORIENTATION_EXEMPT_IS_FAMILIES above: IS26/IS6 builds translocatable
    units with its copies in direct orientation, so applying the generic
    same-orientation test to it would throw away the most clinically important
    AMR architecture there is. We also catch the case where the family column is
    blank but the element was NAMED IS26 by the naming cascade.
    """
    if normalise_family(element["family"]) in ORIENTATION_EXEMPT_IS_FAMILIES:
        return True
    return element["name"].strip().upper().startswith("IS26")


def families_match(left, right):
    """Are these two IS copies the same kind of element? -> (bool, basis).

    Two copies of the SAME element flanking a gene is what makes a composite
    transposon: the transposase of either copy can act on the outer ends and
    move the whole block. Two DIFFERENT IS flanking a gene is just two
    insertions that happen to be nearby.

    Primary test is the IS family. When the family is uninformative ('new',
    'ISNCY', blank) we fall back to ISEScan's `cluster`, which groups elements
    at a finer level - two members of the same cluster are the same element even
    if neither has a family name. If neither test can be applied we return
    False, and the caller writes the reason to the audit.
    """
    left_family = normalise_family(left["family"])
    right_family = normalise_family(right["family"])
    if left_family and right_family:
        return (left_family == right_family), "family"

    left_cluster = left["cluster"].strip().upper()
    right_cluster = right["cluster"].strip().upper()
    if left_cluster and right_cluster and left_cluster not in {"NA", "-"}:
        return (left_cluster == right_cluster), "cluster"

    return False, "unknown"


def relative_orientation(element_strand, gene_strand):
    """Compare an element's strand with the gene's -> 'same'|'opposite'|'unknown'.

    'unknown' is a real and common answer: ISEScan leaves the strand blank for
    an IS with no predicted transposase, and AMRFinderPlus writes 'NA' when it
    has no coordinates. Guessing would invent evidence.
    """
    element = str(element_strand).strip()
    gene = str(gene_strand).strip()
    if element not in {"+", "-"} or gene not in {"+", "-"}:
        return "unknown"
    return "same" if element == gene else "opposite"


def find_disrupting_is(gene, insertion_sequences):
    """Find an IS that has landed INSIDE the AMR coding sequence.

    Biology: when an IS transposes into a resistance gene it breaks the reading
    frame, and the strain is usually SUSCEPTIBLE despite AMRFinderPlus reporting
    a hit (the hit is the surviving fragment). That is the opposite of
    mobilisation, so this is reported in its own `is_inside_amr_cds` flag and
    never folded into the mobility tier.

    Rule: the IS must overlap at least half of the AMR gene
    (MIN_AMR_FRACTION_OVERLAPPED_FOR_DISRUPTION, the spec's `bedtools intersect
    -f 0.5`). The measured overlap in bp is reported either way, so a smaller
    overlap is still visible in the table.

    Returns (disrupting_element_or_None, overlap_bp, overlapping_element_or_None).

    The third value exists because an overlap BELOW the disruption threshold used
    to disappear completely: it is not a disruption call, and find_nearest_is
    skips anything that overlaps the gene at all, so an IS covering (say) 40% of
    the CDS was invisible to every context test and the gene was reported as
    "intrinsic candidate, nothing nearby" at high confidence while an IS sat on
    it. Handing the element back lets the caller flag that honestly.
    """
    gene_length = gene["end"] - gene["start"] + 1
    best_element = None
    best_overlap = 0

    for element in insertion_sequences:
        shared = overlap_bp(gene["start"], gene["end"], element["start"], element["end"])
        if shared > best_overlap:
            best_overlap = shared
            best_element = element

    if best_element is None:
        return None, 0, None
    if gene_length > 0 and (best_overlap / gene_length) >= MIN_AMR_FRACTION_OVERLAPPED_FOR_DISRUPTION:
        return best_element, best_overlap, best_element
    # Too small to call disruption: report the number AND the element, not the call.
    return None, best_overlap, best_element


def find_truncating_abutting_is(gene, insertion_sequences):
    """An IS butted up against a gene that AMRFinderPlus only found PART of.

    Biology, and why the overlap test above is not enough on its own. When an IS
    transposes INTO a resistance gene it splits the coding sequence, and
    AMRFinderPlus then reports only the surviving fragment - the part that still
    matches its reference. So the coordinates it gives END where the IS BEGINS.
    The overlap between the reported gene and the IS is therefore ~0 by
    construction, and the >=50% overlap test can essentially never fire for the
    very architecture it exists to catch. What that architecture looks like on
    paper is a PARTIAL hit sitting flush against an IS.

    This is a second, independent signal, added ALONGSIDE the overlap test rather
    than replacing it - the overlap test still catches an IS called across a
    complete gene, which this one would miss.

    Two guards keep it honest:
      * The hit must be partial for a reason OTHER than the contig running out.
        AMRFinderPlus says which: a method naming CONTIG_END means the assembly
        truncated the gene, not an IS, so those are excluded.
      * The IS must be within ABUTTING_IS_MAX_GAP_BP of the gene's edge. A
        transposition leaves the two flush (often with a short target-site
        duplication), so this window is small on purpose.

    Returns (element_or_None, gap_bp, which_end) - which_end is "start" or "end",
    naming the side of the gene the IS sits against.
    """
    method = str(gene.get("method", "")).upper()
    coverage = to_float(gene.get("pct_coverage", ""))

    # "Partial for a reason other than the contig ending." PARTIALX / PARTIALP etc.
    # mean AMRFinderPlus matched only part of its reference; *CONTIG_END* means the
    # assembly is to blame, which is already reported separately.
    method_says_partial = "PARTIAL" in method and "CONTIG_END" not in method
    coverage_says_partial = coverage is not None and coverage < MIN_COVERAGE_FOR_INTACT_GENE
    if not (method_says_partial or coverage_says_partial):
        return None, None, None

    best_element = None
    best_gap = None
    best_end = None
    for element in insertion_sequences:
        gap_before_gene = gene["start"] - element["end"] - 1   # IS ends just before the gene
        gap_after_gene = element["start"] - gene["end"] - 1    # IS starts just after the gene
        for gap, which_end in ((gap_before_gene, "start"), (gap_after_gene, "end")):
            if 0 <= gap <= ABUTTING_IS_MAX_GAP_BP:
                if best_gap is None or gap < best_gap:
                    best_gap = gap
                    best_element = element
                    best_end = which_end
    return best_element, best_gap, best_end


def is_upstream_of_gene(element, gene):
    """Is this element upstream of the gene, in the gene's own reading direction?

    A gene on the + strand is read left-to-right, so its upstream side is the
    LOWER coordinates; a gene on the - strand is read right-to-left, so its
    upstream side is the HIGHER coordinates. This matters because an outward-
    reading promoter in an IS can only drive a gene that lies downstream of it.

    Returns False when the gene's strand is unknown - we cannot say which side
    is upstream, so we do not guess.
    """
    strand = str(gene["strand"]).strip()
    if strand == "+":
        return element["end"] < gene["start"]
    if strand == "-":
        return element["start"] > gene["end"]
    return False


def find_upstream_promoter_is(gene, insertion_sequences):
    """Find the closest IS that could be driving the gene from an outward promoter.

    This is ladder tier 2. Three things must all be true:
      1. the IS lies upstream in the gene's reading direction (see above);
      2. the IS is on the SAME strand as the gene - the convention we use for
         "pointing at the gene", since ISEScan reports the transposase strand and
         an IS end reading in the gene's direction is the one that can form a
         hybrid promoter;
      3. it is close enough that no other gene is likely to sit in between
         (HYBRID_PROMOTER_MAX_BP).

    Returns (element, gap_bp) for the closest qualifying IS, or (None, None).
    Also returns, as a third value, the closest IS that satisfied 1 and 2 but was
    TOO FAR - so the caller can write an honest audit line instead of silently
    dropping a near miss.
    """
    within_window = []
    beyond_window = []

    for element in insertion_sequences:
        if not is_upstream_of_gene(element, gene):
            continue
        if relative_orientation(element["strand"], gene["strand"]) != "same":
            continue
        distance = gap_bp(gene["start"], gene["end"], element["start"], element["end"])
        if distance <= HYBRID_PROMOTER_MAX_BP:
            within_window.append((distance, element))
        else:
            beyond_window.append((distance, element))

    within_window.sort(key=lambda pair: pair[0])
    beyond_window.sort(key=lambda pair: pair[0])

    best = within_window[0] if within_window else (None, None)
    near_miss = beyond_window[0] if beyond_window else (None, None)
    return best[1], best[0], near_miss[1], near_miss[0]


def find_nearest_is(gene, insertion_sequences):
    """Closest IS on the same contig that does not overlap the gene -> (element, gap).

    Reported for every gene regardless of tier, because the measured distance to
    the nearest mobile element is useful context even when it changes nothing:
    "nearest IS is 40 kb away" and "nearest IS is 300 bp away" are very different
    statements about the same tier-1 call.
    """
    best_element = None
    best_distance = None
    for element in insertion_sequences:
        if overlap_bp(gene["start"], gene["end"], element["start"], element["end"]) > 0:
            continue
        distance = gap_bp(gene["start"], gene["end"], element["start"], element["end"])
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_element = element
    return best_element, best_distance


def find_nearest_element(gene, candidate_elements):
    """Closest element of ANY type -> (element, distance), 0 when it overlaps.

    The sibling of find_nearest_is above, for the large elements (ICE, IME,
    conjugative region, genomic island) rather than insertion sequences. Two
    deliberate differences:

      * It does NOT skip overlapping elements. An element that overlaps the gene
        without covering the 90% needed to say the gene is "inside" it is exactly
        the ambiguous case worth reporting, not hiding; it comes back at distance 0.
      * It is only ever used to add context and cap confidence, never to raise a
        tier, so a loose match here cannot inflate a mobility claim.

    Returns (None, None) when the list is empty.
    """
    best_element = None
    best_distance = None
    for element in candidate_elements:
        if overlap_bp(gene["start"], gene["end"], element["start"], element["end"]) > 0:
            distance = 0
        else:
            distance = gap_bp(gene["start"], gene["end"], element["start"], element["end"])
        if best_distance is None or distance < best_distance:
            best_distance = distance
            best_element = element
    return best_element, best_distance


def find_composite_pair(gene, insertion_sequences, max_span_bp):
    """Look for two copies of the same IS flanking the gene -> a composite transposon.

    This is ladder tier 3. The structure being tested for is:

        ---[ IS_x ]--------- AMR gene ---------[ IS_x ]---
            left                                  right
        |<------------- span (<= max_span_bp) ------------>|

    All of the following must hold (spec §7):
      * both IS on the SAME contig as the gene (guaranteed - the caller only
        passes IS from this contig);
      * the gene lies strictly BETWEEN them;
      * the two copies are the same element (same family, or same ISEScan
        cluster when the family is uninformative);
      * they are in the SAME orientation - EXCEPT for IS26/IS6, which is exempt
        (see ORIENTATION_EXEMPT_IS_FAMILIES: it flanks cargo in direct
        orientation and is the most important architecture in clinical AMR);
      * the whole span is no longer than max_span_bp.

    Where several pairs qualify, the TIGHTEST one (smallest span) is returned:
    the innermost pair is the most plausible boundary of the moving unit.

    Returns (best_pair_or_None, rejection_summary) where a pair is
    (left_element, right_element, span_bp) and rejection_summary is a dict of
    {reason_slug: [detail strings]} for the audit trail - a candidate structure
    we looked at and turned down is exactly what the audit file is for.

    LIMITATION, and it is the important one: on a short-read assembly the two IS
    copies usually collapse into one contig break, so a real composite very often
    CANNOT be seen. A negative here is weak evidence; a positive is strong.
    """
    left_candidates = [
        element for element in insertion_sequences
        if element["end"] < gene["start"]
        and gap_bp(gene["start"], gene["end"], element["start"], element["end"]) <= max_span_bp
    ]
    right_candidates = [
        element for element in insertion_sequences
        if element["start"] > gene["end"]
        and gap_bp(gene["start"], gene["end"], element["start"], element["end"]) <= max_span_bp
    ]

    accepted_pairs = []
    rejections = {}

    def reject(reason, detail):
        rejections.setdefault(reason, []).append(detail)

    for left in left_candidates:
        for right in right_candidates:
            span = right["end"] - left["start"] + 1
            label = (
                f"{left['id']} ({left['family'] or 'no family'}, {left['start']}-{left['end']}) / "
                f"{right['id']} ({right['family'] or 'no family'}, {right['start']}-{right['end']})"
            )

            if span > max_span_bp:
                reject("flanking_is_pair_span_too_long",
                       f"{label}: span {span} bp > {max_span_bp} bp limit")
                continue

            same_element, basis = families_match(left, right)
            if not same_element:
                if basis == "unknown":
                    reject("flanking_is_family_unknown",
                           f"{label}: neither IS family nor cluster is informative, "
                           "so the two copies cannot be shown to be the same element")
                else:
                    reject("flanking_is_different_family",
                           f"{label}: the two flanking IS are different elements (by {basis})")
                continue

            orientation_exempt = is_orientation_exempt(left) or is_orientation_exempt(right)
            left_strand = left["strand"].strip()
            right_strand = right["strand"].strip()
            strands_known = left_strand in {"+", "-"} and right_strand in {"+", "-"}
            same_orientation = strands_known and left_strand == right_strand

            if orientation_exempt:
                # IS26/IS6: accept either arrangement. Recorded so the reader can
                # see the exception was used rather than wonder why it passed.
                accepted_pairs.append((left, right, span, same_orientation, strands_known, True))
                continue

            if not strands_known:
                reject("flanking_is_strand_unknown",
                       f"{label}: at least one copy has no transposase strand, so the "
                       "same-orientation test cannot be applied")
                continue

            if not same_orientation:
                reject("flanking_is_pair_inverted_orientation",
                       f"{label}: the two copies are in inverted orientation; the "
                       "composite rule requires the same orientation for every family "
                       "except IS26/IS6")
                continue

            accepted_pairs.append((left, right, span, True, True, False))

    if not accepted_pairs:
        return None, rejections

    accepted_pairs.sort(key=lambda item: item[2])   # tightest span first
    return accepted_pairs[0], rejections


def find_containing_element(gene, elements, wanted_types):
    """Find a named element (transposon / integron / ICE / IME) carrying the gene.

    Ladder tiers 4 and 6. Unlike the composite search, this does not infer a
    structure from a pattern - the element interval was called by something else
    (the TnCentral naming cascade, or the long-read ICE module) and we only ask
    whether the AMR gene sits inside it.

    "Inside" means at least MIN_AMR_FRACTION_INSIDE_ELEMENT (0.9) of the gene
    lies within the element's interval, so a gene clipped by a few bases at an
    element boundary still counts.

    Returns the SMALLEST containing element of the wanted types (the tightest,
    most specific call), or None.
    """
    gene_length = gene["end"] - gene["start"] + 1
    if gene_length <= 0:
        return None

    containing = []
    for element in elements:
        if element["element_type"] not in wanted_types:
            continue
        shared = overlap_bp(gene["start"], gene["end"], element["start"], element["end"])
        if (shared / gene_length) >= MIN_AMR_FRACTION_INSIDE_ELEMENT:
            containing.append(element)

    if not containing:
        return None
    containing.sort(key=lambda element: element["end"] - element["start"])
    return containing[0]


CONFIDENCE_RANK = {"low": 1, "medium": 2, "high": 3}


def apply_confidence_caps(starting_level, caps):
    """Lower a confidence level to the strictest cap that applies.

    `caps` is a list of (level, reason, detail). Confidence can only ever go
    DOWN here - that is the point. Every cap that fired is also written to the
    audit file, so a "low" in the table always has a stated cause.
    """
    level = starting_level
    for cap_level, _reason, _detail in caps:
        if CONFIDENCE_RANK[cap_level] < CONFIDENCE_RANK[level]:
            level = cap_level
    return level


# ===========================================================================
# Per-gene assessment: pulls the rules above together into one output row.
# ===========================================================================


def assess_gene(sample, amr, elements_on_contig, replicons, contig_lengths,
                multi_contig_element_ids, max_span_bp, is_table_supplied,
                replicon_table_supplied):
    """Work out one AMR gene's mobile-element context, tier and confidence.

    Input:  one parsed AMRFinderPlus hit, plus every mobile element called on the
            SAME contig, plus the replicon call and contig length for that contig.
    Output: (row_dict keyed by OUTPUT_COLUMNS, list_of_audit_row_dicts).

    The steps below run in order and each one only answers its own question; the
    tier is decided at the end by a single top-down cascade, so there is exactly
    one place where the ladder is applied.
    """
    audit_rows = []
    caps = []          # (level, reason, detail) - confidence can only go down

    def add_audit(decision, reason, detail, tier="NA", context="NA"):
        audit_rows.append({
            "sample": sample,
            "contig": amr["contig"],
            "amr_gene": amr["gene"],
            "amr_start": amr["start_raw"],
            "amr_end": amr["end_raw"],
            "mobility_tier": tier,
            "mge_context": context,
            "decision": decision,
            "reason": reason,
            "detail": detail,
        })

    # Start from an all-unknown row and fill in what we can establish.
    row = {column: "NA" for column in OUTPUT_COLUMNS}
    row.update({
        "sample": sample,
        "contig": amr["contig"],
        "amr_gene": amr["gene"],
        "amr_name": amr["name"],
        "amr_element_type": amr["element_type"],
        "amr_class": amr["amr_class"],
        "amr_subclass": amr["subclass"],
        "amr_start": amr["start_raw"],
        "amr_end": amr["end_raw"],
        "amr_strand": amr["strand"],
        "amrfinder_method": amr["method"],
        "amr_pct_identity": amr["pct_identity"],
        "amr_pct_coverage": amr["pct_coverage"],
        "mge_context": "none",
        "is_inside_amr_cds": "no",
        "n_flanking_is": "0",
        "spans_contigs": "no",
    })

    # --- AMRFinderPlus's own fragmentation flag -----------------------------
    # Method values starting with PARTIAL_CONTIG_END mean the hit was cut off BY
    # THE CONTIG BOUNDARY - the gene continues onto DNA we never assembled. Note
    # the distinction from a bare PARTIAL*, which is truncated but internal
    # (a real pseudogene, not an assembly artefact). Only the former is the
    # fragmentation signal, and it caps confidence at low: we do not know what
    # the gene's real neighbourhood is.
    method = str(amr["method"]).strip().upper()
    partial_at_end = method.startswith("PARTIAL_CONTIG_END")
    row["amr_partial_at_contig_end"] = "yes" if partial_at_end else "no"
    if partial_at_end:
        caps.append((
            "low", "amr_hit_partial_at_contig_end",
            f"AMRFinderPlus Method='{amr['method']}': the gene runs off the end of "
            "the contig, so its true flanking context was not assembled",
        ))

    # --- Step 0: can this hit be placed on the assembly at all? -------------
    start = to_int(amr["start_raw"])
    end = to_int(amr["end_raw"])
    contig = amr["contig"]
    if contig in {"", "NA"} or start is None or end is None:
        # Happens when AMRFinderPlus was run protein-only (no -n): there are no
        # genome coordinates, so no co-localisation is possible. Say so rather
        # than default the gene to "chromosomal, no context", which would look
        # like a real tier-1 result.
        row["mobility_tier"] = "NA"
        row["mobility_tier_label"] = MOBILITY_TIER_LABELS["NA"]
        row["confidence"] = "low"
        add_audit("row_skipped", "no_contig_coordinates",
                  "the AMRFinderPlus row has no usable contig/Start/Stop, so this "
                  "gene cannot be placed against any mobile element")
        return row, audit_rows

    # The gene as the biology rules below see it. method and pct_coverage ride
    # along because find_truncating_abutting_is needs to know whether
    # AMRFinderPlus reported the WHOLE gene or only a surviving fragment - a
    # fragment flush against an IS is what an IS-inactivated gene looks like.
    gene = {
        "start": start,
        "end": end,
        "strand": amr["strand"],
        "method": amr["method"],
        "pct_coverage": amr["pct_coverage"],
    }

    # --- Step 1: where is the gene relative to the ends of its contig? ------
    contig_length = contig_lengths.get(contig)
    if contig_length is None:
        row["contig_length"] = "NA"
        row["dist_to_contig_end"] = "NA"
        row["is_at_contig_boundary"] = "NA"
        caps.append((
            "medium", "contig_length_unknown",
            f"no length known for contig '{contig}', so we cannot check whether the "
            "gene sits near a contig end (the main short-read caveat)",
        ))
    else:
        distance_to_end = min(start - 1, contig_length - end)
        distance_to_end = max(distance_to_end, 0)
        row["contig_length"] = str(contig_length)
        row["dist_to_contig_end"] = str(distance_to_end)
        at_boundary = distance_to_end < CONTIG_END_WINDOW_BP
        row["is_at_contig_boundary"] = "yes" if at_boundary else "no"
        if at_boundary:
            caps.append((
                "low", "amr_gene_near_contig_end",
                f"gene is {distance_to_end} bp from the end of contig '{contig}' "
                f"(< {CONTIG_END_WINDOW_BP} bp): a missing flanking element here may "
                "simply mean the contig ran out",
            ))

    # --- Step 2: which replicon is it on? -----------------------------------
    replicon_call = replicons.get(contig)
    if replicon_call is None:
        row["replicon"] = "unknown"
        row["replicon_id"] = "NA"
        row["plasmid_mobility"] = "NA"
        if not replicon_table_supplied:
            cap_detail = ("no replicon table was supplied (Platon not run, or no "
                          "plasmids found), so chromosome vs plasmid is unknown")
        else:
            cap_detail = f"contig '{contig}' is absent from the replicon table"
        caps.append(("medium", "replicon_call_unavailable", cap_detail))
    else:
        row["replicon"] = replicon_call["replicon"]
        row["replicon_id"] = replicon_call["replicon_id"]
        row["plasmid_mobility"] = (
            replicon_call["plasmid_mobility"] if replicon_call["replicon"] == "plasmid" else "NA"
        )
        row["mobility_evidence"] = replicon_call["mobility_evidence"]
        if replicon_call["replicon"] == "unknown":
            caps.append((
                "medium", "replicon_call_unavailable",
                f"contig '{contig}' is listed but not classified as chromosome or plasmid",
            ))

    # --- Step 3: split the elements on this contig by what they are ---------
    insertion_sequences = [
        element for element in elements_on_contig
        if element["element_type"] == "insertion_sequence"
    ]
    unrecognised = [
        element for element in elements_on_contig if element["element_type"] is None
    ]
    if unrecognised:
        # Do not let an element type we do not understand quietly drive a tier.
        add_audit("evidence_rejected", "unknown_element_type",
                  f"{len(unrecognised)} element(s) on this contig have an "
                  f"unrecognised element_type "
                  f"({sorted({element['raw_type'] for element in unrecognised})}) "
                  "and were ignored")

    # --- Step 4: has an IS landed inside the gene? (inactivation, not mobility)
    disrupting_is, disruption_overlap, overlapping_is = \
        find_disrupting_is(gene, insertion_sequences)
    # Second, independent signal for the same biology - see the long note on
    # find_truncating_abutting_is. Only consulted when the overlap test found
    # nothing, so the stronger evidence always wins.
    abutting_is, abutting_gap, abutting_end = \
        find_truncating_abutting_is(gene, insertion_sequences)
    row["is_amr_overlap_bp"] = str(disruption_overlap)
    if disrupting_is is not None:
        row["is_inside_amr_cds"] = "yes"
        add_audit(
            "evidence_rejected", "is_inside_amr_cds_likely_inactivation",
            f"IS {disrupting_is['id']} ({disrupting_is['family'] or 'no family'}) overlaps "
            f"{disruption_overlap} bp of this gene: the coding sequence is probably "
            "disrupted and the gene inactivated. Reported in is_inside_amr_cds and "
            "deliberately NOT counted as mobilisation.",
        )
    elif abutting_is is not None:
        # The IS-into-gene architecture as it actually appears in the data: a
        # PARTIAL AMRFinderPlus hit sitting flush against an IS. Treated exactly
        # like the overlap case - flagged as likely inactivation and deliberately
        # NOT counted as mobilisation (spec §2.5's special case).
        row["is_inside_amr_cds"] = "yes"
        add_audit(
            "evidence_rejected", "is_abuts_partial_amr_hit_likely_inactivation",
            f"IS {abutting_is['id']} sits {abutting_gap} bp from the {abutting_end} of "
            f"this hit, which AMRFinderPlus reported as partial "
            f"(method={gene.get('method', 'NA')}, coverage={gene.get('pct_coverage', 'NA')}%). "
            "That is what an IS landing INSIDE a resistance gene looks like: the "
            "reported gene is the surviving fragment and its coordinates stop where "
            "the IS starts, so the two barely overlap. Reported as likely "
            "inactivation and deliberately NOT counted as mobilisation.",
        )

    elif overlapping_is is not None:
        # An IS physically overlaps the gene, but by too little to call the CDS
        # disrupted. It is also excluded from every other test - find_nearest_is
        # skips overlapping elements, and the flanking/composite tests need the IS
        # to sit clear of the gene - so without this the row would read as
        # "nothing nearby" while an IS is sitting on the gene. Do not change the
        # tier (we cannot say what this means biologically), but do not claim
        # certainty either.
        caps.append((
            "medium", "is_partially_overlaps_amr_gene",
            f"IS {overlapping_is['id']} overlaps {disruption_overlap} bp of this gene "
            "- too little to call the coding sequence disrupted, but the gene's "
            "context is not clean",
        ))

    # --- Step 5: nearest IS, and the tier-2 upstream-promoter test ----------
    nearest_is, nearest_distance = find_nearest_is(gene, insertion_sequences)
    # An IS that has landed IN the gene must not then be read as an upstream
    # promoter for it. The two calls are mutually exclusive biology: one says the
    # reading frame is broken (the strain is probably susceptible), the other says
    # the gene is intact and being over-expressed. Without this exclusion the
    # abutting case above came out as tier 2 "expression modulation" - the
    # opposite conclusion, from the same coordinates. Spec §2.5 keeps inactivation
    # out of the mobilisation ladder entirely.
    inactivating_is_ids = {element["id"] for element in (disrupting_is, abutting_is)
                           if element is not None}
    promoter_candidates = [element for element in insertion_sequences
                           if element["id"] not in inactivating_is_ids]
    promoter_is, promoter_distance, near_miss_is, near_miss_distance = \
        find_upstream_promoter_is(gene, promoter_candidates)

    flanking_left = [element for element in insertion_sequences
                     if element["end"] < gene["start"]
                     and gap_bp(gene["start"], gene["end"], element["start"], element["end"]) <= max_span_bp]
    flanking_right = [element for element in insertion_sequences
                      if element["start"] > gene["end"]
                      and gap_bp(gene["start"], gene["end"], element["start"], element["end"]) <= max_span_bp]
    row["n_flanking_is"] = str(len(flanking_left) + len(flanking_right))

    # Only report a "just missed the promoter window" line when the IS is at
    # least in the neighbourhood. Without this bound, every gene on a contig
    # collects a line about an IS a megabase away, which buries the real near
    # misses in noise.
    if (near_miss_is is not None and promoter_is is None
            and near_miss_distance <= IS_ADJACENT_MAX_BP):
        add_audit(
            "tier_not_raised", "upstream_oriented_is_beyond_promoter_window",
            f"IS {near_miss_is['id']} is upstream and on the same strand as the gene "
            f"but {near_miss_distance} bp away (> {HYBRID_PROMOTER_MAX_BP} bp), too far "
            "to assume it supplies an outward-reading promoter",
        )

    # --- Step 6: composite transposon (tier 3) ------------------------------
    composite_pair, composite_rejections = find_composite_pair(
        gene, insertion_sequences, max_span_bp
    )
    for reason, details in sorted(composite_rejections.items()):
        # One summary line per reason rather than one per candidate pair, so the
        # audit stays readable on IS-rich contigs. The count plus the first two
        # examples is enough to see what was turned down and why.
        examples = "; ".join(details[:2])
        add_audit(
            "tier_not_raised", reason,
            f"{len(details)} candidate flanking pair(s) rejected. Examples: {examples}",
        )

    # --- Step 7: named elements that CONTAIN the gene (tiers 4 and 6) -------
    named_element = find_containing_element(gene, elements_on_contig,
                                            {"unit_transposon", "integron"})
    ice_element = find_containing_element(gene, elements_on_contig, {"ice"})
    ime_element = find_containing_element(gene, elements_on_contig, {"ime"})

    # Elements that describe a real neighbourhood but must not raise the tier
    # (conjugative region, genomic island). Looked up here so the tier-1 branch
    # can report them instead of claiming the gene has no context at all.
    containing_context_element = find_containing_element(
        gene, elements_on_contig, CONTEXT_ONLY_ELEMENT_TYPES)

    # The nearest NON-IS element of any kind on this contig, with its measured
    # distance. The nearest-IS search deliberately looks at insertion sequences
    # only, which left a gene sitting 16 kb from a called IME reported as
    # "intrinsic candidate, nothing nearby" - true of insertion sequences, and
    # badly misleading about the genome.
    nearest_context_element, nearest_context_distance = find_nearest_element(
        gene,
        [element for element in elements_on_contig
         if element["element_type"] in {"ice", "ime"} | CONTEXT_ONLY_ELEMENT_TYPES],
    )

    # --- Step 8: the ladder, applied once, from the top down ----------------
    # Highest applicable rung wins. The order below IS the ladder (spec §2.5);
    # a more specific, interval-level call (ICE) is preferred over a whole-
    # replicon one (conjugative plasmid) when both would give the same tier.
    replicon = row["replicon"]
    plasmid_mobility = row["plasmid_mobility"]
    base_confidence = "medium"
    # Normally the label is just MOBILITY_TIER_LABELS[tier]. A branch sets this
    # when the generic tier wording would misstate its own evidence (see the
    # non-mobilisable plasmid case below).
    tier_label_override = None
    # The element(s) that actually SET the context below. Each branch fills this
    # in, so the contig-end check in step 9 examines the right features and not
    # merely whichever element happened to be found first.
    context_elements = []

    if ice_element is not None:
        # Tier 6: the gene is inside an integrative and conjugative element - the
        # element carries its own conjugation machinery, so a chromosomal
        # location does NOT mean "not transferable".
        row["mobility_tier"] = 6
        row["mge_context"] = "ice"
        row["mge_id"] = ice_element["id"]
        row["mge_name"] = ice_element["name"] or "NA"
        row["distance_bp"] = "0"
        base_confidence = "high"
        context_elements = [ice_element]
        # The ICE step's own verdict on this element (degraded machinery, spans a
        # contig break, low own confidence) must cap the gene's confidence too -
        # otherwise the strongest claim the module can make, "predicted
        # self-transmissible", would be asserted at high confidence on top of an
        # element the module itself flagged as doubtful.
        caps.extend(element_quality_caps(ice_element))

    elif replicon == "plasmid" and plasmid_mobility == "conjugative":
        # Tier 6: a plasmid that encodes its own conjugation machinery.
        row["mobility_tier"] = 6
        row["mge_context"] = "plasmid"
        row["mge_id"] = row["replicon_id"]
        row["distance_bp"] = "0"
        # NOT "high", and the asymmetry this fixes is worth spelling out.
        #
        # This branch's evidence is Platon's "# Conjugation" column, which is a
        # count of HMM hits against conjugation proteins - not an assessment that
        # a working mating apparatus is present. A conjugative system is a
        # multi-protein machine (relaxase + coupling protein + a mating-pair
        # apparatus of a dozen or so genes), so one hit does not evidence one.
        #
        # On the CHROMOSOME, the very same evidence is treated far more carefully:
        # conjscan_to_ice.py refuses to call relaxase-plus-coupling-protein an ICE
        # and downgrades it to an IME marked "machinery incomplete". The plasmid
        # path applied no completeness test at all, so a single Platon hit reached
        # the top of the ladder at high confidence with no audit line - a stronger
        # claim, on weaker evidence, than the chromosome path would ever allow.
        #
        # The TIER stays 6: the plasmid genuinely carries conjugation evidence, and
        # "is this acquired and potentially transferable" is the regulatory
        # question. Only the certainty is capped. Raising it back to high needs the
        # mating-pair apparatus actually verified on the plasmid contig, which
        # would mean running a conjugation-machinery caller over plasmid contigs
        # (today CONJscan is only asked about chromosomal ICEs).
        base_confidence = "high"
        caps.append((
            "medium", "plasmid_conjugation_from_hit_counts_only",
            f"tier 6 here rests on Platon's conjugation gene count "
            f"({row['mobility_evidence']}), which counts HMM hits rather than "
            "checking that a complete mating apparatus is present; the mating-pair "
            "apparatus was not verified on this contig",
        ))

    elif ime_element is not None:
        # Tier 5: an integrative MOBILISABLE element - it has a relaxase but no
        # conjugation machinery of its own, so it needs a helper element.
        row["mobility_tier"] = 5
        row["mge_context"] = "ime"
        row["mge_id"] = ime_element["id"]
        row["mge_name"] = ime_element["name"] or "NA"
        row["distance_bp"] = "0"
        base_confidence = "high"
        context_elements = [ime_element]
        # Same reasoning as the ICE branch above: carry the element's own verdict.
        caps.extend(element_quality_caps(ime_element))

    elif replicon == "plasmid":
        # Tier 5: on a plasmid. Whether it can actually move depends on the
        # mobility typing, which is reported next to the tier - an untyped
        # plasmid is still an ACQUIRED element, which is the regulatory question,
        # but we must not imply it is transferable.
        row["mobility_tier"] = 5
        row["mge_context"] = "plasmid"
        row["mge_id"] = row["replicon_id"]
        row["distance_bp"] = "0"
        if plasmid_mobility == "mobilisable":
            base_confidence = "high"
        elif plasmid_mobility == "non_mobilisable":
            base_confidence = "high"
            # The tier number stays 5 - the gene IS on a plasmid, which is the
            # acquired-vs-intrinsic answer the regulator asks for - but tier 5's
            # generic label reads "mobilisable_needs_helper", and saying that about
            # a plasmid we just typed NON-mobilisable contradicts our own evidence
            # in the one column most people read. Override the wording; the audit
            # line below explains it.
            tier_label_override = "on_plasmid_typed_non_mobilisable"
            add_audit(
                "tier_not_raised", "plasmid_typed_non_mobilisable",
                "the gene is on a plasmid (tier 5 = plasmid context) but the plasmid "
                "was typed non-mobilisable: it is ACQUIRED but not predicted to "
                "transfer by conjugation",
                tier=5, context="plasmid",
            )
        else:
            base_confidence = "medium"
            add_audit(
                "tier_not_raised", "plasmid_mobility_untyped",
                "the gene is on a plasmid but no conjugative/mobilisable typing was "
                "supplied (CONJscan not run, or no machinery found), so tier 6 cannot "
                "be assessed",
                tier=5, context="plasmid",
            )

    elif named_element is not None:
        # Tier 4: inside a curated unit transposon or integron cassette. A curated
        # named hit OVERRIDES the pattern-based composite call below, because the
        # curated entry carries the true architecture (and gets IS26 right for
        # free) - spec §7, "a TnCentral named hit should override".
        row["mobility_tier"] = 4
        row["mge_context"] = named_element["element_type"]
        row["mge_id"] = named_element["id"]
        row["mge_name"] = named_element["name"] or "NA"
        row["distance_bp"] = "0"
        base_confidence = "high"
        context_elements = [named_element]
        if composite_pair is not None:
            add_audit(
                "tier_not_raised", "named_element_overrides_composite_pattern",
                f"the gene also sits between two {composite_pair[0]['family'] or 'unnamed'} "
                f"copies (span {composite_pair[2]} bp), but the curated element "
                f"{named_element['id']} contains it and takes precedence",
                tier=4, context=named_element["element_type"],
            )

    elif composite_pair is not None:
        # Tier 3: flanked by two copies of the same IS -> a composite transposon,
        # able to hop within the cell.
        left, right, span, same_orientation, strands_known, used_exemption = composite_pair
        row["mobility_tier"] = 3
        row["mge_context"] = "composite"
        row["mge_id"] = f"{left['id']}+{right['id']}"
        row["is_family"] = left["family"] or right["family"] or "NA"
        row["is_cluster"] = left["cluster"] or "NA"
        row["flanking_span_bp"] = str(span)
        row["same_orientation"] = (
            "yes" if same_orientation else ("no" if strands_known else "NA")
        )
        row["orientation"] = relative_orientation(left["strand"], gene["strand"])
        row["distance_bp"] = str(min(
            gap_bp(gene["start"], gene["end"], left["start"], left["end"]),
            gap_bp(gene["start"], gene["end"], right["start"], right["end"]),
        ))
        base_confidence = "high"
        context_elements = [left, right]
        if used_exemption:
            add_audit(
                "evidence_rejected", "is26_orientation_exemption_applied",
                f"flanking pair {left['id']}/{right['id']} is family "
                f"{left['family'] or right['family']}: the same-orientation test was "
                "skipped because IS26/IS6 forms translocatable units with its copies "
                "in direct orientation",
                tier=3, context="composite",
            )
        if left["complete"] != "complete" or right["complete"] != "complete":
            caps.append((
                "medium", "flanking_is_partial",
                f"at least one flanking copy is a partial IS call "
                f"({left['id']}={left['complete']}, {right['id']}={right['complete']}), "
                "so the composite boundaries are provisional",
            ))

    elif promoter_is is not None:
        # Tier 2: an IS upstream and pointing at the gene. This raises EXPRESSION
        # (an outward-reading hybrid promoter); it does NOT make the gene mobile.
        row["mobility_tier"] = 2
        row["mge_context"] = "is_adjacent"
        row["mge_id"] = promoter_is["id"]
        row["mge_name"] = promoter_is["name"] or "NA"
        row["is_family"] = promoter_is["family"] or "NA"
        row["is_cluster"] = promoter_is["cluster"] or "NA"
        row["distance_bp"] = str(promoter_distance)
        row["orientation"] = relative_orientation(promoter_is["strand"], gene["strand"])
        context_elements = [promoter_is]
        # Deliberately capped: the mechanism is inferred from coordinates and
        # strand alone. We have not seen a transcript.
        base_confidence = "medium"
        caps.append((
            "medium", "promoter_inferred_from_position_only",
            f"IS {promoter_is['id']} is {promoter_distance} bp upstream on the same "
            "strand; an outward-reading promoter is plausible but was not measured",
        ))

    else:
        # Tier 1: nothing found. On the chromosome with no mobile-element context,
        # this is the intrinsic-resistance candidate - the "not acquired, not
        # transferable" end of the ladder.
        row["mobility_tier"] = 1
        base_confidence = "high" if replicon == "chromosome" else "medium"

        if nearest_is is not None and nearest_distance is not None \
                and nearest_distance <= IS_ADJACENT_MAX_BP:
            # There IS an IS nearby, it just does not meet the tier-2 test (wrong
            # side of the gene, or wrong orientation). Report the context and the
            # measured distance; do not inflate the tier.
            row["mge_context"] = "is_adjacent"
            row["mge_id"] = nearest_is["id"]
            row["is_family"] = nearest_is["family"] or "NA"
            row["is_cluster"] = nearest_is["cluster"] or "NA"
            row["distance_bp"] = str(nearest_distance)
            row["orientation"] = relative_orientation(nearest_is["strand"], gene["strand"])
            base_confidence = "medium"
            context_elements = [nearest_is]
            add_audit(
                "tier_not_raised", "nearby_is_not_upstream_or_not_oriented",
                f"nearest IS {nearest_is['id']} is {nearest_distance} bp away "
                f"(orientation vs gene: "
                f"{relative_orientation(nearest_is['strand'], gene['strand'])}) but is not "
                "both upstream of the gene and on its strand, so it is reported as "
                "context only, not as expression modulation",
                tier=1, context="is_adjacent",
            )
        elif containing_context_element is not None:
            # No IS nearby, but the gene sits INSIDE a conjugative region or a
            # genomic island. These deliberately do not raise the tier - a
            # conjugative region with no integrase has no established boundaries,
            # and an island carries no conjugation machinery - but calling such a
            # gene an "intrinsic candidate" at high confidence is plainly wrong:
            # something mobile-element-shaped is sitting on top of it. Report the
            # context, cap the confidence, and say why.
            row["mge_context"] = containing_context_element["element_type"]
            row["mge_id"] = containing_context_element["id"]
            row["mge_name"] = containing_context_element["name"] or "NA"
            row["distance_bp"] = "0"
            context_elements = [containing_context_element]
            base_confidence = "medium"
            caps.extend(element_quality_caps(containing_context_element))
            add_audit(
                "tier_not_raised", "inside_context_only_element",
                f"the gene lies inside {containing_context_element['id']}, a "
                f"{containing_context_element['element_type'].replace('_', ' ')}. That is "
                "not enough to call it mobilisable - a conjugative region without an "
                "integrase has no established boundaries, and an island carries no "
                "conjugation machinery of its own - but it is not an intrinsic "
                "chromosomal gene either",
                tier=1, context=row["mge_context"],
            )

        elif nearest_context_element is not None and nearest_context_distance is not None \
                and nearest_context_distance <= IS_ADJACENT_MAX_BP:
            # An ICE/IME/region is CLOSE but does not contain the gene. Worth
            # reporting with its measured distance: the intervals these elements
            # carry are machinery spans, not resolved element boundaries (the ICE
            # step writes boundary_method=none when it could not find att sites),
            # so "just outside the span" does not mean "outside the element".
            row["mge_context"] = "near_" + nearest_context_element["element_type"]
            row["mge_id"] = nearest_context_element["id"]
            row["mge_name"] = nearest_context_element["name"] or "NA"
            row["distance_bp"] = str(nearest_context_distance)
            base_confidence = "medium"
            add_audit(
                "tier_not_raised", "near_but_outside_element_machinery_span",
                f"{nearest_context_element['id']} "
                f"({nearest_context_element['element_type']}) lies "
                f"{nearest_context_distance} bp away. The interval recorded for such an "
                "element is the span of its machinery genes, not a resolved boundary, "
                "so a gene just outside it may still be carried by the element - the "
                "tier is not raised, but this is not a clean intrinsic call either",
                tier=1, context=row["mge_context"],
            )

        else:
            # Genuinely no context. Every such gene gets an explicit reason here -
            # this is the "no-call audit" the BacFlux convention requires.
            if not is_table_supplied:
                add_audit(
                    "no_mge_context", "no_is_calls_available",
                    "no mobile-element table was supplied for this sample (ISEScan "
                    "writes no files when it finds no IS), so absence of context is "
                    "absence of DATA, not evidence of a stable chromosomal gene",
                    tier=1, context="none",
                )
            elif nearest_is is None:
                # Say what was actually checked. This branch is reached from
                # nearest_is, which looks at INSERTION SEQUENCES only, so claiming
                # "no mobile element" was a wider statement than the evidence -
                # a reader reconciling this line against {sample}_ice_candidates.tsv
                # would find the two files contradicting each other.
                add_audit(
                    "no_mge_context", "no_insertion_sequence_on_this_contig",
                    f"no insertion sequence was called anywhere on contig '{contig}'"
                    + (
                        f"; the nearest non-IS element on this contig is "
                        f"{nearest_context_element['id']} "
                        f"({nearest_context_element['element_type']}), "
                        f"{nearest_context_distance} bp away"
                        if nearest_context_element is not None else
                        " and no ICE/IME/region was called on it either"
                    ),
                    tier=1, context="none",
                )
            else:
                add_audit(
                    "no_mge_context", "nearest_mobile_element_too_far",
                    f"nearest IS {nearest_is['id']} is {nearest_distance} bp away, beyond "
                    f"the {IS_ADJACENT_MAX_BP} bp context window",
                    tier=1, context="none",
                )

    # A "no context found" call is only as good as the IS calls behind it. If
    # ISEScan produced nothing at all we cannot claim high confidence in absence.
    if row["mge_context"] == "none" and not is_table_supplied:
        caps.append((
            "medium", "no_is_calls_available",
            "no mobile-element calls were available, so 'no MGE context' means "
            "'not looked at', not 'looked at and found nothing'",
        ))

    # --- Step 9: does the winning element itself sit across a contig break? --
    # Anything spanning contigs is capped at low regardless of everything else
    # (spec §8 phase 6): the coordinates of such an element are a guess.
    # "NA" is excluded explicitly: it is our own placeholder for "no element", and
    # a loader that wrote the literal string NA into the id column would otherwise
    # make every unnamed element look like one element spread over the assembly.
    # Two independent ways an element can span contigs, and BOTH are needed:
    #
    #  1. We see the same element id on two different contig rows. This is how a
    #     split IS shows up.
    #  2. The element's own table says so. An ICE/IME row already carries a
    #     spans_contigs flag that conjscan_to_ice.py worked out from where its
    #     anchor genes landed - and detector 1 is structurally BLIND to it,
    #     because an ICE id embeds its contig ("contig_1|ice-100:200") and each
    #     element is written as exactly one row. So before this, spans_contigs
    #     could never read "yes" for an ICE no matter what the ICE step found.
    spanning_ids = {element["id"] for element in context_elements
                    if element["id"] != "NA" and element["id"] in multi_contig_element_ids}
    if row["mge_id"] != "NA" and row["mge_id"] in multi_contig_element_ids:
        spanning_ids.add(row["mge_id"])
    for context_element in context_elements:
        if str(context_element.get("spans_contigs", "")).strip().lower() in {"true", "yes", "1"}:
            spanning_ids.add(context_element["id"])
    if spanning_ids:
        row["spans_contigs"] = "yes"
        caps.append((
            "low", "mge_spans_contigs",
            f"element(s) {sorted(spanning_ids)} are reported on more than one contig, "
            "so their extent is not established by this assembly",
        ))

    # An element whose own edge sits at a contig end was cut off by the assembly:
    # its far boundary is unknown, so any structure resting on it is provisional.
    # For a composite call this covers BOTH flanking copies.
    if contig_length is not None:
        for context_element in context_elements:
            element_end_distance = min(context_element["start"] - 1,
                                       contig_length - context_element["end"])
            if element_end_distance < CONTIG_END_WINDOW_BP:
                caps.append((
                    "low", "context_element_at_contig_end",
                    f"the element setting this context ({context_element['id']}) lies "
                    f"{max(element_end_distance, 0)} bp from a contig end and is probably "
                    "truncated by the assembly",
                ))

    # --- Step 10: final confidence, and an audit line per cap that fired ----
    row["mobility_tier_label"] = tier_label_override or MOBILITY_TIER_LABELS[row["mobility_tier"]]
    row["confidence"] = apply_confidence_caps(base_confidence, caps)
    for cap_level, reason, detail in caps:
        add_audit("confidence_capped", reason, f"capped at {cap_level}: {detail}",
                  tier=row["mobility_tier"], context=row["mge_context"])

    return row, audit_rows


# ===========================================================================
# Whole-sample assembly and output
# ===========================================================================


def build_report(sample, amr_hits, elements, replicons, contig_lengths, max_span_bp,
                 is_table_supplied, replicon_table_supplied):
    """Run assess_gene over every AMR hit and collect the two output tables.

    Input:  the four parsed inputs for ONE sample.
    Output: (report_rows, audit_rows) ready to be written.
    """
    # Index the elements by contig once, so each gene only looks at its own
    # contig. This is also what makes cross-contig pairing impossible: on a
    # short-read assembly we must never claim two IS flank a gene when they are
    # on different contigs.
    elements_by_contig = {}
    contigs_per_element_id = {}
    for element in elements:
        elements_by_contig.setdefault(element["contig"], []).append(element)
        contigs_per_element_id.setdefault(element["id"], set()).add(element["contig"])

    multi_contig_element_ids = {
        element_id for element_id, contigs in contigs_per_element_id.items()
        if len(contigs) > 1
    }

    report_rows = []
    audit_rows = []
    for amr in amr_hits:
        elements_on_contig = elements_by_contig.get(amr["contig"], [])
        row, gene_audit = assess_gene(
            sample, amr, elements_on_contig, replicons, contig_lengths,
            multi_contig_element_ids, max_span_bp, is_table_supplied,
            replicon_table_supplied,
        )
        report_rows.append(row)
        audit_rows.extend(gene_audit)

    return report_rows, audit_rows


def write_tsv(path, columns, rows):
    """Write rows as a tab-separated table with a header.

    The header is written even when there are no rows: an empty table with a
    header says "we looked and found nothing", a missing file says nothing at
    all and breaks the downstream report rule.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv=None):
    # Arguments are supplied by the amr_mge_colocalisation rule in
    # workflow/rules/mobilome.smk (spec §10).
    parser = argparse.ArgumentParser(
        description="Assign each AMR gene its mobile-element context and mobility tier."
    )
    parser.add_argument("--sample", required=True,
                        help="Sample name, written into the first column of both outputs.")
    parser.add_argument("--amrfinder", required=True,
                        help="AMRFinderPlus 4.x TSV. REQUIRED: it defines the report rows.")
    parser.add_argument("--is-table", action="append", default=None,
                        help="Mobile-element table with coordinates. REPEATABLE: pass it "
                             "once per source and the rows are pooled, because the module "
                             "produces elements from two independent places — ISEScan "
                             "(insertion sequences) and CONJscan (ICE/IME candidates). "
                             "Each file is read on its own terms (they have different "
                             "columns), so no pre-merge step is needed. May be absent or "
                             "empty; every file that is missing or has no usable rows is "
                             "recorded in the audit.")
    parser.add_argument("--replicons", default=None,
                        help="Per-contig chromosome/plasmid call plus plasmid mobility "
                             "typing. May be absent or empty.")
    parser.add_argument("--contig-lengths", default=None,
                        help="contig<TAB>length table, or a samtools .fai. May be absent.")
    parser.add_argument("--max-composite-span", type=int,
                        default=DEFAULT_MAX_COMPOSITE_SPAN_BP,
                        help="Longest span (bp) accepted for a composite transposon "
                             f"(default {DEFAULT_MAX_COMPOSITE_SPAN_BP}). A convention, "
                             "not biology - the measured span is always reported.")
    parser.add_argument("--out-table", required=True,
                        help="Destination for the per-AMR-gene mobility table.")
    parser.add_argument("--out-audit", required=True,
                        help="Destination for the decision/no-call audit TSV.")
    args = parser.parse_args(argv)

    # Load each input. Only the AMRFinderPlus table is mandatory; the others are
    # allowed to be missing so the module still runs on an isolate with no IS
    # calls and no plasmids.
    amr_hits = parse_amrfinder(args.amrfinder)
    # --is-table is repeatable (ISEScan's insertion sequences and CONJscan's
    # ICE/IME candidates arrive as separate files with different columns), so
    # parse each one on its own and pool the rows.
    element_table_paths = args.is_table or []
    elements = []
    for element_table_path in element_table_paths:
        elements.extend(parse_mobile_elements(element_table_path))
    replicons = parse_replicons(args.replicons)
    contig_lengths = parse_contig_lengths(args.contig_lengths)

    # Per-SOURCE, not pooled. --is-table is given twice (ISEScan's IS table and
    # CONJscan's ICE/IME table), and this flag drives the "absence of context is
    # absence of DATA" audit line and its confidence cap. Asking merely whether ANY
    # row was parsed let one unrelated ICE row stand in as proof that IS detection
    # had run: a genome where ISEScan legitimately found nothing then reported its
    # AMR genes as "intrinsic candidate" at HIGH confidence, while the very same
    # genome with both tables empty correctly reported medium. Same IS evidence
    # (none), opposite claim. Count only what actually answers the question.
    is_table_supplied = any(element["element_type"] == "insertion_sequence"
                            for element in elements)
    replicon_table_supplied = bool(replicons)

    # Sample-level audit lines for whatever was not supplied, so a reader of the
    # audit file can tell "no evidence" from "no data" without checking the run log.
    startup_audit = []

    def note_missing(reason, detail):
        startup_audit.append({
            "sample": args.sample, "contig": "NA", "amr_gene": "NA",
            "amr_start": "NA", "amr_end": "NA", "mobility_tier": "NA",
            "mge_context": "NA", "decision": "input_missing",
            "reason": reason, "detail": detail,
        })

    if not elements:
        note_missing(
            "is_table_absent_or_empty",
            f"no usable mobile-element rows from '{args.is_table}'. ISEScan writes no "
            "output files at all when it finds no IS, so this is an expected state; "
            "every gene therefore reports mge_context=none with confidence capped at "
            "medium.",
        )
    if not replicons:
        note_missing(
            "replicon_table_absent_or_empty",
            f"no usable replicon rows from '{args.replicons}'. Chromosome vs plasmid is "
            "unknown for every contig, so tiers 5 and 6 cannot be reached and confidence "
            "is capped at medium.",
        )
    if not contig_lengths:
        note_missing(
            "contig_lengths_absent_or_empty",
            f"no contig lengths from '{args.contig_lengths}'. The distance-to-contig-end "
            "signal is unavailable, so we cannot tell whether a missing flanking element "
            "means 'absent' or 'the contig ended'.",
        )

    report_rows, audit_rows = build_report(
        args.sample, amr_hits, elements, replicons, contig_lengths,
        args.max_composite_span, is_table_supplied, replicon_table_supplied,
    )

    write_tsv(args.out_table, OUTPUT_COLUMNS, report_rows)
    write_tsv(args.out_audit, AUDIT_COLUMNS, startup_audit + audit_rows)

    # A short, honest summary for the run log. Never used for filtering, and
    # deliberately worded as "predicted" for tier 6 (spec §2.6).
    tier_counts = {}
    for row in report_rows:
        tier_counts[row["mobility_tier"]] = tier_counts.get(row["mobility_tier"], 0) + 1
    tier_summary = ", ".join(
        f"tier {tier}: {count}" for tier, count in sorted(tier_counts.items(), key=str)
    ) or "no AMR genes"
    predicted_transmissible = tier_counts.get(6, 0)
    disrupted = sum(1 for row in report_rows if row["is_inside_amr_cds"] == "yes")
    print(
        f"Sample {args.sample}: {len(report_rows)} AMR gene(s) placed ({tier_summary}); "
        f"{predicted_transmissible} predicted self-transmissible; "
        f"{disrupted} with an IS inside the coding sequence (likely inactivated); "
        f"{len(audit_rows) + len(startup_audit)} audit line(s)."
    )


if __name__ == "__main__":
    main()
