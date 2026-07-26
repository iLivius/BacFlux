#!/usr/bin/env python3
"""Turn ISEScan's raw output for ONE sample into the single tidy insertion-sequence
(IS) table the rest of the mobilome module consumes, and compute the honest
short-read QC signals that go with it (WP-C in docs/mobilome_module_SPEC.md §6).

WHY THIS EXISTS — the biology
    Insertion sequences are small mobile elements that copy themselves around a
    genome. They matter here because an AMR gene sitting next to, or between, IS
    copies may be mobilisable, while the same gene with no mobile context is more
    likely intrinsic. So the mobilome module needs one clean list of "where are
    the IS elements on this assembly", which is what this script produces.

    The catch, and the reason for half the code below: IS elements are the single
    biggest cause of contig breaks in short-read assemblies. A repeat that occurs
    in several identical copies cannot be resolved by the assembler, so the
    assembly is cut at exactly those copies. Two consequences we must report
    rather than hide:
      * the number of IS we can LOCATE is a FLOOR, never a count — collapsed
        copies are simply absent from the assembly;
      * an IS found at the very end of a contig is, quite literally, the place
        where the assembly fell apart, so its genomic context is unknown.
    Hence the per-IS distance-to-contig-end and at_contig_boundary flag, and the
    per-sample fraction of IS calls sitting at a contig end. Published IS-calling
    false-discovery rates are 8-24% even on curated data, so nothing downstream
    should ever quote a bare IS count without these numbers next to it.

WHERE THE INPUT COMES FROM
    rule isescan runs ISEScan 1.7.3 on the final assembly
    (02.assembly/{sample}/contigs_final.fasta) WITHOUT --removeShortIS, so both
    complete and partial (truncated / single-copy, no perfect terminal repeat)
    copies are reported and we tier them ourselves instead of letting the tool
    silently drop the weak ones.

    ISEScan mirrors the input path into its output directory: for an input
    02.assembly/{sample}/contigs_final.fasta and --output 08.mobilome/{sample}/isescan
    the results land in 08.mobilome/{sample}/isescan/{sample}/contigs_final.fasta.tsv.
    --isescan-out therefore accepts EITHER that .tsv file directly OR the output
    directory, and we find the results file inside it.

WHAT THIS SCRIPT PRODUCES (and what consumes it)
    --out-table    one row per IS copy, normalised column names, 1-based inclusive
                   coordinates. This is the IS side of the AMR x MGE
                   co-localisation step (WP-D), which will convert it to BED and
                   intersect it with the AMRFinderPlus hits.
    --out-summary  one row per sample: totals, complete vs partial, how many IS sit
                   within --boundary-bp of a contig end, and that fraction. This is
                   the QC metric the spec asks for; it travels with the report so a
                   reader can see how fragmented the evidence is.
    --out-audit    one row per IS record we DROPPED or FLAGGED, with an explicit
                   reason column, PLUS one sample-level row when there was nothing
                   to report at all, saying which of the two very different
                   situations applies: ISEScan wrote no results file (so the run
                   itself, not the biology, is why the table is empty) or ISEScan
                   reported zero IS for this genome. BacFlux convention (see
                   contig_taxonomy_decisions.tsv): every filtering decision must be
                   inspectable afterwards; nothing disappears silently, and "the
                   tool produced nothing" must never look like "the genome has
                   nothing".

COORDINATES
    ISEScan reports isBegin/isEnd as 1-based and inclusive of both ends. We keep
    them exactly as they are — NO conversion happens here — because 1-based
    inclusive is what Bakta, AMRFinderPlus and GFF all use, so every table in this
    module joins without an off-by-one trap. Whoever writes the BED file for
    bedtools downstream must subtract 1 from `start` there (BED is 0-based,
    half-open); that conversion belongs in the BED writer, not here, so there is
    exactly one place to check it.

DEFENSIVE PARSING
    ISEScan's .tsv column NAMES were read from the v1.7.3 source, but this parser
    still keys on the header line rather than on column positions, so a future
    version that inserts a column cannot silently shift our coordinates. If the
    header is missing a column we truly need, we STOP with a message naming the
    file and the column — a loud failure beats a table of wrong coordinates.
    A genome with no IS at all, by contrast, is a perfectly normal biological
    result and is handled gracefully (empty but well-formed outputs, exit 0).

TWO WAYS TO RUN IT
    1) main use — normalise ISEScan output for one sample:
         isescan_to_table.py --sample S \
             --isescan-out 08.mobilome/S/isescan \
             --contig-lengths 08.mobilome/S/contig_lengths.tsv \
             --boundary-bp 100 \
             --out-table  08.mobilome/S/S_is_elements.tsv \
             --out-summary 08.mobilome/S/S_is_summary.tsv \
             --out-audit  08.mobilome/S/S_is_discarded.tsv

    2) helper used by the same rule to make that --contig-lengths file from the
       assembly FASTA (plain stdlib FASTA reading, no Biopython):
         isescan_to_table.py --genome-fasta 02.assembly/S/contigs_final.fasta \
             --out-contig-lengths 08.mobilome/S/contig_lengths.tsv

Standard library only — this runs inside the ISEScan conda env, which carries no
pandas. Unit tests: workflow/scripts/mobilome/test_isescan_to_table.py
"""

import argparse
import csv
import os
import sys


# ── What ISEScan gives us ────────────────────────────────────────────────────

# The 24 columns ISEScan 1.7.3 writes into <contigs>.fasta.tsv, in the exact
# order the tool emits them (read from the v1.7.3 source, not from the docs).
# We normally look columns up BY NAME from the file's own header; this list is
# only the fallback for a results file that somehow arrives without a header, and
# it is also the reference a reader can check the mapping against.
ISESCAN_COLUMNS = [
    "seqID",         # contig name = first whitespace token of the FASTA header
    "family",        # IS family, e.g. IS3, IS6/IS26, IS200/IS605
    "cluster",       # ISEScan's within-family cluster label
    "isBegin",       # IS start on the contig, 1-based inclusive
    "isEnd",         # IS end   on the contig, 1-based inclusive
    "isLen",         # isEnd - isBegin + 1
    "ncopy4is",      # copies of this IS that ISEScan located in the assembly
    "start1", "end1",  # left  terminal inverted repeat (TIR) coordinates
    "start2", "end2",  # right terminal inverted repeat coordinates
    "score",         # alignment score of the two TIR halves against each other
    "irId",          # identical bases in that TIR alignment (a COUNT, not a %)
    "irLen",         # TIR length; 0 when no inverted repeat was found
    "nGaps",         # gaps in the TIR alignment
    "orfBegin", "orfEnd",  # the predicted transposase ORF
    "strand",        # strand of that transposase; EMPTY when none was predicted
    "orfLen",
    "E-value",       # best transposase pHMM E-value across all copies
    "E-value4copy",  # transposase pHMM E-value for THIS copy
    "type",          # 'c' = complete copy, 'p' = partial copy
    "ov",            # transposase copy number (NOT the hmmer overlap the README claims)
    "tir",           # the two TIR sequences as "seq1:seq2"
]

# Without these three we cannot say where an IS is, so their absence is a hard
# stop rather than something to work around.
REQUIRED_ISESCAN_COLUMNS = ["seqID", "isBegin", "isEnd"]

# These we use but can live without: if a future ISEScan drops or renames one we
# warn, write NA, and carry on, because losing (say) the family label is much
# less damaging than losing the whole sample's IS inventory.
OPTIONAL_ISESCAN_COLUMNS = [
    "family", "cluster", "strand", "type", "irLen", "irId", "nGaps", "score",
    "ncopy4is", "E-value", "E-value4copy", "orfBegin", "orfEnd", "tir",
]

# ISEScan's completeness letter. 'c' = a complete copy (full length, with a
# proper terminal inverted repeat); 'p' = partial, i.e. truncated or a single
# copy without a perfect TIR. We deliberately KEEP the partials (ISEScan is run
# without --removeShortIS) because on a fragmented short-read assembly a partial
# call is often a real IS cut in half by a contig break — dropping them would
# throw away exactly the evidence the boundary flag exists to quantify.
ISESCAN_TYPE_COMPLETE = "c"
ISESCAN_TYPE_PARTIAL = "p"


# ── What we write ────────────────────────────────────────────────────────────

# Sequence Ontology term SO:0000973 "insertion_sequence" — the vocabulary the
# whole mobilome module labels its element types with, so IS rows, prophage rows
# and (later) ICE rows can be concatenated into one MGE table.
MGE_TYPE_INSERTION_SEQUENCE = "insertion_sequence"

# The deliverable IS table, in the order the columns are written. One row per IS
# copy ISEScan located on the assembly.
OUTPUT_COLUMNS = [
    "sample",
    "mge_id",                     # contig|insertion_sequence-start:end (stable join key)
    "mge_type",                   # always insertion_sequence here (SO:0000973)
    "contig",
    "start",                      # 1-based inclusive (see module docstring)
    "end",                        # 1-based inclusive
    "length_bp",                  # end - start + 1, recomputed from the coordinates
    "strand",                     # + / - of the transposase; NA when none predicted
    "is_family",
    "cluster",
    "is_complete",                # TRUE = ISEScan 'c', FALSE = 'p', NA = unrecognised
    "is_type",                    # ISEScan's raw letter, kept so nothing is lost
    "ir_present",                 # TRUE when a terminal inverted repeat was found
    "ir_len_bp",
    "ir_identical_bp",            # identical bases in the TIR alignment (a count)
    "ir_gaps",
    "ir_score",
    "n_copies_is",                # ISEScan's located copy number (also a floor)
    "evalue",                     # best transposase E-value across copies
    "evalue_this_copy",
    "tpase_orf_start",            # transposase ORF, useful for IS-inside-CDS work
    "tpase_orf_end",
    "contig_length",
    "dist_to_contig_start",       # bp of contig before the IS
    "dist_to_contig_end",         # bp of contig after the IS
    "dist_to_nearest_contig_end",
    "at_contig_boundary",         # TRUE when the nearest distance <= --boundary-bp
]

# One row per sample. Wide (not key/value) so several samples' summaries can be
# concatenated into one table for the run-level report.
SUMMARY_COLUMNS = [
    "sample",
    "isescan_results_file",       # which file was parsed, or NONE — provenance
    "boundary_bp",                # the threshold used, so the fraction is interpretable
    "min_length_bp",
    "n_is_total",                 # a FLOOR on the true IS content, never a count
    "n_is_complete",
    "n_is_partial",
    "n_is_completeness_unknown",
    "n_is_with_ir",
    "n_contigs_with_is",
    "n_is_at_contig_boundary",
    "fraction_at_contig_boundary",   # the key QC metric (spec §6 item 2)
    "n_records_dropped",
    "n_records_flagged",
]

# One row per IS record we dropped or flagged, with the reason spelled out, plus
# (at most) one sample-level row explaining an empty table — see
# startup_audit_rows. contig/start/end are NA on that sample-level row, because it
# is about the whole run rather than about one element.
AUDIT_COLUMNS = [
    "sample",
    "contig",
    "start",
    "end",
    "action",     # dropped | kept_flagged | input_missing | input_empty
    "reason",     # short machine-readable token, see normalise_records and startup_audit_rows
    "detail",     # human-readable explanation + the offending raw line
]

# Reporting convention, NOT biology: when more than this share of the located IS
# sit at a contig end, the assembly broke at its repeats badly enough that the
# located inventory is clearly a floor, and we say so in the log. Changing the
# number changes only how loud we are, never what is written to the table.
BOUNDARY_FRACTION_WARN = 0.5

# Directories ISEScan keeps its intermediates in, under the same --output tree.
# We never look for results inside them.
ISESCAN_INTERMEDIATE_DIRS = {"proteome", "hmm"}

# How much of an offending line to quote in the audit file. Whole lines would
# make the audit unreadable; this is enough to recognise the record.
AUDIT_DETAIL_MAX_CHARS = 200


# ── Small shared helpers ─────────────────────────────────────────────────────

def _warn(message):
    """Write one warning line to stderr, tagged so it is greppable in a Snakemake
    log where many rules interleave their output."""
    sys.stderr.write("[isescan_to_table] WARNING: " + message + "\n")


def _text_or_na(value):
    """Return a stripped field value, or the string 'NA' when it is absent or
    empty. Every table in BacFlux writes a literal NA rather than an empty cell,
    so a missing value is visible instead of looking like a formatting slip."""
    if value is None:
        return "NA"
    text = str(value).strip()
    return text if text else "NA"


def _tsv_bool(value):
    """Render a yes/no/unknown answer for the TSV.

    We write TRUE / FALSE (not Python's True / False) because these tables are
    read in R as often as in Python, and R's read.delim turns TRUE/FALSE straight
    into a logical column. Unknown stays 'NA', which both languages understand.
    """
    if value is None:
        return "NA"
    return "TRUE" if value else "FALSE"


def _int_or_none(value):
    """Parse an integer field, returning None instead of raising when the field
    is missing or not a number. Used for coordinates, where 'not a number' is a
    reason to send the record to the audit file, not a reason to crash."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _split_row(line):
    """Split one ISEScan result line into fields.

    The .tsv is tab-separated, which is what we expect and what preserves empty
    fields (ISEScan leaves `strand` empty for an IS with no predicted
    transposase). Splitting on tabs is therefore the only safe option, and we
    only fall back to whitespace splitting for a line that contains no tab at all
    — that would be ISEScan's space-padded .raw file, and the caller is warned
    about the field count before we trust it.
    """
    if "\t" in line:
        return line.split("\t")
    return line.split()


def write_tsv(path, columns, rows):
    """Write rows (list of dicts) as a tab-separated table with a header.

    The header is ALWAYS written, even for zero rows: an explicit empty table is
    a result ("this genome has no IS"), whereas a missing or headerless file
    looks like a crashed rule to whoever reads the output next.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ── Helper: contig lengths from the assembly FASTA ───────────────────────────

def fasta_contig_lengths(path):
    """Read a genome FASTA and return {contig_id: length_in_bp}.

    Input:  the assembly this sample's ISEScan run was given
            (02.assembly/{sample}/contigs_final.fasta).
    Output: a dict used to work out how close each IS sits to a contig end.

    The contig ID is the FIRST whitespace token of the header, because that is
    what ISEScan (and Bakta, and AMRFinderPlus) put in their `seqID`/`Contig id`
    columns — using the whole description line would break every join.
    Sequence is counted line by line without holding it in memory, and any
    whitespace inside a sequence line is ignored, so wrapped FASTA is fine.
    """
    lengths = {}
    current_contig = None
    current_length = 0

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                # Store the contig we just finished before starting the next one.
                if current_contig is not None:
                    lengths[current_contig] = current_length
                header = line[1:].strip()
                current_contig = header.split()[0] if header else ""
                current_length = 0
            else:
                # Sequence line: count only the residues, not the line break or
                # any stray spaces some tools insert.
                current_length += len(line.strip().replace(" ", ""))

    # The last contig in the file has no following '>' to trigger the store.
    if current_contig is not None:
        lengths[current_contig] = current_length

    return lengths


def write_contig_lengths(path, lengths):
    """Write {contig: length} as the two-column TSV the main mode reads back.

    Sorted by contig name so the file is byte-identical between runs, which keeps
    Snakemake's re-run logic and any diffing honest.
    """
    rows = [{"contig": contig, "length": lengths[contig]} for contig in sorted(lengths)]
    write_tsv(path, ["contig", "length"], rows)


def read_contig_lengths(path):
    """Read the two-column contig-length TSV back into {contig: length}.

    Accepts the file with or without its 'contig<TAB>length' header (a user may
    well hand-make this file with awk), by simply skipping any line whose second
    field is not a number. Lines that are neither a header nor a valid pair are
    reported and skipped rather than guessed at.
    """
    lengths = {}
    if not path or not os.path.exists(path):
        raise ValueError(f"Contig-lengths file not found: {path}")

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = _split_row(line.rstrip("\n"))
            if len(fields) < 2:
                continue
            contig = fields[0].strip()
            length = _int_or_none(fields[1])
            if not contig or length is None:
                # The header row lands here, which is exactly what we want.
                continue
            lengths[contig] = length

    if not lengths:
        raise ValueError(
            f"No usable 'contig<TAB>length' rows in {path}. Expected the two-column "
            "table written by --genome-fasta/--out-contig-lengths."
        )
    return lengths


# ── Finding ISEScan's results file ───────────────────────────────────────────

def find_isescan_results(path):
    """Resolve --isescan-out to the single ISEScan results .tsv, or None.

    Why this is not just a path: ISEScan re-creates the input's parent directory
    inside its --output directory, so the results file sits at
    {output}/{parent-dir-of-the-input-fasta}/{input-fasta-name}.tsv. Rather than
    reproduce that rule (and break when a path changes), the rule can hand us the
    output DIRECTORY and we look inside it. A direct path to the .tsv also works.

    Returns None when nothing is there — which happens routinely, because ISEScan
    writes NO output files at all for a genome in which it found no IS. That is a
    normal biological result, so it is not an error here.
    """
    if not path or not os.path.exists(path):
        return None

    if os.path.isfile(path):
        return path

    # A directory: walk it, ignoring ISEScan's own intermediate trees (proteome/,
    # hmm/) so we cannot accidentally pick up something that is not a result.
    candidates = []
    for directory, subdirectories, filenames in os.walk(path):
        subdirectories[:] = [
            name for name in subdirectories if name not in ISESCAN_INTERMEDIATE_DIRS
        ]
        for filename in filenames:
            if filename.endswith(".tsv"):
                candidates.append(os.path.join(directory, filename))

    if not candidates:
        return None

    candidates.sort()
    if len(candidates) > 1:
        # One assembly should produce one results file. More than one means the
        # output directory was reused for another sample (ISEScan also caches its
        # intermediates, so stale results really can linger). Say so loudly.
        _warn(
            f"found {len(candidates)} .tsv files under {path}; using the first one. "
            f"All: {', '.join(candidates)}"
        )
    return candidates[0]


# ── Reading and normalising ISEScan's table ──────────────────────────────────

def read_isescan_table(path):
    """Read the ISEScan results .tsv into a list of dicts keyed by ISEScan's own
    column names.

    Input:  the file found by find_isescan_results (may be None or empty).
    Output: a list of per-IS dicts, plus two bookkeeping keys per record:
            '_raw' (the original line, for the audit file) and '_n_fields'.
            An empty list means "no IS in this genome", which is normal.

    Two different kinds of problem, treated differently on purpose:
      * no file / empty file / header only  -> return [], the caller writes empty
        but well-formed outputs and exits 0. A genome with no detectable IS is a
        real result, not a failure.
      * a header that is missing a column we need to place the IS on the contig
        -> raise ValueError naming the file and the column. Guessing here would
        produce a table of wrong coordinates that looks perfectly fine, which is
        far worse than stopping.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []

    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return []

    first_row_fields = _split_row(lines[0])

    if "seqID" in first_row_fields:
        # Normal case: trust the file's own header, so an added/moved column in a
        # future ISEScan cannot shift our coordinates.
        column_names = first_row_fields
        data_lines = lines[1:]
    else:
        # No header. Fall back to the documented v1.7.3 column order, but only
        # when the shape matches exactly — the field count must be right AND the
        # three coordinate-ish columns must actually contain integers. Otherwise
        # we would be inventing a mapping.
        # (Not observed in practice; kept because the exact header behaviour of
        # future ISEScan versions is unverified.)
        begin_index = ISESCAN_COLUMNS.index("isBegin")
        end_index = ISESCAN_COLUMNS.index("isEnd")
        shape_ok = (
            len(first_row_fields) == len(ISESCAN_COLUMNS)
            and _int_or_none(first_row_fields[begin_index]) is not None
            and _int_or_none(first_row_fields[end_index]) is not None
        )
        if not shape_ok:
            raise ValueError(
                f"{path}: first line is neither an ISEScan header (no 'seqID' column) "
                f"nor a {len(ISESCAN_COLUMNS)}-field data row with integer isBegin/isEnd, "
                f"so the columns cannot be identified. First line: {lines[0][:AUDIT_DETAIL_MAX_CHARS]!r}"
            )
        _warn(
            f"{path} has no header line; falling back to ISEScan 1.7.3's documented "
            "24-column order. Check the coordinates in the output table."
        )
        column_names = list(ISESCAN_COLUMNS)
        data_lines = lines

    missing_required = [c for c in REQUIRED_ISESCAN_COLUMNS if c not in column_names]
    if missing_required:
        raise ValueError(
            f"{path}: ISEScan results file is missing required column(s) "
            f"{', '.join(missing_required)}. Without them an IS cannot be placed on a "
            f"contig. Header seen: {column_names}"
        )

    missing_optional = [c for c in OPTIONAL_ISESCAN_COLUMNS if c not in column_names]
    if missing_optional:
        _warn(
            f"{path}: column(s) {', '.join(missing_optional)} not present; the "
            "matching output fields will be NA. Check whether the ISEScan version "
            "changed its output format."
        )

    records = []
    for line in data_lines:
        fields = _split_row(line)
        # zip stops at the shorter of the two, so a short (ragged) line simply
        # yields fewer keys; _n_fields is what we actually judge the line by.
        record = dict(zip(column_names, fields))
        record["_raw"] = line
        record["_n_fields"] = len(fields)
        record["_n_columns"] = len(column_names)
        records.append(record)
    return records


def normalise_records(sample, records, contig_lengths, boundary_bp, min_length_bp):
    """Turn ISEScan's records into the tidy output rows, and decide what to drop.

    Input:  records from read_isescan_table, the {contig: length} map from the
            assembly, the boundary threshold in bp, and a minimum IS length.
    Output: (rows, audit_rows) — the deliverable table, and one audit row for
            every record dropped or flagged, each with an explicit reason.

    The drop rules, in the order they are applied, and why each one exists:
      unexpected_field_count       a ragged line means the column->value mapping
                                   for THAT line cannot be trusted, so its
                                   coordinates might be silently wrong.
      missing_contig_id            nothing to join the IS to.
      unparseable_coordinates      isBegin/isEnd not integers.
      invalid_coordinate_range     start < 1, or end before start.
      below_min_length_bp          only when the user asks for it (default 0).
      contig_not_in_contig_lengths the contig-lengths file does not describe the
                                   assembly ISEScan ran on, so the contig-end
                                   distances — the whole honesty mechanism — could
                                   not be computed for this IS.
      coordinates_beyond_contig_length  same mismatch, seen from the other side.

    The one FLAG (kept, not dropped) is an unrecognised completeness letter: we
    would rather report an IS with is_complete=NA than lose a real element over a
    label we did not expect.
    """
    rows = []
    audit_rows = []

    def audit(contig, start, end, action, reason, detail):
        audit_rows.append({
            "sample": sample,
            "contig": _text_or_na(contig),
            "start": _text_or_na(start),
            "end": _text_or_na(end),
            "action": action,
            "reason": reason,
            "detail": detail[:AUDIT_DETAIL_MAX_CHARS],
        })

    for record in records:
        raw_line = record.get("_raw", "")
        contig = (record.get("seqID") or "").strip()
        start_text = record.get("isBegin")
        end_text = record.get("isEnd")

        # --- 1. the line itself has to have the expected shape ---
        if record.get("_n_fields") != record.get("_n_columns"):
            audit(
                contig, start_text, end_text, "dropped", "unexpected_field_count",
                f"line has {record.get('_n_fields')} fields, header has "
                f"{record.get('_n_columns')}; column mapping unreliable: {raw_line}",
            )
            continue

        # --- 2. we need a contig to place the IS on ---
        if not contig:
            audit(
                contig, start_text, end_text, "dropped", "missing_contig_id",
                f"empty seqID: {raw_line}",
            )
            continue

        # --- 3. coordinates must be readable and sane ---
        start = _int_or_none(start_text)
        end = _int_or_none(end_text)
        if start is None or end is None:
            audit(
                contig, start_text, end_text, "dropped", "unparseable_coordinates",
                f"isBegin/isEnd are not integers: {raw_line}",
            )
            continue
        if start < 1 or end < start:
            audit(
                contig, start, end, "dropped", "invalid_coordinate_range",
                f"expected 1 <= isBegin <= isEnd, got {start}..{end}: {raw_line}",
            )
            continue

        # ISEScan's coordinates are 1-based and inclusive of both ends, so the
        # length includes both the first and last base.
        length_bp = end - start + 1

        # --- 4. optional minimum length ---
        # Default 0, i.e. off. The spec's general 500 bp floor for MGE predictions
        # is deliberately NOT applied to IS here: a short call is usually a
        # partial IS cut by a contig break, and those are exactly what we want to
        # keep and tier honestly.
        if min_length_bp > 0 and length_bp < min_length_bp:
            audit(
                contig, start, end, "dropped", "below_min_length_bp",
                f"IS length {length_bp} bp < --min-length-bp {min_length_bp}",
            )
            continue

        # --- 5. the contig must be one we know the length of ---
        if contig not in contig_lengths:
            audit(
                contig, start, end, "dropped", "contig_not_in_contig_lengths",
                "contig is absent from the contig-lengths file, so the distance to "
                "the contig end (the short-read honesty flag) cannot be computed; "
                "check that ISEScan and the length file used the same assembly",
            )
            continue

        contig_length = contig_lengths[contig]
        if end > contig_length:
            audit(
                contig, start, end, "dropped", "coordinates_beyond_contig_length",
                f"IS ends at {end} but contig {contig} is only {contig_length} bp; "
                "the contig-lengths file does not match the assembly ISEScan ran on",
            )
            continue

        # --- 6. completeness: ISEScan's own complete/partial letter ---
        raw_type = (record.get("type") or "").strip().lower()
        if raw_type == ISESCAN_TYPE_COMPLETE:
            is_complete = True
        elif raw_type == ISESCAN_TYPE_PARTIAL:
            is_complete = False
        else:
            is_complete = None
            audit(
                contig, start, end, "kept_flagged", "unknown_isescan_type_value",
                f"'type' was {record.get('type')!r}, expected 'c' or 'p'; the IS is "
                "kept with is_complete=NA",
            )

        # --- 7. terminal inverted repeats ---
        # A pair of inverted repeats at the two ends is the structural signature of
        # an intact IS, so its presence is one of the confidence signals the
        # module tiers on. ISEScan reports irLen = 0 when it found none; if that
        # column is unavailable we fall back to whether it printed the repeat
        # sequences at all.
        ir_len = _int_or_none(record.get("irLen"))
        if ir_len is not None:
            ir_present = ir_len > 0
        else:
            tir_text = (record.get("tir") or "").strip()
            ir_present = bool(tir_text) and tir_text not in {"-", "NA", "."}

        # --- 8. how close is this IS to where the assembly broke? ---
        # An IS sitting at a contig end is the assembler telling us it could not
        # resolve that repeat: whatever was next to the IS is missing from the
        # assembly, so any statement about its genomic context is unsupported.
        distance_to_start = start - 1
        distance_to_end = contig_length - end
        distance_to_nearest_end = min(distance_to_start, distance_to_end)
        at_boundary = distance_to_nearest_end <= boundary_bp

        # A stable identifier, following the module's ID convention
        # contig|mge_type-start:end, so this IS can be referred to from the
        # co-localisation table and the final report without re-deriving it.
        mge_id = f"{contig}|{MGE_TYPE_INSERTION_SEQUENCE}-{start}:{end}"

        rows.append({
            "sample": sample,
            "mge_id": mge_id,
            "mge_type": MGE_TYPE_INSERTION_SEQUENCE,
            "contig": contig,
            "start": start,
            "end": end,
            "length_bp": length_bp,
            "strand": _text_or_na(record.get("strand")),
            "is_family": _text_or_na(record.get("family")),
            "cluster": _text_or_na(record.get("cluster")),
            "is_complete": _tsv_bool(is_complete),
            "is_type": _text_or_na(record.get("type")),
            "ir_present": _tsv_bool(ir_present),
            "ir_len_bp": _text_or_na(record.get("irLen")),
            "ir_identical_bp": _text_or_na(record.get("irId")),
            "ir_gaps": _text_or_na(record.get("nGaps")),
            "ir_score": _text_or_na(record.get("score")),
            "n_copies_is": _text_or_na(record.get("ncopy4is")),
            "evalue": _text_or_na(record.get("E-value")),
            "evalue_this_copy": _text_or_na(record.get("E-value4copy")),
            "tpase_orf_start": _text_or_na(record.get("orfBegin")),
            "tpase_orf_end": _text_or_na(record.get("orfEnd")),
            "contig_length": contig_length,
            "dist_to_contig_start": distance_to_start,
            "dist_to_contig_end": distance_to_end,
            "dist_to_nearest_contig_end": distance_to_nearest_end,
            "at_contig_boundary": _tsv_bool(at_boundary),
        })

    # Deterministic order: by contig, then by position along it. Makes the file
    # diffable between runs and reads like a walk along each contig.
    rows.sort(key=lambda row: (row["contig"], row["start"], row["end"]))
    return rows, audit_rows


def startup_audit_rows(sample, isescan_out, results_file, n_records):
    """Write down WHY this sample has no IS to report, when that is the case.

    Input:  the --isescan-out path as the rule passed it, the results file
            find_isescan_results resolved it to (or None), and how many records
            read_isescan_table got out of that file.
    Output: a list of zero or one audit rows, written to --out-audit ahead of the
            per-record drop/flag rows. Nothing else consumes it; it exists purely
            so a human reading the audit file can interpret an empty IS table.

    Why this matters, in biology terms: an empty {sample}_is_elements.tsv can mean
    two completely different things, and a reader judging a mobility call needs to
    know which.
      * ISEScan reported no insertion sequences  -> a real biological statement
        about this genome (still a FLOOR, see the module docstring, but a result).
      * ISEScan wrote no results file at all     -> we know nothing about this
        genome's IS content. ISEScan does write nothing when it finds no IS, so
        this is the EXPECTED look of a genuinely IS-free genome — but a crashed or
        killed run looks exactly the same from here, which is why we say so
        instead of quietly reporting zero.
    An absent --isescan-out path is a third case and a different kind of problem
    (the rule's wiring, not the tool), so it gets its own reason token.

    BacFlux hard rule (CLAUDE.md): every such decision gets an audit row with a
    reason, so "the tool produced nothing" is never indistinguishable from "the
    biology genuinely had nothing".
    """
    # Records were parsed, so the table is not empty for any of these reasons and
    # there is nothing to explain here.
    if results_file is not None and n_records > 0:
        return []

    if results_file is None and (not isescan_out or not os.path.exists(isescan_out)):
        # The path the rule handed us is not there at all. ISEScan's own rule
        # creates its output directory unconditionally, so this normally means the
        # directory was removed or the path was mis-wired — not a statement about
        # the genome.
        action = "input_missing"
        reason = "isescan_output_path_missing"
        detail = (
            f"--isescan-out '{isescan_out}' does not exist, so ISEScan's results were "
            "never looked at. Nothing about this genome's insertion sequences is known; "
            "the empty IS table is a wiring problem, not a biological result."
        )
    elif results_file is None:
        # The directory is there but holds no results .tsv. This is exactly what a
        # genuinely IS-free genome looks like, and also exactly what a crashed run
        # looks like, so we report both readings rather than pick one.
        action = "input_missing"
        reason = "isescan_wrote_no_results_file"
        detail = (
            f"no ISEScan results .tsv under '{isescan_out}'. ISEScan writes no output at "
            "all when it finds no insertion sequences, so the expected reading is 'no IS "
            "detected in this genome', but an ISEScan run that failed leaves exactly the "
            "same empty directory: check the isescan log for this sample before quoting "
            "zero."
        )
    else:
        # A results file WAS found and read; it simply lists no IS. This is the
        # one case where the empty table is a statement about the genome.
        action = "input_empty"
        reason = "isescan_results_file_has_no_rows"
        detail = (
            f"ISEScan results file '{results_file}' was found but contains no IS rows "
            "(empty, or header only), i.e. ISEScan reported no insertion sequences for "
            "this genome. Unlike a missing results file, this is a statement about the "
            "assembly, not about the run."
        )

    return [{
        "sample": sample,
        "contig": "NA",
        "start": "NA",
        "end": "NA",
        "action": action,
        "reason": reason,
        # Not truncated to AUDIT_DETAIL_MAX_CHARS: unlike the per-record rows this
        # detail quotes no raw ISEScan line, it is a whole sentence a reader needs.
        "detail": detail,
    }]


def summarise(sample, rows, audit_rows, boundary_bp, min_length_bp, results_file):
    """Build the one-row per-sample QC summary.

    Input:  the normalised rows and audit rows, plus the settings used.
    Output: a single dict following SUMMARY_COLUMNS.

    The number that matters is fraction_at_contig_boundary. It says what share of
    the IS we located sit within --boundary-bp of a contig end, i.e. what share of
    the inventory has unknown genomic context because the assembly broke there.
    A high fraction is not a bug in this script — it is the expected behaviour of
    a short-read assembly of an IS-rich genome, and it is the reason the total is
    reported as a floor rather than as a count.
    """
    n_total = len(rows)
    n_complete = sum(1 for row in rows if row["is_complete"] == "TRUE")
    n_partial = sum(1 for row in rows if row["is_complete"] == "FALSE")
    n_unknown = sum(1 for row in rows if row["is_complete"] == "NA")
    n_with_ir = sum(1 for row in rows if row["ir_present"] == "TRUE")
    n_at_boundary = sum(1 for row in rows if row["at_contig_boundary"] == "TRUE")
    contigs_with_is = {row["contig"] for row in rows}

    # With no IS at all there is no fraction to speak of; 0 is the honest,
    # non-crashing answer and is written with the same 4 decimals as any other
    # value so the column type never changes between samples.
    fraction_at_boundary = (n_at_boundary / n_total) if n_total else 0.0

    n_dropped = sum(1 for row in audit_rows if row["action"] == "dropped")
    n_flagged = sum(1 for row in audit_rows if row["action"] == "kept_flagged")

    return {
        "sample": sample,
        "isescan_results_file": results_file if results_file else "NONE",
        "boundary_bp": boundary_bp,
        "min_length_bp": min_length_bp,
        "n_is_total": n_total,
        "n_is_complete": n_complete,
        "n_is_partial": n_partial,
        "n_is_completeness_unknown": n_unknown,
        "n_is_with_ir": n_with_ir,
        "n_contigs_with_is": len(contigs_with_is),
        "n_is_at_contig_boundary": n_at_boundary,
        "fraction_at_contig_boundary": f"{fraction_at_boundary:.4f}",
        "n_records_dropped": n_dropped,
        "n_records_flagged": n_flagged,
    }


# ── Command line ─────────────────────────────────────────────────────────────

def build_parser():
    """The two usages are described in the module docstring. Everything is
    declared optional here and checked in main(), so that one script can serve
    both the small FASTA-length helper and the main normalisation step without
    argparse subcommands (which would change the documented flag layout)."""
    parser = argparse.ArgumentParser(
        description="Normalise ISEScan output into the mobilome module's IS table.",
    )
    # Main mode.
    parser.add_argument("--sample", help="Sample name, written into every output row.")
    parser.add_argument("--isescan-out",
                        help="ISEScan results .tsv, OR the ISEScan --output directory.")
    parser.add_argument("--contig-lengths",
                        help="Two-column TSV (contig, length) for the same assembly.")
    parser.add_argument("--boundary-bp", type=int, default=100,
                        help="An IS within this many bp of a contig end is flagged "
                             "as sitting at a contig boundary (default: 100).")
    parser.add_argument("--min-length-bp", type=int, default=0,
                        help="Drop IS shorter than this, recording them in the audit "
                             "file. Default 0 = keep everything, because short calls "
                             "are usually partial IS cut by a contig break.")
    parser.add_argument("--out-table", help="Destination for the tidy IS table.")
    parser.add_argument("--out-summary", help="Destination for the per-sample QC summary.")
    parser.add_argument("--out-audit", help="Destination for the dropped/flagged audit TSV.")
    # Helper mode.
    parser.add_argument("--genome-fasta",
                        help="Helper mode: assembly FASTA to measure contig lengths from.")
    parser.add_argument("--out-contig-lengths",
                        help="Helper mode: destination for the contig-length TSV.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    # ── Helper mode: assembly FASTA -> contig lengths TSV ────────────────────
    # Run by the same Snakemake rule, just before the main mode, so the main mode
    # always has the lengths of the exact assembly ISEScan was given.
    if args.genome_fasta or args.out_contig_lengths:
        if not (args.genome_fasta and args.out_contig_lengths):
            sys.stderr.write(
                "ERROR: --genome-fasta and --out-contig-lengths must be given together.\n"
            )
            return 1
        if not os.path.exists(args.genome_fasta):
            sys.stderr.write(f"ERROR: genome FASTA not found: {args.genome_fasta}\n")
            return 1
        lengths = fasta_contig_lengths(args.genome_fasta)
        write_contig_lengths(args.out_contig_lengths, lengths)
        total_bp = sum(lengths.values())
        print(f"Wrote lengths for {len(lengths)} contig(s), {total_bp} bp total, "
              f"to {args.out_contig_lengths}")
        return 0

    # ── Main mode ────────────────────────────────────────────────────────────
    required = {
        "--sample": args.sample,
        "--isescan-out": args.isescan_out,
        "--contig-lengths": args.contig_lengths,
        "--out-table": args.out_table,
        "--out-summary": args.out_summary,
        "--out-audit": args.out_audit,
    }
    missing = [flag for flag, value in required.items() if not value]
    if missing:
        sys.stderr.write(
            "ERROR: missing required argument(s): " + ", ".join(missing) + "\n"
        )
        return 1

    # Contig lengths first: without them there is no boundary signal at all, and a
    # missing/garbled file is a rule wiring problem, not a biological result.
    try:
        contig_lengths = read_contig_lengths(args.contig_lengths)
    except ValueError as error:
        sys.stderr.write(f"ERROR: {error}\n")
        return 1

    # Locate ISEScan's results. None is normal: ISEScan writes no files at all for
    # a genome in which it found no IS.
    results_file = find_isescan_results(args.isescan_out)
    if results_file is None:
        _warn(
            f"no ISEScan results file found at {args.isescan_out}. ISEScan writes "
            "nothing when it finds no IS, so this is treated as 'no insertion "
            "sequences detected' and empty (but complete) outputs are written."
        )

    # A results file whose columns we cannot identify is a different matter — that
    # is a format change, and continuing would mean publishing wrong coordinates.
    try:
        records = read_isescan_table(results_file)
    except ValueError as error:
        sys.stderr.write(f"ERROR: {error}\n")
        return 1

    rows, audit_rows = normalise_records(
        args.sample, records, contig_lengths, args.boundary_bp, args.min_length_bp
    )
    # When there is nothing to report, say in the audit file WHY there is nothing:
    # "ISEScan found no IS in this genome" and "ISEScan produced no results file"
    # are different statements and a reader must not have to guess which happened.
    # These rows are deliberately NOT counted as dropped/flagged records in the
    # summary — nothing was filtered out, there was simply nothing to filter.
    startup_audit = startup_audit_rows(
        args.sample, args.isescan_out, results_file, len(records)
    )
    summary = summarise(
        args.sample, rows, audit_rows, args.boundary_bp, args.min_length_bp, results_file
    )

    write_tsv(args.out_table, OUTPUT_COLUMNS, rows)
    write_tsv(args.out_summary, SUMMARY_COLUMNS, [summary])
    write_tsv(args.out_audit, AUDIT_COLUMNS, startup_audit + audit_rows)

    # A short, honest line for the run log. Never used for filtering — it just
    # makes the fragmentation visible without opening the summary file.
    print(
        f"Sample {args.sample}: {summary['n_is_total']} IS located "
        f"({summary['n_is_complete']} complete, {summary['n_is_partial']} partial) "
        f"on {summary['n_contigs_with_is']} contig(s); "
        f"{summary['n_is_at_contig_boundary']} within {args.boundary_bp} bp of a "
        f"contig end (fraction {summary['fraction_at_contig_boundary']}); "
        f"{summary['n_records_dropped']} record(s) dropped, "
        f"{summary['n_records_flagged']} flagged (see {args.out_audit})."
    )
    # Repeat the sample-level explanation in the log, so an empty IS table is
    # explained both in the run log and in the audit file that outlives it.
    for startup_row in startup_audit:
        print(f"  nothing to report - reason: {startup_row['reason']}. "
              f"{startup_row['detail']}")

    if float(summary["fraction_at_contig_boundary"]) >= BOUNDARY_FRACTION_WARN:
        _warn(
            f"{summary['fraction_at_contig_boundary']} of the IS calls for "
            f"{args.sample} sit within {args.boundary_bp} bp of a contig end. The "
            "assembly broke at its repeats, so this IS inventory is a FLOOR, not a "
            "count, and IS-to-gene distances near those ends are unsupported."
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
