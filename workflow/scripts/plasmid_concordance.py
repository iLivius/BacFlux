#!/usr/bin/env python3
"""Join Platon and geNomad plasmid calls into one confidence-tiered table (D9).

BacFlux v2 calls plasmids with TWO independent tools and reports where they
agree, instead of trusting either one alone:

- Platon (rule plasmid_search) classifies every assembly contig as plasmid or
  chromosome from replicon-distribution scores (RDS). BacFlux also keeps a v1
  "supplementary" check: for each Platon-plasmid contig it greps the general
  contamination-screen BLAST hits for the word "plasmid" in the subject title.
  That check is a weak, non-authoritative signal (both it and Platon can be
  fooled by the same mobile element on a chromosomal contig), so it rides along
  as one visible column, never as the decision.
- geNomad (rule genomad_end_to_end) classifies contigs with a completely
  different, gene-content/marker-based method, so its errors are less likely to
  share Platon's specific confound. Genuine Platon/geNomad agreement is stronger
  evidence than Platon agreeing with a re-purposed screening BLAST.

This script reads both tools' outputs for ONE sample and writes a per-contig
concordance TSV. Confidence is driven by Platon/geNomad AGREEMENT:

    both tools call plasmid          -> agreement=both        confidence=high
    only Platon calls plasmid        -> agreement=platon_only  confidence=medium
    only geNomad calls plasmid       -> agreement=genomad_only confidence=medium
    Platon says chromosome, geNomad  -> agreement=conflict     confidence=low
      says plasmid (a real clash)                              (flagged, kept)

Disagreements are FLAGGED, never discarded. The file itself is the audit trail
(the agreement / confidence / platon_blast_hit columns are the per-contig
reason), matching BacFlux's existing contig_taxonomy_decisions.tsv convention.

Inputs (all produced upstream on the SAME contigs_final.fasta, so contig IDs
match directly — the front-end trims FASTA headers to their first token):
- Platon per-plasmid table  <prefix>.tsv                -> plasmid-called contigs + RDS
- Platon chromosome FASTA    <prefix>.chromosome.fasta   -> chromosome-called contig IDs
- verified_plasmids.txt      (the kept v1 BLAST-text check)
- geNomad plasmid summary    <prefix>_plasmid_summary.tsv -> geNomad-called contigs + score/fdr

Output: {sample}_plasmid_concordance.tsv, one row per plasmid CANDIDATE contig
(the union of the two tools' plasmid calls).

This is a standalone CLI (argparse), stdlib only, and is exercised by
workflow/scripts/tests/test_plasmid_concordance.py without any tool or database.
Called by the plasmid_concordance rule in workflow/rules/shared/60_plasmid.smk.
"""

import argparse
import csv
import os
import sys


# The columns of the deliverable TSV, in the exact order they are written. Kept
# as one list so the header and every row use the same field names and order.
OUTPUT_COLUMNS = [
    "sample",
    "contig",
    "platon_call",       # plasmid | chromosome | not_called
    "platon_rds",        # Platon replicon-distribution score, or NA
    "platon_blast_hit",  # hit | no_hit | NA  (the kept v1 check — supplementary)
    "genomad_call",      # plasmid | absent
    "genomad_score",     # geNomad plasmid_score, or NA
    "genomad_fdr",       # geNomad fdr, often NA (see note in parse_genomad_plasmids)
    "agreement",         # both | platon_only | genomad_only | conflict | undetermined
    "confidence",        # high | medium | low
]

# The two fixed sentences the plasmid_search rule writes into verified_plasmids.txt.
# We match on these exact trailing phrases to recover the per-contig BLAST-text
# state. Kept as named constants so the parser and the producing rule can be
# cross-checked at a glance (rule text lives in 60_plasmid.smk).
VERIFIED_SUFFIX_HIT = " is a plasmid."
VERIFIED_SUFFIX_NO_HIT = " was not verified by BLAST search."

# What the plasmid_search rule writes into verified_plasmids.txt when Platon
# itself exited non-zero (60_plasmid.smk): "{sample}: Platon exited with status
# {rc}; see the log." That line is the ONLY reliable signal that Platon crashed.
#
# Note what is deliberately NOT treated as a failure: an empty Platon directory.
# A completed, exit-0 Platon run legitimately writes almost nothing when every
# contig is longer than its 500 kb size filter - which is the NORMAL case for a
# closed genome, and true of several samples in this project's own validation
# set. Reading "no output files" as "Platon failed" would therefore mislabel
# exactly the assemblies we most want to be right about.
PLATON_CRASH_MARKER = "Platon exited with status"


def first_token(header_line):
    """Return the contig ID from a FASTA header line.

    A FASTA header like ">contig_1 length=1234" refers to the contig as
    "contig_1" everywhere else in the pipeline. Splitting on whitespace and
    taking the first token gives the same ID Platon and geNomad use, so all
    three sources join on identical strings.
    """
    without_marker = header_line.lstrip(">").strip()
    if not without_marker:
        return ""
    return without_marker.split()[0]


def _read_table(path):
    """Read a tab-separated file into (header_list, list_of_row_lists).

    Returns ([], []) for a missing or empty file so an isolate with no plasmids
    (Platon or geNomad writing only a header, or nothing) degrades gracefully
    instead of crashing. Any real content is parsed with the csv module so
    embedded quoting is handled the same way every tool writes it.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row]
    if not rows:
        return [], []
    header = [cell.strip() for cell in rows[0]]
    data_rows = rows[1:]
    return header, data_rows


def _column_index(header, name, source):
    """Find one required column by name, or fail with a clear message.

    Matching by NAME (not by fixed position) keeps the parsers robust if a tool
    reorders or adds columns between versions. A missing required column means
    the upstream output format changed and the join can no longer be trusted, so
    we stop loudly rather than silently produce wrong tiers.
    """
    if name not in header:
        raise ValueError(
            f"Expected column '{name}' not found in {source} header: {header}"
        )
    return header.index(name)


def parse_platon_tsv(path):
    """Read Platon's per-plasmid table into {contig_id: rds}.

    Platon writes one row PER PLASMID-CALLED contig (chromosome contigs are not
    listed here — they go to the .chromosome.fasta instead). We keep each
    plasmid contig's ID and its RDS (replicon-distribution score) for the report.
    The RDS is kept as its raw string; it is shown for review, not thresholded
    here.

    Input:  <prefix>.tsv from the plasmid_search rule.
    Output: dict mapping plasmid contig ID -> RDS string. Empty when Platon found
            no plasmids.
    """
    header, data_rows = _read_table(path)
    if not header:
        return {}
    id_index = _column_index(header, "ID", "Platon .tsv")
    rds_index = _column_index(header, "RDS", "Platon .tsv")

    platon_plasmids = {}
    for row in data_rows:
        # Guard against a short/ragged row rather than indexing past its end.
        # Log it, don't drop it silently: a malformed row that skips a real
        # plasmid contig would change that contig's concordance tier invisibly.
        if len(row) <= max(id_index, rds_index):
            sys.stderr.write(f"[plasmid_concordance] WARNING: skipping ragged Platon .tsv row: {row!r}\n")
            continue
        contig_id = row[id_index].strip()
        rds_value = row[rds_index].strip()
        if not contig_id:
            continue
        platon_plasmids[contig_id] = rds_value if rds_value else "NA"
    return platon_plasmids


def parse_fasta_ids(path):
    """Read a FASTA file and return the set of contig IDs in it.

    Used on Platon's <prefix>.chromosome.fasta to learn which contigs Platon
    actively called CHROMOSOME. That set lets the concordance distinguish a real
    Platon/geNomad clash (Platon chromosome vs geNomad plasmid -> "conflict")
    from a contig Platon simply never scored as a plasmid ("not_called").

    Input:  a FASTA path (may be empty).
    Output: set of contig IDs (first header token each).
    """
    chromosome_ids = set()
    if not path or not os.path.exists(path):
        return chromosome_ids
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                contig_id = first_token(line)
                if contig_id:
                    chromosome_ids.add(contig_id)
    return chromosome_ids


def platon_run_failed(verified_plasmids_path):
    """Did Platon actually CRASH for this sample? -> True/False.

    Input:  verified_plasmids.txt, written by the plasmid_search rule, which is
            already an input to this script.
    Does:   look for the one line that rule writes when Platon exits non-zero.
    Output: True only for a real crash.

    Why this signal and not "the Platon directory looks empty": an exit-0 Platon
    run writes almost nothing when every contig exceeds its 500 kb size filter,
    which is the ordinary outcome for a closed genome. Treating that as a failure
    would flag the best assemblies in the set as unassessed. A crash is a
    different thing, and the rule already records it explicitly.
    """
    if not verified_plasmids_path or not os.path.exists(verified_plasmids_path):
        return False
    with open(verified_plasmids_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if PLATON_CRASH_MARKER in line:
                return True
    return False


def parse_verified_plasmids(path):
    """Read the kept v1 BLAST-text check into {contig_id: "hit"|"no_hit"}.

    verified_plasmids.txt is written by the plasmid_search rule with one line per
    Platon-plasmid contig:
        "{sample}: {contig} is a plasmid."                    -> hit
        "{sample}: {contig} was not verified by BLAST search." -> no_hit
    plus a single "Platon found no plasmid ..." line when there were none, which
    carries no per-contig state and is skipped.

    Sample names cannot contain ':' (rejected by BacFlux's BAD_CHARS), so the
    text before the first ": " is always the sample prefix and never part of the
    contig ID. Contig IDs are single tokens (no spaces), so stripping the fixed
    trailing phrase leaves exactly the ID.

    Input:  verified_plasmids.txt.
    Output: dict mapping contig ID -> "hit" | "no_hit". Only Platon-plasmid
            contigs appear (the check is never run for other contigs).
    """
    blast_states = {}
    if not path or not os.path.exists(path):
        return blast_states
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or ": " not in line:
                continue
            # Everything after the first ": " is "{contig} <fixed phrase>".
            after_sample = line.split(": ", 1)[1]
            if after_sample.endswith(VERIFIED_SUFFIX_HIT):
                contig_id = after_sample[: -len(VERIFIED_SUFFIX_HIT)].strip()
                state = "hit"
            elif after_sample.endswith(VERIFIED_SUFFIX_NO_HIT):
                contig_id = after_sample[: -len(VERIFIED_SUFFIX_NO_HIT)].strip()
                state = "no_hit"
            else:
                # e.g. the "Platon found no plasmid ..." line — no per-contig state.
                continue
            if contig_id:
                blast_states[contig_id] = state
    return blast_states


def parse_genomad_plasmids(path):
    """Read geNomad's plasmid summary into {seq_name: (plasmid_score, fdr)}.

    geNomad lists one row per contig it calls a plasmid. Plasmid seq_names are
    whole-contig IDs (geNomad appends a "|provirus_..." suffix only to VIRUS
    fragments, never to plasmids), so seq_name joins directly to the Platon IDs.

    The fdr column is only populated when geNomad ran with score calibration; by
    default it is "NA" (or the column may be absent entirely). We therefore treat
    fdr as OPTIONAL and never filter on it — it is informational only.

    Input:  <prefix>_plasmid_summary.tsv from the genomad_end_to_end rule.
    Output: dict mapping plasmid contig ID -> (plasmid_score, fdr) as strings.
            Empty when geNomad found no plasmids.
    """
    header, data_rows = _read_table(path)
    if not header:
        return {}
    seq_index = _column_index(header, "seq_name", "geNomad plasmid summary")
    score_index = _column_index(header, "plasmid_score", "geNomad plasmid summary")
    # fdr is optional (see docstring): use it only if geNomad emitted the column.
    fdr_index = header.index("fdr") if "fdr" in header else None

    genomad_plasmids = {}
    for row in data_rows:
        if len(row) <= max(seq_index, score_index):
            sys.stderr.write(f"[plasmid_concordance] WARNING: skipping ragged geNomad summary row: {row!r}\n")
            continue
        seq_name = row[seq_index].strip()
        if not seq_name:
            continue
        score_value = row[score_index].strip() or "NA"
        if fdr_index is not None and len(row) > fdr_index:
            fdr_value = row[fdr_index].strip() or "NA"
        else:
            fdr_value = "NA"
        genomad_plasmids[seq_name] = (score_value, fdr_value)
    return genomad_plasmids


def classify(platon_call, genomad_call):
    """Turn a (Platon, geNomad) call pair into (agreement, confidence).

    These are pure rules — the whole point of D9. Confidence follows how much the
    two independent tools agree:
      - both call plasmid                 -> strongest evidence          (high)
      - one calls plasmid, other silent   -> single-tool evidence        (medium)
      - Platon chromosome vs geNomad plasmid -> a real clash, kept & flagged (low)

    Input:  platon_call in {plasmid, chromosome, not_called},
            genomad_call in {plasmid, absent}.
    Output: (agreement, confidence) strings for the report row.

    LIMITATION (D9 open issue — verify/upgrade on the first real geNomad run):
    the geNomad input here is only its plasmid_summary (POSITIVE plasmid calls),
    so genomad_call="absent" conflates "geNomad actively called it chromosome/
    virus" with "geNomad never scored it". The FORWARD conflict (Platon chromosome
    vs geNomad plasmid) IS flagged low; the REVERSE clash (Platon plasmid vs
    geNomad-NOT-plasmid) is currently reported "platon_only"/medium, not
    "conflict"/low. Reading geNomad's per-contig aggregated_classification.tsv
    would make this symmetric — deferred until geNomad can be run to confirm that
    file's exact format. Nothing is discarded either way; only the tier label of
    that reverse case is conservative.
    """
    # Platon crashed: there is no second opinion to agree or disagree with. Saying
    # "genomad_only / medium" here would claim a two-tool comparison that never
    # happened, so the failure is named instead.
    if platon_call == "not_assessed":
        return "platon_unavailable", "low"

    if platon_call == "plasmid" and genomad_call == "plasmid":
        return "both", "high"
    if platon_call == "plasmid" and genomad_call == "absent":
        return "platon_only", "medium"
    if platon_call == "not_called" and genomad_call == "plasmid":
        return "genomad_only", "medium"
    if platon_call == "chromosome" and genomad_call == "plasmid":
        return "conflict", "low"
    # The row set is the union of the two tools' plasmid calls, so no other
    # combination should reach here. Return an explicit sentinel rather than
    # crash, so an unexpected pairing is visible in the table, not fatal.
    return "undetermined", "low"


def build_rows(sample, platon_plasmids, chromosome_ids, blast_states, genomad_plasmids,
               platon_assessed=True):
    """Assemble one output row per plasmid CANDIDATE contig.

    The candidate set is the UNION of Platon-plasmid and geNomad-plasmid contigs.
    Contigs both tools agree are chromosome are intentionally left out — they are
    not plasmid candidates, matching v1's verified_plasmids.txt scope.

    Input:  the four parsed structures plus the sample name.
    Output: a list of dicts keyed by OUTPUT_COLUMNS, sorted by contig ID so the
            table is deterministic across runs.
    """
    plasmid_candidate_ids = set(platon_plasmids) | set(genomad_plasmids)

    rows = []
    for contig in sorted(plasmid_candidate_ids):
        # --- Platon side: plasmid, chromosome, or simply never called ---
        if not platon_assessed:
            # Platon crashed for this sample, so it has no opinion on any contig.
            # Calling this "not_called" would read as "Platon looked and passed
            # over it", which is exactly the false two-tool claim being fixed.
            platon_call = "not_assessed"
            platon_rds = "NA"
        elif contig in platon_plasmids:
            platon_call = "plasmid"
            platon_rds = platon_plasmids[contig]
        elif contig in chromosome_ids:
            platon_call = "chromosome"
            platon_rds = "NA"
        else:
            platon_call = "not_called"
            platon_rds = "NA"

        # The kept v1 BLAST-text signal — only ever set for Platon-plasmid contigs.
        platon_blast_hit = blast_states.get(contig, "NA")

        # --- geNomad side: plasmid or absent ---
        if contig in genomad_plasmids:
            genomad_call = "plasmid"
            genomad_score, genomad_fdr = genomad_plasmids[contig]
        else:
            genomad_call = "absent"
            genomad_score = "NA"
            genomad_fdr = "NA"

        # --- concordance tier from the two independent calls ---
        agreement, confidence = classify(platon_call, genomad_call)

        rows.append({
            "sample": sample,
            "contig": contig,
            "platon_call": platon_call,
            "platon_rds": platon_rds,
            "platon_blast_hit": platon_blast_hit,
            "genomad_call": genomad_call,
            "genomad_score": genomad_score,
            "genomad_fdr": genomad_fdr,
            "agreement": agreement,
            "confidence": confidence,
        })
    return rows


def write_rows(path, rows):
    """Write the concordance rows as a tab-separated table with a header.

    Input:  output path + the list of row dicts from build_rows.
    Output: the TSV file consumed by the user / report. Always has the header,
            even when there are zero plasmid candidates (an explicit "no
            candidates" table is more useful than a missing file).
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    # Arguments are supplied by the plasmid_concordance rule in 60_plasmid.smk.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True,
                        help="Sample name, written into the first column.")
    parser.add_argument("--platon-tsv", required=True,
                        help="Platon <prefix>.tsv (plasmid-called contigs + RDS).")
    parser.add_argument("--platon-chromosome", required=True,
                        help="Platon <prefix>.chromosome.fasta (chromosome contig IDs).")
    parser.add_argument("--verified-plasmids", required=True,
                        help="verified_plasmids.txt (kept v1 BLAST-text check).")
    parser.add_argument("--genomad-plasmid-summary", required=True,
                        help="geNomad <prefix>_plasmid_summary.tsv (score + fdr).")
    parser.add_argument("--output", required=True,
                        help="Destination concordance TSV.")
    args = parser.parse_args()

    # Load each source into a tidy Python structure.
    platon_plasmids = parse_platon_tsv(args.platon_tsv)
    chromosome_ids = parse_fasta_ids(args.platon_chromosome)
    blast_states = parse_verified_plasmids(args.verified_plasmids)
    platon_assessed = not platon_run_failed(args.verified_plasmids)
    genomad_plasmids = parse_genomad_plasmids(args.genomad_plasmid_summary)

    # Join-key audit — D9's whole agreement mechanism rests on byte-identical
    # contig IDs across the three sources. If both tools called something but NONE
    # of the IDs overlap, the join has almost certainly broken (a header not
    # trimmed to one token upstream, a geNomad seq_name suffix, a Platon rename),
    # and every genuine "both/high" would silently split into two single-tool
    # rows — a plausible-looking all-medium table. Warn loudly rather than hide it.
    platon_ids = set(platon_plasmids) | set(chromosome_ids)
    genomad_ids = set(genomad_plasmids)
    matched = len(genomad_ids & platon_ids)
    if genomad_ids and platon_ids and matched == 0:
        sys.stderr.write(
            "[plasmid_concordance] WARNING: geNomad called "
            f"{len(genomad_ids)} plasmid(s) and Platon classified {len(platon_ids)} "
            "contig(s), but NONE of their contig IDs match. The two tools' IDs "
            "likely disagree, so the concordance tiers are UNRELIABLE for this "
            "sample. Check that upstream FASTA headers are trimmed to one token.\n"
        )

    # Join them into the confidence-tiered table and write it.
    if not platon_assessed:
        sys.stderr.write(
            "[plasmid_concordance] WARNING: Platon exited non-zero for sample "
            f"{args.sample} (see verified_plasmids.txt and the Platon log). This is "
            "a FAILED two-tool comparison, NOT a 'Platon found no plasmids' result: "
            "every row is reported platon_call=not_assessed / "
            "agreement=platon_unavailable / confidence=low, and geNomad's calls "
            "stand alone and unconfirmed.\n"
        )

    rows = build_rows(
        args.sample, platon_plasmids, chromosome_ids, blast_states, genomad_plasmids,
        platon_assessed=platon_assessed,
    )
    write_rows(args.output, rows)

    # A short, honest summary for the log (never used for filtering). The match
    # count surfaces the join-key health in the run log even when it is non-zero.
    both = sum(1 for row in rows if row["agreement"] == "both")
    conflicts = sum(1 for row in rows if row["agreement"] == "conflict")
    print(
        f"Sample {args.sample}: {len(rows)} plasmid candidate contig(s)"
        f"{'' if platon_assessed else ' [PLATON FAILED - single-tool evidence only]'}; "
        f"{both} high-confidence (both tools), {conflicts} flagged conflict(s); "
        f"{matched} geNomad-plasmid ID(s) matched a Platon contig ID."
    )


if __name__ == "__main__":
    main()
