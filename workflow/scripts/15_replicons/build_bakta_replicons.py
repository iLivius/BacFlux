#!/usr/bin/env python3
"""Build the Bakta --replicons table for a long-read assembly (nanopore/hybrid).

Why Bakta has to be told which contigs are circular
---------------------------------------------------
Bakta annotates every sequence as a LINEAR CONTIG unless it is told otherwise.
For a long-read assembly that is wrong and it costs real genes: Pyrodigal (the
gene caller Bakta uses) branches on topology, and on a sequence marked circular
it is allowed to call a gene that runs across the origin. On a closed bacterial
chromosome that is typically a handful of genes at position 1 — including, quite
often, dnaA itself, which is exactly the gene dnaapler rotated to the start.

Long-read assemblies are the only ones where we KNOW the topology: Flye reports,
per contig, whether it closed it into a circle. So in nanopore and hybrid mode we
hand Bakta a small table saying "this contig is circular".

What the table says, and what it deliberately does not
------------------------------------------------------
The table has five fields per contig (see write_replicons for the exact format).
Two of them carry information:

  topology  ALWAYS from Flye's assembly_info.txt "circ." column: Y -> circular,
            anything else -> linear. This is a measurement, not an inference.

  type      chromosome | plasmid | contig. Taken from dnaapler ALONE, and only
            when its hit is strong (see classify_type). Everything else falls
            back to the neutral value "contig", which is what Bakta would have
            assumed anyway. Platon is deliberately NOT consulted: Platon runs at
            stage 06 and annotation at stage 04, so requiring its opinion would
            invert the pipeline and serialise annotation behind plasmid calling.

The asymmetry is intentional. A false negative (a real chromosome called
"contig") costs nothing beyond today's behaviour. A false positive bakes a wrong
replicon type into an INSDC-shaped annotation record that people will read as
fact. So we prefer false negatives.

One pin-sensitive detail — re-check on every Bakta version bump
---------------------------------------------------------------
In Bakta 1.12.0's parser (bakta/utils.py) the line that would force a
type="contig" row back to linear is written as a COMPARISON, not an assignment
(`topology == TOPOLOGY_LINEAR`), so it does nothing. That is what lets us keep
the conservative type="contig" AND still get topology="circular" honoured. If a
future Bakta fixes that typo, every "contig" row would silently become linear and
the whole benefit of this table would disappear without an error. The rule that
calls this script (rule build_replicons, shared/15_replicons.smk) carries the
same note.

What it reads, all for one sample
---------------------------------
  --contigs           FINAL_CONTIGS, the assembly Bakta will annotate. This is
                      the JOIN ANCHOR: one output row per record in this file,
                      in file order.
  --flye-info         Flye's assembly_info.txt  -> topology
  --dnaapler-summary  dnaapler's {sample}_all_reorientation_summary.tsv -> type

All three join on the FIRST WHITESPACE TOKEN of the contig name, which is the ID
Bakta itself uses (BacFlux always passes --keep-contig-headers, and the long-read
front ends already trim headers to one token with FASTA_HEAD_CMD).

What it writes
--------------
  --out-replicons  the 5-column, header-less TSV Bakta reads
  --out-audit      one row per contig with the numbers behind every call, so a
                   marginal decision is visible without re-running anything
                   (the project convention: every filtering decision gets an
                   audit TSV with a reason column)

Stdlib only; runs in the environment Snakemake was launched from, like
select_contigs_by_taxonomy.py and plasmid_concordance.py. Exercised by
scripts/15_replicons/test_build_bakta_replicons.py with no tools or databases.
"""

import argparse
import csv
import os
import sys


# Which dnaapler start-gene marker implies which Bakta replicon type.
# dnaapler's `all` mode searches for four markers at once and writes the winner
# into the Gene_Reoriented column:
#   dnaA     the bacterial chromosomal replication initiator
#   cog1474  the archaeal/Cdc6-Orc1 equivalent — same meaning, different domain
#   repA     a plasmid replication initiator
#   terL     a phage large terminase — Bakta has NO phage replicon type, so terL
#            is deliberately left as the neutral "contig". It is still recorded
#            in the audit file, where it doubles as a free cross-check against
#            the prophage stage (07.phages).
MARKER_TO_TYPE = {
    "dnaA": "chromosome",
    "cog1474": "chromosome",
    "repA": "plasmid",
}

# Markers dnaapler can report that we intentionally do not translate into a Bakta
# type. Kept separate from the status strings below so the audit reason can say
# WHICH of the two situations happened.
MARKERS_NOT_TYPED = ("terL", "custom")

# The audit table's columns, in the order they are written.
AUDIT_COLUMNS = [
    "sample",
    "contig",
    "topology",          # circular | linear
    "topology_source",   # flye | default
    "type",              # chromosome | plasmid | contig
    "marker",            # dnaapler Gene_Reoriented value, or NA
    "marker_coverage",   # dnaapler Coverage value, or NA
    "marker_identity",   # dnaapler Identity_Percentage value, or NA
    "reason",            # "<topology reason>;<type reason>" — see build_rows
]


# ── Read the three inputs: contig IDs, Flye topology, dnaapler markers ───────

def first_token(text):
    """Return the first whitespace-separated token of a string, without '>'.

    Every join in this script (and in Bakta itself) keys on this token. dnaapler
    writes the FULL FASTA description into its Contig column, so this is not
    optional there — see parse_dnaapler_summary.
    """
    cleaned = text.lstrip(">").strip()
    if not cleaned:
        return ""
    return cleaned.split()[0]


def read_contig_ids(path):
    """Read FINAL_CONTIGS and return its contig IDs, in file order.

    Order matters only for readability of the two output files; Bakta itself
    looks rows up by ID. Returns [] for a missing or empty FASTA, which the
    caller reports loudly (an empty replicons file makes Bakta hard-exit).
    """
    contig_ids = []
    if not path or not os.path.exists(path):
        return contig_ids
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                contig_id = first_token(line)
                if contig_id:
                    contig_ids.append(contig_id)
    return contig_ids


def parse_flye_info(path):
    """Read Flye's assembly_info.txt and return {contig_id: "circular"|"linear"}.

    The real header on a Flye 2.9 run is:
        #seq_name  length  cov.  circ.  repeat  mult.  alt_group  graph_path
    We resolve the "circ." column BY NAME rather than by position. That is a
    deliberate improvement over the blind positional awk that builds the dnaapler
    ignore list: if a future Flye release reorders its columns, a positional
    parser would mislabel every replicon silently, whereas this one falls back to
    the historical index 3 and says so in the log.

    Returns {} (with a warning) when the file is missing or unreadable — the
    caller then writes an all-linear table rather than crashing.
    """
    topology_by_contig = {}
    if not path or not os.path.exists(path):
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: Flye assembly_info not found: {path!r}. "
            "Every contig will be reported as linear.\n"
        )
        return topology_by_contig

    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: Flye assembly_info is empty: {path!r}. "
            "Every contig will be reported as linear.\n"
        )
        return topology_by_contig

    # Flye always writes a '#'-prefixed header line. If it is there, use it to
    # find the circularity column by name; if it is not, fall back to the
    # historical position and start reading from the very first line.
    circ_index = 3
    data_lines = lines
    if lines[0].startswith("#"):
        header = lines[0].lstrip("#").split("\t")
        data_lines = lines[1:]
        if "circ." in header:
            circ_index = header.index("circ.")
        else:
            sys.stderr.write(
                "[build_bakta_replicons] WARNING: no 'circ.' column in the Flye header "
                f"({header}); falling back to column 4 (index 3).\n"
            )
    else:
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: {path!r} has no '#' header line; "
            "assuming the historical Flye column order and reading from line 1.\n"
        )

    for line in data_lines:
        fields = line.split("\t")
        if len(fields) <= circ_index:
            sys.stderr.write(
                f"[build_bakta_replicons] WARNING: skipping short Flye row: {line!r}\n"
            )
            continue
        contig_id = first_token(fields[0])
        is_circular = fields[circ_index].strip().upper() == "Y"
        topology_by_contig[contig_id] = "circular" if is_circular else "linear"

    return topology_by_contig


def parse_dnaapler_summary(path):
    """Read dnaapler's reorientation summary into {contig_id: marker_record}.

    The file dnaapler writes ({sample}_all_reorientation_summary.tsv) has these
    columns, of which we use four:
        Contig  Gene_Reoriented  Start  Strand  Top_Hit  Top_Hit_Length
        Covered_Length  Coverage  Identical_AAs  Identity_Percentage
        Overlapping_Contig_End

    Two things about this file are easy to get wrong:

    1. The Contig column holds the FULL FASTA description, not the ID (dnaapler
       stores record.description). We therefore key on its first token, which is
       what every other file in the pipeline uses.
    2. When dnaapler cannot reorient a contig it fills EVERY column with the same
       status string — Contig_ignored, No_MMseqs2_hits, Contig_already_reoriented
       or autocomplete_method_<something>. So Coverage and Identity_Percentage
       are not always numbers; classify_type treats a non-numeric value as "not
       strong", never as zero.

    Returns {} (with a warning) when the file is missing, so the caller degrades
    to an all-"contig" table instead of crashing.
    """
    marker_by_contig = {}
    if not path or not os.path.exists(path):
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: dnaapler summary not found: {path!r}. "
            "No contig will be typed as chromosome or plasmid.\n"
        )
        return marker_by_contig

    with open(path, "r", encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    if len(rows) < 2:
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: dnaapler summary has no data rows: {path!r}.\n"
        )
        return marker_by_contig

    header = [cell.strip() for cell in rows[0]]
    required = ("Contig", "Gene_Reoriented", "Coverage", "Identity_Percentage")
    missing = [name for name in required if name not in header]
    if missing:
        sys.stderr.write(
            "[build_bakta_replicons] WARNING: dnaapler summary is missing column(s) "
            f"{missing}; header was {header}. No contig will be typed.\n"
        )
        return marker_by_contig

    contig_index = header.index("Contig")
    marker_index = header.index("Gene_Reoriented")
    coverage_index = header.index("Coverage")
    identity_index = header.index("Identity_Percentage")
    last_index = max(contig_index, marker_index, coverage_index, identity_index)

    for row in rows[1:]:
        if len(row) <= last_index:
            sys.stderr.write(
                f"[build_bakta_replicons] WARNING: skipping ragged dnaapler row: {row!r}\n"
            )
            continue
        contig_id = first_token(row[contig_index])
        if not contig_id:
            continue
        marker_by_contig[contig_id] = {
            "marker": row[marker_index].strip(),
            "coverage": row[coverage_index].strip(),
            "identity": row[identity_index].strip(),
        }

    return marker_by_contig


# ── Turn a dnaapler marker into a Bakta replicon type ────────────────────────

def _as_float(text):
    """Return text as a float, or None when it is not a number.

    dnaapler writes status strings ("Contig_ignored", …) into the numeric columns
    when it could not reorient a contig, so "not a number" is a normal case here,
    not an error.
    """
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def classify_type(marker_record, min_coverage, min_identity):
    """Decide a contig's replicon type from its dnaapler marker.

    Returns (type, reason) where type is chromosome | plasmid | contig.

    The "strong hit" rule: Coverage >= min_coverage AND Identity >= min_identity,
    with the defaults 80.0 and 40.0.

      * Coverage is the alignment length as a percentage of the reference
        protein, so it can exceed 100 when the alignment carries gaps. On real
        data a genuine chromosomal dnaA covers essentially the whole reference
        (observed 100.6 and 101.3 on this project's isolates) while marginal
        repA/terL hits cover 14–61%. Coverage is what separates them.
      * Identity does NOT separate them on its own — an 18%-coverage terminase
        fragment scored 57.5% identity while a real DnaA scored 68.3%. The 40%
        floor is only there to exclude the protein-alignment "twilight zone"
        (~20–35%), where alignment-based homology inference stops being reliable.
        It is deliberately NOT set near the observed 68–79%: dnaapler's reference
        proteins routinely come from another genus, and a high cutoff would
        demote real chromosomes.

    KNOWN WEAKNESS: one threshold pair covers both DnaA and RepA, and plasmid
    replication initiators are far more diverse than chromosomal ones, so the
    repA leg will mostly fall back to "contig". That is the intended failure
    direction, and the audit file carries the raw numbers so it stays visible.

    Known blind spot: for a contig that was ALREADY dnaA-first, dnaapler writes
    "Contig_already_reoriented" into every column and the marker identity is
    lost, so no type can be assigned. Rare on a first pass over raw Flye output;
    the norm if dnaapler is ever re-run on an already-oriented assembly.
    """
    if marker_record is None:
        return "contig", "no_dnaapler_row"

    marker = marker_record["marker"]

    if marker in MARKERS_NOT_TYPED:
        return "contig", f"marker_not_typed({marker})"

    if marker not in MARKER_TO_TYPE:
        # Everything left here is one of dnaapler's status strings
        # (Contig_ignored, No_MMseqs2_hits, Contig_already_reoriented,
        # autocomplete_method_*) or an unexpected value from a newer version.
        return "contig", f"no_marker({marker})"

    coverage = _as_float(marker_record["coverage"])
    identity = _as_float(marker_record["identity"])
    if coverage is None or identity is None:
        return "contig", "marker_below_threshold"
    if coverage >= min_coverage and identity >= min_identity:
        return MARKER_TO_TYPE[marker], f"strong_{marker}"
    return "contig", "marker_below_threshold"


# ── Join the three sources into one row per contig ───────────────────────────

def build_rows(sample, contig_ids, topology_by_contig, marker_by_contig,
               min_coverage, min_identity):
    """Join the three sources into one audit record per contig.

    FINAL_CONTIGS drives the loop, so there is exactly one row per sequence Bakta
    will see, in the same order. The other two tables are LEFT-joined onto it: a
    contig missing from either simply gets the neutral default.

    The `reason` column holds two semicolon-separated tokens — the topology
    decision then the type decision — so one column answers both questions:
        topology: circular_flye | linear_flye | no_flye_row
        type:     strong_dnaA | strong_repA | strong_cog1474
                  | marker_below_threshold | marker_not_typed(terL|custom)
                  | no_marker(<dnaapler status string>) | no_dnaapler_row
    """
    rows = []
    for contig_id in contig_ids:
        topology = topology_by_contig.get(contig_id)
        if topology is None:
            topology = "linear"
            topology_source = "default"
            topology_reason = "no_flye_row"
        else:
            topology_source = "flye"
            topology_reason = f"{topology}_flye"

        marker_record = marker_by_contig.get(contig_id)
        replicon_type, type_reason = classify_type(
            marker_record, min_coverage, min_identity
        )

        rows.append({
            "sample": sample,
            "contig": contig_id,
            "topology": topology,
            "topology_source": topology_source,
            "type": replicon_type,
            "marker": marker_record["marker"] if marker_record else "NA",
            "marker_coverage": marker_record["coverage"] if marker_record else "NA",
            "marker_identity": marker_record["identity"] if marker_record else "NA",
            "reason": f"{topology_reason};{type_reason}",
        })
    return rows


# ── Write the Bakta table and the audit trail ────────────────────────────────

def write_replicons(path, rows):
    """Write the table Bakta reads: 5 tab-separated fields, NO header line.

        1 original sequence id  the contig ID as it appears in FINAL_CONTIGS
        2 new sequence id       always "-" (BacFlux passes --keep-contig-headers,
                                so Bakta ignores renames; "-" states "no rename")
        3 type                  chromosome | plasmid | contig
        4 topology              circular | linear
        5 name                  always "-"

    Example row:  contig_1<TAB>-<TAB>chromosome<TAB>circular<TAB>-

    The field count is not negotiable. Bakta unpacks every row into exactly five
    variables inside a bare try/except and, on any mismatch, exits with
    "ERROR: wrong replicon table file format!" and nothing else. There is also no
    header line: Bakta would try to unpack it as data.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            writer.writerow([row["contig"], "-", row["type"], row["topology"], "-"])


def write_audit(path, rows):
    """Write the per-contig audit table: header plus one row per contig.

    Nothing downstream reads it — it exists so a marginal call can be checked
    later, which is why the raw dnaapler Coverage and Identity_Percentage go in
    unrounded rather than only the verdict they produced. Requested by name in
    _frontend_targets_for() (00_common.smk) or it would never be built.
    """
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(AUDIT_COLUMNS)
        for row in rows:
            writer.writerow([row[column] for column in AUDIT_COLUMNS])


# ── Stop the run when the contig IDs did not line up ─────────────────────────

def report_join_health(contig_ids, topology_by_contig, marker_by_contig):
    """Log how well the two tables joined onto the assembly, and warn if badly.

    This is the detector for the one failure mode that is otherwise SILENT: if
    Medaka (nanopore mode) renamed the contigs, FINAL_CONTIGS would no longer
    share IDs with the Flye and dnaapler tables, every lookup would miss, and the
    script would happily write a perfectly-formatted all-linear, all-"contig"
    table — i.e. exactly today's behaviour, with no error anywhere.

    Returns (flye_matches, dnaapler_matches) for the caller's summary line.
    """
    contig_set = set(contig_ids)
    flye_matches = len(contig_set & set(topology_by_contig))
    dnaapler_matches = len(contig_set & set(marker_by_contig))

    # A TOTAL join failure is FATAL, not a warning. If we only warned, the script
    # would still write a perfectly well-formed all-linear / all-"contig" table,
    # Bakta would run happily with a table that says nothing, the annotation would
    # be silently identical to passing no table at all, and the only trace would be
    # one line in a log nobody reads during a batch run. Zero overlap is never a
    # legitimate biological outcome — it always means an upstream step renamed the
    # contigs — so it must stop the run and say so. (Partial overlap IS legitimate:
    # dnaapler only reports contigs it could reorient, so that stays a warning.)
    if contig_ids and topology_by_contig and flye_matches == 0:
        sys.exit(
            "[build_bakta_replicons] ERROR: the Flye table lists "
            f"{len(topology_by_contig)} contig(s) but NONE of their IDs match the "
            f"{len(contig_ids)} contig(s) in the assembly. Topology would be lost "
            "entirely and Bakta would silently annotate as if no replicon table "
            "existed. A later step has almost certainly renamed the contigs — check "
            "that FASTA headers are still trimmed to their first token."
        )
    if contig_ids and marker_by_contig and dnaapler_matches == 0:
        sys.exit(
            "[build_bakta_replicons] ERROR: the dnaapler summary lists "
            f"{len(marker_by_contig)} contig(s) but NONE of their IDs match the "
            f"{len(contig_ids)} contig(s) in the assembly. Replicon types would be "
            "lost entirely. A later step has almost certainly renamed the contigs."
        )
    return flye_matches, dnaapler_matches


# ── Build both tables for one sample ─────────────────────────────────────────

def main():
    # Every flag is filled in by rule build_replicons in shared/15_replicons.smk,
    # which is defined only when HAS_LONG_READS — so this script never runs in
    # illumina or contigs mode, where the topology is genuinely unknown.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True,
                        help="Sample name, written into the audit table.")
    parser.add_argument("--contigs", required=True,
                        help="FINAL_CONTIGS — the assembly Bakta will annotate "
                             "(the join anchor: one output row per record).")
    parser.add_argument("--flye-info", required=True,
                        help="Flye assembly_info.txt (source of the topology).")
    parser.add_argument("--dnaapler-summary", required=True,
                        help="dnaapler {sample}_all_reorientation_summary.tsv "
                             "(source of the replicon type).")
    parser.add_argument("--out-replicons", required=True,
                        help="Destination 5-column TSV for Bakta --replicons.")
    parser.add_argument("--out-audit", required=True,
                        help="Destination per-contig audit TSV.")
    parser.add_argument("--min-coverage", type=float, default=80.0,
                        help="Minimum dnaapler Coverage for a strong marker hit "
                             "(default: 80.0).")
    parser.add_argument("--min-identity", type=float, default=40.0,
                        help="Minimum dnaapler Identity_Percentage for a strong "
                             "marker hit (default: 40.0).")
    args = parser.parse_args()

    contig_ids = read_contig_ids(args.contigs)
    if not contig_ids:
        # Bakta hard-exits on an EMPTY replicons file, so say so loudly. rule
        # annotation in shared/40_annotation.smk also tests the file with [ -s ]
        # before passing --replicons, so this is guarded on both sides.
        sys.stderr.write(
            f"[build_bakta_replicons] WARNING: no FASTA records in {args.contigs!r}. "
            "Writing an EMPTY replicon table; the annotation rule will skip "
            "--replicons rather than pass an empty file to Bakta.\n"
        )

    topology_by_contig = parse_flye_info(args.flye_info)
    marker_by_contig = parse_dnaapler_summary(args.dnaapler_summary)

    flye_matches, dnaapler_matches = report_join_health(
        contig_ids, topology_by_contig, marker_by_contig
    )

    rows = build_rows(
        args.sample, contig_ids, topology_by_contig, marker_by_contig,
        args.min_coverage, args.min_identity,
    )

    write_replicons(args.out_replicons, rows)
    write_audit(args.out_audit, rows)

    # One summary line into logs/build_replicons_{sample}.log. The two match
    # counts are the part worth reading: report_join_health only kills the run on
    # ZERO overlap, so a partial rename shows up here as a count well below the
    # contig total and nowhere else.
    circular = sum(1 for row in rows if row["topology"] == "circular")
    chromosomes = sum(1 for row in rows if row["type"] == "chromosome")
    plasmids = sum(1 for row in rows if row["type"] == "plasmid")
    print(
        f"Sample {args.sample}: {len(rows)} contig(s); {circular} circular; "
        f"{chromosomes} typed chromosome, {plasmids} typed plasmid, "
        f"{len(rows) - chromosomes - plasmids} left as neutral 'contig'. "
        f"Join: {flye_matches} contig(s) matched the Flye table, "
        f"{dnaapler_matches} matched the dnaapler summary."
    )


if __name__ == "__main__":
    main()
