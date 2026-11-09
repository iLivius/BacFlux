#!/usr/bin/env python3
"""Write the per-strain Excel workbooks: one for AMR genes, one for traits.

What the workbooks are for
--------------------------
Someone screening a strain collection for safety wants one file per strain
that answers, for every antimicrobial-resistance (AMR) gene in the genome:
which gene it is, how close the match is, and whether it could move to another
bacterium. That last question depends on where the gene sits: on a plasmid,
inside an integrative and conjugative or mobilisable element (ICE or IME), or
inside a phage or prophage. With the trait profile on (traits.run), a second
workbook lists the strain's predicted plant-beneficial and bioremediation
traits, with the genes, caveats, suggested assays and references behind each.

Both workbooks are a view of results BacFlux has already produced. This
script makes no new gene call. It reads the tables of earlier stages, groups
the database hits of one gene into one row, and adds each gene's genomic
context by comparing coordinates. The rows it starts from are kept on their
own sheets, so every answer can be traced to the tool that gave it.

The two commands, one per rule in shared/90_report.smk
------------------------------------------------------
    final_reports.py amr    ... --out 10.workbooks/SAMPLE1_amr.xlsx
        rule amr_workbook, in every run
    final_reports.py traits ... --out 10.workbooks/SAMPLE1_traits.xlsx
        rule traits_workbook, only with traits.run: true
The inputs of each command are listed in make_parser() at the end of this file.

Four answers, not two
---------------------
Each context question (plasmid? ICE or IME? phage?) is answered with one of:
    yes           the evidence was found (the rule for each question is in
                  plasmid_answer, ice_ime_answer and phage_answer below);
    no            the analysis ran and found no such evidence;
    uncertain     the tools disagree, or the evidence stops short of a yes;
    not_assessed  the analysis that could answer did not run: the opt-in
                  mobilome module was off, or Platon failed for this strain.
So an analysis that did not run is never written as "no".

Where it runs
-------------
In envs/report.yaml (Python 3.12 and XlsxWriter). XlsxWriter is imported only
inside write_workbook(), so the reading and deciding code can be tested in any
Python 3: see test_final_reports.py next to this file.
"""

import argparse
import collections
import hashlib
import json
import math
import os
import re
import sys

# Two scripts of the mobilome module are imported rather than written a second
# time here, so that the workbook and 08.mobilome can never disagree:
#   platon_replicons.py makes the chromosome-or-plasmid call for every contig,
#     from Platon, geNomad's second opinion when it ran, and the contig lengths;
#   colocalise.py holds the mobilome's rules for when a gene lies "inside" an
#     element and how far "near" reaches.
# Both use only Python's standard library. They sit in the sibling folder
# 80_mobilome/, which is added to the places Python looks for modules.
# Python would otherwise save a compiled copy of each imported script (a .pyc
# file) in 80_mobilome/__pycache__/, inside the workflow folder; this switches
# that off, so a run writes nothing into the installed workflow.
MOBILOME_SCRIPTS = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "80_mobilome"))
sys.dont_write_bytecode = True
sys.path.insert(0, MOBILOME_SCRIPTS)
import colocalise        # noqa: E402  (imported after the folder is added)
import platon_replicons  # noqa: E402


# ── The four answers, and the version of the layout ──────────────────────────

YES = "yes"
NO = "no"
UNCERTAIN = "uncertain"
NOT_ASSESSED = "not_assessed"
ANSWERS = (YES, NO, UNCERTAIN, NOT_ASSESSED)

# Written to Run_QC. Raise it when a sheet or a column changes meaning, so that
# workbooks made by two BacFlux versions can be told apart.
WORKBOOK_LAYOUT = "1.0"

# Cell values that mean "nothing here" in the tables this script reads. "-" is
# deliberately not one of them: it is the minus strand in ABRicate's STRAND
# column and in Bakta's, and must be copied as it is.
MISSING_VALUES = {"", "NA", "N/A", "."}


# ── The AMR databases ────────────────────────────────────────────────────────

# The three ABRicate databases the AMR workbook reports, each with the name
# ABRicate writes in its own DATABASE column. BacFlux runs ABRicate against more
# databases (rule amr_contigs in shared/50_amr.smk); these three are the
# collections of resistance genes, which EFSA's guidance on whole-genome
# sequencing asks to search with at least two maintained databases.
ABRICATE_DATABASES = (
    ("CARD", "card"),
    ("ResFinder", "resfinder"),
    ("ARG-ANNOT", "argannot"),
)
ABRICATE_COLUMNS = (
    "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE", "COVERAGE_MAP",
    "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT",
    "RESISTANCE",
)

# Two coordinate intervals are taken to be the same gene when at least 80 % of
# the shorter one lies inside the longer one. Database references of one gene
# often differ by a few bases at either end, so exact coordinates would split
# one gene into several rows; neighbouring genes rarely overlap by more than a
# few bases, so 80 % keeps them apart. Used to group database hits into loci and
# to match an AMRFinderPlus row of the mobilome module to a locus.
SAME_GENE_FRACTION = 0.80

# The AMRFinderPlus rows matched to the database loci. AMR is a resistance gene;
# STRESS covers biocide and metal resistance, such as the quaternary-ammonium
# efflux gene qacEdelta1, which CARD, ResFinder and ARG-ANNOT list as resistance
# genes too. Point mutations are left out although their element type is AMR:
# their method (amrfinder_method POINTP, POINTX, ...) says AMRFinderPlus found a
# mutation in a gene every strain of the species carries, such as gyrA, and
# their coordinates are the whole gene's, so a database hit on that gene would
# otherwise take the mutation's name and tier. VIRULENCE rows and point
# mutations stay on the Mobility_evidence sheet without a locus.
MATCHED_AMRFINDER_TYPES = {"AMR", "STRESS"}
POINT_MUTATION_METHOD_PREFIX = "POINT"


# ── The mobile elements ──────────────────────────────────────────────────────

# The mobilome module's own rules, read from colocalise.py so the two cannot
# drift apart: a gene lies inside an element when at least 90 % of it lies
# within the element's interval, and an element is "near" when it lies within
# 5,000 bp of the gene. The near window exists because the interval written for
# an ICE or IME is the span of its machinery genes (integrase, relaxase, mating
# apparatus), not the element's resolved ends, so a gene just outside that span
# may still be carried by the element.
INSIDE_FRACTION = colocalise.MIN_AMR_FRACTION_INSIDE_ELEMENT
NEAR_ELEMENT_BP = colocalise.IS_ADJACENT_MAX_BP

# How the element types of the mobilome's ICE table (written by conjscan_to_ice.py)
# answer the question "is this gene inside an ICE or IME?":
#   ice, ime, aice      yes. An ICE moves itself by conjugation; an IME is moved
#                       by a conjugative element in the same cell; an AICE (the
#                       actinomycete kind of ICE) conjugates through its own
#                       FtsK-like translocase.
#   conjugative_region  uncertain. Conjugation machinery whose element is not
#                       resolved: no integrase was found beside it, or it lies on
#                       a plasmid, which makes it a conjugative plasmid region.
#   genomic_island      no, with the island named in the evidence. An integrase
#                       without a relaxase: the island can sit in the genome but
#                       cannot be moved by conjugation on its own.
# A gene near (not inside) an ice, ime or aice is uncertain, for the reason
# given at NEAR_ELEMENT_BP.
ICE_IME_TYPES = ("ice", "ime", "aice")
CONJUGATIVE_REGION = "conjugative_region"
GENOMIC_ISLAND = "genomic_island"

# The line rule plasmid_search (shared/60_plasmid.smk) writes into
# verified_plasmids.txt when Platon exits with an error. It is the only reliable
# sign that Platon failed; a nearly empty Platon folder is normal for a closed
# genome (see PLATON_MAX_CONTIG_LENGTH in platon_replicons.py).
PLATON_CRASH_MARKER = "Platon exited with status"


# ── Excel's limits and the typing of copied values ───────────────────────────

# An Excel cell holds at most 32,767 characters. A longer text is cut with a
# visible mark, never silently.
EXCEL_CELL_LIMIT = 32767
CUT_MARK = " [... cut: an Excel cell holds at most 32,767 characters; see the source table on Run_QC]"

# The sheets that copy an upstream table cell for cell (the mobility table, the
# trait gene table, the two audits) write a value as a number when it is a plain
# number, so it sorts and filters as one. Columns that hold names stay text
# even when a name happens to be all digits (a contig called "1", an EC number).
PLAIN_NUMBER = re.compile(r"-?\d+(\.\d+)?([eE][-+]?\d+)?")
NAME_COLUMNS = {
    "sample", "contig", "replicon_id", "locus_tag", "gene", "ko", "ec",
    "trait_id", "traits", "mge_id", "amr_gene", "ecofolddb_target",
}


# ── Reading tables ───────────────────────────────────────────────────────────

def read_tsv(path, required=(), skip_comment_lines=False):
    """Read a tab-separated table into (header, rows).

    Takes in: the path of a table written by an earlier BacFlux step or tool,
              and the columns this script needs from it.
    Gives back: the header as a list, and one dict per row (column -> text).
              Each row also carries the number of its line in the file under
              "__line__": error messages name it, Database_hits shows it as
              source_line, and copied_sheet() uses it to find a row again.
    A table with a header and no rows is valid: that is how ABRicate and the
    virus callers report a search that found nothing. A row with the wrong
    number of fields stops the script with the file and line named, because
    padding it would put values under the wrong column. Blank lines are
    skipped, and skip_comment_lines also skips the "#" lines at the top of the
    two trait tables; line numbers still count every line of the file.
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            # (line number in the file, text) for every line that is not blank.
            lines = [(number, line.rstrip("\r\n"))
                     for number, line in enumerate(handle, start=1) if line.strip()]
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if skip_comment_lines:
        while lines and lines[0][1].startswith("#"):
            lines.pop(0)
    if not lines:
        raise ValueError(f"{path}: the file is empty; a header line was expected")

    header = lines[0][1].split("\t")
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")

    rows = []
    for line_number, line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"{path}: line {line_number} has {len(fields)} fields, "
                             f"but the header has {len(header)}")
        row = dict(zip(header, fields))
        row["__line__"] = line_number
        rows.append(row)
    return header, rows


def is_missing(value):
    """True for an empty cell or one of the "nothing here" words above."""
    return str(value if value is not None else "").strip().upper() in MISSING_VALUES


def joined(values):
    """Join the distinct, non-missing values with "; ", in first-seen order.

    Used where several database hits or elements give one cell of a summary
    row. Returns "NA" when nothing is left.
    """
    kept = []
    for value in values:
        text = str(value).strip()
        if not is_missing(text) and text not in kept:
            kept.append(text)
    return "; ".join(kept) if kept else "NA"


def whole_number(value, path, column, line):
    """Read a coordinate or a count, naming the file and line when it is not one."""
    try:
        return int(str(value).strip())
    except ValueError:
        raise ValueError(f"{path}: line {line}, column {column} is not a whole "
                         f"number: {value!r}") from None


def decimal_number(value, path, column, line):
    """Read a percentage or a score, naming the file and line when it is not one."""
    try:
        return float(str(value).strip())
    except ValueError:
        raise ValueError(f"{path}: line {line}, column {column} is not a number: "
                         f"{value!r}") from None


def typed(value, column):
    """Return a copied cell as a number when it is a plain number, else as text.

    Used only on the sheets that copy an upstream table unchanged; see
    PLAIN_NUMBER and NAME_COLUMNS for which values become numbers.
    """
    text = str(value).strip()
    if column in NAME_COLUMNS or not PLAIN_NUMBER.fullmatch(text):
        return value
    if "." in text or "e" in text or "E" in text:
        return float(text)
    return int(text)


def shown_path(path, output_root):
    """The path to show in a workbook, without the computer's own folders.

    A file inside the BacFlux output folder is shown relative to it (for
    example 05.amr/abricate/SAMPLE1/card.tsv); any other file, such as the
    trait tables that ship with BacFlux, by its file name alone. A workbook can
    then be passed on without revealing where the analysis ran.
    """
    absolute = os.path.abspath(path)
    root = os.path.abspath(output_root)
    if os.path.commonpath([absolute, root]) == root:
        return os.path.relpath(absolute, root)
    return os.path.basename(absolute)


def check_sample(rows, sample, path):
    """Stop if a table belongs to another sample than the one being reported.

    The tables of the mobilome module and the trait profile carry the sample
    name on every row; a mix-up between two strains' files is caught here.
    """
    for row in rows:
        if row.get("sample") != sample:
            raise ValueError(f"{path}: line {row['__line__']} is for sample "
                             f"{row.get('sample')!r}, not {sample!r}")


# ── Which strain, and which genome ───────────────────────────────────────────

def read_primary_assembly(path):
    """Return the row of the genome BacFlux delivered for this sample.

    Takes in: 02.assembly/{sample}/eval/{sample}_qc_genomes.tsv, which lists the
              assemblies the QC tools were run on. In hybrid mode it holds the
              Illumina draft (role comparator) beside the final long-read genome
              (role primary); in the other modes only the primary.
    Gives back: the primary row (bin_id, technology, source_path, ...). Its
              bin_id is the name GTDB-Tk used for the genome.
    """
    _header, rows = read_tsv(path, ("bin_id", "role", "technology", "source_path"))
    primary = [row for row in rows if row["role"].strip().lower() == "primary"]
    if len(primary) != 1:
        raise ValueError(f"{path}: expected one row with role primary, found {len(primary)}")
    return primary[0]


def read_taxonomy(taxonomy_dir, genome_name):
    """Return the GTDB-Tk classification of the delivered genome.

    Takes in: the GTDB-Tk folder of the sample (03.taxonomy/{sample}/, from rule
              gtdbtk) and the genome's name there (the primary bin_id above).
    Gives back: the summary table's path, the full lineage, its lowest named
              rank (the species when GTDB-Tk assigned one) and the method.
    GTDB-Tk writes gtdbtk.bac120.summary.tsv (gtdbtk.ar53.summary.tsv for an
    archaeon) in classify/ and a copy at the top of its folder. The classify/
    one is read, as the gtdbtk rule documents it; the top one only when
    classify/ has none. Exactly one is expected in the folder used.
    """
    candidates = []
    for folder in (os.path.join(taxonomy_dir, "classify"), taxonomy_dir):
        if os.path.isdir(folder):
            candidates = sorted(
                os.path.join(folder, name) for name in os.listdir(folder)
                if name.startswith("gtdbtk.") and name.endswith(".summary.tsv"))
        if candidates:
            break
    if len(candidates) != 1:
        raise ValueError(f"{taxonomy_dir}: expected one GTDB-Tk summary table, "
                         f"found {len(candidates)}")
    path = candidates[0]
    _header, rows = read_tsv(path, ("user_genome", "classification"))
    matches = [row for row in rows if row["user_genome"] == genome_name]
    if len(matches) != 1:
        raise ValueError(f"{path}: expected one row for genome {genome_name!r}, "
                         f"found {len(matches)}")
    row = matches[0]

    # "d__Bacteria;p__...;g__Klebsiella;s__Klebsiella pneumoniae" -> the last
    # rank that carries a name. An empty rank ("s__") means GTDB-Tk stopped above.
    lineage = row["classification"].strip()
    named_ranks = []
    for rank in lineage.split(";"):
        name = rank.split("__", 1)[-1].strip()
        if name:
            named_ranks.append(name)
    return {
        "path": path,
        "lineage": lineage or "NA",
        "lowest_taxon": named_ranks[-1] if named_ranks else "unclassified",
        "method": row.get("classification_method", "NA") or "NA",
    }


# ── Coordinates ──────────────────────────────────────────────────────────────
# Every tool read here writes 1-based coordinates with both ends included, so a
# gene from base 100 to base 200 is 101 bases long.

def overlap_bp(start_a, end_a, start_b, end_b):
    """How many bases two intervals share; 0 when they do not overlap."""
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def same_gene(start_a, end_a, start_b, end_b):
    """True when at least 80 % of the shorter interval lies inside the longer.

    See SAME_GENE_FRACTION for why 80 %.
    """
    shorter = min(end_a - start_a + 1, end_b - start_b + 1)
    if shorter <= 0:
        return False
    return overlap_bp(start_a, end_a, start_b, end_b) / shorter >= SAME_GENE_FRACTION


def gene_intervals(start, end, contig_length):
    """Split a gene into the interval(s) it covers on its contig.

    Bakta writes a gene that runs across the origin of a circular contig with
    start greater than end: it covers start..contig end, then 1..end. Such a
    gene becomes two intervals, so an overlap test cannot read it backwards.
    ABRicate never does this, so its hits always give one interval.
    """
    if start <= end:
        return [(start, end)]
    return [(start, contig_length), (1, end)]


def counted_answers(rows, column):
    """One summary cell: how many rows gave each answer, e.g. "yes 2, no 31"."""
    counts = collections.Counter(row[column] for row in rows)
    parts = [f"{answer} {counts[answer]}" for answer in ANSWERS if counts[answer]]
    return ", ".join(parts) if parts else "no rows"


# ── AMR genes: database hits and loci ────────────────────────────────────────

def read_abricate(path, database, database_word, output_root):
    """Read one ABRicate result table into a list of hits.

    Takes in: 05.amr/abricate/{sample}/{card,resfinder,argannot}.tsv, from rule
              amr_contigs (shared/50_amr.smk), which runs ABRicate on the final
              contigs at 80 % identity and 70 % coverage of the reference.
    Gives back: one dict per row, with the coordinates as numbers and every
              other value as ABRicate wrote it. The strand is copied unchanged
              ("+" or "-"), and GAPS stays text, so a value such as 6/4 cannot
              be turned into a date by Excel.
    The DATABASE column is checked, so two tables passed in the wrong order are
    caught instead of being reported under each other's name.
    """
    _header, rows = read_tsv(path, ABRICATE_COLUMNS)
    hits = []
    for row_number, row in enumerate(rows, start=1):
        line = row["__line__"]
        if row["DATABASE"].strip().lower() != database_word:
            raise ValueError(f"{path}: line {line} is from database {row['DATABASE']!r}, "
                             f"but the file was given as the {database} table")
        start = whole_number(row["START"], path, "START", line)
        end = whole_number(row["END"], path, "END", line)
        hits.append({
            # The database and the row's number in ABRicate's table, which is
            # also how a hit is found again on Database_hits.
            "hit_id": f"{database}:{row_number:04d}",
            "database": database,
            "contig": row["SEQUENCE"].strip(),
            # ABRicate writes the smaller coordinate first on both strands; the
            # min/max only protects the interval tests from a reversed pair.
            "start": min(start, end),
            "end": max(start, end),
            "strand": row["STRAND"].strip(),
            "gene": row["GENE"].strip(),
            "pct_identity": decimal_number(row["%IDENTITY"], path, "%IDENTITY", line),
            "pct_coverage": decimal_number(row["%COVERAGE"], path, "%COVERAGE", line),
            "coverage": row["COVERAGE"].strip(),
            "coverage_map": row["COVERAGE_MAP"].strip(),
            "gaps": row["GAPS"].strip(),
            "reference_accession": row["ACCESSION"].strip(),
            "product": row["PRODUCT"].strip(),
            "resistance": row["RESISTANCE"].strip(),
            "source_file": shown_path(path, output_root),
            "source_line": line,
        })
    return hits


def group_hits_into_loci(hits, sample):
    """Group the database hits of one gene into one locus.

    Takes in: the hits of all three databases (read_abricate).
    Gives back: one dict per locus, sorted by contig and position, each with
              its hits under "hits"; every hit also gets its locus_id.

    The databases often name the same gene differently (aadA2 in one, ANT(3'')-Ia
    in another), so names cannot be used; coordinates are. The rule: the
    longest hit of a gene defines the locus, and a shorter hit joins it when at
    least 80 % of the shorter hit lies inside the longer one (SAME_GENE_FRACTION).
    Hits are taken longest first, so every locus is anchored by its longest hit
    and a hit is compared with anchors only. That keeps one partial hit from
    chaining two neighbouring genes into one locus. Two cases the rule leaves
    open are settled here: a hit that fits inside two anchors joins the one it
    shares more bases with (the first one made, if they share as many); and of
    two hits of equal length, the one starting first (then by database and row)
    becomes the anchor, which sets the locus's coordinates and strand.
    """
    hits_by_contig = collections.defaultdict(list)
    for hit in hits:
        hits_by_contig[hit["contig"]].append(hit)

    loci = []
    for contig in sorted(hits_by_contig):
        longest_first = sorted(
            hits_by_contig[contig],
            key=lambda hit: (-(hit["end"] - hit["start"] + 1), hit["start"], hit["hit_id"]))
        contig_loci = []
        for hit in longest_first:
            hit_length = hit["end"] - hit["start"] + 1
            home = None
            home_shared = 0
            for locus in contig_loci:
                shared = overlap_bp(hit["start"], hit["end"], locus["start"], locus["end"])
                if shared / hit_length >= SAME_GENE_FRACTION and shared > home_shared:
                    home = locus
                    home_shared = shared
            if home is None:
                # No anchor holds this hit: it is the longest hit of a new locus.
                home = {"contig": contig, "start": hit["start"], "end": hit["end"],
                        "strand": hit["strand"], "hits": []}
                contig_loci.append(home)
            home["hits"].append(hit)
        loci.extend(contig_loci)

    database_order = [name for name, _word in ABRICATE_DATABASES]
    for locus in loci:
        # The ID is built from the anchor's coordinates, so it moves if a new
        # database version changes the longest hit.
        locus["locus_id"] = f"{sample}|{locus['contig']}:{locus['start']}-{locus['end']}"
        for hit in locus["hits"]:
            hit["locus_id"] = locus["locus_id"]
        found_in = {hit["database"] for hit in locus["hits"]}
        locus["databases"] = "; ".join(name for name in database_order if name in found_in)
    loci.sort(key=lambda locus: (locus["contig"], locus["start"], locus["end"]))
    return loci


def read_abricate_databases(path):
    """Read ABRicate's version and its database copies -> {database word: row}.

    Takes in: 05.amr/abricate/abricate_databases.tsv, from rule
              abricate_databases (shared/50_amr.smk): `abricate --list` (DATABASE,
              SEQUENCES, DBTYPE, DATE) with ABRicate's version as a last column.
    Gives back: one row per database, keyed by the name ABRicate uses (card,
              resfinder, argannot, ...). DATE is the date ABRicate gives for its
              copy of the database, the closest thing to a database version it
              records. A result is reported with the database it came from and
              that database's version or date, so the workbook shows both.
    """
    _header, rows = read_tsv(path, ("DATABASE", "SEQUENCES", "DATE", "ABRICATE_VERSION"))
    return {row["DATABASE"].strip().lower(): row for row in rows}


def database_note(database, database_word, listing):
    """One phrase for a database's copy, e.g. "CARD, dated 2026-Apr-3, 6052 sequences"."""
    row = listing.get(database_word)
    if row is None:
        return f"{database}, not in ABRicate's list"
    return f"{database}, dated {row['DATE']}, {row['SEQUENCES']} sequences"


# ── Plasmid context ──────────────────────────────────────────────────────────

def platon_crashed(platon_dir):
    """True when Platon exited with an error for this sample (see PLATON_CRASH_MARKER)."""
    path = os.path.join(platon_dir, "verified_plasmids.txt")
    if not os.path.exists(path):
        return False
    with open(path, encoding="utf-8", errors="replace") as handle:
        return any(PLATON_CRASH_MARKER in line for line in handle)


def read_replicon_calls(platon_dir, prefix, contig_lengths, concordance, replicon_table):
    """Give every contig its chromosome-or-plasmid call -> {contig: row}.

    Takes in: Platon's folder (06.plasmids/{sample}/platon, rule plasmid_search)
              and the file prefix Platon used; the contig lengths of the final
              assembly; geNomad's concordance table when geNomad ran (rule
              plasmid_concordance), or None; and the mobilome module's replicon
              table when that module ran (rule mobilome_replicons), or None.
    Gives back: one row per contig with the columns of the replicon table:
              replicon (chromosome, plasmid or unknown), plasmid_mobility,
              mobility_evidence and replicon_call_source.

    With the mobilome module on, its table is read as it is, so the workbook
    agrees with the mobility tiers that module computed from it. Without it,
    the same call is made here by the same function, platon_replicons.build_rows.
    That function reads Platon's size filter correctly: Platon skips every
    contig of 500 kb or more, so such a contig is decided by geNomad's call or
    by its length (2 Mb or more: chromosome; between: unknown).
    """
    if replicon_table:
        _header, rows = read_tsv(replicon_table, (
            "contig", "replicon", "plasmid_mobility", "mobility_evidence",
            "replicon_call_source"))
    else:
        genomad_calls = platon_replicons.read_genomad_concordance(concordance or "")
        rows = platon_replicons.build_rows(
            platon_dir, prefix, sorted(contig_lengths), contig_lengths=contig_lengths,
            genomad_calls=genomad_calls)
    return {row["contig"]: row for row in rows}


def plasmid_answer(contig, replicon_calls, crashed):
    """Answer "is this gene on a plasmid?" for a gene on one contig.

    The rule, in order:
      yes           the contig is called plasmid: by Platon, by geNomad alone,
                    or by both (plasmid_call_source says which);
      not_assessed  Platon failed for this sample, so nothing short of a
                    geNomad plasmid call can be said;
      uncertain     the contig is called chromosome, by Platon or from a length
                    of 2 Mb or more, but geNomad calls it a plasmid (call source
                    "conflict"; a contig that long is likely a chromid);
      no            the contig is called chromosome: by Platon, or inferred from
                    a length of 2 Mb or more for a contig Platon skipped;
      uncertain     anything else: a contig neither tool classified, such as one
                    between 500 kb and 2 Mb that could be a megaplasmid, or one
                    under Platon's 1 kb minimum.
    """
    if contig not in replicon_calls:
        raise ValueError(f"contig {contig!r} has a gene but no replicon call; the AMR "
                         "tables and the assembly do not belong together")
    call = replicon_calls[contig]
    replicon = call["replicon"].strip().lower()
    source = call["replicon_call_source"].strip().lower()
    evidence = call["mobility_evidence"]

    if replicon == "plasmid":
        answer = YES
    elif crashed:
        answer = NOT_ASSESSED
        evidence = ("Platon exited with an error for this sample (see "
                    "verified_plasmids.txt in its folder), so the contig was not classified")
    elif source == "conflict":
        answer = UNCERTAIN
    elif replicon == "chromosome":
        answer = NO
    else:
        answer = UNCERTAIN
    return {
        "replicon": replicon or "unknown",
        "plasmid": answer,
        "plasmid_call_source": source or "NA",
        "plasmid_mobility": call["plasmid_mobility"],
        "plasmid_evidence": evidence,
    }


# ── Phage context ────────────────────────────────────────────────────────────

def read_viral_regions(path, caller, contig_lengths):
    """Read the virus caller's regions into a list, in contig coordinates.

    Takes in: the caller's region table, from the folder of rule
              genomad_end_to_end or viral_identification_virsorter2
              (shared/70_phage.smk), and the contig lengths.
      geNomad    {prefix}_virus_summary.tsv. A prophage is named
                 contig|provirus_START_END and its coordinates are "START-END";
                 a whole contig called a virus has coordinates NA, and covers
                 the contig from base 1 to its last base.
      VirSorter2 final-viral-boundary.tsv. A part of a contig is named
                 contig||N_partial in seqname_new, with its ends on the contig in
                 trim_bp_start and trim_bp_end. Any other name (contig||full,
                 contig||lt2gene) is a whole contig: rule
                 viral_identification_virsorter2 runs VirSorter2 with
                 --keep-original-seq, which delivers the untrimmed contig for
                 these, so the region, and CheckV's coordinates, run from base 1
                 to the contig's last base, whatever the trim columns say.
    Gives back: one dict per region: region_id (the name CheckV also uses),
              caller, contig, start, end, length, score, group and taxonomy.
    """
    regions = []
    if caller == "genomad":
        _header, rows = read_tsv(path, ("seq_name", "coordinates", "topology",
                                        "virus_score", "taxonomy"))
        for row in rows:
            name = row["seq_name"].strip()
            contig = name.split("|provirus_", 1)[0]
            if is_missing(row["coordinates"]):
                if contig not in contig_lengths:
                    raise ValueError(f"{path}: virus {name!r} is on a contig the "
                                     "assembly does not have")
                start, end = 1, contig_lengths[contig]
            else:
                first, last = row["coordinates"].split("-", 1)
                start = whole_number(first, path, "coordinates", row["__line__"])
                end = whole_number(last, path, "coordinates", row["__line__"])
            regions.append({
                "region_id": name, "caller": "geNomad", "contig": contig,
                "start": start, "end": end, "length": end - start + 1,
                "score": decimal_number(row["virus_score"], path, "virus_score", row["__line__"]),
                "group": row["topology"], "taxonomy": row["taxonomy"],
            })
    else:
        _header, rows = read_tsv(path, ("seqname", "trim_bp_start", "trim_bp_end",
                                        "seqname_new", "final_max_score",
                                        "final_max_score_group"))
        for row in rows:
            line = row["__line__"]
            name = row["seqname_new"].strip()
            contig = row["seqname"].strip()
            if name.endswith("_partial"):
                start = whole_number(row["trim_bp_start"], path, "trim_bp_start", line)
                end = whole_number(row["trim_bp_end"], path, "trim_bp_end", line)
            else:
                if contig not in contig_lengths:
                    raise ValueError(f"{path}: virus {name!r} is on a contig the "
                                     "assembly does not have")
                start, end = 1, contig_lengths[contig]
            regions.append({
                "region_id": name, "caller": "VirSorter2",
                "contig": contig, "start": start, "end": end,
                "length": end - start + 1,
                "score": decimal_number(row["final_max_score"], path, "final_max_score", line),
                "group": row["final_max_score_group"], "taxonomy": "NA",
            })
    regions.sort(key=lambda region: (region["contig"], region["start"], region["end"]))
    return regions


def checkv_parts(row, path):
    """Split one row of CheckV's contamination.tsv into viral and host parts.

    CheckV looks at every region the caller delivered and, for a prophage,
    finds where the viral genes end and the host genes begin. The row gives
    the parts in region_types ("host,viral,host") and their ends in
    region_coords_bp ("1-3603,3604-48893,48894-163418"), counted from the
    first base of the region, not of the contig. A row with provirus "No"
    has no host part: the whole region is viral.
    Gives back: (viral_parts, host_parts), lists of (start, end) in region
    coordinates.
    """
    region_length = whole_number(row["contig_length"], path, "contig_length", row["__line__"])
    if row["provirus"].strip().lower() != "yes" or is_missing(row["region_types"]):
        return [(1, region_length)], []
    types = row["region_types"].split(",")
    spans = row["region_coords_bp"].split(",")
    if len(types) != len(spans):
        raise ValueError(f"{path}: line {row['__line__']} has {len(types)} region types "
                         f"but {len(spans)} coordinate pairs")
    viral_parts, host_parts = [], []
    for part_type, span in zip(types, spans):
        first, last = span.split("-", 1)
        part = (whole_number(first, path, "region_coords_bp", row["__line__"]),
                whole_number(last, path, "region_coords_bp", row["__line__"]))
        if part_type.strip() == "viral":
            viral_parts.append(part)
        else:
            host_parts.append(part)
    return viral_parts, host_parts


def attach_checkv(regions, quality_path, contamination_path):
    """Add CheckV's assessment and the viral and host parts to every region.

    Takes in: the regions (read_viral_regions) and two tables from the folder of
              rule viral_quality: quality_summary.tsv (quality class,
              completeness, contamination) and contamination.tsv (the parts).
              CheckV's contig_id is the caller's region name, so they join on it.
              With nothing to assess, the rule writes a header-only
              quality_summary.tsv and no contamination.tsv. When CheckV fails,
              the rule also writes a header-only quality_summary.tsv, but
              contamination.tsv, CheckV's first output, may already be there.
    Changes: each region gets "checkv" (its quality row, or None) and
             "viral_parts" / "host_parts" in contig coordinates. A region counts
             as assessed only when it has a row in both tables, so parts left by
             a CheckV run that failed later are not used: without a quality row,
             the parts are None. A part at position x of the region is base
             region_start + x - 1 of the contig.
    Stops when CheckV's length for a region differs from the region's length:
    the two would then not describe the same stretch of DNA.
    Gives back: the number of quality rows read, for Run_QC.
    """
    _header, quality_rows = read_tsv(quality_path, (
        "contig_id", "provirus", "checkv_quality", "miuvig_quality", "completeness",
        "completeness_method", "contamination", "warnings"))
    quality_by_id = {row["contig_id"]: row for row in quality_rows}

    parts_by_id = {}
    checkv_length = {}
    if os.path.exists(contamination_path):
        _header, rows = read_tsv(contamination_path, (
            "contig_id", "contig_length", "provirus", "region_types", "region_coords_bp"))
        for row in rows:
            parts_by_id[row["contig_id"]] = checkv_parts(row, contamination_path)
            checkv_length[row["contig_id"]] = whole_number(
                row["contig_length"], contamination_path, "contig_length", row["__line__"])

    for region in regions:
        region["checkv"] = quality_by_id.get(region["region_id"])
        parts = parts_by_id.get(region["region_id"])
        if parts is None or region["checkv"] is None:
            region["viral_parts"] = None
            region["host_parts"] = None
            continue
        if checkv_length[region["region_id"]] != region["length"]:
            raise ValueError(f"{contamination_path}: CheckV gives {region['region_id']} "
                             f"{checkv_length[region['region_id']]} bp, but the caller's "
                             f"region is {region['length']} bp")
        offset = region["start"] - 1
        region["viral_parts"] = [(offset + first, offset + last) for first, last in parts[0]]
        region["host_parts"] = [(offset + first, offset + last) for first, last in parts[1]]
    return len(quality_rows)


def parts_text(parts):
    """Write a list of contig intervals for a cell, e.g. "3604-48893; 50001-60000"."""
    if not parts:
        return "NA"
    return "; ".join(f"{first}-{last}" for first, last in parts)


def phage_answer(contig, intervals, regions, caller):
    """Answer "is this gene inside a phage or prophage region?".

    Takes in: the gene's contig and interval(s) (gene_intervals), the regions
              with CheckV's parts attached (attach_checkv), and the caller's name.
    The rule:
      yes        the gene overlaps a part of a region that CheckV calls viral
                 (a region without host part is viral as a whole);
      uncertain  the gene overlaps a region the caller delivered, but only the
                 part CheckV calls host DNA, or CheckV did not assess that
                 region: the two tools do not agree that the gene is viral;
      no         the gene overlaps no region the caller delivered.
    Gives back: the answer, the region names and the reason as a sentence.
    """
    viral, host_only, unassessed = [], [], []
    for region in regions:
        if region["contig"] != contig:
            continue
        touches_region = any(overlap_bp(first, last, region["start"], region["end"])
                             for first, last in intervals)
        if not touches_region:
            continue
        if region["viral_parts"] is None:
            unassessed.append(region)
        elif any(overlap_bp(first, last, part_first, part_last)
                 for first, last in intervals
                 for part_first, part_last in region["viral_parts"]):
            viral.append(region)
        else:
            host_only.append(region)

    if viral:
        reasons = [f"overlaps the part of {region['region_id']} that CheckV calls viral "
                   f"(contig bases {parts_text(region['viral_parts'])})" for region in viral]
        return YES, joined(region["region_id"] for region in viral), "; ".join(reasons)
    if host_only or unassessed:
        reasons = [f"overlaps {region['region_id']} from {caller}, but only the part "
                   f"CheckV calls host DNA (viral part: contig bases "
                   f"{parts_text(region['viral_parts'])})" for region in host_only]
        reasons += [f"overlaps {region['region_id']} from {caller}, which CheckV did "
                    "not assess" for region in unassessed]
        names = joined(region["region_id"] for region in host_only + unassessed)
        return UNCERTAIN, names, "; ".join(reasons)
    return NO, "NA", f"overlaps no virus or prophage region called by {caller}"


# ── ICE and IME context ──────────────────────────────────────────────────────

def read_ice_elements(path):
    """Read the mobilome module's ICE table into a list of elements.

    Takes in: 08.mobilome/{sample}/{sample}_ice_candidates.tsv, or its _named
              copy when ICEberg naming is on: the table rule
              amr_mge_colocalisation reads (ICE_TABLE_FOR_COLOCALISE in
              00_common.smk). One row per element, from conjscan_to_ice.py.
    Gives back: one dict per element, with its coordinates as numbers and the
              columns the workbook shows.
    """
    _header, rows = read_tsv(path, ("sample", "mge_id", "mge_name", "element_type",
                                    "mobility", "contig", "start", "end"))
    elements = []
    for row in rows:
        line = row["__line__"]
        elements.append({
            "mge_id": row["mge_id"], "mge_name": row["mge_name"],
            "element_type": row["element_type"].strip().lower(),
            "mobility": row["mobility"], "contig": row["contig"].strip(),
            "start": whole_number(row["start"], path, "start", line),
            "end": whole_number(row["end"], path, "end", line),
            "boundary_method": row.get("boundary_method", "NA"),
            "machinery_intact": row.get("machinery_intact", "NA"),
            "evidence_level": row.get("evidence_level", "NA"),
            "sample": row["sample"], "__line__": line,
        })
    return elements


def element_label(element):
    """Name an element in a sentence: type, ID, and its name when ICEberg gave one."""
    label = f"{element['element_type']} {element['mge_id']}"
    if not is_missing(element["mge_name"]):
        label += f" ({element['mge_name']})"
    return label


def ice_ime_answer(contig, start, end, elements, assessed):
    """Answer "is this gene inside an ICE or IME?" for one gene.

    Takes in: the gene's contig and ends, the elements of the mobilome's ICE
              table (read_ice_elements), and whether that module ran.
    The rule, decided by coordinates for every gene (see ICE_IME_TYPES):
      yes           the gene lies inside an ice, ime or aice;
      uncertain     the gene lies inside a conjugative region, or within
                    5,000 bp of an ice, ime or aice without lying inside it;
      no            none of these; a genomic island around the gene, or a
                    conjugative region near it, is named in the evidence;
      not_assessed  the mobilome module did not run.
    "Inside" and "near" are the mobilome's own (INSIDE_FRACTION, NEAR_ELEMENT_BP),
    and distances are the bases strictly between gene and element, as
    colocalise.py measures them.
    Gives back: the answer, the element IDs it rests on, and the reason.
    """
    if not assessed:
        return NOT_ASSESSED, "NA", "the mobilome module did not run (mobilome.run: false)"

    gene_length = end - start + 1
    inside, near = [], []
    for element in elements:
        if element["contig"] != contig:
            continue
        shared = overlap_bp(start, end, element["start"], element["end"])
        if shared / gene_length >= INSIDE_FRACTION:
            inside.append(element)
            continue
        distance = colocalise.gap_bp(start, end, element["start"], element["end"])
        if distance <= NEAR_ELEMENT_BP:
            near.append((distance, element))
    near.sort(key=lambda pair: pair[0])

    inside_ice_ime = [e for e in inside if e["element_type"] in ICE_IME_TYPES]
    if inside_ice_ime:
        reasons = [f"inside {element_label(e)}: {e['mobility']}" for e in inside_ice_ime]
        return YES, joined(e["mge_id"] for e in inside_ice_ime), "; ".join(reasons)

    inside_region = [e for e in inside if e["element_type"] == CONJUGATIVE_REGION]
    near_ice_ime = [(d, e) for d, e in near if e["element_type"] in ICE_IME_TYPES]
    if inside_region or near_ice_ime:
        reasons = [f"inside {element_label(e)}: {e['mobility']}" for e in inside_region]
        reasons += [f"{d:,} bp from {element_label(e)}; the interval of such an element "
                    "is the span of its machinery genes, not its resolved ends, so the "
                    "gene may still lie inside it" for d, e in near_ice_ime]
        elements_used = inside_region + [e for _d, e in near_ice_ime]
        return UNCERTAIN, joined(e["mge_id"] for e in elements_used), "; ".join(reasons)

    # What is left inside is a genomic island (or a type this script does not
    # know, named as it is); what is left near is a conjugative region or an
    # island. Neither makes the answer more than no, but both are named.
    reasons = []
    for element in inside:
        reason = f"inside {element_label(element)}: {element['mobility']}"
        if element["element_type"] == GENOMIC_ISLAND:
            reason += ("; an island with an integrase but no relaxase cannot be moved "
                       "by conjugation on its own")
        reasons.append(reason)
    reasons += [f"{d:,} bp from {element_label(e)}" for d, e in near]
    if not reasons:
        reasons = [f"no ICE, IME, AICE or conjugative region within {NEAR_ELEMENT_BP:,} bp"]
    elements_used = inside + [e for _d, e in near]
    return NO, joined(e["mge_id"] for e in elements_used), "; ".join(reasons)


# ── The mobilome module's mobility tiers ─────────────────────────────────────

def mobility_columns(locus, mobility_rows, assessed):
    """Add the mobilome module's verdict for the gene at this locus.

    Takes in: one locus, the rows of the mobility table
              (08.mobilome/{sample}/{sample}_amr_mobility.tsv, rule
              amr_mge_colocalisation, one row per AMRFinderPlus gene), and
              whether that module ran.
    The mobilome module rates AMRFinderPlus's genes, not ABRicate's, so the
    rows are matched to the locus by coordinates: same contig, and the same
    gene by SAME_GENE_FRACTION. Only AMR and STRESS rows are matched, and no
    point mutation (see MATCHED_AMRFINDER_TYPES). When several rows match, the
    highest tier is shown, with its label and confidence.
    Gives back: the dict of mobility columns, and the matched rows.
    """
    if not assessed:
        columns = {name: NOT_ASSESSED for name in (
            "amrfinderplus_genes", "amrfinderplus_class", "amrfinderplus_subclass",
            "mge_context", "mobility_tier", "mobility_tier_label", "mobility_confidence")}
        return columns, []

    matches = []
    for row in mobility_rows:
        if row["amr_element_type"].strip().upper() not in MATCHED_AMRFINDER_TYPES:
            continue
        if row["amrfinder_method"].strip().upper().startswith(POINT_MUTATION_METHOD_PREFIX):
            continue
        if row["contig"].strip() != locus["contig"]:
            continue
        if same_gene(locus["start"], locus["end"], row["start"], row["end"]):
            matches.append(row)
    if not matches:
        return {
            "amrfinderplus_genes": "NA", "amrfinderplus_class": "NA",
            "amrfinderplus_subclass": "NA", "mge_context": "NA", "mobility_tier": "NA",
            "mobility_tier_label": "no AMRFinderPlus AMR or STRESS gene at this position",
            "mobility_confidence": "NA",
        }, []

    def tier_number(row):
        tier = row["mobility_tier"].strip()
        return int(tier) if tier.isdigit() else 0

    best = max(matches, key=tier_number)
    best_tier = best["mobility_tier"].strip()
    return {
        "amrfinderplus_genes": joined(row["amr_gene"] for row in matches),
        "amrfinderplus_class": joined(row["amr_class"] for row in matches),
        "amrfinderplus_subclass": joined(row["amr_subclass"] for row in matches),
        # Copied as written, "none" included: it is colocalise.py's own word for
        # a gene with no mobile element around it.
        "mge_context": "; ".join(dict.fromkeys(row["mge_context"] for row in matches)),
        "mobility_tier": int(best_tier) if best_tier.isdigit() else best_tier,
        "mobility_tier_label": best["mobility_tier_label"],
        "mobility_confidence": best["confidence"],
    }, matches


def read_mobility_rows(path, sample):
    """Read the mobility table, with each gene's ends as numbers -> (header, rows)."""
    header, rows = read_tsv(path, ("sample", "contig", "amr_gene", "amr_element_type",
                                   "amr_class", "amr_subclass", "amr_start", "amr_end",
                                   "amrfinder_method", "mge_context", "mobility_tier",
                                   "mobility_tier_label", "confidence"))
    check_sample(rows, sample, path)
    for row in rows:
        row["start"] = whole_number(row["amr_start"], path, "amr_start", row["__line__"])
        row["end"] = whole_number(row["amr_end"], path, "amr_end", row["__line__"])
    return header, rows


# ── Sheets ───────────────────────────────────────────────────────────────────
# A sheet is a plain dict: its name, its columns, its rows (lists of values in
# column order) and the columns whose text should wrap. write_workbook() turns
# the list of sheets into the .xlsx file. Each column list below pairs every
# column with its meaning; the Data_dictionary sheet is made from these pairs,
# so a column and its explanation cannot drift apart.

def make_sheet(name, column_meanings, rows, wrap=()):
    """Assemble one sheet from (column, meaning) pairs and rows of dicts."""
    columns = [column for column, _meaning in column_meanings]
    return {
        "name": name,
        "columns": columns,
        "rows": [[row.get(column, "NA") for column in columns] for row in rows],
        "wrap": set(wrap),
        "meanings": column_meanings,
    }


def copied_sheet(name, path, sample, source_note, not_run_note=None, extra_first=None):
    """A sheet that copies an upstream table cell for cell.

    Takes in: the sheet name; the table's path, or None when the step that
              writes it did not run; a sentence naming the table and the rule
              that wrote it, for Data_dictionary; and, optionally, a first
              column to add (name, meaning, {line number: value}).
    Values that are plain numbers are written as numbers (typed()); everything
    else exactly as written, "none" and "-" included.
    """
    if path is None:
        return {"name": name, "columns": ["note"], "rows": [[not_run_note]], "wrap": {"note"},
                "meanings": [("note", "why this sheet is empty")]}
    header, rows = read_tsv(path)
    check_sample(rows, sample, path)
    columns = list(header)
    meanings = [("all columns", source_note)]
    sheet_rows = []
    for row in rows:
        values = [typed(row[column], column) for column in header]
        if extra_first:
            values.insert(0, extra_first[2].get(row["__line__"], "NA"))
        sheet_rows.append(values)
    if extra_first:
        columns.insert(0, extra_first[0])
        meanings.insert(0, (extra_first[0], extra_first[1]))
    return {"name": name, "columns": columns, "rows": sheet_rows, "wrap": set(),
            "meanings": meanings}


def dictionary_sheet(sheets, answer_note):
    """The Data_dictionary sheet: every column of every sheet, with its meaning."""
    rows = [{"sheet": "all", "column": "yes / no / uncertain / not_assessed",
             "meaning": answer_note}]
    for sheet in sheets:
        for column, meaning in sheet["meanings"]:
            rows.append({"sheet": sheet["name"], "column": column, "meaning": meaning})
    return make_sheet("Data_dictionary", [
        ("sheet", "the sheet the column is on"),
        ("column", "the column name"),
        ("meaning", "what the column holds"),
    ], rows, wrap=("meaning",))


ANSWER_NOTE = (
    "yes: the evidence was found. no: the analysis ran and found none. uncertain: the "
    "tools disagree, or the evidence stops short of a yes. not_assessed: the analysis "
    "that could answer did not run. The rule behind each question is given with its "
    "column below.")

SUMMARY_COLUMNS = [
    ("item", "what is summarised"),
    ("value", "the value for this strain"),
    ("note", "how to read it"),
]
RUN_QC_COLUMNS = [
    ("item", "an input, a setting or a check"),
    ("value", "its path relative to the BacFlux output folder (a file name for files "
              "outside it), or its value"),
    ("detail", "what was read from it"),
]


# ── The AMR workbook ─────────────────────────────────────────────────────────

AMR_LOCUS_COLUMNS = [
    ("sample", "BacFlux sample name, from the input file name"),
    ("gtdb_taxon", "lowest rank GTDB-Tk assigned to the genome; the full lineage is on Summary"),
    ("locus_id", "sample|contig:start-end of the longest database hit of the gene; it changes "
                 "when those coordinates change, for example with a new database version"),
    ("contig", "the contig carrying the gene"),
    ("start", "first base of the longest database hit"),
    ("end", "last base of the longest database hit"),
    ("strand", "strand of the longest database hit, + or -"),
    ("genes", "gene names given by the databases, longest hit first; databases often name "
              "one gene differently"),
    ("databases", "which of CARD, ResFinder and ARG-ANNOT found the gene"),
    ("in_resfinder", "yes when ResFinder found the gene. ResFinder collects acquired "
                     "resistance genes; a gene found only by CARD or ARG-ANNOT can be one "
                     "every strain of the species carries, such as an efflux pump, a porin "
                     "or a regulator"),
    ("database_hits", "how many database hits were grouped into this locus (a hit joins the "
                      "longest one when at least 80 % of it lies inside); each is a row on "
                      "Database_hits"),
    ("min_pct_identity", "lowest percent identity to the reference sequence among those hits"),
    ("max_pct_identity", "highest percent identity among those hits"),
    ("min_pct_coverage", "lowest percent of the reference length covered among those hits"),
    ("max_pct_coverage", "highest percent of the reference length covered among those hits"),
    ("products", "product descriptions given by the databases"),
    ("resistance", "drug classes the databases attach to the gene (ABRicate's RESISTANCE column)"),
    ("reference_accessions", "GenBank accessions, with coordinates, of the reference sequences "
                             "the hits were made to; CARD's ARO numbers are not in ABRicate's output"),
    ("replicon", "chromosome, plasmid or unknown: the replicon call for the contig, as the "
                 "mobilome module makes it (platon_replicons.py)"),
    ("plasmid", "is the gene on a plasmid? yes: the contig is called plasmid (by Platon, "
                "geNomad or both); no: called chromosome (by Platon, or 2 Mb or longer when "
                "Platon skipped it); uncertain: geNomad calls a plasmid what Platon or a "
                "length of 2 Mb or more calls chromosome, or neither "
                "classified the contig; not_assessed: Platon failed"),
    ("plasmid_call_source", "the replicon table's call source: genomad (geNomad alone called a "
                            "plasmid), both, conflict (chromosome by Platon or by a length of "
                            "2 Mb or more, plasmid by geNomad), or "
                            "platon for every other call, including a contig Platon skipped whose "
                            "call comes from its length (see plasmid_evidence)"),
    ("plasmid_mobility", "for a plasmid, conjugative, mobilisable or non-mobilisable from "
                         "Platon's counts of conjugation and mobilisation genes; unknown when "
                         "Platon did not classify the contig; NA for a chromosome"),
    ("plasmid_evidence", "the replicon call's own explanation"),
    ("ice_ime", "is the gene inside an ICE or IME? yes: inside an ICE, IME or AICE; uncertain: "
                "inside a conjugative region, or within 5,000 bp of an ICE, IME or AICE; no: "
                "neither; not_assessed: the mobilome module did not run. Inside means at "
                "least 90 % of the gene lies within the element"),
    ("ice_ime_elements", "the mobilome element(s) the answer rests on; see ICE_IME_elements"),
    ("ice_ime_evidence", "why, with the distance in bases when the gene lies outside the element"),
    ("phage", "is the gene inside a phage or prophage region? yes: it overlaps a part CheckV "
              "calls viral; uncertain: it overlaps a region of the virus caller, but only the "
              "part CheckV calls host DNA, or CheckV did not assess the region; no: it "
              "overlaps no region"),
    ("phage_regions", "the virus caller's region(s) the answer rests on; see Viral_regions"),
    ("phage_evidence", "why, with the contig bases CheckV calls viral"),
    ("amrfinderplus_genes", "AMRFinderPlus gene symbol(s) at the same position, from the "
                            "mobilome module (AMR and STRESS genes only)"),
    ("amrfinderplus_class", "AMRFinderPlus drug class"),
    ("amrfinderplus_subclass", "AMRFinderPlus drug subclass"),
    ("mge_context", "the mobilome module's word for what surrounds the gene, as colocalise.py "
                    "writes it (for example plasmid, ice, is_adjacent, none)"),
    ("mobility_tier", "the mobilome module's mobility ladder, from 1 (no mobile context found) "
                      "to 6 (predicted self-transmissible); the highest tier when several "
                      "AMRFinderPlus rows match"),
    ("mobility_tier_label", "the tier in words"),
    ("mobility_confidence", "the mobilome module's confidence in the tier: high, medium or low"),
]

DATABASE_HIT_COLUMNS = [
    ("hit_id", "database and row number of the hit in ABRicate's table"),
    ("locus_id", "the AMR_loci row the hit was grouped into"),
    ("database", "CARD, ResFinder or ARG-ANNOT"),
    ("contig", "the contig of the hit"),
    ("start", "first base of the hit on the contig"),
    ("end", "last base of the hit on the contig"),
    ("strand", "strand of the hit, + or -"),
    ("gene", "the database's name for the gene"),
    ("pct_identity", "percent identity to the reference sequence"),
    ("pct_coverage", "percent of the reference length covered"),
    ("coverage", "bases of the reference covered, and its length (ABRicate's COVERAGE)"),
    ("coverage_map", "ABRicate's picture of the covered part of the reference"),
    ("gaps", "gap openings / gap bases in the alignment (ABRicate's GAPS), kept as text"),
    ("reference_accession", "GenBank accession, with coordinates, of the reference sequence"),
    ("product", "the database's product description"),
    ("resistance", "the drug classes the database attaches to the gene"),
    ("source_file", "the ABRicate table the row was read from"),
    ("source_line", "the line of that table"),
]

ELEMENT_COLUMNS = [
    ("mge_id", "the element's ID in the mobilome module's ICE table"),
    ("mge_name", "the closest ICEberg element, when ICEberg naming is on"),
    ("element_type", "ice, ime, aice, conjugative_region or genomic_island"),
    ("mobility", "the mobilome module's description of how the element can move"),
    ("contig", "the contig of the element"),
    ("start", "first base of the element's interval: the span of its machinery genes, or "
              "the attachment sites when they were found"),
    ("end", "last base of the element's interval"),
    ("boundary_method", "how the ends were set (none: the machinery span)"),
    ("machinery_intact", "whether the element's machinery was found complete"),
    ("evidence_level", "the mobilome module's evidence level for the element"),
    ("loci_inside", "AMR_loci rows lying inside the element"),
    ("loci_near", "AMR_loci rows within 5,000 bp of the element, without lying inside it"),
]

VIRAL_REGION_COLUMNS = [
    ("region_id", "the region's name, shared by the virus caller and CheckV"),
    ("caller", "geNomad or VirSorter2"),
    ("contig", "the contig of the region"),
    ("start", "first base of the region on the contig"),
    ("end", "last base of the region on the contig"),
    ("length", "length of the region in bases"),
    ("caller_score", "the caller's virus score"),
    ("caller_group", "VirSorter2's viral group, or geNomad's topology"),
    ("taxonomy", "geNomad's virus taxonomy (NA for VirSorter2)"),
    ("checkv_provirus", "Yes when CheckV found host DNA at the ends of the region"),
    ("checkv_quality", "CheckV's quality class"),
    ("miuvig_quality", "the MIUViG quality class CheckV gives"),
    ("checkv_completeness", "CheckV's completeness estimate, in percent"),
    ("checkv_contamination", "percent of the region CheckV calls host DNA"),
    ("checkv_warnings", "CheckV's warnings"),
    ("viral_parts", "contig bases CheckV calls viral; NA when CheckV did not assess the region"),
    ("host_parts", "contig bases CheckV calls host DNA"),
    ("loci", "AMR_loci rows overlapping the region"),
]


def build_amr_workbook(args):
    """Read every input of the AMR workbook and return its sheets.

    The steps, in order:
      1. which strain and genome: the QC genome table and GTDB-Tk;
      2. the database hits of CARD, ResFinder and ARG-ANNOT, grouped into loci;
      3. the context of each locus: plasmid, ICE/IME, phage, and the mobilome
         module's mobility tier;
      4. the sheets: Summary, AMR_loci, Database_hits, ICE_IME_elements,
         Viral_regions, Mobility_evidence, Mobility_audit, Run_QC and
         Data_dictionary.
    """
    out_root = args.output_root
    mobilome_ran = args.ice_table is not None

    # 1. Which strain and genome.
    primary = read_primary_assembly(args.qc_genomes)
    taxonomy = read_taxonomy(args.taxonomy_dir, primary["bin_id"])

    # 2. The database hits, then the loci.
    hits = []
    tables = {"CARD": args.card, "ResFinder": args.resfinder, "ARG-ANNOT": args.argannot}
    for database, database_word in ABRICATE_DATABASES:
        hits.extend(read_abricate(tables[database], database, database_word, out_root))
    loci = group_hits_into_loci(hits, args.sample)
    listing = read_abricate_databases(args.abricate_databases)

    # 3. The context of each locus. Each source is read once, here.
    contig_lengths = platon_replicons.read_contig_lengths(args.contigs)
    crashed = platon_crashed(args.platon_dir)
    replicon_calls = read_replicon_calls(args.platon_dir, args.platon_prefix,
                                         contig_lengths, args.plasmid_concordance,
                                         args.replicons)
    regions = read_viral_regions(args.viral_table, args.phage_caller, contig_lengths)
    checkv_rows = attach_checkv(regions, args.checkv_quality, args.checkv_contamination)
    caller_name = "geNomad" if args.phage_caller == "genomad" else "VirSorter2"
    elements = []
    mobility_header, mobility_rows = [], []
    if mobilome_ran:
        elements = read_ice_elements(args.ice_table)
        check_sample(elements, args.sample, args.ice_table)
        mobility_header, mobility_rows = read_mobility_rows(args.mobility, args.sample)

    locus_rows = []
    loci_of_mobility_line = collections.defaultdict(list)
    for locus in loci:
        hits_here = locus["hits"]
        row = {
            "sample": args.sample, "gtdb_taxon": taxonomy["lowest_taxon"],
            "locus_id": locus["locus_id"], "contig": locus["contig"],
            "start": locus["start"], "end": locus["end"], "strand": locus["strand"],
            "genes": joined(hit["gene"] for hit in hits_here),
            "databases": locus["databases"],
            "in_resfinder": YES if "ResFinder" in locus["databases"].split("; ") else NO,
            "database_hits": len(hits_here),
            "min_pct_identity": min(hit["pct_identity"] for hit in hits_here),
            "max_pct_identity": max(hit["pct_identity"] for hit in hits_here),
            "min_pct_coverage": min(hit["pct_coverage"] for hit in hits_here),
            "max_pct_coverage": max(hit["pct_coverage"] for hit in hits_here),
            "products": joined(hit["product"] for hit in hits_here),
            "resistance": joined(hit["resistance"] for hit in hits_here),
            "reference_accessions": joined(hit["reference_accession"] for hit in hits_here),
        }
        row.update(plasmid_answer(locus["contig"], replicon_calls, crashed))
        answer, names, reason = ice_ime_answer(locus["contig"], locus["start"], locus["end"],
                                               elements, mobilome_ran)
        row.update({"ice_ime": answer, "ice_ime_elements": names, "ice_ime_evidence": reason})
        answer, names, reason = phage_answer(locus["contig"], [(locus["start"], locus["end"])],
                                             regions, caller_name)
        row.update({"phage": answer, "phage_regions": names, "phage_evidence": reason})
        mobility, matched = mobility_columns(locus, mobility_rows, mobilome_ran)
        row.update(mobility)
        for mobility_row in matched:
            loci_of_mobility_line[mobility_row["__line__"]].append(locus["locus_id"])
        locus_rows.append(row)

    # 4. The sheets.
    hit_rows = sorted(hits, key=lambda hit: (hit["contig"], hit["start"], hit["end"], hit["hit_id"]))
    element_rows = []
    for element in elements:
        inside, near = [], []
        for locus in loci:
            if locus["contig"] != element["contig"]:
                continue
            length = locus["end"] - locus["start"] + 1
            shared = overlap_bp(locus["start"], locus["end"], element["start"], element["end"])
            if shared / length >= INSIDE_FRACTION:
                inside.append(locus["locus_id"])
            elif colocalise.gap_bp(locus["start"], locus["end"],
                                   element["start"], element["end"]) <= NEAR_ELEMENT_BP:
                near.append(locus["locus_id"])
        element_rows.append(dict(element, loci_inside=joined(inside), loci_near=joined(near)))

    region_rows = []
    for region in regions:
        checkv = region["checkv"] or {}
        overlapping = [locus["locus_id"] for locus in loci
                       if locus["contig"] == region["contig"]
                       and overlap_bp(locus["start"], locus["end"], region["start"], region["end"])]
        completeness = checkv.get("completeness", "NA")
        contamination = checkv.get("contamination", "NA")
        region_rows.append({
            "region_id": region["region_id"], "caller": region["caller"],
            "contig": region["contig"], "start": region["start"], "end": region["end"],
            "length": region["length"], "caller_score": region["score"],
            "caller_group": region["group"], "taxonomy": region["taxonomy"],
            "checkv_provirus": checkv.get("provirus", "NA"),
            "checkv_quality": checkv.get("checkv_quality", "not assessed by CheckV"),
            "miuvig_quality": checkv.get("miuvig_quality", "NA"),
            "checkv_completeness": typed(completeness, "checkv_completeness"),
            "checkv_contamination": typed(contamination, "checkv_contamination"),
            # CheckV leaves the cell empty when it has no warning.
            "checkv_warnings": checkv.get("warnings") or "NA",
            "viral_parts": parts_text(region["viral_parts"]),
            "host_parts": parts_text(region["host_parts"]),
            "loci": joined(overlapping),
        })

    hits_per_database = collections.Counter(hit["database"] for hit in hits)
    hit_counts = ", ".join(f"{name} {hits_per_database[name]}" for name, _w in ABRICATE_DATABASES)
    # One ABRicate ran every database, so any row of the listing gives its version.
    abricate_version = "NA"
    if listing:
        abricate_version = list(listing.values())[0]["ABRICATE_VERSION"]
    summary_rows = [
        {"item": "Sample", "value": args.sample,
         "note": "BacFlux sample name, from the input file name"},
        {"item": "GTDB taxon", "value": taxonomy["lowest_taxon"],
         "note": "lowest rank GTDB-Tk assigned; the full lineage is below"},
        {"item": "GTDB lineage", "value": taxonomy["lineage"],
         "note": f"GTDB-Tk classification method: {taxonomy['method']}"},
        {"item": "Database hits", "value": len(hits),
         "note": f"{hit_counts}; one row each on Database_hits"},
        {"item": "AMR loci", "value": len(loci),
         "note": "the hits of one gene grouped into one row of AMR_loci"},
        {"item": "Found by ResFinder",
         "value": sum(1 for row in locus_rows if row["in_resfinder"] == YES),
         "note": "loci ResFinder found; ResFinder collects acquired resistance genes, "
                 "while CARD and ARG-ANNOT also list genes the species carries anyway"},
        {"item": "Databases",
         "value": f"ABRicate {abricate_version}",
         "note": "; ".join(database_note(name, word, listing)
                           for name, word in ABRICATE_DATABASES)
                 + "; the date is the one ABRicate gives for its copy of each database"},
        {"item": "On a plasmid", "value": counted_answers(locus_rows, "plasmid"),
         "note": "loci by answer in the plasmid column"},
        {"item": "Inside an ICE or IME", "value": counted_answers(locus_rows, "ice_ime"),
         "note": "loci by answer in the ice_ime column"},
        {"item": "Inside a phage region", "value": counted_answers(locus_rows, "phage"),
         "note": "loci by answer in the phage column"},
        {"item": "Mobilome module", "value": "run" if mobilome_ran else "not run",
         "note": ("ICE/IME answers and mobility tiers come from it" if mobilome_ran else
                  "mobilome.run was false, so the ICE/IME and mobility-tier columns read "
                  "not_assessed")},
        {"item": "How to read the answers", "value": "yes / no / uncertain / not_assessed",
         "note": ANSWER_NOTE},
        {"item": "What this is", "value": "screening evidence from the genome sequence",
         "note": "a resistance gene is not a susceptibility test result, and a mobile "
                 "context is a prediction, not proof that the gene transfers"},
    ]

    run_qc_rows = [
        {"item": "workbook layout", "value": WORKBOOK_LAYOUT,
         "detail": "raised when a sheet or column changes meaning"},
        {"item": "sample", "value": args.sample, "detail": "BacFlux sample name"},
        {"item": "primary assembly", "value": primary["bin_id"],
         "detail": f"technology {primary['technology']}; the genome every table describes"},
        {"item": "QC genome table", "value": shown_path(args.qc_genomes, out_root),
         "detail": "which assembly is the delivered one"},
        {"item": "GTDB-Tk summary", "value": shown_path(taxonomy["path"], out_root),
         "detail": "the taxonomy on Summary"},
        {"item": "final contigs", "value": shown_path(args.contigs, out_root),
         "detail": f"{len(contig_lengths)} contigs; their lengths decide contigs Platon skipped"},
    ]
    for database, database_word in ABRICATE_DATABASES:
        run_qc_rows.append({"item": f"ABRicate {database}",
                            "value": shown_path(tables[database], out_root),
                            "detail": f"{hits_per_database[database]} hits; "
                                      + database_note(database, database_word, listing)})
    run_qc_rows.append({"item": "ABRicate databases",
                        "value": shown_path(args.abricate_databases, out_root),
                        "detail": "ABRicate's version and the date of its copy of each database"})
    run_qc_rows += [
        {"item": "Platon", "value": shown_path(args.platon_dir, out_root),
         "detail": "exited with an error" if crashed else "finished"},
        {"item": "replicon calls",
         "value": shown_path(args.replicons, out_root) if args.replicons else "made by this script",
         "detail": ("the mobilome module's table" if args.replicons else
                    "platon_replicons.build_rows, as the mobilome module would")},
        {"item": "geNomad plasmid concordance",
         "value": (shown_path(args.plasmid_concordance, out_root)
                   if args.plasmid_concordance else "not used"),
         "detail": "geNomad's plasmid calls, a second opinion on Platon's (phage.caller: genomad)"},
        {"item": "virus regions", "value": shown_path(args.viral_table, out_root),
         "detail": f"{caller_name}: {len(regions)} regions"},
        {"item": "CheckV quality", "value": shown_path(args.checkv_quality, out_root),
         "detail": f"{checkv_rows} rows"},
        {"item": "CheckV parts",
         "value": (shown_path(args.checkv_contamination, out_root)
                   if os.path.exists(args.checkv_contamination) else "absent"),
         "detail": "viral and host parts of each region; absent when CheckV had nothing to "
                   "assess or failed"},
    ]
    if mobilome_ran:
        run_qc_rows += [
            {"item": "mobilome ICE table", "value": shown_path(args.ice_table, out_root),
             "detail": f"{len(elements)} elements"},
            {"item": "mobilome mobility table", "value": shown_path(args.mobility, out_root),
             "detail": f"{len(mobility_rows)} AMRFinderPlus rows"},
            {"item": "mobilome audit", "value": shown_path(args.mobility_audit, out_root),
             "detail": "copied to Mobility_audit"},
        ]
    else:
        run_qc_rows.append({"item": "mobilome module", "value": "not run",
                            "detail": "mobilome.run: false"})

    not_run = "The mobilome module did not run (mobilome.run: false), so there is no such table."
    sheets = [
        make_sheet("Summary", SUMMARY_COLUMNS, summary_rows, wrap=("note",)),
        make_sheet("AMR_loci", AMR_LOCUS_COLUMNS, locus_rows,
                   wrap=("products", "plasmid_evidence", "ice_ime_evidence", "phage_evidence")),
        make_sheet("Database_hits", DATABASE_HIT_COLUMNS, hit_rows, wrap=("product",)),
    ]
    if mobilome_ran:
        sheets.append(make_sheet("ICE_IME_elements", ELEMENT_COLUMNS, element_rows))
    else:
        sheets.append(copied_sheet("ICE_IME_elements", None, args.sample, "", not_run))
    sheets.append(make_sheet("Viral_regions", VIRAL_REGION_COLUMNS, region_rows))
    sheets.append(copied_sheet(
        "Mobility_evidence", args.mobility, args.sample,
        "copied from 08.mobilome/{sample}/{sample}_amr_mobility.tsv (rule "
        "amr_mge_colocalisation); every column is explained in docs/mobilome/output.md",
        not_run,
        extra_first=("locus_ids", "the AMR_loci row(s) this AMRFinderPlus gene was matched to; "
                                  "NA for genes no database hit matches, such as point "
                                  "mutations and virulence genes",
                     {line: joined(ids) for line, ids in loci_of_mobility_line.items()})))
    sheets.append(copied_sheet(
        "Mobility_audit", args.mobility_audit, args.sample,
        "copied from 08.mobilome/{sample}/{sample}_amr_mobility_audit.tsv (rule "
        "amr_mge_colocalisation): the reason behind each gene's tier; see docs/mobilome/output.md",
        not_run))
    sheets.append(make_sheet("Run_QC", RUN_QC_COLUMNS, run_qc_rows, wrap=("detail",)))
    sheets.append(dictionary_sheet(sheets, ANSWER_NOTE))
    return sheets, {"database hits": len(hits), "loci": len(loci)}


# ── The traits workbook ──────────────────────────────────────────────────────

RULE_COLUMNS = ("trait_id", "trait", "application", "required", "supporting",
                "colocalised", "key_genes", "caveat", "assay", "reference")
PROFILE_COLUMNS = ("sample", "trait_id", "trait", "application", "decision", "genes",
                   "missing", "supporting", "borderline", "structure_only", "caveat")
GENE_COLUMNS = ("sample", "locus_tag", "contig", "start", "stop", "strand", "gene",
                "traits", "counted")

TRAIT_PREDICTION_COLUMNS = [
    ("sample", "BacFlux sample name"),
    ("gtdb_taxon", "lowest rank GTDB-Tk assigned to the genome"),
    ("trait_id", "the trait's ID in the catalogue (trait_rules_v1.tsv)"),
    ("trait", "the trait's name"),
    ("application", "what the trait is relevant to; context traits describe the strain's "
                    "metabolism, not a benefit (see docs/analysis/traits.md)"),
    ("decision", "called: every gene of a required set counts, and the genes that must lie "
                 "together do. partial: the set is complete but its genes are broken "
                 "(pseudogenes) or not together, or at least half of a set of two or more "
                 "genes counts. not called: anything else; it does not show the trait is absent"),
    ("genes", "the genes that count, each as gene:locus_tag(evidence)"),
    ("missing", "for a partial call, the genes of the set that did not count"),
    ("supporting", "supporting genes found; they never make the call"),
    ("borderline", "genes that count from a structure hit close to the filter's limits"),
    ("structure_only", "genes that count from a structure hit alone"),
    ("caveat", "why a genomic call is not a demonstrated phenotype"),
    ("assay", "a laboratory test for the trait, from the catalogue"),
    ("reference", "the literature behind the trait's rule"),
]

TRAIT_RULE_COLUMNS = [
    ("trait_id", "the trait's ID"),
    ("trait", "the trait's name"),
    ("application", "what the trait is relevant to"),
    ("required", "the gene sets that can make the call, separated by |; + joins the genes "
                 "of one set"),
    ("supporting", "genes reported beside the call, never needed for it"),
    ("colocalised", "which genes must lie together on one contig"),
    ("key_genes", "genes without which a partial call is withheld"),
    ("assay", "a laboratory test for the trait"),
    ("reference", "the literature behind the rule"),
]

GENE_CONTEXT_COLUMNS = [
    ("replicon", "chromosome, plasmid or unknown, for the gene's contig"),
    ("plasmid", "is the gene on a plasmid? yes / no / uncertain / not_assessed, by the rule "
                "given for the AMR workbook's plasmid column"),
    ("plasmid_call_source", "genomad, both, conflict, or platon for every other call (see the "
                            "AMR workbook's plasmid_call_source)"),
    ("plasmid_mobility", "conjugative, mobilisable or non-mobilisable for a plasmid"),
    ("plasmid_evidence", "the replicon call's own explanation"),
    ("phage", "is the gene inside a phage or prophage region? yes / no / uncertain, by the "
              "rule given for the AMR workbook's phage column"),
    ("phage_regions", "the virus caller's region(s) the answer rests on"),
    ("phage_evidence", "why, with the contig bases CheckV calls viral"),
]


def sha256_of(path):
    """The SHA-256 fingerprint of a file's contents, read in 1 MB pieces."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for piece in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(piece)
    return digest.hexdigest()


def verify_trait_provenance(path, sample, inputs, outputs):
    """Check that the trait calls were made with the catalogue quoted here.

    Takes in: {sample}_trait_provenance.json, written by trait_profile.py with
              the three trait tables (rule trait_profile), and the files this
              workbook reads: the two catalogue tables and EcoFoldDB's
              descriptions (inputs), and the profile, gene and audit tables
              (outputs). Each argument maps a name to a path.
    The workbook quotes caveats, assays and references from trait_rules_v1.tsv
    beside calls that were made earlier. If the catalogue has been edited since,
    those texts could describe rules the calls were not made with. The
    provenance file holds the SHA-256 fingerprint of every one of these files
    as they were when the calls were made; any difference stops the workbook,
    with the trait profile to be made again. Only contents are compared, so a
    moved output folder still passes.
    Gives back: the provenance record, whose parameters go to Run_QC.
    """
    with open(path, encoding="utf-8") as handle:
        record = json.load(handle)
    if record.get("schema_version") != "1.0":
        raise ValueError(f"{path}: unknown provenance format "
                         f"{record.get('schema_version')!r}")
    if record.get("sample") != sample:
        raise ValueError(f"{path}: written for sample {record.get('sample')!r}, not {sample!r}")
    for section, files in (("inputs", inputs), ("outputs", outputs)):
        recorded = record.get(section) or {}
        for name, file_path in files.items():
            expected = (recorded.get(name) or {}).get("sha256")
            if expected != sha256_of(file_path):
                raise ValueError(
                    f"{path}: {os.path.basename(file_path)} is not the file the trait calls "
                    f"were made with ({section}: {name}); run rule trait_profile again "
                    "for this sample")
    return record


def one_row_per_trait(rows, trait_ids, path):
    """Index a trait table by trait_id, checking it has every catalogue trait once."""
    by_id = {}
    for row in rows:
        if row["trait_id"] in by_id:
            raise ValueError(f"{path}: trait {row['trait_id']!r} appears twice")
        by_id[row["trait_id"]] = row
    if set(by_id) != set(trait_ids):
        raise ValueError(f"{path}: its traits differ from the catalogue's: missing "
                         f"{sorted(set(trait_ids) - set(by_id))}, extra "
                         f"{sorted(set(by_id) - set(trait_ids))}")
    return by_id


def build_traits_workbook(args):
    """Read every input of the traits workbook and return its sheets.

    The steps, in order:
      1. which strain and genome, as for the AMR workbook;
      2. the catalogue and the three trait tables (rule trait_profile), checked
         against the provenance file written with them;
      3. the plasmid and phage context of every gene with trait evidence;
      4. the sheets: Summary, Trait_predictions, Gene_evidence,
         Prediction_audit, Trait_rules, Run_QC and Data_dictionary.
    """
    out_root = args.output_root

    # 1. Which strain and genome.
    primary = read_primary_assembly(args.qc_genomes)
    taxonomy = read_taxonomy(args.taxonomy_dir, primary["bin_id"])

    # 2. The catalogue, the three trait tables and their provenance.
    _header, rules = read_tsv(args.rules, RULE_COLUMNS, skip_comment_lines=True)
    trait_ids = [rule["trait_id"] for rule in rules]
    rules_by_id = one_row_per_trait(rules, trait_ids, args.rules)
    _header, profiles = read_tsv(args.profile, PROFILE_COLUMNS)
    check_sample(profiles, args.sample, args.profile)
    profiles_by_id = one_row_per_trait(profiles, trait_ids, args.profile)
    gene_header, genes = read_tsv(args.genes, GENE_COLUMNS)
    check_sample(genes, args.sample, args.genes)
    provenance = verify_trait_provenance(
        args.provenance, args.sample,
        inputs={"rules": args.rules, "gene_dictionary": args.gene_dictionary,
                "catalogue": args.catalogue},
        outputs={"profile": args.profile, "genes": args.genes, "audit": args.audit})

    # 3. The context of every gene with trait evidence (column traits not NA).
    contig_lengths = platon_replicons.read_contig_lengths(args.contigs)
    crashed = platon_crashed(args.platon_dir)
    replicon_calls = read_replicon_calls(args.platon_dir, args.platon_prefix,
                                         contig_lengths, args.plasmid_concordance,
                                         args.replicons)
    regions = read_viral_regions(args.viral_table, args.phage_caller, contig_lengths)
    attach_checkv(regions, args.checkv_quality, args.checkv_contamination)
    caller_name = "geNomad" if args.phage_caller == "genomad" else "VirSorter2"

    gene_columns = [column for column in gene_header if column != "sample"]
    gene_rows = []
    for gene in genes:
        if is_missing(gene["traits"]):
            continue
        line = gene["__line__"]
        start = whole_number(gene["start"], args.genes, "start", line)
        stop = whole_number(gene["stop"], args.genes, "stop", line)
        contig = gene["contig"]
        if contig not in contig_lengths:
            raise ValueError(f"{args.genes}: line {line} is on contig {contig!r}, which the "
                             "assembly does not have")
        row = {column: typed(gene[column], column) for column in gene_columns}
        row.update(plasmid_answer(contig, replicon_calls, crashed))
        answer, names, reason = phage_answer(
            contig, gene_intervals(start, stop, contig_lengths[contig]), regions, caller_name)
        row.update({"phage": answer, "phage_regions": names, "phage_evidence": reason})
        gene_rows.append(row)

    # 4. The sheets.
    prediction_rows = []
    for trait_id in trait_ids:
        profile = profiles_by_id[trait_id]
        rule = rules_by_id[trait_id]
        row = {"sample": args.sample, "gtdb_taxon": taxonomy["lowest_taxon"]}
        row.update({column: profile[column] for column in PROFILE_COLUMNS if column != "sample"})
        row.update({"assay": rule["assay"], "reference": rule["reference"]})
        prediction_rows.append(row)
    decisions = collections.Counter(row["decision"] for row in prediction_rows)

    summary_rows = [
        {"item": "Sample", "value": args.sample,
         "note": "BacFlux sample name, from the input file name"},
        {"item": "GTDB taxon", "value": taxonomy["lowest_taxon"],
         "note": "lowest rank GTDB-Tk assigned; the full lineage is below"},
        {"item": "GTDB lineage", "value": taxonomy["lineage"],
         "note": f"GTDB-Tk classification method: {taxonomy['method']}"},
        {"item": "Traits called", "value": decisions["called"],
         "note": "every gene of a required set counts, and the genes that must lie together do"},
        {"item": "Traits partial", "value": decisions["partial"],
         "note": "a lead, not half a trait: the column missing names what is lacking"},
        {"item": "Traits not called", "value": decisions["not called"],
         "note": "this does not show that the trait is absent"},
        {"item": "Genes with trait evidence", "value": len(gene_rows),
         "note": f"of {len(genes)} rows in the trait gene table; one row each on Gene_evidence"},
        {"item": "Laboratory assays", "value": "not compared",
         "note": "these are predictions from the genome; the assay column suggests how to "
                 "test each trait"},
        {"item": "Safety", "value": "not assessed here",
         "note": "a trait call says nothing about whether the strain is safe to use: see the "
                 "GTDB taxon above and the strain's AMR workbook"},
        {"item": "How to read the context answers", "value": "yes / no / uncertain / not_assessed",
         "note": ANSWER_NOTE},
    ]

    parameters = provenance.get("parameters") or {}
    run_qc_rows = [
        {"item": "workbook layout", "value": WORKBOOK_LAYOUT,
         "detail": "raised when a sheet or column changes meaning"},
        {"item": "sample", "value": args.sample, "detail": "BacFlux sample name"},
        {"item": "primary assembly", "value": primary["bin_id"],
         "detail": f"technology {primary['technology']}"},
        {"item": "GTDB-Tk summary", "value": shown_path(taxonomy["path"], out_root),
         "detail": "the taxonomy on Summary"},
        {"item": "trait profile", "value": shown_path(args.profile, out_root),
         "detail": f"{len(profiles)} traits"},
        {"item": "trait gene table", "value": shown_path(args.genes, out_root),
         "detail": f"{len(gene_rows)} of {len(genes)} rows have trait evidence"},
        {"item": "trait audit", "value": shown_path(args.audit, out_root),
         "detail": "copied to Prediction_audit"},
        {"item": "trait rules", "value": shown_path(args.rules, out_root),
         "detail": f"{len(rules)} traits; SHA-256 {sha256_of(args.rules)}"},
        {"item": "gene dictionary", "value": shown_path(args.gene_dictionary, out_root),
         "detail": f"SHA-256 {sha256_of(args.gene_dictionary)}"},
        {"item": "EcoFoldDB descriptions", "value": shown_path(args.catalogue, out_root),
         "detail": f"SHA-256 {sha256_of(args.catalogue)}"},
        {"item": "provenance", "value": shown_path(args.provenance, out_root),
         "detail": "every fingerprint above matches the one recorded with the calls"},
        {"item": "calling parameters", "value": json.dumps(parameters, sort_keys=True),
         "detail": "as recorded when the calls were made"},
        {"item": "Platon", "value": shown_path(args.platon_dir, out_root),
         "detail": "exited with an error" if crashed else "finished"},
        {"item": "virus regions", "value": shown_path(args.viral_table, out_root),
         "detail": f"{caller_name}: {len(regions)} regions"},
    ]

    sheets = [
        make_sheet("Summary", SUMMARY_COLUMNS, summary_rows, wrap=("note",)),
        make_sheet("Trait_predictions", TRAIT_PREDICTION_COLUMNS, prediction_rows,
                   wrap=("caveat", "assay", "reference")),
        make_sheet("Gene_evidence",
                   [(column, "copied from the trait gene table; see docs/analysis/traits.md")
                    for column in gene_columns] + GENE_CONTEXT_COLUMNS,
                   gene_rows),
        copied_sheet("Prediction_audit", args.audit, args.sample,
                     "copied from {sample}_trait_profile_audit.tsv (rule trait_profile): the "
                     "reason for each decision; see docs/analysis/traits.md"),
        make_sheet("Trait_rules", TRAIT_RULE_COLUMNS, rules, wrap=("reference",)),
        make_sheet("Run_QC", RUN_QC_COLUMNS, run_qc_rows, wrap=("detail",)),
    ]
    sheets.append(dictionary_sheet(sheets, ANSWER_NOTE))
    return sheets, {"traits": len(prediction_rows), "genes with evidence": len(gene_rows)}


# ── Writing the .xlsx file ───────────────────────────────────────────────────

def cell_text(value):
    """Turn a value into the text of one cell, cut with a mark when too long."""
    text = "NA" if value is None else str(value)
    if len(text) > EXCEL_CELL_LIMIT:
        text = text[:EXCEL_CELL_LIMIT - len(CUT_MARK)] + CUT_MARK
    return text


def column_width(column, values, wrapped):
    """A readable width for one column: as wide as its header and first 200
    values, at most 60 characters (50 for a column whose text wraps)."""
    longest = len(column)
    for value in values[:200]:
        longest = max(longest, len(cell_text(value)))
    return min(longest + 2, 50 if wrapped else 60)


def write_workbook(path, sheets):
    """Write the sheets to one .xlsx file with XlsxWriter.

    Every cell is written with an explicit type: a Python int or float as a
    number, everything else as text. The workbook options switch off
    XlsxWriter's guessing as well, so no text is ever read as a formula, a
    link or a number. Each sheet gets a bold header row that stays in view
    while scrolling, and a filter on every column.
    """
    import xlsxwriter  # here, not at the top: see the module docstring

    workbook = xlsxwriter.Workbook(path, {
        "strings_to_numbers": False,
        "strings_to_formulas": False,
        "strings_to_urls": False,
    })
    header_format = workbook.add_format({"bold": True, "bg_color": "#D9D9D9",
                                         "bottom": 1, "valign": "top"})
    plain_format = workbook.add_format({"valign": "top"})
    wrap_format = workbook.add_format({"valign": "top", "text_wrap": True})

    for sheet in sheets:
        worksheet = workbook.add_worksheet(sheet["name"])
        columns = sheet["columns"]
        for column_index, column in enumerate(columns):
            worksheet.write_string(0, column_index, column, header_format)
            column_values = [row[column_index] for row in sheet["rows"]]
            wrapped = column in sheet["wrap"]
            worksheet.set_column(column_index, column_index,
                                 column_width(column, column_values, wrapped),
                                 wrap_format if wrapped else plain_format)
        for row_index, row in enumerate(sheet["rows"], start=1):
            for column_index, value in enumerate(row):
                cell_format = wrap_format if columns[column_index] in sheet["wrap"] else plain_format
                is_number = isinstance(value, (int, float)) and not isinstance(value, bool)
                if is_number and math.isfinite(value):
                    worksheet.write_number(row_index, column_index, value, cell_format)
                else:
                    worksheet.write_string(row_index, column_index, cell_text(value), cell_format)
        worksheet.freeze_panes(1, 0)
        worksheet.autofilter(0, 0, max(len(sheet["rows"]), 1), len(columns) - 1)
    workbook.close()


# ── Command line ─────────────────────────────────────────────────────────────

def make_parser():
    """The two commands and their arguments, as rules amr_workbook and
    traits_workbook (shared/90_report.smk) pass them."""
    parser = argparse.ArgumentParser(description="Write BacFlux's per-strain Excel workbooks.")
    commands = parser.add_subparsers(dest="command", required=True)

    def add_common(command):
        # What both workbooks need: the strain, the genome, and the context sources.
        command.add_argument("--sample", required=True, help="BacFlux sample name.")
        command.add_argument("--output-root", required=True,
                             help="The BacFlux output folder; paths shown in the workbook "
                                  "are relative to it.")
        command.add_argument("--qc-genomes", required=True,
                             help="02.assembly/{sample}/eval/{sample}_qc_genomes.tsv")
        command.add_argument("--taxonomy-dir", required=True,
                             help="GTDB-Tk folder of the sample (03.taxonomy/{sample}).")
        command.add_argument("--contigs", required=True,
                             help="The final contigs (02.assembly/{sample}/contigs_final.fasta).")
        command.add_argument("--platon-dir", required=True,
                             help="Platon's folder (06.plasmids/{sample}/platon).")
        command.add_argument("--platon-prefix", required=True,
                             help="The prefix of Platon's files (the contigs file's name "
                                  "without its extension).")
        command.add_argument("--plasmid-concordance", default=None,
                             help="geNomad's plasmid concordance table, when geNomad ran.")
        command.add_argument("--replicons", default=None,
                             help="The mobilome module's replicon table, when it ran.")
        command.add_argument("--phage-caller", required=True, choices=("genomad", "virsorter2"))
        command.add_argument("--viral-table", required=True,
                             help="The caller's region table (geNomad virus summary or "
                                  "VirSorter2 final-viral-boundary.tsv).")
        command.add_argument("--checkv-quality", required=True,
                             help="CheckV's quality_summary.tsv.")
        command.add_argument("--checkv-contamination", required=True,
                             help="CheckV's contamination.tsv (may be absent).")
        command.add_argument("--out", required=True, help="The .xlsx file to write.")

    amr = commands.add_parser("amr", help="the AMR workbook")
    add_common(amr)
    amr.add_argument("--card", required=True, help="ABRicate's CARD table.")
    amr.add_argument("--resfinder", required=True, help="ABRicate's ResFinder table.")
    amr.add_argument("--argannot", required=True, help="ABRicate's ARG-ANNOT table.")
    amr.add_argument("--abricate-databases", required=True,
                     help="05.amr/abricate/abricate_databases.tsv: ABRicate's version and "
                          "the date of its copy of each database.")
    amr.add_argument("--ice-table", default=None,
                     help="The mobilome module's ICE table, when it ran.")
    amr.add_argument("--mobility", default=None,
                     help="The mobilome module's mobility table, when it ran.")
    amr.add_argument("--mobility-audit", default=None,
                     help="The mobilome module's mobility audit, when it ran.")

    traits = commands.add_parser("traits", help="the traits workbook")
    add_common(traits)
    traits.add_argument("--profile", required=True, help="{sample}_trait_profile.tsv")
    traits.add_argument("--genes", required=True, help="{sample}_trait_genes.tsv")
    traits.add_argument("--audit", required=True, help="{sample}_trait_profile_audit.tsv")
    traits.add_argument("--rules", required=True, help="trait_rules_v1.tsv")
    traits.add_argument("--gene-dictionary", required=True, help="trait_genes_v1.tsv")
    traits.add_argument("--catalogue", required=True, help="EcoFoldDB_descriptions.txt")
    traits.add_argument("--provenance", required=True, help="{sample}_trait_provenance.json")
    return parser


def main(argv=None):
    """Parse the command line, build the workbook's sheets and write the file."""
    args = make_parser().parse_args(argv)
    if args.command == "amr":
        # The three mobilome tables come together or not at all.
        given = [args.ice_table, args.mobility, args.mobility_audit]
        if any(given) and not all(given):
            raise ValueError("--ice-table, --mobility and --mobility-audit go together")
        sheets, counts = build_amr_workbook(args)
    else:
        sheets, counts = build_traits_workbook(args)
    write_workbook(args.out, sheets)
    print(f"[BacFlux] wrote {args.out}: "
          + ", ".join(f"{count} {name}" for name, count in counts.items()))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        sys.exit(f"[BacFlux] workbook: {error}")
