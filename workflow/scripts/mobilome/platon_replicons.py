#!/usr/bin/env python3
"""Turn Platon's plasmid call into the per-contig replicon table the mobility
ladder needs.

WHY THIS EXISTS
    The mobility ladder's top two rungs are about PLASMIDS:
        tier 5  on a mobilizable plasmid   -> transferable, but needs a helper
        tier 6  on a conjugative plasmid   -> PREDICTED self-transmissible
    Deciding between them needs two facts per contig: (a) is it a plasmid at all,
    and (b) does that plasmid carry the machinery to move itself.

    Platon already answers BOTH, which is easy to miss. It splits the assembly
    into `*.chromosome.fasta` and `*.plasmid.fasta`, and its per-contig table
    carries `# Conjugation`, `# Mobilization` and `# OriT` columns. So no extra
    tool is needed for the plasmid case — this script just reads what Platon
    already worked out and writes it in the shape colocalise.py expects.

    (The CHROMOSOMAL case — conjugation machinery sitting on the chromosome, i.e.
    an ICE — is a different question that Platon structurally cannot answer: it
    excludes long contigs from consideration entirely. That is CONJscan's job,
    handled by conjscan_to_ice.py.)

WHAT IT TAKES IN
    A Platon output directory, as produced by rule plasmid_search:
        <prefix>.tsv                 per-contig table (may be header-only)
        <prefix>.chromosome.fasta    contigs Platon called chromosomal
        <prefix>.plasmid.fasta       contigs Platon called plasmid

WHAT IT PRODUCES
    A TSV with one row per contig:
        contig, replicon, replicon_id, plasmid_mobility, mobility_evidence
    consumed by colocalise.py (--replicons).

    IMPORTANT: a contig MISSING from this table is treated downstream as
    `unknown`, not as chromosome. That is deliberate — calling something
    "chromosomal, intrinsic candidate" when we simply never looked would be the
    most misleading mistake this module could make. So this script lists EVERY
    contig it can see, and says `unknown` honestly when Platon did not classify
    one.

WHY THERE IS NO --out-audit HERE, WHEN EVERY OTHER SCRIPT IN THIS MODULE HAS ONE
    The project rule is that every filtering decision is auditable with a stated
    reason. This script honours that rule in a COLUMN rather than in a separate
    file, and the reason is that it filters nothing: it writes exactly one row per
    contig it can see, and drops none. There is no discard list to explain.

    What it does instead is REASON about each contig, and that reasoning is
    written out per row in `mobility_evidence` (the Platon gene counts behind a
    conjugative/mobilisable call) and `replicon_call_source` (which of Platon's
    outputs the chromosome/plasmid call came from, including the
    --min-chromosome-bp size rule). So the audit is the table itself. If this
    script ever starts DROPPING contigs, it needs a real audit file at that point.
"""

import argparse
import os
import sys


def read_fasta_ids(path):
    """Collect the sequence IDs from a FASTA file.

    Takes in: a path that may not exist (Platon writes no plasmid FASTA at all
              for a genome with no plasmid — a perfectly normal result).
    Does:     read only the header lines, taking the first whitespace-delimited
              token as the ID, which is the same convention every other join in
              BacFlux uses.
    Returns:  a set of contig IDs, empty when the file is absent or empty.
    """
    ids = set()
    if not path or not os.path.exists(path):
        return ids
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def read_platon_table(path):
    """Read Platon's per-contig table into {contig: {column: value}}.

    Takes in: <prefix>.tsv from Platon. It is legitimately HEADER-ONLY when
              Platon found no plasmid, so an empty result is not an error.
    Does:     a plain TSV read keyed on Platon's own `ID` column.
    Returns:  a dict; empty when there is nothing to read.
    """
    rows = {}
    if not path or not os.path.exists(path):
        return rows

    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if len(lines) < 2:
        # Header only (or completely empty) = Platon ran and found no plasmid.
        return rows

    header = lines[0].split("\t")
    for line in lines[1:]:
        values = line.split("\t")
        record = dict(zip(header, values))
        contig = record.get("ID", "").strip()
        if contig:
            rows[contig] = record
    return rows


def read_genomad_concordance(path):
    """Read geNomad's second opinion about which contigs are plasmids.

    WHERE IT COMES FROM
        rule plasmid_concordance (60_plasmid.smk) compares Platon's plasmid call
        with geNomad's and writes `{sample}_plasmid_concordance.tsv`. Its rows are
        the UNION of the two tools' plasmid calls, so a contig appears here if
        EITHER tool thought it was a plasmid — which is exactly the set of contigs
        where a second opinion can change our mind.

    WHY THE MOBILITY LADDER CARES
        Platon alone decides tiers 5 and 6, and it can miss a plasmid. When it
        does, every AMR gene on that contig is reported as "chromosomal, intrinsic
        candidate" (tier 1) — the module's worst possible error, because it is the
        direction that hides transferability from the reader.

    NOTE ON OPTIONALITY
        geNomad is opt-in (it is academic/non-commercial licensed, so it cannot be
        BacFlux's default — see the spec §3.1 licensing note). When it was not run
        this file does not exist, and that is a normal, silent no-op: the caller
        simply gets an empty dict and Platon decides alone, exactly as before.

    Returns: {contig: {"call": ..., "agreement": ..., "score": ...}}, empty when
             the file is absent, empty or header-only.
    """
    calls = {}
    if not path or not os.path.exists(path):
        return calls

    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if len(lines) < 2:
        # Header-only is normal and means "no plasmid candidate from either tool".
        return calls

    header = lines[0].split("\t")
    for line in lines[1:]:
        record = dict(zip(header, line.split("\t")))
        contig = (record.get("contig") or "").strip()
        if not contig:
            continue
        calls[contig] = {
            "call": (record.get("genomad_call") or "").strip().lower(),
            "agreement": (record.get("agreement") or "").strip().lower(),
            "score": (record.get("genomad_score") or "NA").strip(),
        }
    return calls


def count_from(record, column):
    """Read one of Platon's count columns as an integer, tolerantly.

    Platon writes plain integers here, but a missing column or a non-numeric
    placeholder (it uses 'NA' for coverage, for instance) must not crash the
    module — an unreadable count simply means "no evidence of this feature".
    """
    raw = (record.get(column) or "").strip()
    try:
        return int(raw)
    except (TypeError, ValueError):
        return 0


def classify_plasmid_mobility(record):
    """Decide how transferable one PLASMID contig is, from Platon's own counts.

    The biology, and why the order of these tests matters:
      * conjugation genes present  -> the plasmid encodes its own mating
        apparatus, so it can move itself. Ladder tier 6, "PREDICTED
        self-transmissible" (never bare "transmissible" — the confirmatory
        experiment is a mating assay, not software).
      * no conjugation, but a relaxase (mobilization) or an origin of transfer
        (OriT) -> the plasmid can be moved, but only if some OTHER element in
        the same cell supplies the machinery. Ladder tier 5, "mobilisable".
      * neither -> nothing suggests it moves at all.

    Returns (mobility_class, evidence_text). The evidence text is passed
    straight through to the report so a reader can see WHY, not just WHAT.
    """
    n_conjugation = count_from(record, "# Conjugation")
    n_mobilization = count_from(record, "# Mobilization")
    n_orit = count_from(record, "# OriT")
    inc_types = (record.get("Inc Type(s)") or "").strip()

    evidence_parts = []
    if n_conjugation:
        evidence_parts.append(f"conjugation={n_conjugation}")
    if n_mobilization:
        evidence_parts.append(f"mobilization={n_mobilization}")
    if n_orit:
        evidence_parts.append(f"oriT={n_orit}")
    # An Inc (incompatibility) type is a replicon-family label, not mobility
    # evidence, but it is useful context so it rides along when Platon found one.
    if inc_types and inc_types not in {"0", "-", "NA"}:
        evidence_parts.append(f"inc={inc_types}")

    evidence = ";".join(evidence_parts) if evidence_parts else "no mobility genes detected by Platon"

    if n_conjugation > 0:
        return "conjugative", evidence
    if n_mobilization > 0 or n_orit > 0:
        return "mobilisable", evidence
    return "non-mobilisable", evidence


# Platon refuses to even look at a contig longer than this — read from its own
# source (platon/constants.py: MAX_CONTIG_LENGTH = 500000). Anything above it is
# reported as "too long" and never classified, so it appears in NEITHER the
# chromosome nor the plasmid FASTA. On a closed genome that means the chromosome
# itself is missing from Platon's output entirely.
PLATON_MAX_CONTIG_LENGTH = 500000

# Above THIS length we are willing to call an unclassified contig the chromosome.
# Why not just use Platon's 500 kb cut-off? Because megaplasmids exist, and some
# run to 1–2 Mb — a 600 kb contig Platon skipped could genuinely be one. Above
# 2 Mb that is no longer a realistic worry for a bacterial isolate, so this is
# the point where "Platon skipped it because it is far too big to be a plasmid"
# becomes safe to state. Between the two thresholds we say `unknown` and let the
# confidence cap downstream reflect that we really do not know.
DEFAULT_MIN_CHROMOSOME_BP = 2000000


def build_rows(platon_dir, prefix, extra_contigs, contig_lengths=None,
               min_chromosome_bp=DEFAULT_MIN_CHROMOSOME_BP, genomad_calls=None):
    """Join Platon's three outputs into one row per contig.

    Takes in: the Platon output directory and the file prefix Platon used (which
              is the input genome's basename, so BacFlux passes it in rather
              than guessing), plus any contig IDs seen elsewhere (from the
              assembly) so contigs Platon never mentioned are still listed, plus
              geNomad's second opinion when it was run.
    Produces: a list of dicts ready to write, one per contig, sorted by name so
              the output is stable between runs.
    """
    chromosome_ids = read_fasta_ids(os.path.join(platon_dir, f"{prefix}.chromosome.fasta"))
    plasmid_ids = read_fasta_ids(os.path.join(platon_dir, f"{prefix}.plasmid.fasta"))
    table = read_platon_table(os.path.join(platon_dir, f"{prefix}.tsv"))
    genomad_calls = genomad_calls or {}

    all_contigs = set(chromosome_ids) | set(plasmid_ids) | set(table) | set(extra_contigs)

    rows = []
    for contig in sorted(all_contigs):
        genomad = genomad_calls.get(contig, {})
        genomad_says_plasmid = genomad.get("call") == "plasmid"

        if contig in plasmid_ids:
            replicon = "plasmid"
            mobility, evidence = classify_plasmid_mobility(table.get(contig, {}))
        elif contig in chromosome_ids and genomad_says_plasmid:
            # BOTH tools looked and DISAGREED. Platon positively called this
            # chromosome; geNomad positively called it plasmid.
            #
            # The call is NOT flipped. Platon is BacFlux's default replicon
            # caller and made an active call, and quietly overriding it on a
            # disagreement would be the same overclaim this module exists to
            # avoid. But the disagreement is recorded in the evidence text and
            # the source is marked 'conflict', which caps confidence downstream -
            # because reporting "chromosomal, intrinsic candidate" at full
            # confidence when a second tool says plasmid would be the most
            # misleading thing this table could say.
            replicon = "chromosome"
            mobility = "NA"
            evidence = (
                f"CONFLICT: Platon called this contig chromosomal but geNomad "
                f"called it a plasmid (score {genomad.get('score', 'NA')}). The "
                "Platon call is kept because it is the default caller, but the "
                "two tools disagree - check this contig by hand before treating "
                "genes on it as intrinsic."
            )
        elif contig in chromosome_ids:
            replicon = "chromosome"
            # Platon says nothing about chromosomal mobility, and pretending it
            # did would be wrong. Whether this chromosome carries an ICE is
            # CONJscan's question, answered separately.
            mobility = "NA"
            evidence = "chromosomal contig; Platon does not assess chromosomal mobility"
        elif genomad_says_plasmid:
            # Platon made NO call about this contig, and geNomad calls it a
            # plasmid. Nothing is being overridden here - this is the only
            # opinion there is, so using it is not a promotion but simply reading
            # the available evidence.
            #
            # WHY THIS CASE EXISTS AND WHY IT MATTERS: without it, such a contig
            # falls through to 'unknown' below, and every AMR gene on it is
            # reported as tier 1, "chromosomal, intrinsic candidate". That is the
            # worst error the ladder can make, because it is the direction that
            # HIDES transferability: an acquired, mobile resistance gene reported
            # as an intrinsic species trait.
            #
            # The mobility is left unknown on purpose. geNomad's plasmid summary
            # counts conjugation genes, but Platon's per-contig mobility columns
            # are what classify_plasmid_mobility reads and they are absent here.
            # CONJscan runs over the whole proteome independently, so
            # colocalise.py can still find a typed conjugative system on this
            # contig and resolve tier 5 vs 6 from that - which is better evidence
            # than a gene count anyway.
            replicon = "plasmid"
            mobility = "unknown"
            evidence = (
                f"geNomad called this contig a plasmid (score "
                f"{genomad.get('score', 'NA')}); Platon did not classify it at "
                "all. Mobility is not typed - Platon's conjugation/mobilization "
                "counts are what types it and they are absent, so tier 5 vs 6 "
                "rests on whether CONJscan finds machinery on this contig."
            )
        else:
            # Seen in the assembly but absent from BOTH Platon FASTAs. The usual
            # cause is Platon's size filter: it skips anything over 500 kb, so on
            # a closed genome the chromosome itself never appears in its output.
            length = (contig_lengths or {}).get(contig)
            if length is not None and length >= min_chromosome_bp:
                # Far larger than any realistic plasmid, so this is the
                # chromosome. Say so — but record that it is INFERRED from
                # Platon's size exclusion, not a positive call Platon made.
                replicon = "chromosome"
                mobility = "NA"
                evidence = (
                    f"inferred chromosomal: {length} bp, above Platon's "
                    f"{PLATON_MAX_CONTIG_LENGTH} bp size filter and above the "
                    f"{min_chromosome_bp} bp megaplasmid ceiling, so Platon never "
                    "classified it"
                )
            elif length is not None and length > PLATON_MAX_CONTIG_LENGTH:
                # In the awkward band: too long for Platon to assess, but small
                # enough that a megaplasmid is still possible. Do not guess.
                replicon = "unknown"
                mobility = "NA"
                evidence = (
                    f"not classified: {length} bp exceeds Platon's "
                    f"{PLATON_MAX_CONTIG_LENGTH} bp size filter, but is small "
                    "enough that a megaplasmid cannot be ruled out"
                )
            else:
                replicon = "unknown"
                mobility = "NA"
                evidence = "contig not classified by Platon"

        # WHO decided this, so downstream can weight it. A one-tool call is real
        # evidence but weaker than two tools agreeing, and a conflict is weaker
        # still - colocalise.py caps confidence accordingly rather than treating
        # every replicon call as equally certain.
        if not genomad_calls:
            call_source = "platon"          # geNomad was not run at all
        elif replicon == "plasmid" and contig in plasmid_ids and genomad_says_plasmid:
            call_source = "both"
        elif replicon == "plasmid" and contig in plasmid_ids:
            call_source = "platon"
        elif replicon == "plasmid":
            call_source = "genomad"         # the rescue case above
        elif contig in chromosome_ids and genomad_says_plasmid:
            call_source = "conflict"
        else:
            call_source = "platon"

        rows.append({
            "contig": contig,
            "replicon": replicon,
            "replicon_id": contig,
            "plasmid_mobility": mobility,
            "mobility_evidence": evidence,
            "replicon_call_source": call_source,
        })
    return rows


def read_contig_lengths(path):
    """Contig ID -> length, read straight from the assembly FASTA.

    Needed because Platon's output alone cannot tell us how long a contig it
    SKIPPED was, and that length is what decides whether an unclassified contig
    is safely the chromosome or genuinely ambiguous.
    """
    lengths = {}
    if not path or not os.path.exists(path):
        return lengths
    current = None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())
    return lengths


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Build the per-contig replicon table from Platon output."
    )
    parser.add_argument("--sample", required=True, help="Sample name (for messages).")
    parser.add_argument("--platon-dir", required=True,
                        help="Platon output directory for this sample.")
    # NAMING NOTE, so nobody "tidies" these and breaks the rule that calls them.
    # This script says --contigs where isescan_to_table.py says --genome-fasta and
    # conjscan_to_ice.py/att_search.py say --genome, and it says --out where the
    # others say --out-table. All of them are handed the same assembly FASTA and
    # all write one table; the names simply drifted. They are left alone because
    # renaming a flag silently breaks any command a colleague has saved, and the
    # matching change in 80_mobilome.smk would have to land in the same edit.
    parser.add_argument("--prefix", required=True,
                        help="Platon's file prefix. Both Platon and geNomad name "
                             "every output file after their input's basename, so "
                             "this is that basename - see the reconciliation in "
                             "00_common.smk, which is where the value comes from.")
    parser.add_argument("--contigs", default="",
                        help="Assembly FASTA, so contigs Platon never mentioned "
                             "are still listed as 'unknown' rather than omitted.")
    parser.add_argument("--min-chromosome-bp", type=int, default=DEFAULT_MIN_CHROMOSOME_BP,
                        help="An unclassified contig at least this long is taken to be "
                             "the chromosome (Platon skips anything over 500 kb, and "
                             "megaplasmids can reach ~2 Mb). Default "
                             f"{DEFAULT_MIN_CHROMOSOME_BP}.")
    parser.add_argument("--genomad-concordance", default="",
                        help="Optional {sample}_plasmid_concordance.tsv from the "
                             "plasmid_concordance rule. When geNomad was run, its "
                             "plasmid calls are used as a second opinion: a contig "
                             "Platon never classified but geNomad calls a plasmid "
                             "is reported as a plasmid rather than falling through "
                             "to 'unknown' (and thence to a wrong 'intrinsic' "
                             "call), and a straight disagreement is flagged. "
                             "Absent file = geNomad was not run = no-op.")
    parser.add_argument("--out", required=True, help="Destination TSV.")
    args = parser.parse_args(argv)

    contig_lengths = read_contig_lengths(args.contigs) if args.contigs else {}
    genomad_calls = read_genomad_concordance(args.genomad_concordance)
    rows = build_rows(args.platon_dir, args.prefix, sorted(contig_lengths),
                      contig_lengths=contig_lengths,
                      min_chromosome_bp=args.min_chromosome_bp,
                      genomad_calls=genomad_calls)

    columns = ["contig", "replicon", "replicon_id", "plasmid_mobility",
               "mobility_evidence", "replicon_call_source"]
    with open(args.out, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    n_plasmid = sum(1 for row in rows if row["replicon"] == "plasmid")
    n_chromosome = sum(1 for row in rows if row["replicon"] == "chromosome")
    n_unknown = sum(1 for row in rows if row["replicon"] == "unknown")
    n_conjugative = sum(1 for row in rows if row["plasmid_mobility"] == "conjugative")
    n_mobilisable = sum(1 for row in rows if row["plasmid_mobility"] == "mobilisable")

    print(
        f"Sample {args.sample}: {len(rows)} contig(s) — {n_chromosome} chromosome, "
        f"{n_plasmid} plasmid ({n_conjugative} conjugative, {n_mobilisable} mobilisable), "
        f"{n_unknown} unclassified."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
