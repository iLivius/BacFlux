#!/usr/bin/env python3
"""Put curated names on the ICE/IME candidates, using ICEberg.

conjscan_to_ice.py works out that an element IS an ICE — it has an integrase, a
relaxase and a mating-pair apparatus — but has no way to say WHICH ICE. Every row
came out with mge_name = NA, so a report could say "predicted self-transmissible
element" and never "ICEKp1". A name is what lets a reader look the thing up,
compare it with the literature, and recognise it in the next isolate.

Reads {sample}_ice_candidates.tsv (rule conjscan_ice) and a blastn of the whole
genome against ICEberg (rule iceberg_blast). Writes the same table with mge_name
and the naming evidence filled in, plus an audit TSV explaining every candidate
that could NOT be named. It creates, moves and drops nothing: the only things
that change are the name column and the evidence behind it, so no AMR gene's
mobility tier can move because of this step.

Why it joins by overlap instead of blasting each element: blasting the whole
genome once costs less than extracting every candidate to its own FASTA, and it
buys something — a curated ICE that OVERHANGS our interval still shows up, and
the overhang gets reported. That matters because our interval is a FLOOR, the
machinery span rather than the element's true ends, unless a tRNA-anchored att
pair was found. On the K. pneumoniae positive control (ATCC BAA-2146, CP006659.2)
the ICE call runs 4,604,073-4,644,558 while ICEberg's record for the same element
(ICEKpnATCCBAA-2146-1) runs 4,603,840-4,661,887: the right end is under-called by
about 17 kb, and this can say so rather than leave the reader to guess.

The caveat worth knowing: ICEberg entries are derived from published genomes, so
when the sample IS one of those genomes (or a close relative) the match comes out
at 100% and the name is exact rather than approximate. That is a real result, not
a bug — but "100% identity to ICEKpnATCCBAA-2146-1" is not independent
confirmation when the sample is ATCC BAA-2146.

The whole layer is off unless mobilome.iceberg.urls or mobilome.iceberg.dir is
set (MOBILOME_NAME_ICE in shared/00_common.smk, which also decides whether
colocalise.py reads the named table or the raw one). ICEberg publishes no licence
or terms of use, so BacFlux ships a URL and never the data.
"""

import argparse
import os
import sys


# ── What BLAST gives us, and the three naming thresholds ─────────────────────

# Columns requested from `blastn -outfmt 6`, in this order. rule iceberg_blast
# spells the same twelve out in its -outfmt string, so the rule and this list
# must be changed together — a mismatch shifts every field silently.
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length",
    "qstart", "qend", "sstart", "send",
    "evalue", "bitscore", "slen", "qlen",
]

# How much of OUR candidate the curated element must cover before we accept that
# the two are the same thing. Deliberately lenient compared with the transposon
# naming cascade (which wants 80% of the reference): an ICE is mosaic and its
# cargo varies between strains, so demanding near-completeness would refuse to
# name exactly the divergent elements a name would help most with.
# Overridden by mobilome.iceberg.min_overlap_fraction, which rule
# name_ice_elements always passes.
DEFAULT_MIN_OVERLAP_FRACTION = 0.50

# Identity floor. Below this the elements are related but not the same, and the
# name would mislead. Overridden by mobilome.iceberg.min_identity.
DEFAULT_MIN_IDENTITY = 80.0

# At or above this fraction of the REFERENCE element present, the name is used
# bare; below it the name is suffixed "-like", because what we have is clearly
# related to the curated element but is not the whole of it.
#
# This is the ONLY thing that controls the "-like" suffix, and it is deliberately
# NOT exposed on the command line — there is no --exact-name-coverage flag, so a
# user who wants to change when a name is hedged has to edit this line. It is
# named here so that grepping for "-like" or for this constant finds the one
# place that decides it (applied in name_elements, below).
EXACT_NAME_REFERENCE_COVERAGE = 0.80

# Element types that can carry a curated ICEberg name. conjscan_to_ice.py writes
# exactly one of five — ice, ime, aice, genomic_island, conjugative_region — and
# every one has to appear either here or in NOT_NAMEABLE_ELEMENT_TYPES below;
# test_name_ice_elements.py fails if one is in neither, or if something is listed
# here that the classifier cannot emit.
#
# `aice` is in the set because there is something to match. The database this
# layer searches (ICEberg's ICE_seq_all + IME_seq_all, 1,774 records) holds 26
# whose name begins AICE — AICEFraal5456, AICESare1562, … — and the classic
# actinomycete elements are catalogued alongside them under their historical
# names instead: SLP1 and pMEA100 are both there, and both are AICEs (te Poele
# et al. 2008, PMID 18523858). The class was added to conjscan_to_ice.py after
# this file was written and never added here, so AICE calls were refused a name
# for no stated reason.
#
# A name changes no tier: `aice` stays in colocalise.py's
# CONTEXT_ONLY_ELEMENT_TYPES, so an AMR gene inside one keeps tier 1 and its
# capped confidence either way. The name is for the reader, who can then look the
# element up.
#
# Untested in practice, and say so: an AICE call needs the ICEscan model set
# (mobilome.icescan.run, off by default), a name needs this layer
# (mobilome.iceberg.urls, empty by default), and AICEs live in actinomycetes —
# Streptomyces, Frankia, Salinispora, Mycobacterium. No benchmark run has taken
# this path. It is written from the catalogue's contents, not from a result.
NAMEABLE_ELEMENT_TYPES = {"ice", "ime", "aice", "genomic_island"}

# The one type that must NEVER carry an ICEberg name, and the sentence the audit
# uses to say why. The reason lives next to the type rather than in the audit
# call, because the old hard-coded sentence explained integrases to every refused
# element — which is how an AICE, which has an integrase by definition, was told
# it did not have one.
NOT_NAMEABLE_ELEMENT_TYPES = {
    "conjugative_region":
        "a conjugative region has a relaxase and a mating apparatus but no "
        "integrase, so it is deliberately not called an ICE (spec §8 phase 4, "
        "\"report it, do not call it an ICE\") and an ICE name would undo that",
}


# ── Reading ICEberg deflines and the two input files ─────────────────────────

def parse_iceberg_name(subject_id):
    """Pull the element name out of an ICEberg defline.

    ICEberg deflines are pipe-delimited and regular:

        ICEberg|1174|ICEKpnATCCBAA-2146-1|GenBank|CP006659.2|4603840..4661887 ...
        |0     |1   |2                   |3      |4         |5

    Field 2 is the element name, field 4 the source accession. Names themselves
    contain hyphens, which is why this splits on the pipe and not on punctuation.

    Returns (name, accession). An unexpected shape yields the raw string as the
    name rather than dropping the hit - a malformed defline should cost us the
    accession, not the element.
    """
    text = subject_id.strip()
    fields = text.split("|")
    if len(fields) >= 5 and fields[0].lower().startswith("iceberg"):
        return fields[2], fields[4]
    return text, ""


def read_blast_hits(path):
    """Read the BLAST tabular output. A missing or empty file is not an error.

    A line with fewer than the twelve expected fields is skipped rather than
    padded, because a short line means the columns no longer line up with
    BLAST_COLUMNS and every value read from it would be attributed to the wrong
    field. An empty result is normal — a genome with no ICEberg match at all.
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
                continue
            hits.append(dict(zip(BLAST_COLUMNS, fields)))
    return hits


def read_tsv(path):
    """Read a TSV into (header, rows-as-dicts), preserving column order.

    Column order is kept because the output is the SAME table with extra columns
    appended, and a reader comparing the two files side by side should not have
    to hunt for a column that moved.
    """
    if not path or not os.path.isfile(path):
        return [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return header, rows


# ── Small shared helpers ─────────────────────────────────────────────────────
# Coordinates throughout are 1-based and inclusive, the convention BLAST, Bakta
# and the rest of the mobilome module all use, which is why overlap_bp adds 1.

def to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def to_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def overlap_bp(a_start, a_end, b_start, b_end):
    """Bases shared by two closed intervals; 0 when they do not meet."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def group_hits_by_contig(hits):
    """Index BLAST hits by contig so each candidate only scans its own.

    An ICE on contig 3 can only be named by a curated element that matched
    contig 3 — an interval on one contig says nothing about coordinates on
    another, and joining across contigs would produce nonsense overlaps.
    """
    by_contig = {}
    for hit in hits:
        by_contig.setdefault(hit["qseqid"], []).append(hit)
    return by_contig


# ── Which curated element best explains this candidate ───────────────────────

def best_match_for_element(element_start, element_end, contig_hits,
                           min_overlap_fraction, min_identity):
    """Pick the curated element that best explains this candidate.

    Takes in: our candidate's coordinates and every ICEberg hit on its contig.
    Does:     keeps hits that overlap the candidate by at least
              min_overlap_fraction OF THE CANDIDATE and reach the identity floor,
              then takes the highest-scoring one.
    Returns:  (best_hit_dict_or_None, n_comparable), where n_comparable counts the
              distinct OTHER curated elements that fit about as well — within half
              a percentage point of the winner's identity and within 10% of its
              overlap, the bands applied below.

    Why n_comparable matters: ICEs of one species are often near-identical across
    strains, so a single genuine element routinely matches a dozen ICEberg entries
    at ~100%. Reporting one name without saying that would imply a precision the
    data does not have. On the K. pneumoniae positive control the winning name
    beats six others that are all within 0.01% identity.

    The pool it counts over is whatever rule iceberg_blast returned, which is
    capped at -max_target_seqs 50, so the number saturates rather than growing
    without limit.
    """
    element_length = element_end - element_start + 1
    if element_length <= 0:
        return None, 0

    candidates = []
    for hit in contig_hits:
        query_start = min(to_int(hit["qstart"]), to_int(hit["qend"]))
        query_end = max(to_int(hit["qstart"]), to_int(hit["qend"]))
        shared = overlap_bp(element_start, element_end, query_start, query_end)
        if shared <= 0:
            continue
        if shared / element_length < min_overlap_fraction:
            continue
        if to_float(hit["pident"]) < min_identity:
            continue
        candidates.append((to_float(hit["bitscore"]), shared, hit))

    if not candidates:
        return None, 0

    candidates.sort(key=lambda item: item[0], reverse=True)
    best_score, best_shared, best_hit = candidates[0]

    # How many OTHER curated elements explain this candidate about as well.
    #
    # Banded on IDENTITY and overlap, not on bitscore. Bitscore scales with
    # alignment LENGTH, so a 1% band on it is a different question depending on
    # how long the hit is: a short, perfect match to a small element scores far
    # below a long, slightly worse match to a big one, and would not be counted
    # as comparable even though it explains our element just as well. Since the
    # point of this number is to warn that the name is one of a near-identical
    # group, the comparison has to be length-independent.
    #
    # Counted per distinct element NAME so several HSPs of one reference do not
    # inflate it.
    winner_name, _ = parse_iceberg_name(best_hit["sseqid"])
    winner_identity = to_float(best_hit["pident"])
    winner_overlap = best_shared
    comparable_names = set()
    for _score, shared, hit in candidates[1:]:
        name, _ = parse_iceberg_name(hit["sseqid"])
        if name == winner_name:
            continue
        # Within half a percentage point of identity, and covering a comparable
        # amount of our element (within 10%). Both bands are conventions chosen to
        # catch the near-identical group, not biological boundaries.
        if abs(to_float(hit["pident"]) - winner_identity) > 0.5:
            continue
        if winner_overlap > 0 and abs(shared - winner_overlap) / winner_overlap > 0.10:
            continue
        comparable_names.add(name)

    best_hit = dict(best_hit)
    best_hit["_overlap_bp"] = best_shared
    best_hit["_overlap_fraction"] = best_shared / element_length
    return best_hit, len(comparable_names)


# ── What we write ────────────────────────────────────────────────────────────

def audit_row(sample, action, reason, detail, contig="NA", start="NA", end="NA"):
    """One line of the decision trail, same shape as the other mobilome audits."""
    return {
        "sample": sample, "contig": contig,
        "start": str(start), "end": str(end),
        "action": action, "reason": reason, "detail": detail,
    }


AUDIT_COLUMNS = ["sample", "contig", "start", "end", "action", "reason", "detail"]

# EVERY `action` / `reason` PAIR THIS SCRIPT CAN WRITE. Nothing is ever dropped
# from the ICE table here — only a name column is added — so there is no
# "discarded" action. What the audit records is why an element did or did not get
# a curated name.
#
#   action=named           a curated ICEberg name was applied to the element
#     iceberg_match                 the detail line carries the identity, the
#                                   overlap fraction and any comparable runners-up
#
#   action=kept_flagged    the element stays in the table with mge_name=NA
#     no_curated_name_found         ICEberg was searched and nothing cleared the
#                                   identity/overlap thresholds. Normal for a
#                                   genuinely novel element, or a host genus
#                                   ICEberg covers thinly.
#
#   action=not_applicable  naming was never attempted, for a stated reason
#     no_ice_candidates             the ICE table was empty; nothing to name
#     no_iceberg_hits               BLAST returned nothing against ICEberg
#     element_type_not_nameable     the element is a type that must NOT carry an
#                                   ICE name. Today that means exactly one thing,
#                                   a conjugative_region: no integrase, so not an
#                                   ICE, so no ICE name. The detail sentence is
#                                   looked up per type from
#                                   NOT_NAMEABLE_ELEMENT_TYPES, so a type added
#                                   later gets its own reason instead of
#                                   inheriting this one.

# Columns appended to the ICE table. mge_name already exists there (as NA); these
# are the evidence behind whatever it now says, and they let a reader judge a name
# rather than take it on trust.
ADDED_COLUMNS = [
    "iceberg_accession",       # the GenBank record the curated element came from
    "iceberg_identity",        # percent identity of the best hit
    "iceberg_overlap_fraction",  # how much of OUR element the hit covers
    "iceberg_reference_coverage",  # how much of the CURATED element is present
    "iceberg_alternatives",    # other curated elements that fit about as well
    "iceberg_reference_span",  # despite the name, the curated element's full
                               # LENGTH in bp (BLAST slen), not a start-end span
]


# ── Filling in the name, and the evidence behind it ──────────────────────────

def name_elements(sample, element_rows, hits,
                  min_overlap_fraction=DEFAULT_MIN_OVERLAP_FRACTION,
                  min_identity=DEFAULT_MIN_IDENTITY):
    """Fill in mge_name on each ICE/IME row, auditing every one left unnamed.

    Takes in: the rows of {sample}_ice_candidates.tsv, and BLAST hits of the
              genome against ICEberg.
    Does:     for each nameable element, finds the best overlapping curated hit.
    Returns:  (rows, audit_rows). Rows are modified in place; nothing is dropped.

    Every candidate gets every ADDED_COLUMN filled, with NA where naming did not
    happen, so the table has one shape whether ICEberg matched or not. Downstream
    only mge_name is read (colocalise.py carries it onto the AMR row); the
    iceberg_* columns are there for a person deciding whether to trust a name.
    """
    audit_rows = []
    by_contig = group_hits_by_contig(hits)

    if not hits:
        audit_rows.append(audit_row(
            sample, "not_applicable", "no_iceberg_hits",
            "BLAST returned no hit against ICEberg, so no candidate could be "
            "named. Normal for an element that is genuinely novel, or for a host "
            "genus ICEberg covers thinly."))

    for row in element_rows:
        element_type = (row.get("element_type") or "").strip().lower()
        contig = row.get("contig", "NA")
        start = to_int(row.get("start"))
        end = to_int(row.get("end"))

        for column in ADDED_COLUMNS:
            row.setdefault(column, "NA")

        if element_type not in NAMEABLE_ELEMENT_TYPES:
            # Why this type is refused, in its own words. An unknown type — one a
            # future class forgot to declare — falls back to a sentence that
            # states only what we actually know, rather than borrowing another
            # type's biology.
            why = NOT_NAMEABLE_ELEMENT_TYPES.get(
                element_type,
                "it is not one of the types this layer names ("
                + ", ".join(sorted(NAMEABLE_ELEMENT_TYPES)) + ")")
            audit_rows.append(audit_row(
                sample, "not_applicable", "element_type_not_nameable",
                f"{row.get('mge_id', '?')}: type '{element_type}' does not take an "
                f"ICEberg name - {why}.",
                contig=contig, start=start, end=end))
            continue

        best, n_alternatives = best_match_for_element(
            start, end, by_contig.get(contig, []), min_overlap_fraction, min_identity)

        if best is None:
            audit_rows.append(audit_row(
                sample, "kept_flagged", "no_curated_name_found",
                f"{row.get('mge_id', '?')}: no ICEberg element overlaps this "
                f"candidate by at least {100 * min_overlap_fraction:.0f}% of its "
                f"length at >={min_identity:.0f}% identity. The element stands on "
                "its own machinery evidence; it simply is not one ICEberg has "
                "catalogued.",
                contig=contig, start=start, end=end))
            continue

        name, accession = parse_iceberg_name(best["sseqid"])
        reference_length = to_int(best["slen"])

        # Measure the reference coverage over the part of the hit that is
        # actually INSIDE our element, not over the whole genome-wide HSP.
        #
        # The whole genome is blasted, so a curated element can match far beyond
        # our interval — on the K. pneumoniae positive control the ICEberg record
        # runs 17.5 kb past our call.
        # Judging "-like" on the full HSP therefore answers "how much of the
        # curated element exists anywhere on this contig?" when the question is
        # "how much of it is in the thing we are naming?". The first reading
        # awards a bare, confident name to an element we have only partly found.
        #
        # The clip is proportional: alignments here are near-collinear, so the
        # fraction of the HSP inside the element is a fair proxy for the fraction
        # of reference bases inside it.
        hit_start = min(to_int(best["qstart"]), to_int(best["qend"]))
        hit_end = max(to_int(best["qstart"]), to_int(best["qend"]))
        hit_span = max(1, hit_end - hit_start + 1)
        inside_bp = overlap_bp(start, end, hit_start, hit_end)
        aligned_inside = to_int(best["length"]) * (inside_bp / hit_span)
        reference_covered = (aligned_inside / reference_length) if reference_length else 0.0

        # "-like" when only part of the curated element is present. The element is
        # clearly related, but calling a 55% match by the bare name would claim an
        # identity the sequence does not support. EXACT_NAME_REFERENCE_COVERAGE is
        # the only thing that decides this, and it is not a command-line flag.
        display_name = name if reference_covered >= EXACT_NAME_REFERENCE_COVERAGE else f"{name}-like"

        row["mge_name"] = display_name
        row["iceberg_accession"] = accession or "NA"
        row["iceberg_identity"] = f"{to_float(best['pident']):.2f}"
        row["iceberg_overlap_fraction"] = f"{best['_overlap_fraction']:.3f}"
        row["iceberg_reference_coverage"] = f"{reference_covered:.3f}"
        row["iceberg_alternatives"] = str(n_alternatives)
        reference_start = min(to_int(best["sstart"]), to_int(best["send"]))
        reference_end = max(to_int(best["sstart"]), to_int(best["send"]))
        row["iceberg_reference_span"] = f"{reference_length}"

        audit_rows.append(audit_row(
            sample, "named", "iceberg_match",
            f"{row.get('mge_id', '?')}: named {display_name} "
            f"({to_float(best['pident']):.2f}% identity, covering "
            f"{100 * best['_overlap_fraction']:.0f}% of our interval and "
            f"{100 * reference_covered:.0f}% of the {reference_length} bp curated "
            f"element, accession {accession or 'NA'})."
            + (f" {n_alternatives} other curated element(s) fit about as well, so "
               "treat the exact name as one of a near-identical group."
               if n_alternatives else "")
            + (f" NOTE the curated element is {reference_length} bp while our "
               f"interval is {end - start + 1} bp: our boundaries are a floor, "
               "not the element's true ends."
               if reference_length > (end - start + 1) * 1.1 else ""),
            contig=contig, start=start, end=end))

        # reference_start / reference_end are where the hit falls in the CURATED
        # element's own coordinates. Nothing reports them: iceberg_reference_span
        # carries the reference's full length (slen) instead, which is the number
        # the audit line above sets our interval against. The `del` says so out
        # loud, so nobody hunts for a use that is not there.
        del reference_start, reference_end

    return element_rows, audit_rows


# ── Writing the two outputs ──────────────────────────────────────────────────

def write_tsv(path, columns, rows):
    """Write a TSV with a header, creating the parent directory if needed."""
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "NA")) for column in columns) + "\n")


# ── Command line ─────────────────────────────────────────────────────────────
# rule name_ice_elements passes all six flags, the two thresholds coming from
# mobilome.iceberg.min_identity and mobilome.iceberg.min_overlap_fraction. The
# defaults below are for a hand run.

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Name ICE/IME candidates from a BLAST of the genome against ICEberg. "
            "Adds a curated name and its evidence; never creates, moves or drops "
            "an element."
        )
    )
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--ice-table", required=True,
                        help="{sample}_ice_candidates.tsv from conjscan_to_ice.py.")
    parser.add_argument("--blast", required=True,
                        help="BLAST tabular output (outfmt 6) of the genome vs ICEberg.")
    parser.add_argument("--min-overlap-fraction", type=float,
                        default=DEFAULT_MIN_OVERLAP_FRACTION,
                        help="Fraction of OUR candidate a curated element must "
                             "cover before the two are treated as the same thing. "
                             f"Default {DEFAULT_MIN_OVERLAP_FRACTION}.")
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY,
                        help="Percent identity floor for a name. Default "
                             f"{DEFAULT_MIN_IDENTITY}.")
    parser.add_argument("--out-table", required=True,
                        help="The ICE table with names filled in.")
    parser.add_argument("--out-audit", required=True,
                        help="Why each candidate was or was not named.")
    args = parser.parse_args(argv)

    header, element_rows = read_tsv(args.ice_table)
    if not header:
        # conjscan_to_ice.py always writes a header, so an empty file means the
        # sample had no candidate at all. Write an empty table so the join
        # downstream still works, rather than failing the run.
        write_tsv(args.out_table, ["sample", "mge_id", "mge_name"], [])
        write_tsv(args.out_audit, AUDIT_COLUMNS, [audit_row(
            args.sample, "not_applicable", "no_ice_candidates",
            "The ICE table is empty, so there is nothing to name.")])
        print(f"Sample {args.sample}: no ICE/IME candidates to name.")
        return 0

    hits = read_blast_hits(args.blast)
    rows, audit = name_elements(
        args.sample, element_rows, hits,
        min_overlap_fraction=args.min_overlap_fraction,
        min_identity=args.min_identity)

    out_header = list(header) + [c for c in ADDED_COLUMNS if c not in header]
    write_tsv(args.out_table, out_header, rows)
    write_tsv(args.out_audit, AUDIT_COLUMNS, audit)

    n_named = sum(1 for row in rows
                  if row.get("mge_name") not in ("", "NA", None))
    print(
        f"Sample {args.sample}: {n_named} of {len(rows)} ICE/IME candidate(s) "
        f"named from ICEberg (reasons in {args.out_audit})."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
