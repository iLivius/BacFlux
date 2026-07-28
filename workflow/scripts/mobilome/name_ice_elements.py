#!/usr/bin/env python3
"""Put curated names on the ICE/IME candidates, using ICEberg.

WHY THIS EXISTS
    conjscan_to_ice.py works out that an element IS an ICE - it has an integrase,
    a relaxase and a mating-pair apparatus - but it has no way to say WHICH ICE.
    Every row came out with mge_name = NA, so a report could say "predicted
    self-transmissible element" and never "ICEKp1". A name is what lets a reader
    look the thing up, compare it with the literature, and recognise it in the
    next isolate.

WHAT IT DOES
    Takes the ICE/IME table and a BLAST of the whole genome against ICEberg, and
    for each candidate finds the best curated element overlapping it. It does not
    create, move or drop any element: the only thing that changes is the name
    column and the evidence behind it.

WHY IT JOINS BY OVERLAP RATHER THAN BLASTING EACH ELEMENT
    Blasting the whole genome once is cheaper than extracting every candidate to
    its own FASTA, and it has a useful side effect: a curated ICE that overhangs
    our interval still shows up, and the overhang is reported. That matters
    because our interval is a FLOOR - the machinery span, not the element's true
    ends, unless a tRNA-anchored att pair was found. On the KPNIH1 positive
    control the ICE call runs 4,604,073-4,644,558 while ICEberg's record for the
    same element runs 4,603,840-4,661,887: we under-call the right end by about
    17 kb, and this script can say so instead of leaving the reader to guess.

A CAVEAT WORTH KNOWING
    ICEberg entries are derived from published genomes, so if the sample IS one
    of those genomes (or a close relative) the match will be 100% and the name is
    exact rather than approximate. That is a real result, not a bug - but do not
    read "100% identity to ICEKpnATCCBAA-2146-1" as independent confirmation when
    the sample is ATCC BAA-2146.

WHAT IT PRODUCES
    A copy of the ICE table with mge_name and the naming evidence filled in, plus
    an audit TSV explaining every candidate that could NOT be named.
"""

import argparse
import os
import sys


# Columns requested from `blastn -outfmt 6`, in this order. The rule and this
# list must be changed together.
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
DEFAULT_MIN_OVERLAP_FRACTION = 0.50

# Identity floor. Below this the elements are related but not the same, and the
# name would mislead.
DEFAULT_MIN_IDENTITY = 80.0

# At or above this fraction of the REFERENCE element present, the name is used
# bare; below it the name is suffixed "-like", because what we have is clearly
# related to the curated element but is not the whole of it.
EXACT_NAME_REFERENCE_COVERAGE = 0.80

# Element types that can carry a curated ICEberg name. A conjugative_region has
# no integrase and is explicitly NOT an ICE (spec §8 phase 4, "report it, do not
# call it an ICE"), so giving it an ICE name would undo that distinction.
NAMEABLE_ELEMENT_TYPES = {"ice", "ime", "cime", "genomic_island"}


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
    """Read the BLAST tabular output. A missing or empty file is not an error."""
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
    """Read a TSV into (header, rows-as-dicts), preserving column order."""
    if not path or not os.path.isfile(path):
        return [], []
    with open(path, encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return header, rows


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
    """Index BLAST hits by contig so each candidate only scans its own."""
    by_contig = {}
    for hit in hits:
        by_contig.setdefault(hit["qseqid"], []).append(hit)
    return by_contig


def best_match_for_element(element_start, element_end, contig_hits,
                           min_overlap_fraction, min_identity):
    """Pick the curated element that best explains this candidate.

    Takes in: our candidate's coordinates and every ICEberg hit on its contig.
    Does:     keeps hits that overlap the candidate by at least
              min_overlap_fraction OF THE CANDIDATE and reach the identity floor,
              then takes the highest-scoring one.
    Returns:  (best_hit_dict_or_None, n_comparable) where n_comparable counts the
              other references that scored within 1% of the winner.

    Why n_comparable matters: ICEs of one species are often near-identical across
    strains, so a single genuine element routinely matches a dozen ICEberg entries
    at ~100%. Reporting one name without saying that would imply a precision the
    data does not have. On KPNIH1 the winning name beats six others that are all
    within 0.01% identity.
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
    # Compared on bitscore, within 1%, and counted per distinct element name so
    # several HSPs of one reference do not inflate it.
    winner_name, _ = parse_iceberg_name(best_hit["sseqid"])
    comparable_names = set()
    for score, _shared, hit in candidates[1:]:
        if best_score > 0 and score < best_score * 0.99:
            continue
        name, _ = parse_iceberg_name(hit["sseqid"])
        if name != winner_name:
            comparable_names.add(name)

    best_hit = dict(best_hit)
    best_hit["_overlap_bp"] = best_shared
    best_hit["_overlap_fraction"] = best_shared / element_length
    return best_hit, len(comparable_names)


def audit_row(sample, action, reason, detail, contig="NA", start="NA", end="NA"):
    """One line of the decision trail, same shape as the other mobilome audits."""
    return {
        "sample": sample, "contig": contig,
        "start": str(start), "end": str(end),
        "action": action, "reason": reason, "detail": detail,
    }


AUDIT_COLUMNS = ["sample", "contig", "start", "end", "action", "reason", "detail"]

# Columns this script adds to the ICE table. mge_name already exists there (as
# NA); the rest are the evidence behind whatever it now says.
ADDED_COLUMNS = [
    "iceberg_accession",       # the GenBank record the curated element came from
    "iceberg_identity",        # percent identity of the best hit
    "iceberg_overlap_fraction",  # how much of OUR element the hit covers
    "iceberg_reference_coverage",  # how much of the CURATED element is present
    "iceberg_alternatives",    # other curated elements that fit about as well
    "iceberg_reference_span",  # the curated element's own extent, for comparison
]


def name_elements(sample, element_rows, hits,
                  min_overlap_fraction=DEFAULT_MIN_OVERLAP_FRACTION,
                  min_identity=DEFAULT_MIN_IDENTITY):
    """Fill in mge_name on each ICE/IME row, auditing every one left unnamed.

    Takes in: the rows of {sample}_ice_candidates.tsv, and BLAST hits of the
              genome against ICEberg.
    Does:     for each nameable element, finds the best overlapping curated hit.
    Returns:  (rows, audit_rows). Rows are modified in place; nothing is dropped.
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
            # A conjugative_region has no integrase and is deliberately NOT called
            # an ICE. Hanging an ICE name on it would quietly undo that.
            audit_rows.append(audit_row(
                sample, "not_applicable", "element_type_not_nameable",
                f"{row.get('mge_id', '?')}: type '{element_type}' does not take an "
                "ICEberg name. A conjugative region has no integrase, so it is not "
                "an ICE, and naming it as one would contradict the classification.",
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
        reference_covered = (to_int(best["length"]) / reference_length) if reference_length else 0.0

        # "-like" when only part of the curated element is present. The element is
        # clearly related, but calling a 55% match by the bare name would claim an
        # identity the sequence does not support.
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

        # Not used for the name, but computed above and worth keeping honest.
        del reference_start, reference_end

    return element_rows, audit_rows


def write_tsv(path, columns, rows):
    """Write a TSV with a header, creating the parent directory if needed."""
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
