#!/usr/bin/env python3
"""Turn the two CARD read-mapping passes into one annotated, tiered report.

Reads the strict and relaxed BBMap coverage tables written by rule map_amr_db
(shared/50_amr.smk) — the same trimmed reads mapped onto CARD's protein homolog
model twice, at read identity 0.99 and 0.95 — plus CARD's own aro_index.tsv,
and writes one row per CARD reference sequence that reached the coverage
threshold at either setting. rule card_mapping_report runs it and puts the
result in 05.amr/mapping/{sample}/{sample}_CARD_report.tsv. That file is the
one to open; the two covstats tables and v1's AMR_legend stay beside it as the
raw evidence.

The two covstats must be the SAME sample against the SAME reference, differing
only in the identity filter — the rule guarantees that by having the relaxed
pass reuse the index the strict pass built. --min-covered (70 by default, from
CARD_MIN_COVERED in 00_common.smk) is the fraction of a reference sequence's
length that reads must cover, at one setting or the other, for it to reach the
report at all.

Why the same reads are mapped twice
-----------------------------------
The CARD leg maps trimmed reads straight onto CARD's reference sequences, so it
sees genes the assembly collapsed or lost. It has always mapped at a read
identity of 0.99 — near-exact matching. That gives excellent specificity and
poor sensitivity: a real but divergent member of a resistance family simply
goes missing. The rule now maps twice, strict and relaxed, and this script joins
the two passes so a gene found only at the relaxed setting is REPORTED AS SUCH
rather than lost. A `divergent` call is information, not merely extra
sensitivity: it says a variant of this family is present but the reference
allele is not. How much that buys was measured rather than assumed, and it is
less than one would expect — see rule map_amr_db in shared/50_amr.smk.

Why a raw CARD hit list overstates the resistome
------------------------------------------------
Not every CARD entry is an acquired resistance gene. CARD's protein homolog
model also contains multi-drug efflux components, the transcriptional regulators
of those components, accessory/structural proteins, and — most awkwardly —
entries whose resistance mechanism is "resistance by absence". Counting all of
those as "AMR genes found" badly overstates what a genome carries, especially
for environmental Gram-negatives whose chromosomes encode whole RND efflux
repertoires. CARD already classifies every entry; the workflow already
downloads its index; nothing used the classification. This script surfaces it.

What the report deliberately does not say
-----------------------------------------
It does NOT decide whether a gene is intrinsic or acquired. That is a
population-level judgement about a species (is the gene shared by the vast
majority of wild-type strains?), it needs a comparison across many genomes, and
nothing here does that — see docs/methods_amr_intrinsic_acquired.md. The
mechanism category is a HINT, not a determination: efflux is enriched for
chromosomal core systems but includes plenty of mobile ones (tet(A), tet(L),
qacA, mef(A)), and "antibiotic inactivation" spans both the intrinsic
chromosomal AmpC and the acquired CTX-M. The script reports the facts CARD
already states and leaves the judgement to the reader, exactly as the mobility
ladder reports `intrinsic_candidate` rather than `intrinsic`.

The table it writes, column by column
-------------------------------------
  aro_accession         CARD's ARO accession — the key that joins a mapped
                        sequence to CARD's own classification
                        source: the ARO:NNNN field of the covstats defline
  aro_name              the gene name, e.g. OXA-58
                        source: aro_index.tsv "ARO Name", falling back to the
                                last pipe-field of the defline when the
                                accession is missing from the index
  detection             exact | divergent
                        source: which pass cleared --min-covered
  covered_percent_id99  how much of the reference sequence's LENGTH the reads
  covered_percent_id95  covered, one column per identity filter. The numbers in
                        the names come from --strict-id / --relaxed-id, so a
                        figure always says which filter produced it.
                        source: column 5 (Covered_percent) of each covstats
  category              resistance_determinant | efflux_other |
                        efflux_component | regulator |
                        presence_indicates_susceptibility
                        source: classify()
  resistance_mechanism  CARD's mechanism text, semicolon-separated when an entry
                        carries several
                        source: aro_index.tsv "Resistance Mechanism"
  amr_gene_family       source: aro_index.tsv "AMR Gene Family"
  drug_class            source: aro_index.tsv "Drug Class"
  reference_organism    the organism CARD's REFERENCE sequence came from, never
                        a claim about this sample
                        source: the trailing [brackets] of the covstats defline
  note                  the plain-language warning classify() attached, empty
                        for the rows that need none
                        source: classify()

Consumed by: whoever reads 05.amr/mapping/{sample}/. It is a terminal product —
rule all asks for it (00_common.smk) and no other rule parses it.
"""

import argparse
import csv
import re
import sys


# ── What a covstats row and a CARD defline look like ────────────────────────

# BBMap covstats: column 1 is the reference defline, column 5 is Covered_percent.
# Kept as named constants because the file has no stable header beyond "#ID".
COVSTATS_ID_COLUMN = 0
COVSTATS_COVERED_PERCENT_COLUMN = 4

# A CARD defline looks like:
#   gb|AE004091.2|+|2810008-2813197|ARO:3000804|MexF [Pseudomonas aeruginosa PAO1]
# The ARO accession is the stable join key to aro_index.tsv. It is pulled by
# pattern rather than by field position because the number of pipe-separated
# fields differs between CARD's models (the homolog model carries a coordinate
# field that some other models omit), and a positional split silently picks up
# the wrong field when that happens.
ARO_PATTERN = re.compile(r"(ARO:\d+)")

# The organism CARD attributes the reference sequence to, in trailing brackets.
# Reported for context only: it is the organism the REFERENCE came from, NOT a
# claim about the sample. Someone reading "MexF [Pseudomonas aeruginosa]" in a
# report on a soil isolate needs that distinction spelled out.
ORGANISM_PATTERN = re.compile(r"\[([^\]]+)\]\s*$")


# ── CARD's own words for the things this report groups by ───────────────────

# CARD's controlled vocabulary for Resistance Mechanism. An entry may carry
# several, semicolon-separated.
MECHANISM_EFFLUX = "antibiotic efflux"
MECHANISM_ABSENCE = "resistance by absence"

# Gene families that are multi-component chromosomal transport machinery rather
# than acquired resistance determinants. Matched as substrings of CARD's own
# "AMR Gene Family" text, so a family naming several systems still matches.
STRUCTURAL_EFFLUX_FAMILIES = (
    "resistance-nodulation-cell division",
    "ATP-binding cassette",
    "major facilitator superfamily",
    "small multidrug resistance",
    "multidrug and toxic compound extrusion",
    "outer membrane porin",
)

# Names ending in R are, by a long-standing bacterial-genetics convention, the
# REGULATOR of the operon named without it: MexR regulates mexAB, ArmR regulates
# armZ/mexR, BltR regulates blt. A regulator modulates expression of a resistance
# system; it does not itself confer resistance, so counting it as a resistance
# gene double-counts the system it controls. The convention is not a guarantee,
# which is why this produces a FLAG for a human to weigh and never a filter.
#
# The convention is ALSO only applied where CARD's own mechanism for the entry is
# efflux, and that guard is not cosmetic. CARD's regulator entries are almost all
# efflux-pump regulators (MexR, MexZ, AcrR, MarR, RamR, SoxR, CpxR, BltR), filed
# under "antibiotic efflux". Without the guard the suffix rule misfires on real
# determinants that merely end in R: measured on a Bacillus isolate from this
# project's own test set, vmlR — an ABC-F ribosomal protection protein, CARD
# mechanism "antibiotic target protection", a genuine resistance gene — was
# demoted to `regulator`, which is exactly the wrong direction to be wrong in.
REGULATOR_SUFFIX_PATTERN = re.compile(r"(?:^|\s)([A-Za-z][A-Za-z0-9]*R)$")


# ── Reading the two inputs: BBMap's coverage, then CARD's own index ─────────


def parse_covstats(path):
    """Read a BBMap covstats file into {defline: covered_percent}.

    Rows whose Covered_percent cannot be parsed are skipped rather than crashing:
    BBMap emits a trailing blank line, and a truncated run can leave a partial
    final row. A dropped row can only ever cost a hit, never invent one, so a
    malformed line is preferable to a failed report.

    A missing header, though, is fatal. rule map_amr_db (shared/50_amr.smk) runs
    without `set -e` (a deliberate v1 decision — see the rule's comment), so a
    BBMap pass that dies leaves an empty covstats behind and the rule still exits
    0. Read leniently, that empty file would become a clean report saying the
    genome carries no AMR gene, which is the worst possible failure mode: wrong,
    plausible, and silent. BBMap always writes the '#ID' header and then one row
    per reference sequence, so no header means the pass did not finish.
    """
    coverage = {}
    saw_header = False
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#ID"):
                saw_header = True
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) <= COVSTATS_COVERED_PERCENT_COLUMN:
                continue
            try:
                covered = float(fields[COVSTATS_COVERED_PERCENT_COLUMN])
            except ValueError:
                continue
            coverage[fields[COVSTATS_ID_COLUMN].strip()] = covered

    if not saw_header:
        raise ValueError(
            f"'{path}' has no '#ID' header line, so BBMap did not finish writing "
            "it. Refusing to report an empty coverage table as 'no AMR genes "
            "found' — check the map_amr log for this sample."
        )
    return coverage


def parse_aro_index(path):
    """Read CARD's aro_index.tsv into {ARO accession: {family, drug_class, mechanism}}.

    Only the four columns this report needs are kept — ARO Name, AMR Gene
    Family, Drug Class, Resistance Mechanism — and they are looked up BY HEADER
    NAME rather than by position, because CARD has added columns between releases
    and a positional read would silently start reporting the wrong field.
    """
    index = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "ARO Accession" not in reader.fieldnames:
            raise ValueError(
                f"'{path}' does not look like CARD's aro_index.tsv "
                "(no 'ARO Accession' column)."
            )
        for row in reader:
            accession = (row.get("ARO Accession") or "").strip()
            if not accession:
                continue
            index[accession] = {
                "gene_family": (row.get("AMR Gene Family") or "").strip(),
                "drug_class": (row.get("Drug Class") or "").strip(),
                "mechanism": (row.get("Resistance Mechanism") or "").strip(),
                "aro_name": (row.get("ARO Name") or "").strip(),
            }
    return index


# ── Label one CARD entry: determinant, pump part, regulator, or absence ─────

# Every label below is a statement about the CARD ENTRY, not about the sample,
# and none of them is a filter: nothing is dropped, only labelled. They exist so
# a reader can group the report instead of reading one flat count that mixes
# acquired beta-lactamases with the chromosomal efflux machinery every
# Pseudomonad carries.
def classify(aro_name, gene_family, mechanism):
    """Return (category, notes) for one CARD entry — the single most useful
    label, plus any plain-language observations for the note column. The three
    arguments are the ARO Name, AMR Gene Family and Resistance Mechanism strings
    from aro_index.tsv, each empty when the accession is missing from the index.
    """
    mechanisms = [m.strip() for m in mechanism.split(";") if m.strip()]
    family_lower = gene_family.lower()
    notes = []

    # Checked FIRST, because it inverts the meaning of the whole row. CARD holds
    # entries — mgrB, OmpK35, OmpK36, OprD, LamB, carO — where resistance arises
    # from the gene being ABSENT. Detecting such a gene by presence therefore
    # indicates the OPPOSITE of resistance. Reporting it as a hit is not a
    # borderline call, it is a sign error, so it gets its own category and a note
    # saying what the presence actually means.
    if MECHANISM_ABSENCE in mechanisms:
        notes.append(
            "resistance arises from ABSENCE of this gene; its presence here "
            "indicates the susceptible state, not resistance"
        )
        return "presence_indicates_susceptibility", notes

    # A regulator of an efflux system rather than a determinant of resistance.
    # Checked before the efflux branches below so that CpxR is reported as the
    # regulator it is rather than as one more pump subunit, but gated on the
    # mechanism so that a determinant whose name happens to end in R is not
    # demoted — see the note on REGULATOR_SUFFIX_PATTERN above.
    if MECHANISM_EFFLUX in mechanisms and REGULATOR_SUFFIX_PATTERN.search(aro_name.strip()):
        notes.append(
            "name follows the -R regulator convention; regulates a resistance "
            "system rather than conferring resistance itself"
        )
        return "regulator", notes

    # Multi-component chromosomal transport machinery. These are the entries that
    # inflate counts on Gram-negative genomes: a single RND system contributes an
    # inner-membrane transporter, a periplasmic adaptor and an outer-membrane
    # channel, each a separate CARD entry.
    if MECHANISM_EFFLUX in mechanisms and any(
        family in family_lower for family in STRUCTURAL_EFFLUX_FAMILIES
    ):
        notes.append(
            "component of a multi-subunit efflux system; such systems are "
            "commonly chromosomal core machinery, and one system contributes "
            "several separate entries to this table"
        )
        return "efflux_component", notes

    # Efflux, but not one of the multi-subunit families listed above: a
    # single-protein pump, or a family this list does not name. Still efflux, so
    # it stays out of the determinant count rather than being folded into it.
    if MECHANISM_EFFLUX in mechanisms:
        return "efflux_other", notes

    # Everything else — inactivating enzymes, target protection, target
    # alteration, target replacement: the entries a reader normally means by "an
    # AMR gene". What the label does NOT claim is that the gene was acquired;
    # that judgement is not made anywhere in this script.
    return "resistance_determinant", notes


# ── Join the two passes and annotate what cleared either one ────────────────


def build_rows(strict, relaxed, aro_index, min_covered, strict_id, relaxed_id):
    """Join the two coverage maps and annotate every sequence that passed either.

    A sequence qualifies if it reached min_covered at EITHER identity setting.
    Detection tier:
        exact     - reached the threshold at the strict identity filter
        divergent - reached it only at the relaxed one, i.e. a variant of this
                    family is present but the reference allele is not
    """
    rows = []
    for defline in sorted(set(strict) | set(relaxed)):
        covered_strict = strict.get(defline, 0.0)
        covered_relaxed = relaxed.get(defline, 0.0)
        if covered_strict < min_covered and covered_relaxed < min_covered:
            continue

        aro_match = ARO_PATTERN.search(defline)
        accession = aro_match.group(1) if aro_match else ""
        organism_match = ORGANISM_PATTERN.search(defline)
        reference_organism = organism_match.group(1).strip() if organism_match else ""

        entry = aro_index.get(accession, {})
        gene_family = entry.get("gene_family", "")
        drug_class = entry.get("drug_class", "")
        mechanism = entry.get("mechanism", "")
        # Prefer CARD's curated ARO Name; fall back to the defline field when the
        # accession is missing from the index (a database/version mismatch).
        aro_name = entry.get("aro_name", "")
        if not aro_name:
            without_organism = ORGANISM_PATTERN.sub("", defline).strip()
            aro_name = without_organism.split("|")[-1].strip()

        category, notes = classify(aro_name, gene_family, mechanism)
        detection = "exact" if covered_strict >= min_covered else "divergent"

        rows.append({
            "aro_accession": accession or "NA",
            "aro_name": aro_name or "NA",
            "detection": detection,
            f"covered_percent_id{strict_id}": f"{covered_strict:.2f}",
            f"covered_percent_id{relaxed_id}": f"{covered_relaxed:.2f}",
            "category": category,
            "resistance_mechanism": mechanism or "NA",
            "amr_gene_family": gene_family or "NA",
            "drug_class": drug_class or "NA",
            "reference_organism": reference_organism or "NA",
            "note": "; ".join(notes) if notes else "",
        })
    return rows


# ── Command line, ranking, and the per-category summary ─────────────────────

# Arguments come from rule card_mapping_report in shared/50_amr.smk.
# --strict-id and --relaxed-id only NAME the two coverage columns
# (covered_percent_id99 / covered_percent_id95). The filtering they refer to
# happened earlier, inside BBMap in rule map_amr_db, and nothing here can change
# it: they are passed so that each figure says which pass produced it.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--covstats-strict", required=True)
    parser.add_argument("--covstats-relaxed", required=True)
    parser.add_argument("--aro-index", required=True)
    parser.add_argument("--strict-id", default="99",
                        help="strict identity filter, for the column name (e.g. 99)")
    parser.add_argument("--relaxed-id", default="95",
                        help="relaxed identity filter, for the column name (e.g. 95)")
    parser.add_argument("--min-covered", type=float, default=70.0,
                        help="minimum Covered_percent for a sequence to be reported")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    strict = parse_covstats(args.covstats_strict)
    relaxed = parse_covstats(args.covstats_relaxed)
    aro_index = parse_aro_index(args.aro_index)

    rows = build_rows(strict, relaxed, aro_index, args.min_covered,
                      args.strict_id, args.relaxed_id)

    # Sort so the rows a reader should look at first come first: real
    # determinants above efflux machinery, then by how well covered they are.
    #
    # Coverage is taken as the HIGHER of the two figures, not the strict one. A
    # divergent hit has a strict coverage of 0 by definition, so ranking on the
    # strict column alone would push every divergent gene to the bottom of its
    # category — burying exactly the rows the relaxed pass was added to surface.
    category_order = {
        "resistance_determinant": 0,
        "efflux_other": 1,
        "efflux_component": 2,
        "regulator": 3,
        "presence_indicates_susceptibility": 4,
    }
    strict_column = f"covered_percent_id{args.strict_id}"
    relaxed_column = f"covered_percent_id{args.relaxed_id}"
    rows.sort(key=lambda r: (category_order.get(r["category"], 9),
                             -max(float(r[strict_column]),
                                  float(r[relaxed_column]))))

    # Written column order — keep it in step with the column-by-column list in
    # the module docstring, which is where a reader looks up what each one means.
    columns = ["aro_accession", "aro_name", "detection", strict_column,
               f"covered_percent_id{args.relaxed_id}", "category",
               "resistance_mechanism", "amr_gene_family", "drug_class",
               "reference_organism", "note"]
    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # A one-line summary per category, so the log says what the table contains
    # without anyone opening it. Counting by category is the whole point: a bare
    # total is the number this script exists to stop people quoting.
    counts = {}
    for row in rows:
        counts[row["category"]] = counts.get(row["category"], 0) + 1
    divergent = sum(1 for row in rows if row["detection"] == "divergent")
    print(f"Sample {args.sample}: {len(rows)} CARD sequences at >= "
          f"{args.min_covered:.0f}% covered length.")
    for category in sorted(counts, key=lambda c: category_order.get(c, 9)):
        print(f"  {category}: {counts[category]}")
    print(f"  detected only at the relaxed identity filter (divergent): {divergent}")
    if counts.get("presence_indicates_susceptibility"):
        print("  NOTE: rows flagged 'presence_indicates_susceptibility' are CARD "
              "entries whose resistance mechanism is loss of the gene. Their "
              "presence is NOT evidence of resistance.")


if __name__ == "__main__":
    sys.exit(main())
