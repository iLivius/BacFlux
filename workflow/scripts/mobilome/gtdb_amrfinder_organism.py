#!/usr/bin/env python3
"""Translate a GTDB-Tk species call into the matching AMRFinderPlus
`--organism` value — or into nothing at all, which is the common case here.

WHY THIS EXISTS
    AMRFinderPlus only reports POINT MUTATIONS (gyrA/rpoB-style chromosomal
    substitutions) when it is told which organism it is looking at, because
    every curated organism carries its own hand-checked list of
    resistance-conferring changes. Those point mutations are precisely the
    "intrinsic, chromosomal, not transferable" category that a homology screen
    like ABRicate structurally cannot see, and tier 1 of the mobilome module's
    mobility ladder (chromosomal / intrinsic candidate) is exactly that
    category. So when the isolate IS one of the curated organisms we want the
    flag; when it is not, we must run without it.

    AMRFinderPlus 4.2.7 curates only 31 organisms, and it names them
    NCBI-style. GTDB names do not always agree: GTDB splits genera and marks
    the split-off lineages with an alphanumeric suffix (`Pseudomonas_E`,
    `Klebsiella_A`), and gives genomes with no named species a placeholder
    epithet (`sp024807945`). A naive substring match would happily turn
    `s__Pseudomonas_E sp010095445` into `--organism Pseudomonas_aeruginosa`,
    and AMRFinderPlus would then score that genome's proteins against
    P. aeruginosa's curated mutation list. That is worse than no call at all,
    so every mismatch below is handled explicitly.

WHERE IT RUNS
    In the mobilome module's WP-A step, between GTDB-Tk (stage 03.taxonomy) and
    AMRFinderPlus (stage 05.amr / 08.mobilome). It is pure text handling —
    standard library only, no AMRFinderPlus call — so it unit-tests without any
    tool or database present.

INPUT (give exactly one of the first two)
    --gtdbtk-summary  the GTDB-Tk `gtdbtk.*.summary.tsv` written for this sample
                      (BacFlux writes one directory per sample, so every row in
                      the file belongs to this isolate). Only two of its columns
                      are read: `user_genome` and `classification`, the latter
                      looking like
                      `d__Bacteria;p__Pseudomonadota;...;g__Pseudomonas_E;s__Pseudomonas_E sp010095445`.
    --gtdbtk-dir      the GTDB-Tk output DIRECTORY for this sample
                      (03.taxonomy/{sample}), from which the summary is located
                      the same way the existing report rule does it: top level
                      first, then the classify/ subdirectory. Convenient for a
                      Snakemake rule that takes the whole directory as input.
    --sample          the BacFlux sample name, used to pick this sample's rows
                      (hybrid mode classifies both assemblies and writes
                      `<sample>_illumina` and `<sample>_ont`) and to label the
                      audit row.

OUTPUT
    --out-organism    one line with the AMRFinderPlus organism name on a match,
                      or a completely EMPTY file on no match. The Snakemake rule
                      turns it into the flag, e.g.
                      `org=$(cat {input.organism}); [ -n "$org" ] && set -- --organism "$org"`.
                      An empty file therefore means "run AMRFinderPlus without
                      --organism", which is the graceful, silent fallback the
                      spec asks for.
    --out-audit       a four-column TSV (sample, gtdb_classification,
                      matched_organism, reason) recording the decision and why
                      it was taken — the BacFlux convention that every
                      filtering/selection decision is auditable, as in
                      contig_taxonomy_decisions.tsv.

EXIT STATUS
    0 in every ordinary situation, INCLUDING no match, no row for the sample,
    an unreadable summary and a missing summary file. Environmental isolates
    normally have no curated organism (none of the six isolates screened in this
    project do), so "no match" is the normal path, not an error.
"""

import argparse
import os
import re
import sys


# ── The 31 curated AMRFinderPlus organisms ───────────────────────────────────
# Copied verbatim from `amrfinder --list_organisms -d {bakta_db}/amrfinderplus-db/latest`
# (version 4.2.7, recorded in docs/mobilome_wpA_ground_truth.md §"Spec §12 Q2").
# They are split by the RANK at which AMRFinderPlus curates them, because that
# decides how a GTDB name is allowed to match.

# Curated per SPECIES: the genome must be that species, not merely that genus.
# Key format is AMRFinderPlus's own: Genus_epithet, underscore-joined.
SPECIES_ORGANISMS = {
    "Acinetobacter_baumannii",
    "Bordetella_pertussis",
    "Burkholderia_cepacia",
    "Burkholderia_mallei",
    "Burkholderia_pseudomallei",
    "Citrobacter_freundii",
    "Clostridioides_difficile",
    "Corynebacterium_diphtheriae",
    "Enterobacter_asburiae",
    "Enterobacter_cloacae",
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Haemophilus_influenzae",
    "Helicobacter_pylori",
    "Klebsiella_oxytoca",
    "Klebsiella_pneumoniae",
    "Neisseria_gonorrhoeae",
    "Neisseria_meningitidis",
    "Pseudomonas_aeruginosa",
    "Serratia_marcescens",
    "Staphylococcus_aureus",
    "Staphylococcus_pseudintermedius",
    "Streptococcus_agalactiae",
    "Streptococcus_pneumoniae",
    "Streptococcus_pyogenes",
    "Vibrio_cholerae",
    "Vibrio_parahaemolyticus",
    "Vibrio_vulnificus",
}

# Curated per GENUS: AMRFinderPlus's mutation list for these covers the whole
# genus, so any species in it (named or not) is allowed to match.
#   Campylobacter — the curated set is the thermophilic food-chain species
#                   (C. jejuni / C. coli); GTDB keeps those in the unsuffixed
#                   g__Campylobacter and puts other lineages in Campylobacter_A,
#                   _B, _D ... which we deliberately do NOT match (see below).
#   Escherichia   — AMRFinderPlus's "Escherichia" covers E. coli AND Shigella.
#                   GTDB has already folded Shigella into g__Escherichia, so the
#                   genus-level rule handles both without a special case.
#   Salmonella    — one genus, effectively S. enterica.
GENUS_ORGANISMS = {
    "Campylobacter",
    "Escherichia",
    "Salmonella",
}

# Every genus AMRFinderPlus curates ANYTHING under. Used only to write a more
# honest "no match" reason: for Klebsiella we can say "the species is a GTDB
# placeholder", for Arthrobacter the truth is simply "nothing curated here".
CURATED_GENERA = {name.split("_")[0] for name in SPECIES_ORGANISMS} | GENUS_ORGANISMS

# All 31 valid --organism strings, for the caller (and the tests) to check against.
ALL_ORGANISMS = SPECIES_ORGANISMS | GENUS_ORGANISMS

# Column names read from the GTDB-Tk summary. Stable across GTDB-Tk 2.x.
GENOME_COLUMN = "user_genome"
CLASSIFICATION_COLUMN = "classification"

# Header of the audit TSV, in order.
AUDIT_HEADER = ["sample", "gtdb_classification", "matched_organism", "reason"]

# GTDB marks a lineage it has split off from the NCBI name with an underscore
# plus capital letters at the END of the word: Pseudomonas_E, Klebsiella_A,
# Escherichia coli_D. The suffix means "related to, but NOT the same taxon as"
# the unsuffixed name, so it must always block a match.
GTDB_SUFFIX = re.compile(r"_[A-Z]+$")

# GTDB gives genomes with no validly published species name a placeholder
# epithet: the literal "sp" followed by the GTDB accession digits, e.g.
# "s__Pseudomonas_E sp010095445". It is not a species name, so it can never
# match a curated species.
GTDB_PLACEHOLDER_SPECIES = re.compile(r"^sp\d*$")


# ── Reading the GTDB-Tk summary ──────────────────────────────────────────────

def read_gtdbtk_summary(summary_path):
    """Pull the (genome name, classification) pairs out of a GTDB-Tk summary TSV.

    Input:  path to `gtdbtk.bac120.summary.tsv` (or the ar53 equivalent) as
            written by rule gtdbtk into 03.taxonomy/{sample}/.
    Output: (rows, problem) where rows is a list of (user_genome,
            classification) tuples in file order, and problem is an empty
            string when all was well or a short human-readable message when the
            file could not be used.

    Nothing here raises: a missing or malformed summary must degrade to "no
    organism" rather than kill the run (spec §1.4), so the problem message is
    returned as data and ends up in the audit TSV.
    """
    if not os.path.isfile(summary_path):
        return [], "GTDB-Tk summary not found at %s" % summary_path

    with open(summary_path) as handle:
        # Strip \r as well as \n so a file that has been through Windows still
        # parses; GTDB fields never contain tabs, so a plain split is enough
        # (no csv quoting rules to honour).
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]

    if not lines:
        return [], "GTDB-Tk summary is empty: %s" % summary_path

    header = lines[0].split("\t")
    if GENOME_COLUMN not in header or CLASSIFICATION_COLUMN not in header:
        return [], "GTDB-Tk summary has no %s/%s columns: %s" % (
            GENOME_COLUMN, CLASSIFICATION_COLUMN, summary_path)

    genome_index = header.index(GENOME_COLUMN)
    classification_index = header.index(CLASSIFICATION_COLUMN)

    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        # A truncated line (interrupted GTDB-Tk run) is skipped rather than
        # crashing the parse; if that leaves no rows the caller reports it.
        if len(fields) <= max(genome_index, classification_index):
            continue
        rows.append((fields[genome_index].strip(), fields[classification_index].strip()))

    if not rows:
        return [], "GTDB-Tk summary has a header but no usable rows: %s" % summary_path
    return rows, ""


def find_summary_in_dir(gtdbtk_dir):
    """Locate the GTDB-Tk summary inside a per-sample GTDB-Tk output directory.

    Input:  03.taxonomy/{sample} as produced by rule gtdbtk.
    Output: (summary path, problem) — the path of the summary to read, or an
            empty path plus a message when there is none.

    GTDB-Tk writes `gtdbtk.bac120.summary.tsv` (bacteria) and/or
    `gtdbtk.ar53.summary.tsv` (archaea) into classify/, and classify_wf also
    leaves a symlink to it at the top level. We look at the top level first and
    then in classify/, which is the same order the existing report rule in the
    Snakefile uses. BacFlux is a bacterial workflow, so the bac120 file wins
    when both are present; the ar53 file is still accepted as a fallback so an
    archaeal isolate does not silently produce nothing at all.
    """
    if not os.path.isdir(gtdbtk_dir):
        return "", "GTDB-Tk output directory not found at %s" % gtdbtk_dir

    for directory in (gtdbtk_dir, os.path.join(gtdbtk_dir, "classify")):
        for basename in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
            candidate = os.path.join(directory, basename)
            if os.path.isfile(candidate):
                return candidate, ""
    return "", "no gtdbtk.*.summary.tsv found under %s" % gtdbtk_dir


def rows_for_sample(rows, sample):
    """Keep only the rows that belong to this sample.

    Input:  rows from read_gtdbtk_summary, plus the BacFlux sample name.
    Output: (selected_rows, used_every_row) — the second value is True when we
            had to fall back to using the whole file.

    In hybrid mode BacFlux classifies both assemblies of one isolate, so the
    summary holds two rows named `<sample>_illumina` and `<sample>_ont`; in
    single-technology mode there is one row named `<sample>`. Both are matched
    by "equals the sample name, or starts with the sample name and an
    underscore". If nothing matches by name we still use every row, because
    GTDB-Tk is run per sample into its own directory, so the file cannot
    contain another isolate — the genome may just have been named after the
    contigs file instead of the sample.
    """
    selected = [row for row in rows
                if row[0] == sample or row[0].startswith(sample + "_")]
    if selected:
        return selected, False
    return rows, True


# ── Reading a single GTDB classification string ──────────────────────────────

def split_classification(classification):
    """Break a GTDB lineage string into a {rank letter: name} dictionary.

    Input:  `d__Bacteria;p__Pseudomonadota;...;g__Pseudomonas_E;s__Pseudomonas_E sp010095445`
    Output: {"d": "Bacteria", ..., "g": "Pseudomonas_E", "s": "Pseudomonas_E sp010095445"}

    Ranks GTDB-Tk left blank (a genome placed only to genus ends in `;s__`) are
    simply absent from the dictionary, which the caller reads as "no species".
    A string with no `x__` fields at all (GTDB-Tk writes `Unclassified Bacteria`
    or `N/A` when placement fails) yields an empty dictionary.
    """
    ranks = {}
    for field in classification.split(";"):
        field = field.strip()
        # Every real field looks like a single rank letter, two underscores,
        # then the name. Anything shorter is an empty rank (`s__`) and anything
        # else is not a GTDB field at all.
        if len(field) > 3 and field[1:3] == "__":
            ranks[field[0]] = field[3:].strip()
    return ranks


def strip_gtdb_suffix(name):
    """Remove a GTDB lineage suffix and say whether there was one.

    Input:  a single GTDB name token, e.g. "Pseudomonas_E", "Klebsiella", "coli_D".
    Output: (base name, had a suffix) e.g. ("Pseudomonas", True), ("Klebsiella", False).

    The flag is what matters: a suffixed name is a DIFFERENT taxon from the
    unsuffixed NCBI one, so it must never be matched to a curated organism.
    """
    if not name:
        return "", False
    base = GTDB_SUFFIX.sub("", name)
    return base, base != name


def is_placeholder_species(epithet):
    """True for GTDB's stand-in species epithets ("sp024807945", "sp").

    These mark a genome with no validly published species name, so there is by
    definition no curated AMRFinderPlus organism for it.
    """
    return bool(GTDB_PLACEHOLDER_SPECIES.match(epithet))


def map_classification(classification):
    """Decide which AMRFinderPlus organism (if any) a GTDB lineage justifies.

    Input:  one `classification` string from the GTDB-Tk summary.
    Output: (organism, reason) — organism is one of the 31 curated names, or an
            empty string for no match; reason always explains the decision in
            plain words and is written straight into the audit TSV.

    The order of the checks is the whole point of this function, so it is spelt
    out step by step rather than compressed.
    """
    ranks = split_classification(classification)
    genus_field = ranks.get("g", "")
    species_field = ranks.get("s", "")

    # The species field is a binomial: "Genus epithet", and its genus token
    # carries the same GTDB suffix as the g__ field ("Pseudomonas_E sp010095445").
    # Prefer the explicit g__ field, but fall back to the species field's first
    # word if GTDB-Tk only filled the species in.
    species_words = species_field.split()
    genus = genus_field or (species_words[0] if species_words else "")
    epithet = species_words[1] if len(species_words) > 1 else ""

    if not genus and not epithet:
        return "", ("no GTDB classification (the summary has no genus or species "
                    "for this genome)")

    genus_base, genus_is_suffixed = strip_gtdb_suffix(genus)
    epithet_base, epithet_is_suffixed = strip_gtdb_suffix(epithet)

    # Step 1 — exact species match.
    # Only allowed when NEITHER token carries a GTDB suffix and the epithet is a
    # real name, because "Pseudomonas_E aeruginosa" or "Enterobacter cloacae_A"
    # would otherwise collapse onto the NCBI species they are explicitly not.
    if (epithet
            and not genus_is_suffixed
            and not epithet_is_suffixed
            and not is_placeholder_species(epithet)):
        species_key = genus + "_" + epithet
        if species_key in SPECIES_ORGANISMS:
            return species_key, "exact species match (GTDB s__%s %s)" % (genus, epithet)

    # Step 2 — genus-level match.
    # Three AMRFinderPlus organisms are curated for a whole genus, so any
    # species in it qualifies — including unnamed ones, and including Shigella,
    # which GTDB has already folded into g__Escherichia. A GTDB-suffixed genus
    # is still excluded: Campylobacter_A is not the curated Campylobacter.
    if genus and not genus_is_suffixed and genus in GENUS_ORGANISMS:
        return genus, ("genus-level match (AMRFinderPlus curates '%s' at genus "
                       "level; GTDB g__%s)" % (genus, genus))

    # Step 3 — no match. Say which of the GTDB/NCBI mismatches caused it, so the
    # audit row is useful six months later. The generic reason comes first
    # because for most environmental isolates (Paenibacillus, Arthrobacter) the
    # genus simply has nothing curated, and the suffix/placeholder detail would
    # be a red herring.
    if genus_base not in CURATED_GENERA:
        return "", ("no curated organism for this taxon (AMRFinderPlus curates 31 "
                    "organisms; genus '%s' is not one of them)" % genus)

    if genus_is_suffixed:
        return "", ("no curated organism for this taxon: GTDB-suffixed genus '%s' "
                    "is a lineage split off from NCBI '%s' and is not the same "
                    "taxon" % (genus, genus_base))

    if is_placeholder_species(epithet):
        return "", ("GTDB placeholder species (sp<digits>: '%s %s') - not a curated "
                    "organism" % (genus, epithet))

    if epithet_is_suffixed:
        return "", ("no curated organism for this taxon: GTDB-suffixed species "
                    "'%s %s' is a lineage split off from NCBI '%s %s' and is not "
                    "the same taxon" % (genus, epithet, genus, epithet_base))

    if not epithet:
        return "", ("no curated organism for this taxon: GTDB placed this genome "
                    "only to genus '%s', which AMRFinderPlus curates per species "
                    "rather than per genus" % genus)

    return "", ("no curated organism for this taxon: '%s %s' is not one of the "
                "curated species in genus '%s'" % (genus, epithet, genus))


# ── Turning a whole summary file into one decision ───────────────────────────

def decide_organism(summary_path, sample):
    """Read the summary, pick this sample's rows, and settle on one organism.

    Input:  path to the GTDB-Tk summary and the BacFlux sample name.
    Output: (organism, classification_text, reason), all strings; organism is
            empty when no curated organism applies. classification_text is what
            goes in the audit's gtdb_classification column.

    Hybrid runs give two rows for one isolate. They are two assemblies of the
    SAME DNA, so they should agree; the deterministic rule is:
      * every row is mapped independently;
      * if they all land on the same answer (the normal case, including "no
        organism"), that answer is used;
      * if they disagree, NO organism is emitted. A wrong --organism means
        point mutations are scored against the wrong curated list, so when the
        evidence is contradictory the safe move is to run AMRFinderPlus without
        it and say so in the audit.
    """
    rows, problem = read_gtdbtk_summary(summary_path)
    if problem:
        # Missing / empty / malformed file: no organism, but not a failure.
        return "", "", "no curated organism for this taxon: %s" % problem

    selected, used_every_row = rows_for_sample(rows, sample)

    # Map each row on its own, keeping the row name so the reason can point at
    # the offending assembly when the two disagree.
    decisions = []
    for genome_name, classification in selected:
        organism, reason = map_classification(classification)
        decisions.append((genome_name, classification, organism, reason))

    distinct_organisms = sorted({d[2] for d in decisions})
    distinct_classifications = []
    for decision in decisions:
        if decision[1] not in distinct_classifications:
            distinct_classifications.append(decision[1])
    classification_text = " | ".join(distinct_classifications)

    # The rows disagree about the organism -> refuse to guess.
    if len(distinct_organisms) > 1:
        per_row = ", ".join("%s=%s" % (d[0], d[2] or "none") for d in decisions)
        return "", classification_text, (
            "no curated organism for this taxon: the assemblies of this sample "
            "were classified differently (%s), so no --organism is used rather "
            "than risk the wrong point-mutation list" % per_row)

    organism = decisions[0][2]
    reason = decisions[0][3]

    # Book-keeping notes appended to the reason, so the audit row explains
    # exactly which rows the decision rests on.
    if len(decisions) > 1:
        reason += " [%d assemblies of this sample agree: %s]" % (
            len(decisions), ", ".join(d[0] for d in decisions))
    if used_every_row:
        reason += (" [no row named after sample '%s'; used every row in the "
                   "per-sample summary]" % sample)
    return organism, classification_text, reason


# ── Writing the two outputs ──────────────────────────────────────────────────

def write_outputs(out_organism, out_audit, sample, classification, organism, reason):
    """Write the organism file consumed by the AMRFinderPlus rule and the audit
    TSV kept for the record.

    The organism file is deliberately EMPTY (zero bytes) when there is no match,
    so the shell test `[ -n "$(cat ...)" ]` in the rule is all that is needed to
    decide whether to pass --organism.

    In the audit, an empty classification or organism is written as the literal
    "NA" — the same convention AMRFinderPlus itself uses for missing values, and
    it keeps every column of the TSV visibly filled.
    """
    for path in (out_organism, out_audit):
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)

    with open(out_organism, "w") as handle:
        if organism:
            handle.write(organism + "\n")

    with open(out_audit, "w") as handle:
        handle.write("\t".join(AUDIT_HEADER) + "\n")
        handle.write("\t".join([
            sample,
            classification or "NA",
            organism or "NA",
            reason,
        ]) + "\n")


def main(argv=None):
    """Command-line entry point used by the Snakemake rule.

    Always returns 0 in ordinary use, including every no-match path: an
    environmental isolate with no curated organism is the expected result, not
    an error, and the AMRFinderPlus rule must still run.
    """
    parser = argparse.ArgumentParser(
        description=("Map a GTDB-Tk species call to an AMRFinderPlus --organism "
                     "value; emit nothing when there is no curated organism."))
    # Two ways in, because a Snakemake rule may declare either the summary file
    # or the whole GTDB-Tk directory as its input. Exactly one must be given.
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--gtdbtk-summary", default="",
                        help="GTDB-Tk gtdbtk.*.summary.tsv for this sample")
    source.add_argument("--gtdbtk-dir", default="",
                        help="GTDB-Tk output directory for this sample "
                             "(03.taxonomy/{sample}); the summary is found inside it")
    parser.add_argument("--sample", required=True,
                        help="BacFlux sample name (selects the summary rows)")
    parser.add_argument("--out-organism", required=True,
                        help="write the organism name here, or nothing on no match")
    parser.add_argument("--out-audit", required=True,
                        help="write the four-column decision TSV here")
    args = parser.parse_args(argv)

    # Resolve whichever input form was given to a single summary file. A
    # directory with no summary in it is treated exactly like a missing file:
    # no organism, no failure.
    if args.gtdbtk_dir:
        summary_path, problem = find_summary_in_dir(args.gtdbtk_dir)
    else:
        summary_path, problem = args.gtdbtk_summary, ""

    if problem:
        organism = ""
        classification = ""
        reason = "no curated organism for this taxon: %s" % problem
    else:
        organism, classification, reason = decide_organism(summary_path, args.sample)

    write_outputs(args.out_organism, args.out_audit, args.sample,
                  classification, organism, reason)

    # One line to the Snakemake log so the choice is visible without opening the
    # audit file.
    if organism:
        sys.stdout.write("%s: AMRFinderPlus --organism %s (%s)\n"
                         % (args.sample, organism, reason))
    else:
        sys.stdout.write("%s: no AMRFinderPlus --organism (%s)\n"
                         % (args.sample, reason))
    return 0


if __name__ == "__main__":
    sys.exit(main())
