#!/usr/bin/env python3
"""Translate a GTDB-Tk species call into the matching AMRFinderPlus
`--organism` value — or into nothing at all, which is the common case here.

Why AMRFinderPlus has to be told which organism it is looking at
----------------------------------------------------------------
AMRFinderPlus only reports POINT MUTATIONS (gyrA/rpoB-style chromosomal
substitutions) when it is told which organism it is looking at, because every
curated organism carries its own hand-checked list of resistance-conferring
changes. Those point mutations are precisely the "intrinsic, chromosomal, not
transferable" category that a homology screen like ABRicate structurally cannot
see, and tier 1 of the mobilome module's mobility ladder (chromosomal /
intrinsic candidate) is exactly that category. So when the isolate IS one of the
curated organisms we want the flag; when it is not, we must run without it.

Why a GTDB name cannot be matched to that list naively
------------------------------------------------------
AMRFinderPlus 4.2.7 curates only 31 organisms, and it names them NCBI-style.
GTDB names do not always agree: GTDB splits genera and marks the split-off
lineages with an alphanumeric suffix (`Pseudomonas_E`, `Klebsiella_A`), and
gives genomes with no named species a placeholder epithet (`sp024807945`). A
naive substring match would happily turn `s__Pseudomonas_E sp010095445` into
`--organism Pseudomonas_aeruginosa`, and AMRFinderPlus would then score that
genome's proteins against P. aeruginosa's curated mutation list. That is worse
than no call at all, so every mismatch below is handled explicitly.

When a GTDB suffix means a different taxon, and when it does not
----------------------------------------------------------------
The suffix rule is not symmetric. When GTDB breaks up a genus, the species that
move out keep their epithet, so `s__Enterococcus_B faecium` really is NCBI's
Enterococcus faecium and does deserve the flag. Those cases live in
GTDB_SPECIES_EQUIVALENCES below — not hand-typed one bug report at a time, but
GENERATED from GTDB's own metadata by generate_gtdb_organism_table.py, which
applies the exact rule that originally justified the first few entries (see the
comment above that dict) uniformly to every GTDB species cluster, for every one
of AMRFinderPlus's 31 curated organisms. Likewise GENUS_ORGANISMS /
NAME_ONLY_ORGANISMS below are generated, not hand-picked: whether matching an
organism from its bare GTDB genus is safe is a question about genome counts,
answered once per GTDB release rather than guessed.

What goes in, what comes out, and where it sits in the workflow
---------------------------------------------------------------
Runs in the mobilome module's WP-A step, between GTDB-Tk (stage 03.taxonomy) and
AMRFinderPlus. Both the caller and AMRFinderPlus itself live in the mobilome
stage — rules amrfinder_organism and amrfinderplus in shared/80_mobilome.smk,
writing into 08.mobilome, not into the 05.amr stage where ABRicate and the CARD
read mapping sit. Pure text handling, standard library only and no AMRFinderPlus
call, so it unit-tests without any tool or database present.

Give exactly one of the first two:
  --gtdbtk-summary : the GTDB-Tk gtdbtk.*.summary.tsv written for this sample.
                     BacFlux runs GTDB-Tk into one directory per sample, so
                     every row in the file belongs to this isolate. Only two of
                     its columns are read, user_genome and classification, the
                     latter looking like
                     d__Bacteria;p__Pseudomonadota;...;g__Pseudomonas_E;s__Pseudomonas_E sp010095445
  --gtdbtk-dir     : the GTDB-Tk output DIRECTORY for this sample
                     (03.taxonomy/{sample}), inside which find_summary_in_dir
                     locates the summary. Convenient for a Snakemake rule that
                     takes the whole directory as input, which is what rule
                     amrfinder_organism does.
  --sample         : the BacFlux sample name, used to pick this sample's rows
                     (hybrid mode classifies both assemblies and writes
                     <sample>_illumina and <sample>_ont) and to label the audit
                     row.

And it writes two files:
  --out-organism   : one line holding the AMRFinderPlus organism name on a
                     match, or a completely EMPTY file on no match. Rule
                     amrfinderplus reads it, strips the whitespace, and passes
                     `--organism` (plus `--mutation_all`, which is meaningless
                     without it) only when something is left — so an
                     empty file means "run AMRFinderPlus without --organism",
                     the graceful, silent fallback the spec asks for.
  --out-audit      : a four-column TSV (sample, gtdb_classification,
                     matched_organism, reason) recording the decision and why it
                     was taken — the BacFlux convention that every
                     filtering/selection decision is auditable, as in
                     contig_taxonomy_decisions.tsv.

Why nothing here is ever an error
---------------------------------
The exit status is 0 in every ordinary situation, INCLUDING no match, no row for
the sample, an unreadable summary and a missing summary file. Environmental
isolates normally have no curated organism (none of the six isolates screened in
this project do), so "no match" is the normal path, not an error.
"""

import argparse
import os
import re
import sys


# ── The 31 curated AMRFinderPlus organisms ───────────────────────────────────
# Copied verbatim from `amrfinder --list_organisms -d {bakta_db}/amrfinderplus-db/latest`
# (version 4.2.7, recorded in docs/mobilome_wpA_ground_truth.md §"Spec §12 Q2").
# This is the ONE hand-verified fact in this whole file — it comes from running
# the tool, not from GTDB, so no metadata file can tell us this list. Every
# other set below is DERIVED from it, either by simple string shape (does the
# name have an underscore?) or by loading a small generated data file.
ALL_ORGANISMS = frozenset({
    "Acinetobacter_baumannii",
    "Bordetella_pertussis",
    "Burkholderia_cepacia",
    "Burkholderia_mallei",
    "Burkholderia_pseudomallei",
    "Campylobacter",
    "Citrobacter_freundii",
    "Clostridioides_difficile",
    "Corynebacterium_diphtheriae",
    "Enterobacter_asburiae",
    "Enterobacter_cloacae",
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Escherichia",
    "Haemophilus_influenzae",
    "Helicobacter_pylori",
    "Klebsiella_oxytoca",
    "Klebsiella_pneumoniae",
    "Neisseria_gonorrhoeae",
    "Neisseria_meningitidis",
    "Pseudomonas_aeruginosa",
    "Salmonella",
    "Serratia_marcescens",
    "Staphylococcus_aureus",
    "Staphylococcus_pseudintermedius",
    "Streptococcus_agalactiae",
    "Streptococcus_pneumoniae",
    "Streptococcus_pyogenes",
    "Vibrio_cholerae",
    "Vibrio_parahaemolyticus",
    "Vibrio_vulnificus",
})

# Curated per SPECIES: the genome must be that exact species, not merely that
# genus. AMRFinderPlus spells these Genus_epithet, underscore-joined, and that
# shape IS the fact that separates them from the genus-level ones below — a
# species-level organism always has an underscore, a genus-level one never does.
SPECIES_ORGANISMS = frozenset(name for name in ALL_ORGANISMS if "_" in name)

# The three organisms AMRFinderPlus curates at GENUS level (no underscore in the
# name at all: Campylobacter, Escherichia, Salmonella — always exactly these
# three, because that is a fact about AMRFinderPlus's own list, not about any
# one GTDB release). Whether matching one of them from its BARE GTDB genus is
# safe is a genome-counting question, answered by the generated genus-rules file
# below, not assumed here.
_GENUS_LEVEL_CANDIDATE_ORGANISMS = frozenset(name for name in ALL_ORGANISMS if "_" not in name)

# Where the two generated data files live: next to this script, so the loaders
# below find them regardless of the caller's working directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
GENUS_RULES_PATH = os.path.join(_HERE, "gtdb_organism_genus_rules.tsv")
EQUIVALENCES_PATH = os.path.join(_HERE, "gtdb_organism_equivalences.tsv")


def _read_generated_tsv(path):
    """Read one of the two generated data files into (header, list-of-dicts).

    Input:  a path written by generate_gtdb_organism_table.py: one leading line
            starting with '#' (provenance — which release, when, from which
            metadata file), then a normal tab-separated header and data rows.
    Output: (header, rows) where rows is a list of {column: value} dicts, in
            file order.

    Returns ([], []) for anything that stops us trusting the file — missing,
    empty, or a header that does not even parse. Both files are committed next to
    this script, so a missing one means the checkout is incomplete rather than
    fresh; either way the loaders below fall back to their own safe defaults and
    the AMRFinderPlus rule still runs, with the audit row saying no organism was
    matched. Losing point mutations is recoverable; killing the run for every
    sample over a data file is not.
    """
    if not os.path.isfile(path):
        return [], []
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    lines = [line for line in lines if not line.startswith("#")]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            continue  # a truncated/corrupted row is skipped, not fatal
        rows.append(dict(zip(header, fields)))
    return header, rows


def _load_genus_rules(path, candidates):
    """Load which of the genus-level-candidate organisms are safe to match from
    their bare GTDB genus alone.

    Input:  gtdb_organism_genus_rules.tsv (organism, genus_safe, ...), plus the
            fixed 3-member set of organisms this question is even asked about.
    Output: (genus_organisms, name_only_organisms) — two sets partitioning
            `candidates`. An organism the file does not mention, or a file that
            cannot be read at all, defaults to NAME_ONLY (unsafe): the module's
            whole stance is that a wrong --organism is worse than none, so a
            missing verdict must fail towards the more cautious answer, never
            towards silently allowing a genus-wide shortcut nobody has checked.
    """
    _header, rows = _read_generated_tsv(path)
    verdict_by_organism = {row.get("organism"): row.get("genus_safe") for row in rows}

    genus_organisms = set()
    name_only_organisms = set()
    for organism in candidates:
        if verdict_by_organism.get(organism) == "yes":
            genus_organisms.add(organism)
        else:
            name_only_organisms.add(organism)
    return frozenset(genus_organisms), frozenset(name_only_organisms)


def _load_species_equivalences(path, valid_organisms):
    """Load the per-species override table: {GTDB species name -> organism}.

    Input:  gtdb_organism_equivalences.tsv (gtdb_species, amrfinder_organism,
            ...), plus the full set of legal --organism strings.
    Output: a plain dict. Any row whose organism is NOT one of the 31 legal
            values is dropped with a warning rather than trusted — AMRFinderPlus
            exits with an error on an unrecognised --organism, so a stale or
            corrupted generated file must never be allowed to reach the point of
            emitting one. A missing file yields an empty dict: no overrides
            known, which is the safe, always-valid starting state.
    """
    _header, rows = _read_generated_tsv(path)
    equivalences = {}
    for row in rows:
        gtdb_species = row.get("gtdb_species", "")
        organism = row.get("amrfinder_organism", "")
        if not gtdb_species or organism not in valid_organisms:
            sys.stderr.write(
                "[gtdb_amrfinder_organism] WARNING: ignoring a row in %s with an "
                "unrecognised organism '%s' for GTDB species '%s'\n"
                % (path, organism, gtdb_species))
            continue
        equivalences[gtdb_species] = organism
    return equivalences


# ── The two generated tables, loaded once at import time ─────────────────────
# GENUS_ORGANISMS: safe to match from the bare GTDB genus alone (currently
#     Escherichia and Salmonella — see gtdb_organism_genus_rules.tsv for the
#     genome counts behind each verdict).
# NAME_ONLY_ORGANISMS: legal --organism values that must NEVER be picked from a
#     genus alone (currently just Campylobacter: its unsuffixed GTDB genus holds
#     none of the species it is curated for — see GTDB_SPECIES_EQUIVALENCES,
#     which is how its actual species reach it instead).
# GTDB_SPECIES_EQUIVALENCES: the per-species override table itself — see the
#     long comment below for what it means and how each entry earns its place.
GENUS_ORGANISMS, NAME_ONLY_ORGANISMS = _load_genus_rules(
    GENUS_RULES_PATH, _GENUS_LEVEL_CANDIDATE_ORGANISMS)

# Every genus AMRFinderPlus curates ANYTHING under. Used only to write a more
# honest "no match" reason: for Klebsiella we can say "the species is a GTDB
# placeholder", for Arthrobacter the truth is simply "nothing curated here".
CURATED_GENERA = ({name.split("_")[0] for name in SPECIES_ORGANISMS}
                  | GENUS_ORGANISMS | NAME_ONLY_ORGANISMS)

# Column names read from the GTDB-Tk summary. Stable across GTDB-Tk 2.x.
GENOME_COLUMN = "user_genome"
CLASSIFICATION_COLUMN = "classification"

# Header of the audit TSV, in order.
AUDIT_HEADER = ["sample", "gtdb_classification", "matched_organism", "reason"]

# GTDB marks a lineage it has split off from the NCBI name with an underscore
# plus capital letters at the END of the word: Pseudomonas_E, Klebsiella_A,
# Escherichia coli_D. The suffix means "related to, but NOT the same taxon as"
# the unsuffixed name, so by default it blocks a match.
#
# The one thing allowed to override it is an exact GTDB species name listed in
# the generated GTDB_SPECIES_EQUIVALENCES table below, and a suffix on the GENUS
# earns that override far more readily than one on the EPITHET. Which names, and
# on what evidence: see the comment above that table.
GTDB_SUFFIX = re.compile(r"_[A-Z]+$")

# GTDB gives genomes with no validly published species name a placeholder
# epithet: the literal "sp" followed by the GTDB accession digits, e.g.
# "s__Pseudomonas_E sp010095445". It is not a species name, so it can never
# match a curated species.
GTDB_PLACEHOLDER_SPECIES = re.compile(r"^sp\d*$")

# ── GTDB/NCBI name differences: GTDB moved the name, but it is the same species ─
# The suffix rule above is right about lineages GTDB genuinely split off
# (Pseudomonas_E is not NCBI Pseudomonas), but it is too strict in some real
# cases. When GTDB breaks up a genus that NCBI keeps whole, the species that
# move out keep their own epithet, so "s__Enterococcus_B faecium" is exactly
# NCBI's Enterococcus faecium, just filed under a different GTDB genus. GTDB
# assigns that genus suffix by where the genus TYPE SPECIES landed — E. faecalis
# is the type species of Enterococcus, so it keeps the unsuffixed genus and
# faecium was pushed into _B. The same mechanism explains Campylobacter: C.
# fetus is the type species, so C. jejuni and C. coli — the food-chain
# pathogens AMRFinderPlus actually curates — sit in g__Campylobacter_D, not in
# the unsuffixed genus at all. That is why Campylobacter is in
# NAME_ONLY_ORGANISMS: its bare GTDB genus contains none of what it is curated
# for, and the species reach it through the equivalences table instead.
#
# A suffix on the EPITHET means something different and stronger: GTDB itself
# is saying it cannot safely attach that name here (usually because the type
# strain was never sequenced), so an epithet suffix is only overridden on
# overwhelming genome-count evidence, never assumed — Helicobacter pylori_C, for
# instance: 740 genomes, 99.9% of them NCBI Helicobacter pylori.
#
# GTDB_SPECIES_EQUIVALENCES is therefore NOT a hand-typed list. It is GENERATED
# by generate_gtdb_organism_table.py, which reads a GTDB release's own metadata
# (every genome's GTDB name AND its NCBI name, side by side) and applies exactly
# the two rules above — mechanically, to every GTDB species cluster, for all 31
# curated organisms — rather than one bug report at a time. Each row carries the
# rule that produced it (unsuffixed_epithet or suffixed_epithet_strong_evidence),
# the NCBI species it agrees with, and the genome count behind the agreement. See
# that script's docstring for the full rule, and check_gtdb_organism_table.py for
# how to tell whether the generated table still holds after a GTDB release bump,
# before spending the time to regenerate it.
#
# The data itself lives in gtdb_organism_equivalences.tsv, next to this script,
# committed to the repo so BacFlux never needs the multi-hundred-MB GTDB
# metadata file at runtime — only when regenerating the table.
GTDB_SPECIES_EQUIVALENCES = _load_species_equivalences(EQUIVALENCES_PATH, ALL_ORGANISMS)


# ── Reading the GTDB-Tk summary ──────────────────────────────────────────────

def read_gtdbtk_summary(summary_path):
    """Pull the (genome name, classification) pairs out of a GTDB-Tk summary TSV.

    Input:  path to `gtdbtk.bac120.summary.tsv` (or the ar53 equivalent) as
            written by rule taxonomic_assignment (shared/30_taxonomy.smk) into
            03.taxonomy/{sample}/.
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

    Input:  03.taxonomy/{sample} as produced by rule taxonomic_assignment
            (shared/30_taxonomy.smk).
    Output: (summary path, problem) — the path of the summary to read, or an
            empty path plus a message when there is none.

    GTDB-Tk writes `gtdbtk.bac120.summary.tsv` (bacteria) and/or
    `gtdbtk.ar53.summary.tsv` (archaea) into classify/, and classify_wf also
    leaves a symlink to it at the top level, so both places are worth trying.
    This looks at the top level first and then in classify/; the MultiQC staging
    step in shared/90_report.smk does the same search in the OPPOSITE order.
    Either way it lands on the same summary, since the top-level entry points at
    the classify/ one. BacFlux is a bacterial workflow, so within a directory the
    bac120 file is tried first; the ar53 file is still accepted as a fallback so
    an archaeal isolate does not silently produce nothing at all.
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
    unsuffixed NCBI one, so it blocks the plain species and genus matches in
    map_classification. The only thing that gets past it is an exact name listed
    in the generated GTDB_SPECIES_EQUIVALENCES table, which is checked first.
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

    The order of the four steps carries the logic and cannot be shuffled. The
    generated equivalences table is consulted FIRST, because the names in it are
    exactly the ones the suffix rule in step 2 would otherwise refuse; only then
    does the plain species match run, then the genus match, and last the "no
    match" branch, which exists purely to write a reason worth reading later.
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

    # The species name spelled exactly as GTDB writes it ("Enterococcus_B
    # faecium"). Used for the generated-table lookup in step 1 and again in the
    # "why not" reasons at the bottom, so the audit quotes the real GTDB name.
    gtdb_binomial = (genus + " " + epithet).strip()

    # Step 1 — a known GTDB/NCBI naming artefact for the SAME species.
    # Checked before anything else, because these are exactly the names where
    # step 2's suffix rule would give the wrong answer: GTDB has moved the
    # species into a split-off genus (or split the epithet itself), but it is
    # still the same organism AMRFinderPlus curates. See
    # GTDB_SPECIES_EQUIVALENCES for how each entry is generated and verified.
    if gtdb_binomial in GTDB_SPECIES_EQUIVALENCES:
        equivalent_organism = GTDB_SPECIES_EQUIVALENCES[gtdb_binomial]
        return equivalent_organism, (
            "known GTDB/NCBI naming artefact (GTDB s__%s is the same species as "
            "AMRFinderPlus organism '%s'; see gtdb_organism_equivalences.tsv)"
            % (gtdb_binomial, equivalent_organism))

    # Step 2 — exact species match.
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

    # Step 3 — genus-level match.
    # AMRFinderPlus curates three organisms for a whole genus, but only the ones
    # in GENUS_ORGANISMS reach this branch: the generated genus-rules file counts
    # how much of each GTDB genus really is that organism, and Campylobacter fails
    # that count outright (its unsuffixed GTDB genus holds none of the curated
    # species, which arrive through step 1 instead). For the genera that do pass,
    # any species qualifies — including unnamed ones, and including Shigella,
    # which GTDB has already folded into g__Escherichia. A GTDB-suffixed genus is
    # still excluded here: a lineage GTDB split off is not the curated organism,
    # whatever the base name says.
    if genus and not genus_is_suffixed and genus in GENUS_ORGANISMS:
        return genus, ("genus-level match (AMRFinderPlus curates '%s' at genus "
                       "level; GTDB g__%s)" % (genus, genus))

    # Step 4 — no match. Say which of the GTDB/NCBI mismatches caused it, so the
    # audit row is useful six months later. The generic reason comes first
    # because for most environmental isolates (Paenibacillus, Arthrobacter) the
    # genus simply has nothing curated, and the suffix/placeholder detail would
    # be a red herring.
    if genus_base not in CURATED_GENERA:
        return "", ("no curated organism for this taxon (AMRFinderPlus curates 31 "
                    "organisms; genus '%s' is not one of them)" % genus)

    # A suffixed genus is normally a lineage GTDB split off from the NCBI genus,
    # so it is not the same taxon. Since step 1 accepts the generated exceptions
    # to that, the reason names where they are recorded, so a reader who sees
    # one Enterococcus_B accepted and another refused can tell why.
    if genus_is_suffixed:
        return "", ("no curated organism for this taxon: GTDB-suffixed genus '%s' "
                    "is a lineage split off from NCBI '%s', and '%s' is not one of "
                    "the known GTDB/NCBI naming artefacts (see "
                    "gtdb_organism_equivalences.tsv)" % (genus, genus_base, gtdb_binomial))

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

    # The rows disagree about the organism — refuse to guess.
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

    The organism file is deliberately EMPTY (zero bytes) when there is no match.
    That is the whole signalling mechanism: rule amrfinderplus reads the file,
    strips whitespace (a match is written with a trailing newline) and tests
    whether anything is left, so no exit code, marker word or second file is
    needed to say "no organism".

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

    # One line summarising the choice. Rule amrfinder_organism redirects stdout
    # into logs/mobilome_amrfinder_organism_{sample}.log, so the decision is
    # visible there without opening the audit TSV.
    if organism:
        sys.stdout.write("%s: AMRFinderPlus --organism %s (%s)\n"
                         % (args.sample, organism, reason))
    else:
        sys.stdout.write("%s: no AMRFinderPlus --organism (%s)\n"
                         % (args.sample, reason))
    return 0


if __name__ == "__main__":
    sys.exit(main())
