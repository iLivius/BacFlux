#!/usr/bin/env python3
"""Build the complete GTDB -> AMRFinderPlus organism table for one GTDB release.

Why this exists
    gtdb_amrfinder_organism.py used to carry a short, HAND-TYPED list of GTDB
    names needing special handling (5 entries, covering Campylobacter and
    Enterococcus faecium, found one bug report at a time). Running
    check_gtdb_organism_table.py showed that hand-typing covers only a
    fraction of the real cases: roughly half of AMRFinderPlus's 31 curated
    organisms have at least one GTDB naming quirk in a real release, from
    Helicobacter pylori (dozens of suffixed clusters, thousands of genomes) down
    to single-genome edge cases nobody would ever notice by eye.

    This script replaces "type entries in by hand" with "read GTDB's own
    metadata once, and mechanically apply the exact rule that was already used,
    by hand, to justify the first five entries." It writes two small,
    git-tracked data files that gtdb_amrfinder_organism.py loads at import time:
        gtdb_organism_equivalences.tsv   - the per-species override table
        gtdb_organism_genus_rules.tsv    - which genus-level organisms are safe
                                            to match by genus alone
    Nothing here changes gtdb_amrfinder_organism.py's DECISION LOGIC (the order
    of checks in map_classification is untouched) - only the DATA it reads.
    Re-run this whenever gtdbtk_db is pointed at a new GTDB release; see
    docs/README_notes.md item 11. Both scripts read the same metadata file, so
    run check_gtdb_organism_table.py against it first: when it reports every
    entry OK and no missing candidates, there is nothing here to regenerate.

THE RULE, IN ONE PARAGRAPH
    A GTDB species name has the shape "Genus[_suffix] epithet[_suffix]". A
    suffix on the GENUS means GTDB reshuffled the genus tree - it is a
    naming artefact, not a claim that the SPECIES is different, so once the
    genomes agree with an NCBI species name we accept it outright (RULE A). A
    suffix on the EPITHET means something stronger: GTDB itself is saying it
    cannot safely attach that name to this cluster (usually because the type
    strain was never sequenced), so it is accepted only when the genome-count
    evidence is overwhelming - at least 99% of a cluster of at least 20 genomes
    agreeing with one NCBI name (RULE B). This is exactly how the first five
    entries in the old hand-typed table were justified; here it is applied
    uniformly to every GTDB species cluster in the release, for all 31 curated
    organisms, instead of to five hand-picked ones.

    A cluster the plain rules in gtdb_amrfinder_organism.py would already match
    WITHOUT an override (neither token suffixed, or the genus is unsuffixed AND
    already judged safe to match on its own) is left OUT of the table - the
    table's whole job is to cover the exceptions, not restate the obvious cases.

REUSES, RATHER THAN DUPLICATES
    The fact-gathering - reading the metadata, counting NCBI names inside each
    GTDB cluster, working out which NCBI species an AMRFinderPlus organism
    covers - is exactly what check_gtdb_organism_table.py already does, so this
    script imports those functions rather than re-implementing them. Think of
    the checker as "does the shipped table still hold" and this script as
    "build the whole table from first principles."

INPUT
    --metadata   bac120_metadata_r<N>.tsv.gz, the same file
                 check_gtdb_organism_table.py takes. A few hundred MB; needed
                 only to REGENERATE the table, never at BacFlux runtime.

OUTPUT
    --out-equivalences   gtdb_organism_equivalences.tsv (the per-species table)
    --out-genus-rules    gtdb_organism_genus_rules.tsv (the 3-row genus verdicts)
    Both default to the files gtdb_amrfinder_organism.py actually reads, so the
    ordinary invocation only needs --metadata.

RUN
    python workflow/scripts/80_mobilome/generate_gtdb_organism_table.py \
        --metadata /path/to/bac120_metadata_r232.tsv.gz
    Then re-run check_gtdb_organism_table.py against the SAME metadata file to
    confirm the freshly written table reports everything OK, then commit both
    files. (R232 in the example is the release the committed tables were built
    from; the first run of this script was against R226.)
"""

import argparse
import csv
import datetime
import os
import re
import sys

# Both sibling scripts live next to this one. They are plain files run by path,
# not an installed package, so their directory has to go on the import path by
# hand before either can be imported.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_gtdb_organism_table as gtdb_check
import gtdb_amrfinder_organism as gao


# ── The evidence bars ────────────────────────────────────────────────────────
# Same evidence bar the hand-built table used for a suffixed epithet, imported
# BY NAME from the checker (not re-typed) so the two scripts can never quietly
# drift onto different numbers.
MIN_AGREEMENT_FOR_SUFFIXED_EPITHET = gtdb_check.MIN_AGREEMENT_FOR_SUFFIXED_EPITHET
MIN_GENOMES_FOR_SUFFIXED_EPITHET = gtdb_check.MIN_GENOMES_FOR_SUFFIXED_EPITHET

# How much of an unsuffixed genus has to agree with the curated species before
# genus-alone matching is judged safe. Same bar check_gtdb_organism_table.py
# already used (95.0, inline there); named here so this script's own logic is
# self-contained and readable without cross-referencing the checker's source.
# This is the one number NOT shared by import, so changing it here alone makes
# the checker start calling a freshly generated table unsafe.
MIN_GENUS_SAFETY_PERCENT = 95.0


# ── The two files this script writes ─────────────────────────────────────────
# gtdb_amrfinder_organism.py loads both at import time. Generation writes here
# by default so the ordinary workflow is just "run this script, commit the two
# files it touches."
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EQUIVALENCES_PATH = os.path.join(_HERE, "gtdb_organism_equivalences.tsv")
DEFAULT_GENUS_RULES_PATH = os.path.join(_HERE, "gtdb_organism_genus_rules.tsv")

# gtdb_organism_equivalences.tsv, one row per GTDB species that needs an
# override. Only the first two columns are read back by the workflow; the rest
# is the evidence, kept so a human can see what the row rests on:
#   gtdb_species       : GTDB name exactly as GTDB-Tk writes it ("Campylobacter_D coli")
#                        source: the s__ field of the release's gtdb_taxonomy
#   amrfinder_organism : the --organism value to pass for it
#                        source: build_ncbi_to_organism_map, keyed on ncbi_species
#   ncbi_species       : what the cluster's genomes are called at NCBI
#                        source: majority vote over ncbi_taxonomy (majority_ncbi_species)
#   percent_agreement  : share of the cluster carrying that NCBI name
#   n_genomes          : genomes in the cluster
#   rule               : unsuffixed_epithet (RULE A) | suffixed_epithet_strong_evidence (RULE B)
EQUIVALENCES_HEADER = ["gtdb_species", "amrfinder_organism", "ncbi_species",
                       "percent_agreement", "n_genomes", "rule"]

# gtdb_organism_genus_rules.tsv, exactly three rows - one per AMRFinderPlus
# organism curated at genus level. Again only the first two columns are read
# back:
#   organism            : Campylobacter | Escherichia | Salmonella
#   genus_safe          : yes | no - may this organism be matched from the bare
#                         g__ field, with no species check at all?
#   percent_appropriate : share of the genomes in that unsuffixed genus that
#                         really are one of its curated species
#   n_appropriate       : genomes behind that share
#   n_total             : genomes GTDB files under the unsuffixed genus
#                         source (all four): decide_genus_safety
GENUS_RULES_HEADER = ["organism", "genus_safe", "percent_appropriate",
                     "n_appropriate", "n_total"]

# Matches the release label in a metadata filename: bac120_metadata_r232.tsv.gz
# -> "r232". Used for the provenance line only, never for a decision.
RELEASE_FROM_FILENAME = re.compile(r"_(r\d+)[._]")


# ── Which release did this table come from? ──────────────────────────────────

def guess_release(metadata_path):
    """Pull a release label ("r226") out of the metadata filename for the
    provenance line, so a human reading the generated file later knows, without
    opening this script, which GTDB release it came from. Falls back to
    "unknown release" rather than guessing wrong."""
    match = RELEASE_FROM_FILENAME.search(os.path.basename(metadata_path))
    return match.group(1) if match else "unknown release"


# ── Which AMRFinderPlus organism does an NCBI species belong to? ─────────────

def build_ncbi_to_organism_map():
    """{NCBI species name -> AMRFinderPlus organism}, for all 31 curated
    organisms.

    Reuses target_species_for from the checker, which already knows the
    special cases (Campylobacter -> jejuni/coli only; Escherichia -> coli plus
    the four Shigella species; everything else derived straight from the
    organism string). Built once and used for every GTDB cluster below, so a
    cluster's NCBI majority name is turned into "which organism, if any" by a
    single dictionary lookup rather than a name-by-name comparison.
    """
    ncbi_to_organism = {}
    for organism in gao.ALL_ORGANISMS:
        for ncbi_species in gtdb_check.target_species_for(organism):
            ncbi_to_organism[ncbi_species] = organism
    return ncbi_to_organism


# ── Genus-level verdicts: may this organism be matched from the genus alone? ─

def decide_genus_safety(organism, counts, genus_of_species):
    """Is matching AMRFinderPlus organism `organism` from its bare GTDB genus,
    with no further checking, safe?

    Input:  the organism name (only meaningful for one with NO underscore -
            AMRFinderPlus's three genus-level taxgroups), plus the per-cluster
            NCBI-name counts and each cluster's literal g__ genus field, both
            from check_gtdb_organism_table.read_gtdb_metadata.
    Does:   sums genome counts across every GTDB species cluster whose g__
            field is EXACTLY this organism's name (no suffix at all), splitting
            them into "resolves to one of this organism's curated species" and
            "resolves to something else". Only the LITERAL unsuffixed spelling
            counts here - a genus-suffixed cluster is judged separately, per
            species, in build_equivalence_rows.
    Output: (is_safe, percent_appropriate, n_appropriate, n_total). is_safe is
            True only when at least MIN_GENUS_SAFETY_PERCENT of the genomes
            GTDB placed in that bare genus are the curated species - Escherichia
            and Salmonella clear this bar; Campylobacter does not, because its
            curated species (jejuni, coli) do not even live in the unsuffixed
            genus in the first place. Campylobacter reaches its organism through
            the per-species overrides instead (g__Campylobacter_D), which is why
            a "no" here is not a gap.
    """
    wanted = gtdb_check.target_species_for(organism)
    appropriate = inappropriate = 0
    for gtdb_species, genus in genus_of_species.items():
        if genus != organism:
            continue
        ncbi_name, _percent, total = gtdb_check.majority_ncbi_species(counts[gtdb_species])
        if ncbi_name in wanted:
            appropriate += total
        else:
            inappropriate += total
    total_seen = appropriate + inappropriate
    percent = (100.0 * appropriate / total_seen) if total_seen else 0.0
    return percent >= MIN_GENUS_SAFETY_PERCENT, percent, appropriate, total_seen


def build_genus_rules(counts, genus_of_species):
    """One row per genus-level-candidate organism - the ones AMRFinderPlus
    names with no underscore (Campylobacter, Escherichia, Salmonella; there are
    only ever these three, since that is a fact about AMRFinderPlus's own
    curated list, not about any one GTDB release).

    Rows carry the GENUS_RULES_HEADER columns and are written to
    gtdb_organism_genus_rules.tsv. Every candidate gets a row, "no" verdicts
    included: an organism missing from the file is read back as unsafe anyway,
    so a "no" row is there for the reader, to show the question was asked and
    with what numbers.
    """
    candidates = sorted(organism for organism in gao.ALL_ORGANISMS if "_" not in organism)
    rows = []
    for organism in candidates:
        safe, percent, appropriate, total = decide_genus_safety(
            organism, counts, genus_of_species)
        rows.append({
            "organism": organism,
            "genus_safe": "yes" if safe else "no",
            "percent_appropriate": "%.1f" % percent,
            "n_appropriate": str(appropriate),
            "n_total": str(total),
        })
    return rows


# ── Per-species overrides: RULE A and RULE B, applied to every cluster ───────

def build_equivalence_rows(counts, genus_safe_organisms):
    """One row per GTDB species cluster that needs an EXPLICIT override entry.

    Input:  the per-cluster NCBI-name counts (from read_gtdb_metadata), plus
            the set of organisms build_genus_rules just judged safe to match by
            genus alone.
    Does, for every GTDB species cluster in the release:
      1. Find its NCBI majority name (>=50% of the cluster, else no name at all)
         and look up which curated organism (if any) that NCBI name belongs to.
         A cluster whose majority is not one of the 31 curated species is
         skipped immediately - most GTDB clusters are environmental bacteria
         with nothing curated for them at all.
      2. Reject it if the GTDB epithet, once its own suffix is stripped, is NOT
         literally the epithet of the matched NCBI name. This is the guard
         against a genuinely DIFFERENT species that merely happens to carry the
         curated NCBI label on some of its genomes. It is not hypothetical: the
         first run of this script, against R226, wrote four wrong rows in two
         shapes; adding the guard and re-running removed exactly those four
         rows and nothing else.
           One mislabelled deposit is enough. Arthrobacter_D sp009728235 (n=1)
             -> NCBI says "Vibrio cholerae" for that ONE deposited genome - a
             submission mislabelled at NCBI, not a real Vibrio. The same thing
             happened to Terrimonas_A sp003243445 (-> Citrobacter_freundii) and
             Phascolarctobacterium_A sp963603005 (-> Escherichia). A placeholder
             epithet ("sp<digits>") is never a real species name, so
             is_placeholder_species is checked directly on top of the
             epithet-agreement test, for clarity as much as defence: an epithet
             that IS a placeholder always fails the agreement test too (it
             cannot equal a real curated epithet), but saying so explicitly
             makes the intent obvious to a reader.
           A real species NCBI has not caught up with. Enterococcus_B lactis
             (n=692) -> only 55.1% of that cluster is labelled "Enterococcus
             faecium" at NCBI, evidently many historical submissions predating
             E. lactis being recognised as its own species. E. lactis IS its own
             GTDB species (a real, clinically distinct one) and must never be
             reported as E. faecium, however many old NCBI records blur the two.
         check_gtdb_organism_table.py has carried this same guard since the
         Serratia sarumanii episode; it was just not carried across when this
         script was written fresh, which is how all four rows got through.
      3. Skip it again if the PLAIN rules in gtdb_amrfinder_organism.py would
         already match it with no override: neither the genus nor the epithet
         token carries a GTDB suffix (the ordinary exact-species rule handles
         it), or the genus is unsuffixed and already judged genus-safe (the
         genus-level rule handles it regardless of the epithet).
      4. Otherwise decide RULE A vs RULE B by where the remaining suffix sits:
         epithet unsuffixed -> RULE A, accepted on agreement alone (GTDB's own
         type-strain convention is the real justification; the vote is
         corroboration). Epithet suffixed -> RULE B, accepted only above the
         percent/genome-count bar.
    Output: a list of dict rows, one per accepted override, unsorted (main
            sorts by GTDB name before writing them to
            gtdb_organism_equivalences.tsv). Every row here becomes a lookup
            that map_classification consults FIRST, ahead of its own suffix
            rules - which is exactly why the guards above matter.
    """
    ncbi_to_organism = build_ncbi_to_organism_map()
    rows = []

    for gtdb_species, ncbi_counter in counts.items():
        if " " not in gtdb_species:
            continue  # GTDB-Tk left the species field blank - nothing to key an override on

        ncbi_name, percent, total = gtdb_check.majority_ncbi_species(ncbi_counter)
        organism = ncbi_to_organism.get(ncbi_name)
        if organism is None:
            continue  # this cluster's NCBI name is not one of the 31 we curate for

        genus_token, epithet_token = gtdb_species.split(" ", 1)
        if gao.is_placeholder_species(epithet_token):
            continue  # "sp<digits>" is never a real species name (see note above)

        _genus_base, genus_is_suffixed = gao.strip_gtdb_suffix(genus_token)
        epithet_base, epithet_is_suffixed = gao.strip_gtdb_suffix(epithet_token)

        # The GTDB epithet (suffix stripped) must be the SAME epithet as the
        # matched NCBI name - not merely a cluster whose majority vote happens
        # to land on a curated name. See the note above for the two shapes of
        # real counter-example this catches.
        ncbi_epithet = ncbi_name.split(" ", 1)[1] if " " in ncbi_name else ncbi_name
        if epithet_base != ncbi_epithet:
            continue

        covered_for_free = (
            (not genus_is_suffixed and not epithet_is_suffixed
             and organism in gao.SPECIES_ORGANISMS)
            or (not genus_is_suffixed and organism in genus_safe_organisms)
        )
        if covered_for_free:
            continue

        if epithet_is_suffixed:
            # RULE B - GTDB itself says this name's application is uncertain,
            # so only overwhelming genome-count evidence overrides that caution.
            if percent < MIN_AGREEMENT_FOR_SUFFIXED_EPITHET or total < MIN_GENOMES_FOR_SUFFIXED_EPITHET:
                continue
            rule = "suffixed_epithet_strong_evidence"
        else:
            # RULE A - the epithet is unsuffixed, so by GTDB's own nomenclature
            # rule this cluster holds the TYPE material for that name: the
            # genus reshuffle is a naming artefact, not a different species.
            # We only reached here because ncbi_name already matched the
            # organism's target species, so the majority-vote agreement is
            # corroboration on top of that, not the primary justification.
            rule = "unsuffixed_epithet"

        rows.append({
            "gtdb_species": gtdb_species,
            "amrfinder_organism": organism,
            "ncbi_species": ncbi_name,
            "percent_agreement": "%.1f" % percent,
            "n_genomes": str(total),
            "rule": rule,
        })

    return rows


# ── Write the two generated files ────────────────────────────────────────────

def write_tsv(path, header, rows, provenance):
    """Write one generated TSV: a single '#'-prefixed provenance line, then a
    normal header row, then one row per dict in `rows` (columns in `header`
    order). The loader in gtdb_amrfinder_organism.py knows to skip the leading
    '#' line; nothing else needs to.

    The provenance line is the only record of which GTDB release a committed
    table came from - the rows themselves carry no release stamp - so it names
    the metadata file, the release and the date, and says not to hand-edit.
    """
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("# %s\n" % provenance)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row[column] for column in header])


# ── Read the metadata, rebuild both tables, say what was written ─────────────

def main(argv=None):
    """Read one release's metadata and rewrite both data files.

    The order of the two builds matters: the genus verdicts come first, because
    a cluster sitting in a genus already judged safe needs no per-species
    override, and build_equivalence_rows can only skip those if it is handed the
    set of safe organisms.

    Returns 0 in every case; the one thing it refuses to continue past is a
    missing metadata file. Nothing here checks its own output - that is
    check_gtdb_organism_table.py's job, which the closing message asks for.
    """
    parser = argparse.ArgumentParser(
        description="Build the full GTDB->AMRFinderPlus organism table from a "
                    "GTDB release's own metadata. Overwrites the two data files "
                    "gtdb_amrfinder_organism.py loads.")
    parser.add_argument("--metadata", required=True,
                        help="bac120_metadata_r<N>.tsv.gz for the release to "
                             "generate the table from.")
    parser.add_argument("--out-equivalences", default=DEFAULT_EQUIVALENCES_PATH,
                        help="where to write the per-species override table "
                             "(default: the file the runtime script reads)")
    parser.add_argument("--out-genus-rules", default=DEFAULT_GENUS_RULES_PATH,
                        help="where to write the genus-safety verdicts "
                             "(default: the file the runtime script reads)")
    parser.add_argument("--gtdb-release", default="",
                        help="release label for the provenance line (e.g. "
                             "r226); guessed from the metadata filename if omitted")
    args = parser.parse_args(argv)

    if not os.path.isfile(args.metadata):
        sys.exit("[generate_gtdb_organism_table] no such file: %s" % args.metadata)

    release = args.gtdb_release or guess_release(args.metadata)
    print("Reading %s (GTDB release %s)..." % (args.metadata, release))
    counts, genus_of_species = gtdb_check.read_gtdb_metadata(args.metadata)
    print("Read %d GTDB species clusters." % len(counts))

    genus_rules = build_genus_rules(counts, genus_of_species)
    genus_safe_organisms = {row["organism"] for row in genus_rules
                            if row["genus_safe"] == "yes"}
    print("Genus-level organisms judged safe to match by genus alone: %s"
          % (sorted(genus_safe_organisms) or "none"))

    equivalence_rows = build_equivalence_rows(counts, genus_safe_organisms)
    equivalence_rows.sort(key=lambda row: row["gtdb_species"])
    print("Built %d species-equivalence entries." % len(equivalence_rows))

    today = datetime.date.today().isoformat()
    genus_provenance = (
        "Generated by generate_gtdb_organism_table.py from %s (GTDB release %s) "
        "on %s. Do not hand-edit; re-run the generator instead."
        % (os.path.basename(args.metadata), release, today))
    equivalences_provenance = (
        "Generated by generate_gtdb_organism_table.py from %s (GTDB release %s) "
        "on %s. Do not hand-edit; re-run the generator instead."
        % (os.path.basename(args.metadata), release, today))

    write_tsv(args.out_genus_rules, GENUS_RULES_HEADER, genus_rules, genus_provenance)
    write_tsv(args.out_equivalences, EQUIVALENCES_HEADER, equivalence_rows,
              equivalences_provenance)

    print("\nWrote %s" % args.out_genus_rules)
    print("Wrote %s" % args.out_equivalences)
    print("\nRun check_gtdb_organism_table.py again against the SAME metadata "
          "file - it should now report every entry OK and no missing candidates.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
