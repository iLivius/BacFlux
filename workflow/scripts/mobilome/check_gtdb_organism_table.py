#!/usr/bin/env python3
"""Re-check the GTDB -> AMRFinderPlus organism table against a GTDB release.

WHY THIS EXISTS
    gtdb_amrfinder_organism.py carries a small hand-checked table saying which
    GTDB species names mean which AMRFinderPlus `--organism`. Those names are NOT
    stable: GTDB explicitly states that the alphabetic suffixes are best-effort
    and may change between releases, and it moves species between genera when the
    tree says so. C. jejuni sits in g__Campylobacter_D today; nothing promises it
    will tomorrow.

    So every time you point `gtdbtk_db` at a NEW GTDB release, the table needs
    re-checking. This script does the checking for you. It does NOT edit the
    table - a wrong --organism is worse than none, so the decision stays with a
    human - it tells you exactly which entries still hold, which have moved, and
    which new ones you could now add.

WHAT IT CHECKS, AND WHAT "CORRECT" MEANS
    GTDB's own rule is that the species cluster holding the nomenclatural TYPE
    keeps the unsuffixed name; suffixed names are placeholders. That is why a
    genus suffix (Campylobacter_D jejuni) says nothing about species identity
    while an epithet suffix (jejuni_C) means "we cannot safely attach this name
    here". The evidence we use for both is the same: GTDB's metadata records, for
    every genome, both its GTDB name and its NCBI name. Counting the NCBI names
    inside one GTDB species cluster tells you what that cluster IS in NCBI terms.
    That is the same majority-vote idea as GTDB-Tk's own
    gtdb_to_ncbi_majority_vote.py, done here on just the handful of species we
    care about.

    Three questions are asked:
      1. Does every entry in GTDB_SPECIES_EQUIVALENCES still point at the right
         organism in this release?
      2. Are the genus-level rules (GENUS_ORGANISMS) still safe - i.e. is nearly
         everything GTDB files under that unsuffixed genus really that organism?
      3. Are there NEW GTDB species that now resolve to a curated organism but
         are missing from the table?

INPUT
    --metadata   bac120_metadata_r<release>.tsv.gz, downloaded from
                 https://data.gtdb.ecogenomic.org/releases/release<N>/<N>.0/
                 It is a few hundred MB and is NOT needed to run BacFlux - only
                 to run this check when upgrading GTDB.

OUTPUT
    A report on stdout. Exit status 0 when nothing needs attention, 1 when
    something does, so it can be wired into CI if you ever want that.

RUN
    python workflow/scripts/mobilome/check_gtdb_organism_table.py \
        --metadata /path/to/bac120_metadata_r232.tsv.gz
"""

import argparse
import collections
import csv
import gzip
import os
import sys

# The table under test lives next to this script.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gtdb_amrfinder_organism as gao


# A suffixed EPITHET is admitted only on strong evidence: GTDB is saying it
# cannot safely attach the name, so we need the genome counts to say otherwise
# very clearly. These two numbers are the bar used when the table was built.
MIN_AGREEMENT_FOR_SUFFIXED_EPITHET = 99.0   # percent of the cluster
MIN_GENOMES_FOR_SUFFIXED_EPITHET = 20       # below this, one bad label flips it

# Which NCBI species each AMRFinderPlus organism is meant for. Needed because the
# organism string is not always an NCBI binomial: "Campylobacter" is a genus-level
# taxgroup that NCBI documents as covering C. jejuni and C. coli ONLY, and
# "Escherichia" covers E. coli plus Shigella. Everything else is derived from the
# organism name itself, so this map only lists the awkward ones.
ORGANISM_TARGET_SPECIES = {
    "Campylobacter": {"Campylobacter jejuni", "Campylobacter coli"},
    "Escherichia": {"Escherichia coli", "Shigella sonnei", "Shigella flexneri",
                    "Shigella dysenteriae", "Shigella boydii"},
    "Salmonella": {"Salmonella enterica", "Salmonella bongori"},
}


def target_species_for(organism):
    """NCBI species names an AMRFinderPlus organism is curated for.

    For a species-level organism the name IS the binomial with an underscore
    ("Enterococcus_faecium" -> "Enterococcus faecium"). The genus-level ones are
    listed explicitly above because their scope is not derivable from the string.
    """
    if organism in ORGANISM_TARGET_SPECIES:
        return ORGANISM_TARGET_SPECIES[organism]
    return {organism.replace("_", " ")}


def read_gtdb_metadata(metadata_path):
    """Count NCBI species names inside each GTDB species cluster.

    Input:  bac120_metadata_r<N>.tsv.gz - one row per genome, with both a
            gtdb_taxonomy and an ncbi_taxonomy lineage string.
    Output: {gtdb_species: Counter({ncbi_species: n_genomes})}, plus
            {gtdb_species: gtdb_genus} so callers can ask about genera too.

    Both names are taken as the last (s__) field of their lineage. A genome with
    no NCBI name contributes an empty string, which is counted rather than
    dropped: a cluster that is mostly unnamed at NCBI is itself a useful signal.
    """
    ncbi_names_by_gtdb_species = collections.defaultdict(collections.Counter)
    genus_of_species = {}

    opener = gzip.open if metadata_path.endswith(".gz") else open
    with opener(metadata_path, "rt", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gtdb_index = header.index("gtdb_taxonomy")
        ncbi_index = header.index("ncbi_taxonomy")

        for row in reader:
            if len(row) <= max(gtdb_index, ncbi_index):
                continue
            gtdb_lineage = row[gtdb_index].split(";")
            gtdb_species = gtdb_lineage[-1].replace("s__", "")
            if not gtdb_species:
                continue
            gtdb_genus = ""
            for rank in gtdb_lineage:
                if rank.startswith("g__"):
                    gtdb_genus = rank[3:]
            ncbi_species = row[ncbi_index].split(";")[-1].replace("s__", "")
            ncbi_names_by_gtdb_species[gtdb_species][ncbi_species] += 1
            genus_of_species[gtdb_species] = gtdb_genus

    return ncbi_names_by_gtdb_species, genus_of_species


def majority_ncbi_species(ncbi_counter):
    """The NCBI name for a GTDB cluster -> (name, percent, n_genomes).

    Returns ("", 0.0, n) when no single name reaches half the cluster, which is
    itself an answer: GTDB's cluster does not correspond to one NCBI species.
    """
    total = sum(ncbi_counter.values())
    if not total:
        return "", 0.0, 0
    name, count = ncbi_counter.most_common(1)[0]
    percent = 100.0 * count / total
    if percent < 50.0 or not name:
        return "", percent, total
    return name, percent, total


def check_equivalence_entries(counts, problems):
    """Question 1: does each hand-checked entry still hold in this release?"""
    print("\n1. Hand-checked species equivalences")
    print("   (GTDB species -> AMRFinderPlus --organism)\n")

    for gtdb_species, organism in sorted(gao.GTDB_SPECIES_EQUIVALENCES.items()):
        wanted = target_species_for(organism)

        if gtdb_species not in counts:
            print(f"   GONE     {gtdb_species:34s} -> {organism}")
            print(f"            This GTDB species no longer exists in this release.")
            problems.append(f"{gtdb_species}: no longer in GTDB")
            continue

        ncbi_name, percent, total = majority_ncbi_species(counts[gtdb_species])
        epithet = gtdb_species.split(" ", 1)[1] if " " in gtdb_species else ""
        _base, epithet_is_suffixed = gao.strip_gtdb_suffix(epithet)

        if ncbi_name not in wanted:
            print(f"   CHANGED  {gtdb_species:34s} -> {organism}")
            print(f"            now resolves to '{ncbi_name or 'no majority'}' "
                  f"({percent:.1f}% of {total} genomes), not one of "
                  f"{sorted(wanted)}. REMOVE or re-justify this entry.")
            problems.append(f"{gtdb_species}: now resolves to {ncbi_name or 'no majority'}")
            continue

        # An entry whose epithet carries a suffix has to clear the higher bar,
        # because GTDB is saying the name's application is uncertain.
        if epithet_is_suffixed and (percent < MIN_AGREEMENT_FOR_SUFFIXED_EPITHET
                                    or total < MIN_GENOMES_FOR_SUFFIXED_EPITHET):
            print(f"   WEAK     {gtdb_species:34s} -> {organism}")
            print(f"            suffixed epithet, {percent:.1f}% of {total} genomes; "
                  f"the bar is >={MIN_AGREEMENT_FOR_SUFFIXED_EPITHET}% and "
                  f">={MIN_GENOMES_FOR_SUFFIXED_EPITHET} genomes. Consider removing.")
            problems.append(f"{gtdb_species}: evidence weakened ({percent:.1f}%, n={total})")
            continue

        print(f"   OK       {gtdb_species:34s} -> {organism:22s} "
              f"{ncbi_name} {percent:.1f}% of {total}")


def check_genus_rules(counts, genus_of_species, problems):
    """Question 2: is matching on the unsuffixed genus still safe?"""
    print("\n2. Genus-level rules")
    print("   (matching on the GTDB genus alone)\n")

    species_by_genus = collections.defaultdict(list)
    for gtdb_species, genus in genus_of_species.items():
        species_by_genus[genus].append(gtdb_species)

    for organism in sorted(gao.GENUS_ORGANISMS):
        wanted = target_species_for(organism)
        appropriate = inappropriate = 0
        offenders = []

        for gtdb_species in species_by_genus.get(organism, []):
            ncbi_name, _percent, total = majority_ncbi_species(counts[gtdb_species])
            if ncbi_name in wanted:
                appropriate += total
            else:
                inappropriate += total
                offenders.append((total, gtdb_species, ncbi_name))

        seen = appropriate + inappropriate
        share = (100.0 * appropriate / seen) if seen else 0.0
        verdict = "OK      " if share >= 95.0 else "UNSAFE  "
        print(f"   {verdict} g__{organism:20s} {share:5.1f}% appropriate "
              f"({appropriate} of {seen} genomes)")
        for total, gtdb_species, ncbi_name in sorted(offenders, reverse=True)[:5]:
            print(f"              - {gtdb_species:32s} n={total:<6d} "
                  f"is NCBI '{ncbi_name or 'unnamed'}'")
        if share < 95.0:
            problems.append(f"genus rule {organism}: only {share:.1f}% appropriate")


def check_for_missing_entries(counts, problems):
    """Question 3: which GTDB species SHOULD be in the table but are not?

    This is the half that catches the Campylobacter-style bug, where the curated
    species quietly moved into a suffixed genus and stopped matching. It scans
    every GTDB species, works out its NCBI name, and reports any that resolve to
    an organism we curate but that no rule in the module would currently match.
    """
    print("\n3. Candidate entries that are missing")
    print("   (GTDB species that resolve to a curated organism but do not match)\n")

    ncbi_to_organism = {}
    for organism in gao.ALL_ORGANISMS:
        for ncbi_species in target_species_for(organism):
            ncbi_to_organism[ncbi_species] = organism

    candidates = []
    for gtdb_species, ncbi_counter in counts.items():
        ncbi_name, percent, total = majority_ncbi_species(ncbi_counter)
        organism = ncbi_to_organism.get(ncbi_name)
        if not organism:
            continue
        # Would the module already emit something for this name? Build the same
        # lineage shape map_classification expects, so this asks the real code
        # rather than re-implementing its rules here.
        genus = gtdb_species.split(" ", 1)[0]
        lineage = f"d__Bacteria;g__{genus};s__{gtdb_species}"
        already, _reason = gao.map_classification(lineage)
        if already:
            continue
        candidates.append((total, percent, gtdb_species, ncbi_name, organism))

    if not candidates:
        print("   none - every GTDB species resolving to a curated organism is matched.")
        return

    for total, percent, gtdb_species, ncbi_name, organism in sorted(candidates, reverse=True):
        epithet = gtdb_species.split(" ", 1)[1] if " " in gtdb_species else ""
        epithet_base, _suffixed = gao.strip_gtdb_suffix(epithet)
        ncbi_epithet = ncbi_name.split(" ", 1)[1] if " " in ncbi_name else ""

        # The distinction that decides whether a candidate is safe.
        #
        # SAME epithet (Helicobacter pylori_C -> H. pylori): GTDB has split the
        # species or moved the genus, but it still calls the organism by the same
        # name. That is the Campylobacter_D jejuni situation - a naming artefact.
        #
        # DIFFERENT epithet (Serratia sarumanii -> S. marcescens): GTDB has
        # recognised a SEPARATE species that NCBI has not caught up with. The
        # genomes are labelled marcescens at NCBI because that is where they were
        # deposited, not because GTDB thinks they are marcescens. Mapping it would
        # assert a different taxon - exactly what this module exists to prevent -
        # so it is never proposed, whatever the percentage.
        names_agree = epithet_base and epithet_base == ncbi_epithet
        strong = (percent >= MIN_AGREEMENT_FOR_SUFFIXED_EPITHET
                  and total >= MIN_GENOMES_FOR_SUFFIXED_EPITHET)

        if not names_agree:
            verdict = "DIFFER  "   # GTDB says this is its own species; do not map
        elif not strong:
            verdict = "weak    "   # right name, evidence below the bar
        else:
            verdict = "ADD?    "
            problems.append(f"missing: {gtdb_species} -> {organism}")
        print(f"   {verdict} {gtdb_species:34s} -> {organism:22s} "
              f"{ncbi_name} {percent:.1f}% of {total}")

    print("\n   ADD?   = same species name, strong evidence -> a real candidate.")
    print("   weak   = same species name, but too few genomes or too little agreement.")
    print("   DIFFER = GTDB calls this a DIFFERENT species from the curated one;")
    print("            the NCBI labels are just where those genomes were deposited.")
    print("            Never map these - it would assert a different taxon.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Re-check the GTDB->AMRFinderPlus organism table against a "
                    "GTDB release. Reports; never edits.")
    parser.add_argument("--metadata", required=True,
                        help="bac120_metadata_r<N>.tsv.gz from the GTDB release "
                             "you are moving to.")
    args = parser.parse_args(argv)

    if not os.path.isfile(args.metadata):
        sys.exit(f"[check_gtdb_organism_table] no such file: {args.metadata}")

    print(f"Checking the organism table against {os.path.basename(args.metadata)}")
    counts, genus_of_species = read_gtdb_metadata(args.metadata)
    print(f"Read {len(counts):,} GTDB species clusters.")

    problems = []
    check_equivalence_entries(counts, problems)
    check_genus_rules(counts, genus_of_species, problems)
    check_for_missing_entries(counts, problems)

    print("\n" + "=" * 70)
    if problems:
        print(f"{len(problems)} thing(s) need a human decision:")
        for problem in problems:
            print(f"   - {problem}")
        print("\nEdit GTDB_SPECIES_EQUIVALENCES / GENUS_ORGANISMS in "
              "gtdb_amrfinder_organism.py, then re-run this check.")
        return 1

    print("The organism table is still correct for this GTDB release.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
