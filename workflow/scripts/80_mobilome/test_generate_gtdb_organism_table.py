"""Unit tests for generate_gtdb_organism_table.py.

The generator reads a GTDB release's own metadata — every genome's GTDB name and
its NCBI name, side by side — and writes the two small tables committed next to
it, gtdb_organism_equivalences.tsv and gtdb_organism_genus_rules.tsv. Those are
what gtdb_amrfinder_organism.py loads at runtime to decide whether an isolate
earns an AMRFinderPlus --organism flag, so a wrong row here ends up scoring a
genome against another organism's curated point-mutation list. That is the error
worth testing against, and it is why most of these tests are about REFUSING.

Everything runs on a tiny SYNTHETIC bac120_metadata file built in pytest's
tmp_path — no 225 MB GTDB download needed. Each test builds a handful of genomes
shaped like one real situation (a naming artefact to accept, a look-alike species
to reject, weak evidence to reject) and checks what the generator does with it.

Run: pytest workflow/scripts/80_mobilome/test_generate_gtdb_organism_table.py -q
"""

import csv
import os

import generate_gtdb_organism_table as gen
import gtdb_amrfinder_organism as gao


# ── Building a synthetic bac120_metadata.tsv ─────────────────────────────────

# read_gtdb_metadata only reads these two columns by NAME (gtdb_taxonomy,
# ncbi_taxonomy), so the fixture only needs to carry them — a real GTDB
# metadata file has ~110 other columns nothing here touches.
METADATA_HEADER = ["accession", "gtdb_taxonomy", "ncbi_taxonomy"]


def genome(accession, gtdb_species, ncbi_species):
    """One metadata row: an accession, its GTDB species, and its NCBI species.

    Both names are given WITHOUT the "s__" prefix and wrapped into a full
    lineage string here, matching how a real metadata file's columns look.
    """
    gtdb_lineage = ("d__Bacteria;p__Pseudomonadota;c__C;o__O;f__F;"
                    "g__%s;s__%s" % (gtdb_species.split(" ")[0], gtdb_species))
    ncbi_lineage = ("d__Bacteria;p__P;c__C;o__O;f__F;g__%s;s__%s"
                    % (ncbi_species.split(" ")[0] if ncbi_species else "",
                       ncbi_species))
    return [accession, gtdb_lineage, ncbi_lineage]


def cluster(gtdb_species, ncbi_species, n, agreeing_fraction=1.0):
    """n genomes in one GTDB species cluster, a given fraction of them
    carrying `ncbi_species` and the rest an unrelated filler NCBI name.

    This is what lets a test build, e.g., "70% of this cluster agrees" without
    writing 100 individual genome rows by hand.
    """
    agreeing = round(n * agreeing_fraction)
    rows = []
    for i in range(n):
        ncbi_name = ncbi_species if i < agreeing else "Filler unrelated"
        rows.append(genome("%s_%d" % (gtdb_species.replace(" ", "_"), i),
                           gtdb_species, ncbi_name))
    return rows


def write_metadata(tmp_path, rows):
    path = tmp_path / "bac120_metadata_test.tsv"
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(METADATA_HEADER)
        for row in rows:
            writer.writerow(row)
    return str(path)


def load(tmp_path, rows):
    """Write the fixture and read it back through the real loader."""
    path = write_metadata(tmp_path, rows)
    import check_gtdb_organism_table as gtdb_check
    return gtdb_check.read_gtdb_metadata(path)


def equivalence_for(rows, gtdb_species):
    """The generated row for one GTDB species, or None if it produced nothing."""
    for row in rows:
        if row["gtdb_species"] == gtdb_species:
            return row
    return None


# ── RULE A: unsuffixed epithet, genus reshuffled ─────────────────────────────
# GTDB broke up a genus and the species moved out with its epithet intact, so
# "Campylobacter_D jejuni" IS NCBI's C. jejuni. GTDB's own type-strain convention
# is what anchors that, so agreement alone carries the entry and no minimum
# genome count is demanded. Contrast RULE B below, where the suffix sits on the
# epithet and GTDB is saying the opposite.

def test_unsuffixed_epithet_naming_artefact_is_accepted(tmp_path):
    """The Campylobacter_D jejuni shape: genus suffixed, epithet plain, strong
    agreement — GTDB's own type-strain convention is doing the real work here,
    so no minimum genome count is demanded beyond a plain majority."""
    counts, genus_of_species = load(tmp_path, cluster(
        "Campylobacter_D jejuni", "Campylobacter jejuni", 50, agreeing_fraction=0.97))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    row = equivalence_for(rows, "Campylobacter_D jejuni")
    assert row is not None
    assert row["amrfinder_organism"] == "Campylobacter"
    assert row["rule"] == "unsuffixed_epithet"


def test_a_cluster_matching_no_curated_organism_produces_nothing(tmp_path):
    """Most GTDB clusters are environmental bacteria — the ordinary case."""
    counts, _genus = load(tmp_path, cluster(
        "Arthrobacter_D sp123456", "Arthrobacter globiformis", 10))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Arthrobacter_D sp123456") is None


# ── The two real bugs the first generation run against R226 found ────────────

def test_a_placeholder_epithet_is_never_accepted_however_it_votes(tmp_path):
    """The Arthrobacter_D sp009728235 → "Vibrio cholerae" case: ONE deposited
    genome, mislabelled at NCBI. A placeholder epithet ("sp<digits>") is not a
    real species name, so it must be refused regardless of what a single
    (mis-)labelled genome says at NCBI."""
    counts, _genus = load(tmp_path, cluster(
        "Arthrobacter_D sp009728235", "Vibrio cholerae", 1, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Arthrobacter_D sp009728235") is None


def test_a_look_alike_species_is_refused_even_at_majority_agreement(tmp_path):
    """The Enterococcus_B lactis → "Enterococcus faecium" case: a REAL, named
    GTDB species whose genomes are majority-labelled with a DIFFERENT curated
    species' name at NCBI (historical submissions predating the newer species
    being recognised). The GTDB epithet itself (lactis) disagrees with the
    matched NCBI epithet (faecium), so this must be refused even though the
    vote clears 50%.
    """
    counts, _genus = load(tmp_path, cluster(
        "Enterococcus_B lactis", "Enterococcus faecium", 100, agreeing_fraction=0.55))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Enterococcus_B lactis") is None


# ── RULE B: suffixed epithet, needs strong evidence ──────────────────────────
# A suffix on the EPITHET is GTDB saying it cannot safely attach that name here,
# usually because the type strain was never sequenced. Overriding that takes
# ≥99% agreement across ≥20 genomes (MIN_AGREEMENT_FOR_SUFFIXED_EPITHET and
# MIN_GENOMES_FOR_SUFFIXED_EPITHET in check_gtdb_organism_table.py). Both halves
# of the bar are tested below, because each was breached by a real Campylobacter
# cluster: jejuni_C clears the genome count but not the agreement, jejuni_A/_B/_D
# clear the agreement but sit on one or two genomes.

def test_suffixed_epithet_with_strong_evidence_is_accepted(tmp_path):
    """The Helicobacter pylori_C shape: genus plain, epithet itself suffixed.
    GTDB is saying the name's application here is uncertain, so only
    overwhelming genome-count evidence overrides that — here, 100% of 40."""
    counts, _genus = load(tmp_path, cluster(
        "Helicobacter pylori_C", "Helicobacter pylori", 40, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    row = equivalence_for(rows, "Helicobacter pylori_C")
    assert row is not None
    assert row["amrfinder_organism"] == "Helicobacter_pylori"
    assert row["rule"] == "suffixed_epithet_strong_evidence"


def test_suffixed_epithet_below_the_genome_count_bar_is_refused(tmp_path):
    """100% agreement, but on too few genomes for the number to mean anything -
    the jejuni_A/_B/_D situation (n=1-2 in the real data)."""
    counts, _genus = load(tmp_path, cluster(
        "Helicobacter pylori_Z", "Helicobacter pylori", 2, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Helicobacter pylori_Z") is None


def test_suffixed_epithet_below_the_agreement_bar_is_refused(tmp_path):
    """Plenty of genomes, but the agreement itself is too weak — the jejuni_C
    situation (55.6% real C. jejuni, majority actually C. lari)."""
    counts, _genus = load(tmp_path, cluster(
        "Helicobacter pylori_Y", "Helicobacter pylori", 50, agreeing_fraction=0.60))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Helicobacter pylori_Y") is None


# ── Skipping what the plain runtime rules already cover for free ─────────────
# Every row in the equivalences table is an override, and an override that merely
# restates what gtdb_amrfinder_organism.py would have matched anyway is noise a
# reader has to check. So the generator emits a row only where a plain exact-
# species or genus-level match would have failed.

def test_a_fully_unsuffixed_species_match_gets_no_override_row(tmp_path):
    """Enterococcus faecalis: neither genus nor epithet carries a suffix, and
    Enterococcus_faecalis is a species-level organism — the ordinary exact-
    species rule in gtdb_amrfinder_organism.py already matches this, so the
    table must not restate it."""
    counts, _genus = load(tmp_path, cluster(
        "Enterococcus faecalis", "Enterococcus faecalis", 30, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms=set())
    assert equivalence_for(rows, "Enterococcus faecalis") is None


def test_an_unsuffixed_genus_safe_organism_gets_no_override_row_regardless_of_epithet(tmp_path):
    """Escherichia coli_D: the GENUS is unsuffixed and Escherichia is genus-safe,
    so gtdb_amrfinder_organism.py's genus-level rule already matches this -
    REGARDLESS of the epithet's own suffix — so no override entry is needed."""
    counts, _genus = load(tmp_path, cluster(
        "Escherichia coli_D", "Escherichia coli", 30, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms={"Escherichia"})
    assert equivalence_for(rows, "Escherichia coli_D") is None


def test_a_suffixed_genus_still_needs_an_override_even_for_a_genus_safe_organism(tmp_path):
    """The counterpart: if the GENUS itself carries a suffix, the genus-safe
    shortcut does NOT apply (it only fires on the bare, unsuffixed genus), so a
    real override entry is still required even for Escherichia."""
    counts, _genus = load(tmp_path, cluster(
        "Escherichia_X coli", "Escherichia coli", 30, agreeing_fraction=1.0))
    rows = gen.build_equivalence_rows(counts, genus_safe_organisms={"Escherichia"})
    row = equivalence_for(rows, "Escherichia_X coli")
    assert row is not None
    assert row["amrfinder_organism"] == "Escherichia"


# ── Genus-safety verdicts ────────────────────────────────────────────────────
# Three of AMRFinderPlus's organisms are curated per GENUS rather than per
# species. Matching one of those from the bare GTDB genus is only safe when
# nearly everything GTDB files under that genus really is the curated organism,
# which is a question about genome counts, answered once per GTDB release rather
# than assumed. The bar is MIN_GENUS_SAFETY_PERCENT = 95.0.

def test_genus_safety_above_the_bar_is_judged_safe(tmp_path):
    # 95 S. enterica plus 5 S. bongori, both in the bare g__Salmonella and both
    # on the organism's curated species list (ORGANISM_TARGET_SPECIES in
    # check_gtdb_organism_table.py), so all 100 count as appropriate.
    rows = (cluster("Salmonella enterica", "Salmonella enterica", 95, agreeing_fraction=1.0)
            + cluster("Salmonella bongori", "Salmonella bongori", 5, agreeing_fraction=1.0))
    counts, genus_of_species = load(tmp_path, rows)
    safe, percent, appropriate, total = gen.decide_genus_safety(
        "Salmonella", counts, genus_of_species)
    assert safe is True
    assert appropriate == 100 and total == 100
    assert percent == 100.0


def test_genus_safety_below_the_bar_is_judged_unsafe(tmp_path):
    # Mimics the real Campylobacter situation: the unsuffixed genus is mostly
    # species AMRFinderPlus does not curate for this organism.
    rows = (cluster("Campylobacter fetus", "Campylobacter fetus", 80, agreeing_fraction=1.0)
            + cluster("Campylobacter hyointestinalis", "Campylobacter hyointestinalis",
                     20, agreeing_fraction=1.0))
    counts, genus_of_species = load(tmp_path, rows)
    safe, percent, appropriate, total = gen.decide_genus_safety(
        "Campylobacter", counts, genus_of_species)
    assert safe is False
    assert appropriate == 0
    assert percent == 0.0


def test_build_genus_rules_covers_exactly_the_three_candidates(tmp_path):
    # A verdict is written for all three genus-level organisms whatever the
    # metadata contains, including the two with no genomes in this fixture. A
    # silently missing row would read downstream as "genus matching is unsafe"
    # and quietly cost every E. coli isolate its --organism flag.
    counts, genus_of_species = load(tmp_path, cluster(
        "Escherichia coli", "Escherichia coli", 10, agreeing_fraction=1.0))
    rows = gen.build_genus_rules(counts, genus_of_species)
    organisms = {row["organism"] for row in rows}
    assert organisms == {"Campylobacter", "Escherichia", "Salmonella"}


# ── Reading and writing the generated files ──────────────────────────────────
# The generator writes the two TSVs and gtdb_amrfinder_organism.py reads them
# back, so the two halves have to agree on the file shape — in particular on the
# leading "#" provenance line that records which GTDB release and which date the
# table came from. A loader that choked on that line, or a writer that dropped
# it, would either break the run or leave a table nobody can date.

def test_write_tsv_then_load_species_equivalences_round_trips(tmp_path):
    path = str(tmp_path / "equivalences.tsv")
    rows = [{"gtdb_species": "Campylobacter_D jejuni", "amrfinder_organism": "Campylobacter",
            "ncbi_species": "Campylobacter jejuni", "percent_agreement": "97.3",
            "n_genomes": "3419", "rule": "unsuffixed_epithet"}]
    gen.write_tsv(path, gen.EQUIVALENCES_HEADER, rows, "test provenance line")

    with open(path) as handle:
        first_line = handle.readline()
    assert first_line.startswith("#")   # the loader must be able to skip this

    loaded = gao._load_species_equivalences(path, gao.ALL_ORGANISMS)
    assert loaded == {"Campylobacter_D jejuni": "Campylobacter"}


def test_write_tsv_then_load_genus_rules_round_trips(tmp_path):
    path = str(tmp_path / "genus_rules.tsv")
    rows = [
        {"organism": "Campylobacter", "genus_safe": "no", "percent_appropriate": "0.0",
         "n_appropriate": "0", "n_total": "403"},
        {"organism": "Escherichia", "genus_safe": "yes", "percent_appropriate": "98.3",
         "n_appropriate": "44776", "n_total": "45533"},
        {"organism": "Salmonella", "genus_safe": "yes", "percent_appropriate": "100.0",
         "n_appropriate": "17457", "n_total": "17457"},
    ]
    gen.write_tsv(path, gen.GENUS_RULES_HEADER, rows, "test provenance line")

    genus_organisms, name_only_organisms = gao._load_genus_rules(
        path, {"Campylobacter", "Escherichia", "Salmonella"})
    assert genus_organisms == {"Escherichia", "Salmonella"}
    assert name_only_organisms == {"Campylobacter"}


def test_guess_release_reads_the_r_number_out_of_the_filename():
    # The release number goes into the provenance line, and it is the only record
    # of which GTDB the shipped tables were built from. An unrecognisable
    # filename must say "unknown release" rather than guess a number.
    assert gen.guess_release("bac120_metadata_r226.tsv.gz") == "r226"
    assert gen.guess_release("/some/path/bac120_metadata_r232.tsv.gz") == "r232"
    assert gen.guess_release("no_release_number.tsv.gz") == "unknown release"
