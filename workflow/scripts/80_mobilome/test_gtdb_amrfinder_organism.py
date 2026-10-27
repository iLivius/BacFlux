"""Unit tests for gtdb_amrfinder_organism.py.

The mapper sits between GTDB-Tk and AMRFinderPlus. It reads the summary GTDB-Tk
wrote into 03.taxonomy/{sample} and emits one line naming an AMRFinderPlus
--organism, or a genuinely empty file. Rule amrfinderplus passes --organism (and
--mutation_all, which is meaningless without it) only when that file has
something in it, so what these tests really pin down is when the flag is given
and when it is withheld.

What matters biologically: a WRONG --organism makes AMRFinderPlus score the
genome against another organism's curated point-mutation list, which is worse
than no --organism at all. So the bulk of these tests are about refusing to
match, not about matching.

Everything is pure text handling, so no GTDB-Tk, no AMRFinderPlus and no database
are needed. The summary files are built on the fly in pytest's tmp_path with the
real GTDB-Tk 2.x column layout, so the tests also prove the column lookup works
when `classification` is not the first column. The two regression anchors at the
end of the file are the exception — they read a real GTDB-Tk summary and the real
GTDB taxonomy from disk, and skip themselves wherever those are absent.

Run: pytest workflow/scripts/80_mobilome/test_gtdb_amrfinder_organism.py -q
"""

import os

import gtdb_amrfinder_organism as gao


# ── Fixtures: build realistic GTDB-Tk summary files ──────────────────────────

# The real GTDB-Tk 2.x bac120 summary header, exactly as GTDB-Tk writes it into
# 03.taxonomy/{sample}/gtdbtk.bac120.summary.tsv.
GTDBTK_HEADER = [
    "user_genome", "classification", "closest_genome_reference",
    "closest_genome_reference_radius", "closest_genome_taxonomy",
    "closest_genome_ani", "closest_genome_af", "closest_placement_reference",
    "closest_placement_radius", "closest_placement_taxonomy",
    "closest_placement_ani", "closest_placement_af", "pplacer_taxonomy",
    "classification_method", "note",
    "other_related_references(genome_id,species_name,radius,ANI,AF)",
    "msa_percent", "translation_table", "red_value", "warnings",
]


def lineage(genus, species=""):
    """Build a GTDB classification string with the ranks the mapper reads.

    `genus` is the g__ value ("Pseudomonas_E"), `species` the full s__ binomial
    ("Pseudomonas_E sp010095445"); pass species="" for a genome GTDB placed only
    to genus, which really does come out of GTDB-Tk ending in a bare ";s__".
    """
    return ("d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Enterobacterales;f__Enterobacteriaceae;g__%s;s__%s"
            % (genus, species))


def write_summary(tmp_path, rows, name="gtdbtk.bac120.summary.tsv"):
    """Write a GTDB-Tk-shaped summary TSV containing the given rows.

    rows = [(user_genome, classification), ...]. Every other column is filled
    with plausible filler so the file has the real width.
    """
    lines = ["\t".join(GTDBTK_HEADER)]
    for genome, classification in rows:
        fields = [""] * len(GTDBTK_HEADER)
        fields[0] = genome
        fields[1] = classification
        fields[13] = "ani_screen"
        lines.append("\t".join(fields))
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def read_audit(path):
    """Return the audit TSV as (header list, single data row as a dict)."""
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header = lines[0].split("\t")
    assert len(lines) == 2, "the audit must hold exactly one decision row"
    values = lines[1].split("\t")
    return header, dict(zip(header, values))


# ── The curated list itself ──────────────────────────────────────────────────

def test_curated_list_has_31_organisms():
    # AMRFinderPlus 4.2.7 curates exactly 31 organisms (verified with
    # `amrfinder --list_organisms`).
    assert len(gao.ALL_ORGANISMS) == 31
    assert len(gao.SPECIES_ORGANISMS) == 28
    # Two of the three genus-level organisms may be picked FROM THE GENUS alone.
    # Campylobacter is the third and is deliberately not among them: GTDB's
    # unsuffixed g__Campylobacter holds none of the species AMRFinderPlus curates
    # for it, so it is reached only by name via GTDB_SPECIES_EQUIVALENCES. It is
    # still a legal --organism value, hence NAME_ONLY_ORGANISMS.
    assert len(gao.GENUS_ORGANISMS) == 2
    assert gao.NAME_ONLY_ORGANISMS == {"Campylobacter"}
    assert "Campylobacter" in gao.ALL_ORGANISMS


def test_curated_genera_are_derived_from_the_species_list():
    # The helper set used for the "why not" reasons must contain the genus of
    # every curated species, plus the three genus-level organisms.
    assert "Pseudomonas" in gao.CURATED_GENERA      # from Pseudomonas_aeruginosa
    assert "Staphylococcus" in gao.CURATED_GENERA   # two curated species
    assert "Escherichia" in gao.CURATED_GENERA      # genus-level organism
    assert "Arthrobacter" not in gao.CURATED_GENERA


# ── Matching: the cases where an --organism IS justified ─────────────────────

def test_exact_species_match():
    organism, reason = gao.map_classification(
        lineage("Klebsiella", "Klebsiella pneumoniae"))
    assert organism == "Klebsiella_pneumoniae"
    assert reason.startswith("exact species match")


def test_exact_species_match_for_every_genus_with_two_curated_species():
    # Genera where AMRFinderPlus curates more than one species must resolve to
    # the right one, not just to "something in this genus".
    for species, expected in [
        ("Klebsiella oxytoca", "Klebsiella_oxytoca"),
        ("Klebsiella pneumoniae", "Klebsiella_pneumoniae"),
        ("Enterococcus faecalis", "Enterococcus_faecalis"),
        ("Enterococcus faecium", "Enterococcus_faecium"),
        ("Staphylococcus aureus", "Staphylococcus_aureus"),
        ("Staphylococcus pseudintermedius", "Staphylococcus_pseudintermedius"),
        ("Burkholderia cepacia", "Burkholderia_cepacia"),
        ("Burkholderia pseudomallei", "Burkholderia_pseudomallei"),
        ("Pseudomonas aeruginosa", "Pseudomonas_aeruginosa"),
    ]:
        genus = species.split()[0]
        organism, _ = gao.map_classification(lineage(genus, species))
        assert organism == expected, species


def test_genus_level_match_for_the_genus_organisms():
    # AMRFinderPlus curates these per genus AND nearly everything GTDB files
    # under the unsuffixed genus really is that organism — 97.8% of g__Escherichia
    # and 100% of g__Salmonella in the shipped R232 rules (98.3% / 100% in R226),
    # against a 95% bar. Which genera clear it is not assumed here: the verdicts
    # are generated per GTDB release into gtdb_organism_genus_rules.tsv.
    for genus, species in [
        ("Escherichia", "Escherichia coli"),
        ("Salmonella", "Salmonella enterica"),
    ]:
        organism, reason = gao.map_classification(lineage(genus, species))
        assert organism == genus
        assert reason.startswith("genus-level match")


def test_genus_level_match_survives_a_placeholder_species():
    # An unnamed Salmonella is still a Salmonella, and the curated mutation list
    # is genus-wide, so the flag is still correct.
    organism, _ = gao.map_classification(
        lineage("Salmonella", "Salmonella sp002504165"))
    assert organism == "Salmonella"


def test_genus_level_match_survives_a_suffixed_species_epithet():
    # GTDB splits E. coli into coli, coli_D, coli_E ... and folds Shigella into
    # g__Escherichia. All of them are covered by AMRFinderPlus's genus-level
    # "Escherichia" organism, so the split epithet must not block the match.
    organism, _ = gao.map_classification(
        lineage("Escherichia", "Escherichia coli_D"))
    assert organism == "Escherichia"

    organism, _ = gao.map_classification(
        lineage("Escherichia", "Escherichia flexneri"))   # Shigella in GTDB
    assert organism == "Escherichia"


def test_every_returned_organism_is_a_real_amrfinder_value():
    # Guard against a typo in the table producing a flag AMRFinderPlus rejects.
    for genus, species in [
        ("Klebsiella", "Klebsiella pneumoniae"),
        ("Escherichia", "Escherichia coli"),
        ("Vibrio", "Vibrio cholerae"),
        ("Clostridioides", "Clostridioides difficile"),
    ]:
        organism, _ = gao.map_classification(lineage(genus, species))
        assert organism in gao.ALL_ORGANISMS


# ── Refusing to match: the GTDB-vs-NCBI mismatches ───────────────────────────

def test_gtdb_suffixed_genus_does_not_match_the_unsuffixed_organism():
    # The headline case: Pseudomonas_E is a lineage GTDB split off from
    # Pseudomonas. It is NOT P. aeruginosa, so no --organism.
    organism, reason = gao.map_classification(
        lineage("Pseudomonas_E", "Pseudomonas_E sp010095445"))
    assert organism == ""
    assert "Pseudomonas_E" in reason


def test_gtdb_suffixed_genus_blocks_even_a_named_curated_species():
    # Even spelled with the curated epithet, a suffixed genus is a different
    # taxon and must not collapse onto the NCBI species.
    organism, reason = gao.map_classification(
        lineage("Klebsiella_A", "Klebsiella_A pneumoniae"))
    assert organism == ""
    assert "suffixed genus" in reason


def test_gtdb_suffixed_genus_blocks_a_genus_level_organism_too():
    # Campylobacter_A / _B / _D are separate GTDB genera; only the unsuffixed
    # Campylobacter is the curated one.
    organism, reason = gao.map_classification(
        lineage("Campylobacter_A", "Campylobacter_A jejuni"))
    assert organism == ""
    assert "suffixed genus" in reason


def test_gtdb_suffixed_species_epithet_blocks_a_species_level_organism():
    # Enterobacter cloacae_A is a GTDB split of the NCBI species, so it is not
    # the curated Enterobacter_cloacae.
    organism, reason = gao.map_classification(
        lineage("Enterobacter", "Enterobacter cloacae_A"))
    assert organism == ""
    assert "suffixed species" in reason


# ── GTDB moved the genus but kept the species: the equivalence table ─────────
#
# Why this section exists. GTDB uses the same "_A/_B/_E" suffix for two very
# different things, and the difference decides whether AMRFinderPlus may be told
# the organism:
#   * suffix on the GENUS  — GTDB broke up a genus NCBI keeps whole. The species
#     moves to a new genus name but keeps its epithet, so
#     "s__Enterococcus_B faecium" IS NCBI's Enterococcus faecium.
#   * suffix on the EPITHET — GTDB is saying it cannot safely attach that name
#     here, so "faecalis_A" is presumed a different taxon from NCBI's
#     E. faecalis.
# The first kind is accepted on agreement alone. The second is accepted only on
# overwhelming evidence (≥99% of the cluster agreeing, across ≥20 genomes), which
# most of the shipped table's rows in fact are — every Helicobacter pylori_* row
# is an epithet suffix overridden that way.
#
# Nothing on that list is hand-typed. generate_gtdb_organism_table.py applies
# both rules to every GTDB species cluster in a release's own metadata and writes
# gtdb_organism_equivalences.tsv, which the script loads into
# GTDB_SPECIES_EQUIVALENCES at import.

def test_gtdb_split_genus_still_matches_enterococcus_faecium():
    # The case that motivated this exception: GTDB R226 has no unsuffixed
    # "Enterococcus faecium" at all (nor does R232 — the regression anchor at the
    # end of this file re-checks that against the installed release), so before
    # the exception table every E. faecium isolate silently ran AMRFinderPlus
    # without --organism and none of its curated point mutations (daptomycin
    # liaR/liaS/cls, linezolid 23S) could ever be reported.
    organism, reason = gao.map_classification(
        lineage("Enterococcus_B", "Enterococcus_B faecium"))
    assert organism == "Enterococcus_faecium"
    assert reason.startswith("known GTDB/NCBI naming artefact")
    # The audit must quote the GTDB name as GTDB writes it, so the decision can
    # be re-checked against the taxonomy file later.
    assert "Enterococcus_B faecium" in reason


def test_enterococcus_faecalis_still_matches_by_the_ordinary_species_rule():
    # The sister case, and the reason faecalis is NOT in the exception table:
    # GTDB keeps it in the unsuffixed g__Enterococcus, so nothing special is
    # needed and nothing special must be added.
    organism, reason = gao.map_classification(
        lineage("Enterococcus", "Enterococcus faecalis"))
    assert organism == "Enterococcus_faecalis"
    assert reason.startswith("exact species match")


def test_split_species_inside_the_split_genus_is_still_blocked():
    # "Enterococcus_B faecalis_A" is a real GTDB R226 name and is NOT NCBI's
    # E. faecalis. The exception must not open a door for it.
    organism, reason = gao.map_classification(
        lineage("Enterococcus_B", "Enterococcus_B faecalis_A"))
    assert organism == ""
    assert reason.startswith("no curated organism for this taxon")


def test_placeholder_species_inside_the_split_genus_is_still_blocked():
    # An unnamed Enterococcus_B is not E. faecium, whatever the genus says.
    organism, _ = gao.map_classification(
        lineage("Enterococcus_B", "Enterococcus_B sp900766935"))
    assert organism == ""


def test_suffixed_genus_not_in_the_exception_table_is_still_blocked():
    # The general rule is unchanged: only names the generated table earned a row
    # for get through, and the audit names that file rather than leaving the
    # reader guessing why one Enterococcus_B was accepted and this was not.
    organism, reason = gao.map_classification(
        lineage("Klebsiella_A", "Klebsiella_A pneumoniae"))
    assert organism == ""
    assert "gtdb_organism_equivalences.tsv" in reason


def test_exception_table_maps_only_onto_real_amrfinder_organisms():
    # A typo in the table would produce an --organism value AMRFinderPlus
    # rejects, which fails the run instead of degrading quietly.
    for gtdb_name, organism in gao.GTDB_SPECIES_EQUIVALENCES.items():
        assert organism in gao.ALL_ORGANISMS, gtdb_name


def test_exception_table_holds_only_same_species_naming_artefacts():
    """Structural guard on the GENERATED table: every entry must be a NAMING
    ARTEFACT for the SAME species AMRFinderPlus curates — GTDB moved a suffix
    onto the genus, the epithet, or both — never a genuinely different species
    that merely happens to share an NCBI label.

    The table is built by generate_gtdb_organism_table.py from a GTDB release's
    own metadata, and covers both patterns it finds real examples of: genus
    suffixed / epithet plain (Campylobacter_D jejuni is C. jejuni) AND genus plain
    / epithet suffixed (Helicobacter pylori_C is still H. pylori). This test does
    not re-derive the evidence behind either pattern — that needs the
    multi-hundred-MB metadata file, which a fast unit test must not require — it
    only checks the SHAPE every entry must have: once every GTDB suffix is
    stripped from both tokens, the result must name the exact organism the entry
    claims, not a look-alike.
    """
    for gtdb_name, organism in gao.GTDB_SPECIES_EQUIVALENCES.items():
        gtdb_genus, gtdb_epithet = gtdb_name.split(" ")
        genus_base, genus_is_suffixed = gao.strip_gtdb_suffix(gtdb_genus)
        epithet_base, epithet_is_suffixed = gao.strip_gtdb_suffix(gtdb_epithet)

        # The whole reason an entry is IN this table is that the plain rules
        # (exact species match, or a safe genus rule) would not already have
        # matched it without help — so it must carry a suffix somewhere.
        assert genus_is_suffixed or epithet_is_suffixed, gtdb_name
        assert not gao.is_placeholder_species(gtdb_epithet), gtdb_name

        # The organism may be genus-level ("Campylobacter") or species-level
        # ("Enterococcus_faecium"); only the latter carries an epithet to check.
        if "_" in organism:
            ncbi_genus, ncbi_epithet = organism.split("_", 1)
            assert epithet_base == ncbi_epithet, gtdb_name
        else:
            ncbi_genus = organism
        # Once stripped, the GTDB genus must equal the organism's own genus —
        # not merely start with it, which would also accept a wrong relative.
        assert genus_base == ncbi_genus, gtdb_name


def test_placeholder_species_in_a_curated_genus():
    # Klebsiella has curated species but no genus-level entry, so an unnamed
    # Klebsiella gets nothing — and the audit says exactly why.
    organism, reason = gao.map_classification(
        lineage("Klebsiella", "Klebsiella sp900493575"))
    assert organism == ""
    assert reason.startswith("GTDB placeholder species")


def test_placeholder_species_in_an_uncurated_genus_gets_the_plain_reason():
    # Arthrobacter, a common environmental genus, has nothing curated at all;
    # blaming the placeholder epithet would be a red herring.
    organism, reason = gao.map_classification(
        lineage("Arthrobacter", "Arthrobacter sp026393275"))
    assert organism == ""
    assert reason.startswith("no curated organism for this taxon")
    assert "Arthrobacter" in reason


def test_named_species_in_a_curated_genus_but_not_a_curated_species():
    organism, reason = gao.map_classification(
        lineage("Klebsiella", "Klebsiella variicola"))
    assert organism == ""
    assert "not one of the curated species" in reason


def test_classified_only_to_genus():
    # GTDB-Tk writes a bare ";s__" when it cannot place a genome to species.
    organism, reason = gao.map_classification(lineage("Klebsiella", ""))
    assert organism == ""
    assert "only to genus" in reason


def test_classified_only_to_genus_still_matches_a_genus_level_organism():
    organism, _ = gao.map_classification(lineage("Salmonella", ""))
    assert organism == "Salmonella"


def test_species_field_alone_is_enough():
    # Defensive: if the g__ rank were missing, the genus is taken from the
    # species binomial rather than losing the match.
    organism, _ = gao.map_classification(
        "d__Bacteria;p__Pseudomonadota;s__Pseudomonas aeruginosa")
    assert organism == "Pseudomonas_aeruginosa"


def test_unclassified_and_empty_classifications():
    # GTDB-Tk writes "Unclassified Bacteria" or "N/A" when placement fails.
    for text in ["Unclassified Bacteria", "N/A", "", "d__Bacteria;p__;c__;o__;f__;g__;s__"]:
        organism, reason = gao.map_classification(text)
        assert organism == ""
        assert reason  # never silent about why


def test_uncurated_environmental_genus():
    # Paenibacillus: another environmental genus with nothing curated.
    organism, reason = gao.map_classification(
        lineage("Paenibacillus", "Paenibacillus sp012647845"))
    assert organism == ""
    assert "Paenibacillus" in reason


# ── Reading the summary file ─────────────────────────────────────────────────

def test_single_row_summary(tmp_path):
    summary = write_summary(tmp_path, [
        ("SAMPLE1", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == "Pseudomonas_aeruginosa"
    assert "s__Pseudomonas aeruginosa" in classification
    assert reason.startswith("exact species match")


def test_multi_row_summary_from_hybrid_mode(tmp_path):
    # Hybrid mode classifies both assemblies of one isolate; they agree, so the
    # single shared answer is used and the audit records that they agreed.
    summary = write_summary(tmp_path, [
        ("SAMPLE1_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("SAMPLE1_ont", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == "Klebsiella_pneumoniae"
    assert "2 assemblies of this sample agree" in reason
    assert "SAMPLE1_illumina" in reason and "SAMPLE1_ont" in reason
    # Identical classifications are reported once, not twice.
    assert classification.count("d__Bacteria") == 1


def test_multi_row_summary_agreeing_on_no_organism(tmp_path):
    # The usual outcome for an environmental isolate: both assemblies agree on a
    # lineage with no curated organism, here Pseudomonas_E.
    summary = write_summary(tmp_path, [
        ("SAMPLE2_illumina", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
        ("SAMPLE2_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
    ])
    organism, _, reason = gao.decide_organism(summary, "SAMPLE2")
    assert organism == ""
    assert "Pseudomonas_E" in reason


def test_multi_row_summary_that_disagrees_emits_no_organism(tmp_path):
    # Two assemblies of one isolate landing on different organisms means the
    # taxonomy is not solid enough to pick a point-mutation list.
    summary = write_summary(tmp_path, [
        ("SAMPLE1_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("SAMPLE1_ont", lineage("Klebsiella", "Klebsiella oxytoca")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == ""
    assert "classified differently" in reason
    assert "SAMPLE1_illumina=Klebsiella_pneumoniae" in reason
    assert "SAMPLE1_ont=Klebsiella_oxytoca" in reason
    # Both lineages are kept in the audit so the disagreement is visible.
    assert classification.count("d__Bacteria") == 2


def test_disagreement_between_a_match_and_no_match(tmp_path):
    summary = write_summary(tmp_path, [
        ("SAMPLE3_illumina", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
        ("SAMPLE3_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp010095445")),
    ])
    organism, _, reason = gao.decide_organism(summary, "SAMPLE3")
    assert organism == ""
    assert "SAMPLE3_ont=none" in reason


def test_rows_of_other_genomes_are_ignored(tmp_path):
    # Belt and braces: if a summary ever held more than one isolate, only the
    # requested sample's rows may drive the decision.
    summary = write_summary(tmp_path, [
        ("SAMPLE2_illumina", lineage("Klebsiella", "Klebsiella oxytoca")),
        ("SAMPLE1_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("SAMPLE1_ont", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, _, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == "Klebsiella_pneumoniae"
    assert "SAMPLE2_illumina" not in reason


def test_sample_name_prefix_does_not_leak_into_another_sample(tmp_path):
    # "SAMPLE1" must not swallow a row named "SAMPLE10": matching is on the
    # exact name or the name plus an underscore.
    summary = write_summary(tmp_path, [
        ("SAMPLE10", lineage("Klebsiella", "Klebsiella oxytoca")),
        ("SAMPLE1", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, _, _ = gao.decide_organism(summary, "SAMPLE1")
    assert organism == "Klebsiella_pneumoniae"


def test_row_not_named_after_the_sample_is_still_used(tmp_path):
    # GTDB-Tk is run per sample into its own directory, so a row named after the
    # contigs file (not the sample) still belongs to this isolate — use it, and
    # note that we did.
    summary = write_summary(tmp_path, [
        ("contigs_final", lineage("Vibrio", "Vibrio cholerae")),
    ])
    organism, _, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == "Vibrio_cholerae"
    assert "used every row" in reason


def test_truncated_line_is_skipped_not_fatal(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    good = [""] * len(GTDBTK_HEADER)
    good[0] = "SAMPLE1"
    good[1] = lineage("Vibrio", "Vibrio vulnificus")
    summary.write_text("\n".join([
        "\t".join(GTDBTK_HEADER),
        "SAMPLE1_broken",              # a single field: an interrupted GTDB-Tk write
        "\t".join(good),
    ]) + "\n")
    organism, _, _ = gao.decide_organism(str(summary), "SAMPLE1")
    assert organism == "Vibrio_vulnificus"


def test_windows_line_endings_are_handled(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    row = [""] * len(GTDBTK_HEADER)
    row[0] = "SAMPLE1"
    row[1] = lineage("Helicobacter", "Helicobacter pylori")
    summary.write_text("\r\n".join([
        "\t".join(GTDBTK_HEADER), "\t".join(row)]) + "\r\n")
    organism, _, _ = gao.decide_organism(str(summary), "SAMPLE1")
    assert organism == "Helicobacter_pylori"


# ── Graceful degradation: nothing here may raise or exit non-zero ────────────

def test_missing_summary_file(tmp_path):
    organism, classification, reason = gao.decide_organism(
        str(tmp_path / "does_not_exist.tsv"), "SAMPLE1")
    assert organism == ""
    assert classification == ""
    assert "not found" in reason


def test_empty_summary_file(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    summary.write_text("")
    organism, _, reason = gao.decide_organism(str(summary), "SAMPLE1")
    assert organism == ""
    assert "empty" in reason


def test_header_only_summary_file(tmp_path):
    summary = write_summary(tmp_path, [])
    organism, _, reason = gao.decide_organism(summary, "SAMPLE1")
    assert organism == ""
    assert "no usable rows" in reason


def test_summary_without_the_expected_columns(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    summary.write_text("genome\tlineage\nfoo\tbar\n")
    organism, _, reason = gao.decide_organism(str(summary), "SAMPLE1")
    assert organism == ""
    assert "no user_genome/classification columns" in reason


# ── Finding the summary inside a GTDB-Tk output directory ────────────────────

def test_find_summary_at_the_top_level(tmp_path):
    # classify_wf leaves a copy/symlink of the summary at the top level.
    write_summary(tmp_path, [("SAMPLE1", lineage("Vibrio", "Vibrio cholerae"))])
    found, problem = gao.find_summary_in_dir(str(tmp_path))
    assert problem == ""
    assert found.endswith("gtdbtk.bac120.summary.tsv")


def test_find_summary_in_the_classify_subdirectory(tmp_path):
    # The real file lives in classify/; only the symlink is at the top level,
    # and a rule may be pointed at a directory where it was never made.
    classify = tmp_path / "classify"
    classify.mkdir()
    write_summary(classify, [("SAMPLE1", lineage("Vibrio", "Vibrio cholerae"))])
    found, problem = gao.find_summary_in_dir(str(tmp_path))
    assert problem == ""
    assert found.endswith(os.path.join("classify", "gtdbtk.bac120.summary.tsv"))


def test_bacterial_summary_wins_over_the_archaeal_one(tmp_path):
    # BacFlux is a bacterial workflow; when GTDB-Tk wrote both, read bac120.
    write_summary(tmp_path, [("SAMPLE1", lineage("Vibrio", "Vibrio cholerae"))],
                  name="gtdbtk.ar53.summary.tsv")
    write_summary(tmp_path, [("SAMPLE1", lineage("Klebsiella", "Klebsiella oxytoca"))])
    found, _ = gao.find_summary_in_dir(str(tmp_path))
    assert found.endswith("gtdbtk.bac120.summary.tsv")


def test_directory_without_a_summary_is_not_fatal(tmp_path):
    empty_dir = tmp_path / "03.taxonomy" / "SAMPLE1"
    empty_dir.mkdir(parents=True)
    found, problem = gao.find_summary_in_dir(str(empty_dir))
    assert found == ""
    assert "no gtdbtk.*.summary.tsv" in problem


def test_missing_directory_is_not_fatal(tmp_path):
    found, problem = gao.find_summary_in_dir(str(tmp_path / "nope"))
    assert found == ""
    assert "not found" in problem


def test_cli_accepts_the_directory_form(tmp_path):
    # The Snakemake rule declares the whole 03.taxonomy/{sample} directory as
    # its input, so it passes --gtdbtk-dir instead of --gtdbtk-summary.
    taxonomy_dir = tmp_path / "03.taxonomy" / "SAMPLE1"
    classify = taxonomy_dir / "classify"
    classify.mkdir(parents=True)
    write_summary(classify, [
        ("SAMPLE1_illumina", lineage("Escherichia", "Escherichia coli")),
        ("SAMPLE1_ont", lineage("Escherichia", "Escherichia coli")),
    ])
    organism_file = tmp_path / "org.txt"
    audit_file = tmp_path / "audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-dir", str(taxonomy_dir),
        "--sample", "SAMPLE1",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.read_text().strip() == "Escherichia"


def test_cli_requires_exactly_one_input_form(tmp_path):
    # Bad wiring of the rule (both forms, or neither) should be a loud argparse
    # error, not a silent "no organism" that looks like a real result.
    import pytest

    base = ["--sample", "SAMPLE1",
            "--out-organism", str(tmp_path / "org.txt"),
            "--out-audit", str(tmp_path / "audit.tsv")]
    with pytest.raises(SystemExit):
        gao.main(base)
    with pytest.raises(SystemExit):
        gao.main(["--gtdbtk-summary", "a.tsv", "--gtdbtk-dir", "b"] + base)


def test_cli_directory_form_degrades_gracefully(tmp_path):
    taxonomy_dir = tmp_path / "03.taxonomy" / "SAMPLE1"
    taxonomy_dir.mkdir(parents=True)
    organism_file = tmp_path / "org.txt"
    audit_file = tmp_path / "audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-dir", str(taxonomy_dir),
        "--sample", "SAMPLE1",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.read_text() == ""
    _, row = read_audit(str(audit_file))
    assert "no gtdbtk.*.summary.tsv" in row["reason"]


# ── The two output files (what Snakemake actually consumes) ──────────────────

def test_cli_writes_the_organism_and_the_audit_on_a_match(tmp_path):
    summary = write_summary(tmp_path, [
        ("SAMPLE1", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
    ])
    organism_file = tmp_path / "out" / "SAMPLE1_amrfinder_organism.txt"
    audit_file = tmp_path / "out" / "SAMPLE1_amrfinder_organism_decision.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "SAMPLE1",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    # The rule reads this with $(cat ...), so the trailing newline is harmless.
    assert organism_file.read_text().strip() == "Pseudomonas_aeruginosa"

    header, row = read_audit(str(audit_file))
    assert header == ["sample", "gtdb_classification", "matched_organism", "reason"]
    assert row["sample"] == "SAMPLE1"
    assert row["matched_organism"] == "Pseudomonas_aeruginosa"
    assert "s__Pseudomonas aeruginosa" in row["gtdb_classification"]
    assert row["reason"].startswith("exact species match")


def test_cli_writes_the_organism_for_a_gtdb_split_genus(tmp_path):
    # End to end for the exception: a hybrid E. faecium sample (two assemblies,
    # both classified as Enterococcus_B) must end up with a NON-empty organism
    # file, because that file is the only thing the Snakemake rule looks at when
    # deciding whether to pass --organism to AMRFinderPlus.
    summary = write_summary(tmp_path, [
        ("EF01_illumina", lineage("Enterococcus_B", "Enterococcus_B faecium")),
        ("EF01_ont", lineage("Enterococcus_B", "Enterococcus_B faecium")),
    ])
    organism_file = tmp_path / "EF01_amrfinder_organism.txt"
    audit_file = tmp_path / "EF01_amrfinder_organism_decision.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "EF01",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.read_text().strip() == "Enterococcus_faecium"

    _, row = read_audit(str(audit_file))
    assert row["matched_organism"] == "Enterococcus_faecium"
    assert row["reason"].startswith("known GTDB/NCBI naming artefact")
    assert "2 assemblies of this sample agree" in row["reason"]


def test_cli_writes_an_empty_file_and_exits_zero_on_no_match(tmp_path):
    # The common case for environmental isolates: the file must be genuinely
    # empty so `[ -n "$(cat ...)" ]` in the rule is false and AMRFinderPlus runs
    # without --organism.
    summary = write_summary(tmp_path, [
        ("SAMPLE2_illumina", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
        ("SAMPLE2_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
    ])
    organism_file = tmp_path / "SAMPLE2_amrfinder_organism.txt"
    audit_file = tmp_path / "SAMPLE2_amrfinder_organism_decision.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "SAMPLE2",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.exists()
    assert organism_file.read_text() == ""
    assert os.path.getsize(str(organism_file)) == 0

    _, row = read_audit(str(audit_file))
    assert row["matched_organism"] == "NA"
    assert "Pseudomonas_E" in row["reason"]


def test_cli_still_succeeds_when_the_summary_is_missing(tmp_path):
    organism_file = tmp_path / "org.txt"
    audit_file = tmp_path / "audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", str(tmp_path / "nope.tsv"),
        "--sample", "SAMPLE1",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.read_text() == ""
    _, row = read_audit(str(audit_file))
    assert row["gtdb_classification"] == "NA"
    assert row["matched_organism"] == "NA"
    assert "not found" in row["reason"]


def test_audit_stays_a_valid_four_column_tsv(tmp_path):
    # A reason string with a stray tab or newline would break every downstream
    # reader of the audit, so check the shape rather than the wording.
    summary = write_summary(tmp_path, [
        ("SAMPLE1_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("SAMPLE1_ont", lineage("Klebsiella", "Klebsiella oxytoca")),
    ])
    audit_file = tmp_path / "audit.tsv"
    gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "SAMPLE1",
        "--out-organism", str(tmp_path / "org.txt"),
        "--out-audit", str(audit_file),
    ])
    with open(str(audit_file)) as handle:
        lines = handle.read().splitlines()
    assert len(lines) == 2
    for line in lines:
        assert len(line.split("\t")) == 4


def test_output_directories_are_created(tmp_path):
    # Snakemake normally makes them, but the script must not fall over if a
    # nested output path is handed to it directly.
    summary = write_summary(tmp_path, [
        ("SAMPLE1", lineage("Salmonella", "Salmonella enterica")),
    ])
    organism_file = tmp_path / "08.mobilome" / "SAMPLE1" / "organism.txt"
    audit_file = tmp_path / "08.mobilome" / "SAMPLE1" / "audit.tsv"
    gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "SAMPLE1",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])
    assert organism_file.read_text().strip() == "Salmonella"


# ── Regression anchor against a real GTDB-Tk file ────────────────────────────

# Every test above builds its summary on the fly. This one reads a summary that
# GTDB-Tk really wrote, to prove the mapper copes with the genuine file. Give it
# a summary from a hybrid run (two rows, <sample>_illumina and <sample>_ont) of
# an isolate in a genus AMRFinderPlus does not curate, together with that
# sample's name and the genus the summary reports (the g__ value, without the
# "g__" prefix):
#
#     BACFLUX_TEST_GTDBTK_SUMMARY=/your/path/gtdbtk.bac120.summary.tsv \
#     BACFLUX_TEST_GTDBTK_SAMPLE=<sample> \
#     BACFLUX_TEST_GTDBTK_GENUS=<genus> pytest ...
#
# No such file ships with the repository, so the test skips unless all three
# variables are set.
REAL_SUMMARY = os.environ.get("BACFLUX_TEST_GTDBTK_SUMMARY", "")
REAL_SAMPLE = os.environ.get("BACFLUX_TEST_GTDBTK_SAMPLE", "")
REAL_GENUS = os.environ.get("BACFLUX_TEST_GTDBTK_GENUS", "")


def test_real_gtdbtk_summary_gives_no_organism(tmp_path):
    if not os.path.isfile(REAL_SUMMARY) or not REAL_SAMPLE or not REAL_GENUS:
        import pytest
        pytest.skip("set BACFLUX_TEST_GTDBTK_SUMMARY, BACFLUX_TEST_GTDBTK_SAMPLE "
                    "and BACFLUX_TEST_GTDBTK_GENUS to run this")

    organism, classification, reason = gao.decide_organism(REAL_SUMMARY, REAL_SAMPLE)
    assert organism == ""
    # The genus to expect is passed in rather than written here because GTDB
    # renames genera between releases: a split-off lineage with a placeholder
    # suffix (the _E of Pseudomonas_E) can receive a validly published genus
    # name in a later release, and GTDB's own FAQ warns that suffix letters are
    # not guaranteed to survive. So the right genus is whatever the release that
    # wrote this summary calls it. Nothing is curated for it, so the organism
    # stays empty, and the reason must name the genus that was refused.
    assert "g__" + REAL_GENUS in classification
    assert REAL_GENUS in reason
    # Both hybrid rows were seen and agreed.
    assert REAL_SAMPLE + "_illumina" in reason and REAL_SAMPLE + "_ont" in reason


# ── Regression anchor against the real GTDB taxonomy ─────────────────────────

# The exception table is only correct for as long as GTDB keeps these names, so
# this checks it against the taxonomy file that ships inside the GTDB release
# (the directory the config key directories.gtdbtk_db points at). If a future
# release renames things, this test fails and the table gets revisited instead
# of quietly mapping the wrong taxon. The file is far too large to vendor, so,
# like the summary above, it comes from the environment and the test skips
# without it:
#
#     BACFLUX_TEST_GTDB_TAXONOMY=$GTDBTK_DB/taxonomy/gtdb_taxonomy.tsv pytest ...
GTDB_TAXONOMY = os.environ.get("BACFLUX_TEST_GTDB_TAXONOMY", "")


def gtdb_species_names_matching(prefix):
    """Every distinct GTDB species name (the s__ value) starting with `prefix`.

    Reads the release's own taxonomy TSV — one genome per line, the lineage in
    the second column — and returns the species field only, e.g.
    {"Enterococcus faecalis", "Enterococcus_B faecium", ...}.
    """
    names = set()
    with open(GTDB_TAXONOMY) as handle:
        for line in handle:
            lineage_field = line.rstrip("\n").split("\t")[-1]
            species_field = lineage_field.split(";")[-1]
            if not species_field.startswith("s__"):
                continue
            species_name = species_field[3:]
            if species_name.startswith(prefix):
                names.add(species_name)
    return names


def test_gtdb_really_has_no_unsuffixed_enterococcus_faecium():
    if not os.path.isfile(GTDB_TAXONOMY):
        import pytest
        pytest.skip("set BACFLUX_TEST_GTDB_TAXONOMY to run this")

    enterococci = gtdb_species_names_matching("Enterococcus")
    faecium_names = {n for n in enterococci if n.split()[1] == "faecium"}
    # The whole reason for the exception: the plain NCBI name does not exist in
    # GTDB, only the split-genus form does.
    assert faecium_names == {"Enterococcus_B faecium"}
    # And the reason faecalis needs no exception: it is there unsuffixed.
    assert "Enterococcus faecalis" in enterococci
    # Every key of the exception table must still be a name GTDB actually uses.
    for gtdb_name in gao.GTDB_SPECIES_EQUIVALENCES:
        genus = gtdb_name.split()[0]
        assert gtdb_name in gtdb_species_names_matching(genus), gtdb_name


# ── Campylobacter: the curated species live in a SUFFIXED genus ──────────────
# The unsuffixed genus holds only species AMRFinderPlus does not curate.
#
# GTDB assigns the genus suffix by where the genus TYPE SPECIES landed. The type
# species of Campylobacter is C. fetus, which keeps g__Campylobacter, so jejuni
# and coli were pushed into g__Campylobacter_D. Matching on the unsuffixed genus
# therefore hits only genomes of the wrong species and none of the right ones:
# gtdb_organism_genus_rules.tsv scores it 0 appropriate out of 450 in R232 (403
# in R226). That is why Campylobacter is the one genus-level organism reachable
# by name alone, through the equivalence table, and never from its bare genus.


def test_the_curated_campylobacter_species_are_matched_despite_the_genus_suffix():
    """Verified against NCBI: both GTDB clusters are anchored on type material
    (C. jejuni NCTC 11351 / C. coli LMG 9860), so they ARE the curated species."""
    for species in ("Campylobacter_D jejuni", "Campylobacter_D coli"):
        organism, reason = gao.map_classification(
            lineage("Campylobacter_D", species))
        assert organism == "Campylobacter", species
        assert reason


def test_the_unsuffixed_campylobacter_genus_no_longer_matches():
    """C. fetus and friends must NOT get --organism Campylobacter.

    AMRFinderPlus documents that taxgroup as covering C. jejuni / C. coli only,
    and its mutation references are C. jejuni sequences. Screening C. fetus gyrA
    against them is exactly the wrong-organism call the module docstring warns is
    worse than no call at all.
    """
    for species in ("Campylobacter fetus", "Campylobacter hyointestinalis",
                    "Campylobacter testudinum", "Campylobacter lanienae"):
        organism, reason = gao.map_classification(lineage("Campylobacter", species))
        assert organism == "", species
        assert reason


def test_well_evidenced_suffixed_epithets_are_admitted():
    """coli_A and coli_B clear the RULE B bar: both 100% NCBI Campylobacter coli,
    on n=92 and n=72 genomes in the shipped R232 table (90 and 70 in R226)."""
    for species in ("Campylobacter_D coli_A", "Campylobacter_D coli_B"):
        organism, _reason = gao.map_classification(lineage("Campylobacter_D", species))
        assert organism == "Campylobacter", species


def test_thinly_evidenced_or_contradicted_suffixed_epithets_are_refused():
    """The cases the ≥99%-and-≥20-genomes bar exists to exclude, one for each
    half of it.

    jejuni_C is the important one: only 55.6% of its cluster is NCBI C. jejuni
    and the majority is C. LARI, so naive suffix-stripping would have called it
    C. jejuni. jejuni_A/_B/_D are 100% but on 1–2 genomes, where a single
    mislabelled genome would flip the vote. None of the four has a row in
    gtdb_organism_equivalences.tsv, which is what these assertions rest on.
    """
    for species in ("Campylobacter_D jejuni_A", "Campylobacter_D jejuni_B",
                    "Campylobacter_D jejuni_C", "Campylobacter_D jejuni_D"):
        organism, reason = gao.map_classification(lineage("Campylobacter_D", species))
        assert organism == "", species
        assert reason


def test_every_organism_the_equivalence_table_emits_is_a_legal_amrfinder_value():
    """AMRFinderPlus exits with an error on an unrecognised --organism, so an
    entry pointing at a non-existent name would kill the rule, not degrade."""
    for gtdb_name, organism in gao.GTDB_SPECIES_EQUIVALENCES.items():
        assert organism in gao.ALL_ORGANISMS, gtdb_name


# ── Burkholderia: a GTDB species wider than the curated organism ─────────────
# GTDB merges Burkholderia pseudomallei into its s__Burkholderia mallei (there is
# no s__Burkholderia pseudomallei in GTDB), while AMRFinderPlus curates the two as
# separate organisms with different point-mutation lists: none for
# Burkholderia_mallei, 15 for Burkholderia_pseudomallei in AMRFinderPlus database
# 2026-01-21.1. The plain exact-species rule would therefore give every
# B. pseudomallei isolate the B. mallei list, find no mutations, and write an
# audit row that reads like a confident match. So GTDB's B. mallei is on the
# hand-curated refusal list, gtdb_organism_refusals.tsv: no --organism at all,
# and an audit reason that says why.

def burkholderia_lineage(species):
    """A real-shaped GTDB lineage for a Burkholderia species (the shared
    lineage() helper fills in Enterobacterales ranks, which would be wrong here,
    although the mapper reads only g__ and s__)."""
    return ("d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
            "o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__%s"
            % species)


def test_gtdb_burkholderia_mallei_gets_no_organism():
    # The shipped refusal list must be present and loaded; without it this
    # lineage would fall through to the exact-species match again.
    assert "Burkholderia mallei" in gao.GTDB_SPECIES_REFUSALS, (
        "gtdb_organism_refusals.tsv is missing or no longer lists it")

    organism, reason = gao.map_classification(
        burkholderia_lineage("Burkholderia mallei"))
    assert organism == ""
    # The audit must say that this was a deliberate refusal, of which GTDB
    # name, where it is recorded, and the biology behind it, in plain words.
    assert reason.startswith("organism withheld on purpose")
    assert "s__Burkholderia mallei" in reason
    assert "gtdb_organism_refusals.tsv" in reason
    assert "pseudomallei" in reason
    # And it must not be mistaken for the ordinary "nothing curated" outcome.
    assert "exact species match" not in reason
    assert not reason.startswith("no curated organism")


def test_cli_writes_an_empty_organism_file_for_a_refused_species(tmp_path):
    # End to end, as the Snakemake rule sees it: a hybrid sample whose two
    # assemblies are both GTDB B. mallei. The organism file must be genuinely
    # empty (rule amrfinderplus then runs without --organism), and the audit
    # row must carry the refusal reason with matched_organism = NA.
    summary = write_summary(tmp_path, [
        ("BM01_illumina", burkholderia_lineage("Burkholderia mallei")),
        ("BM01_ont", burkholderia_lineage("Burkholderia mallei")),
    ])
    organism_file = tmp_path / "BM01_amrfinder_organism.txt"
    audit_file = tmp_path / "BM01_amrfinder_organism_audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "BM01",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert os.path.getsize(str(organism_file)) == 0

    header, row = read_audit(str(audit_file))
    assert header == ["sample", "gtdb_classification", "matched_organism", "reason"]
    assert row["matched_organism"] == "NA"
    assert "s__Burkholderia mallei" in row["gtdb_classification"]
    assert row["reason"].startswith("organism withheld on purpose")
    assert "2 assemblies of this sample agree" in row["reason"]
    # The reason comes from a hand-edited file, so check it did not break the
    # four-column shape of the audit.
    with open(str(audit_file)) as handle:
        for line in handle.read().splitlines():
            assert len(line.split("\t")) == 4


def test_neighbouring_burkholderia_species_are_unchanged():
    # The refusal is keyed on one exact GTDB name and must not spill over onto
    # the rest of the genus.
    organism, reason = gao.map_classification(
        burkholderia_lineage("Burkholderia cepacia"))
    assert organism == "Burkholderia_cepacia"
    assert reason.startswith("exact species match")

    # Should a GTDB release ever name B. pseudomallei separately again, that
    # name is unambiguous and still earns its own organism.
    organism, reason = gao.map_classification(
        burkholderia_lineage("Burkholderia pseudomallei"))
    assert organism == "Burkholderia_pseudomallei"
    assert reason.startswith("exact species match")

    # An uncurated Burkholderia keeps the ordinary "not curated" reason, not
    # the refusal one.
    organism, reason = gao.map_classification(
        burkholderia_lineage("Burkholderia thailandensis"))
    assert organism == ""
    assert "not one of the curated species" in reason


def test_refusal_wins_over_an_equivalence_row(monkeypatch):
    # The refusal is checked before the generated equivalences table, so a name
    # that somehow sits in both is still refused. Enterococcus_B faecium is used
    # only because it has a real equivalence row; it is not actually refused.
    assert "Enterococcus_B faecium" in gao.GTDB_SPECIES_EQUIVALENCES
    monkeypatch.setitem(gao.GTDB_SPECIES_REFUSALS, "Enterococcus_B faecium",
                        "test reason")
    organism, reason = gao.map_classification(
        lineage("Enterococcus_B", "Enterococcus_B faecium"))
    assert organism == ""
    assert reason.endswith("test reason")


def test_shipped_refusals_and_equivalences_do_not_overlap():
    # A name in both files would be a contradiction a reader has to untangle;
    # the generator is meant to make it impossible.
    overlap = set(gao.GTDB_SPECIES_REFUSALS) & set(gao.GTDB_SPECIES_EQUIVALENCES)
    assert overlap == set()


def test_refusal_loader_keeps_a_row_with_a_blank_reason(tmp_path):
    # On this list, dropping a row is the unsafe direction (the name would be
    # mapped again), so a row with an empty reason still refuses.
    path = tmp_path / "refusals.tsv"
    path.write_text("# provenance\n"
                    "gtdb_species\treason\n"
                    "Burkholderia mallei\t\n"
                    "\tno species on this row\n")
    refusals = gao._load_species_refusals(str(path))
    assert list(refusals) == ["Burkholderia mallei"]
    assert refusals["Burkholderia mallei"].startswith("no reason recorded")


def test_refusal_loader_keeps_a_row_with_no_tab_at_all(tmp_path):
    # A hand-typed row with the reason left out and no tab: still a refusal.
    path = tmp_path / "refusals.tsv"
    path.write_text("# provenance\n"
                    "gtdb_species\treason\n"
                    "Burkholderia mallei\n")
    refusals = gao._load_species_refusals(str(path))
    assert list(refusals) == ["Burkholderia mallei"]
    assert refusals["Burkholderia mallei"].startswith("no reason recorded")


def test_refusal_loader_keeps_a_row_with_an_extra_column(tmp_path):
    # One tab too many must not drop the row; the pieces become one reason, so
    # the audit TSV keeps its four columns.
    path = tmp_path / "refusals.tsv"
    path.write_text("# provenance\n"
                    "gtdb_species\treason\n"
                    "Burkholderia mallei\tfirst part\tsecond part\n")
    refusals = gao._load_species_refusals(str(path))
    assert refusals == {"Burkholderia mallei": "first part; second part"}
    assert "\t" not in refusals["Burkholderia mallei"]


def test_refusal_loader_accepts_a_pasted_gtdb_prefix(tmp_path):
    # A name copied straight from a GTDB-Tk lineage keeps its "s__"; it must
    # still key on the bare name that map_classification looks up.
    path = tmp_path / "refusals.tsv"
    path.write_text("gtdb_species\treason\n"
                    "s__Burkholderia  mallei \tsome reason\n")
    refusals = gao._load_species_refusals(str(path))
    assert refusals == {"Burkholderia mallei": "some reason"}


def test_refusal_loader_warns_about_a_row_it_cannot_key(tmp_path, capsys):
    # A row typed with spaces instead of a tab gives a "name" of many words,
    # which can never match a GTDB species: it is reported, not silently lost.
    path = tmp_path / "refusals.tsv"
    path.write_text("gtdb_species\treason\n"
                    "Burkholderia mallei  reason typed after spaces\n"
                    "\tno species on this row\n")
    refusals = gao._load_species_refusals(str(path))
    assert refusals == {}
    warnings = capsys.readouterr().err
    assert warnings.count("WARNING") == 2
    assert "reason typed after spaces" in warnings


def test_refusal_loader_warns_when_the_file_is_missing(tmp_path, capsys):
    # A missing list cannot fail towards caution, so it must at least be loud
    # (stderr ends up in the amrfinder_organism log).
    refusals = gao._load_species_refusals(str(tmp_path / "absent.tsv"))
    assert refusals == {}
    assert "WARNING" in capsys.readouterr().err


def test_gtdb_really_merges_pseudomallei_into_mallei():
    # Regression anchor for the refusal, against the installed GTDB taxonomy
    # (same BACFLUX_TEST_GTDB_TAXONOMY variable as the anchor above). If a
    # future release gives B. pseudomallei its own species again, this fails and
    # the refusal gets revisited instead of withholding an organism for nothing.
    if not os.path.isfile(GTDB_TAXONOMY):
        import pytest
        pytest.skip("set BACFLUX_TEST_GTDB_TAXONOMY to run this")

    burkholderia = gtdb_species_names_matching("Burkholderia ")
    assert "Burkholderia mallei" in burkholderia
    assert "Burkholderia pseudomallei" not in burkholderia
    # Every refused name must still be a name GTDB actually uses.
    for gtdb_name in gao.GTDB_SPECIES_REFUSALS:
        genus = gtdb_name.split()[0]
        assert gtdb_name in gtdb_species_names_matching(genus), gtdb_name
