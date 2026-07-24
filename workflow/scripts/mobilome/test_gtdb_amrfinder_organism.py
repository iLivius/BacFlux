"""Unit tests for gtdb_amrfinder_organism.py.

Everything here is pure text handling, so no GTDB-Tk, no AMRFinderPlus and no
database are needed. The summary files are built on the fly in pytest's tmp_path
with the real GTDB-Tk 2.x column layout, so the tests also prove the column
lookup works when `classification` is not the first column.

What matters biologically: a WRONG --organism makes AMRFinderPlus score the
genome against another organism's curated point-mutation list, which is worse
than no --organism at all. So the bulk of these tests are about refusing to
match, not about matching.

Run: pytest workflow/scripts/mobilome/test_gtdb_amrfinder_organism.py -q
"""

import os

import gtdb_amrfinder_organism as gao


# ── Fixtures: build realistic GTDB-Tk summary files ──────────────────────────

# The real GTDB-Tk 2.x bac120 summary header, copied from
# 03.taxonomy/{sample}/gtdbtk.bac120.summary.tsv of a v2 validation run.
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
    # `amrfinder --list_organisms`, docs/mobilome_wpA_ground_truth.md).
    assert len(gao.ALL_ORGANISMS) == 31
    assert len(gao.SPECIES_ORGANISMS) == 28
    assert len(gao.GENUS_ORGANISMS) == 3


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


def test_genus_level_match_for_the_three_genus_organisms():
    # AMRFinderPlus curates these per genus, so any species in them qualifies.
    for genus, species in [
        ("Escherichia", "Escherichia coli"),
        ("Salmonella", "Salmonella enterica"),
        ("Campylobacter", "Campylobacter jejuni"),
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
    # The headline case for this project: Pseudomonas_E is a lineage GTDB split
    # off from Pseudomonas. It is NOT P. aeruginosa, so no --organism.
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


def test_placeholder_species_in_a_curated_genus():
    # Klebsiella has curated species but no genus-level entry, so an unnamed
    # Klebsiella gets nothing — and the audit says exactly why.
    organism, reason = gao.map_classification(
        lineage("Klebsiella", "Klebsiella sp900493575"))
    assert organism == ""
    assert reason.startswith("GTDB placeholder species")


def test_placeholder_species_in_an_uncurated_genus_gets_the_plain_reason():
    # Arthrobacter (a real isolate in this project) has nothing curated at all;
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
    # Paenibacillus: the third genus in this project's screening set.
    organism, reason = gao.map_classification(
        lineage("Paenibacillus", "Paenibacillus sp012647845"))
    assert organism == ""
    assert "Paenibacillus" in reason


# ── Reading the summary file ─────────────────────────────────────────────────

def test_single_row_summary(tmp_path):
    summary = write_summary(tmp_path, [
        ("006", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "006")
    assert organism == "Pseudomonas_aeruginosa"
    assert "s__Pseudomonas aeruginosa" in classification
    assert reason.startswith("exact species match")


def test_multi_row_summary_from_hybrid_mode(tmp_path):
    # Hybrid mode classifies both assemblies of one isolate; they agree, so the
    # single shared answer is used and the audit records that they agreed.
    summary = write_summary(tmp_path, [
        ("006_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("006_ont", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "006")
    assert organism == "Klebsiella_pneumoniae"
    assert "2 assemblies of this sample agree" in reason
    assert "006_illumina" in reason and "006_ont" in reason
    # Identical classifications are reported once, not twice.
    assert classification.count("d__Bacteria") == 1


def test_multi_row_summary_agreeing_on_no_organism(tmp_path):
    # The common case in this project: both assemblies say Pseudomonas_E.
    summary = write_summary(tmp_path, [
        ("015_illumina", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
        ("015_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
    ])
    organism, _, reason = gao.decide_organism(summary, "015")
    assert organism == ""
    assert "Pseudomonas_E" in reason


def test_multi_row_summary_that_disagrees_emits_no_organism(tmp_path):
    # Two assemblies of one isolate landing on different organisms means the
    # taxonomy is not solid enough to pick a point-mutation list.
    summary = write_summary(tmp_path, [
        ("006_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("006_ont", lineage("Klebsiella", "Klebsiella oxytoca")),
    ])
    organism, classification, reason = gao.decide_organism(summary, "006")
    assert organism == ""
    assert "classified differently" in reason
    assert "006_illumina=Klebsiella_pneumoniae" in reason
    assert "006_ont=Klebsiella_oxytoca" in reason
    # Both lineages are kept in the audit so the disagreement is visible.
    assert classification.count("d__Bacteria") == 2


def test_disagreement_between_a_match_and_no_match(tmp_path):
    summary = write_summary(tmp_path, [
        ("386_illumina", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
        ("386_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp010095445")),
    ])
    organism, _, reason = gao.decide_organism(summary, "386")
    assert organism == ""
    assert "386_ont=none" in reason


def test_rows_of_other_genomes_are_ignored(tmp_path):
    # Belt and braces: if a summary ever held more than one isolate, only the
    # requested sample's rows may drive the decision.
    summary = write_summary(tmp_path, [
        ("015_illumina", lineage("Klebsiella", "Klebsiella oxytoca")),
        ("006_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("006_ont", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, _, reason = gao.decide_organism(summary, "006")
    assert organism == "Klebsiella_pneumoniae"
    assert "015_illumina" not in reason


def test_sample_name_prefix_does_not_leak_into_another_sample(tmp_path):
    # "006" must not swallow a row named "0061": matching is on the exact name
    # or the name plus an underscore.
    summary = write_summary(tmp_path, [
        ("0061", lineage("Klebsiella", "Klebsiella oxytoca")),
        ("006", lineage("Klebsiella", "Klebsiella pneumoniae")),
    ])
    organism, _, _ = gao.decide_organism(summary, "006")
    assert organism == "Klebsiella_pneumoniae"


def test_row_not_named_after_the_sample_is_still_used(tmp_path):
    # GTDB-Tk is run per sample into its own directory, so a row named after the
    # contigs file (not the sample) still belongs to this isolate — use it, and
    # note that we did.
    summary = write_summary(tmp_path, [
        ("contigs_final", lineage("Vibrio", "Vibrio cholerae")),
    ])
    organism, _, reason = gao.decide_organism(summary, "006")
    assert organism == "Vibrio_cholerae"
    assert "used every row" in reason


def test_truncated_line_is_skipped_not_fatal(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    good = [""] * len(GTDBTK_HEADER)
    good[0] = "006"
    good[1] = lineage("Vibrio", "Vibrio vulnificus")
    summary.write_text("\n".join([
        "\t".join(GTDBTK_HEADER),
        "006_broken",              # a single field: an interrupted GTDB-Tk write
        "\t".join(good),
    ]) + "\n")
    organism, _, _ = gao.decide_organism(str(summary), "006")
    assert organism == "Vibrio_vulnificus"


def test_windows_line_endings_are_handled(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    row = [""] * len(GTDBTK_HEADER)
    row[0] = "006"
    row[1] = lineage("Helicobacter", "Helicobacter pylori")
    summary.write_text("\r\n".join([
        "\t".join(GTDBTK_HEADER), "\t".join(row)]) + "\r\n")
    organism, _, _ = gao.decide_organism(str(summary), "006")
    assert organism == "Helicobacter_pylori"


# ── Graceful degradation: nothing here may raise or exit non-zero ────────────

def test_missing_summary_file(tmp_path):
    organism, classification, reason = gao.decide_organism(
        str(tmp_path / "does_not_exist.tsv"), "006")
    assert organism == ""
    assert classification == ""
    assert "not found" in reason


def test_empty_summary_file(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    summary.write_text("")
    organism, _, reason = gao.decide_organism(str(summary), "006")
    assert organism == ""
    assert "empty" in reason


def test_header_only_summary_file(tmp_path):
    summary = write_summary(tmp_path, [])
    organism, _, reason = gao.decide_organism(summary, "006")
    assert organism == ""
    assert "no usable rows" in reason


def test_summary_without_the_expected_columns(tmp_path):
    summary = tmp_path / "gtdbtk.bac120.summary.tsv"
    summary.write_text("genome\tlineage\nfoo\tbar\n")
    organism, _, reason = gao.decide_organism(str(summary), "006")
    assert organism == ""
    assert "no user_genome/classification columns" in reason


# ── Finding the summary inside a GTDB-Tk output directory ───────────────────

def test_find_summary_at_the_top_level(tmp_path):
    # classify_wf leaves a copy/symlink of the summary at the top level.
    write_summary(tmp_path, [("006", lineage("Vibrio", "Vibrio cholerae"))])
    found, problem = gao.find_summary_in_dir(str(tmp_path))
    assert problem == ""
    assert found.endswith("gtdbtk.bac120.summary.tsv")


def test_find_summary_in_the_classify_subdirectory(tmp_path):
    # The real file lives in classify/; only the symlink is at the top level,
    # and a rule may be pointed at a directory where it was never made.
    classify = tmp_path / "classify"
    classify.mkdir()
    write_summary(classify, [("006", lineage("Vibrio", "Vibrio cholerae"))])
    found, problem = gao.find_summary_in_dir(str(tmp_path))
    assert problem == ""
    assert found.endswith(os.path.join("classify", "gtdbtk.bac120.summary.tsv"))


def test_bacterial_summary_wins_over_the_archaeal_one(tmp_path):
    # BacFlux is a bacterial workflow; when GTDB-Tk wrote both, read bac120.
    write_summary(tmp_path, [("006", lineage("Vibrio", "Vibrio cholerae"))],
                  name="gtdbtk.ar53.summary.tsv")
    write_summary(tmp_path, [("006", lineage("Klebsiella", "Klebsiella oxytoca"))])
    found, _ = gao.find_summary_in_dir(str(tmp_path))
    assert found.endswith("gtdbtk.bac120.summary.tsv")


def test_directory_without_a_summary_is_not_fatal(tmp_path):
    empty_dir = tmp_path / "03.taxonomy" / "006"
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
    taxonomy_dir = tmp_path / "03.taxonomy" / "006"
    classify = taxonomy_dir / "classify"
    classify.mkdir(parents=True)
    write_summary(classify, [
        ("006_illumina", lineage("Escherichia", "Escherichia coli")),
        ("006_ont", lineage("Escherichia", "Escherichia coli")),
    ])
    organism_file = tmp_path / "org.txt"
    audit_file = tmp_path / "audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-dir", str(taxonomy_dir),
        "--sample", "006",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    assert organism_file.read_text().strip() == "Escherichia"


def test_cli_requires_exactly_one_input_form(tmp_path):
    # Bad wiring of the rule (both forms, or neither) should be a loud argparse
    # error, not a silent "no organism" that looks like a real result.
    import pytest

    base = ["--sample", "006",
            "--out-organism", str(tmp_path / "org.txt"),
            "--out-audit", str(tmp_path / "audit.tsv")]
    with pytest.raises(SystemExit):
        gao.main(base)
    with pytest.raises(SystemExit):
        gao.main(["--gtdbtk-summary", "a.tsv", "--gtdbtk-dir", "b"] + base)


def test_cli_directory_form_degrades_gracefully(tmp_path):
    taxonomy_dir = tmp_path / "03.taxonomy" / "006"
    taxonomy_dir.mkdir(parents=True)
    organism_file = tmp_path / "org.txt"
    audit_file = tmp_path / "audit.tsv"

    exit_code = gao.main([
        "--gtdbtk-dir", str(taxonomy_dir),
        "--sample", "006",
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
        ("006", lineage("Pseudomonas", "Pseudomonas aeruginosa")),
    ])
    organism_file = tmp_path / "out" / "006_amrfinder_organism.txt"
    audit_file = tmp_path / "out" / "006_amrfinder_organism_decision.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "006",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])

    assert exit_code == 0
    # The rule reads this with $(cat ...), so the trailing newline is harmless.
    assert organism_file.read_text().strip() == "Pseudomonas_aeruginosa"

    header, row = read_audit(str(audit_file))
    assert header == ["sample", "gtdb_classification", "matched_organism", "reason"]
    assert row["sample"] == "006"
    assert row["matched_organism"] == "Pseudomonas_aeruginosa"
    assert "s__Pseudomonas aeruginosa" in row["gtdb_classification"]
    assert row["reason"].startswith("exact species match")


def test_cli_writes_an_empty_file_and_exits_zero_on_no_match(tmp_path):
    # The common case for environmental isolates: the file must be genuinely
    # empty so `[ -n "$(cat ...)" ]` in the rule is false and AMRFinderPlus runs
    # without --organism.
    summary = write_summary(tmp_path, [
        ("015_illumina", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
        ("015_ont", lineage("Pseudomonas_E", "Pseudomonas_E sp024807945")),
    ])
    organism_file = tmp_path / "015_amrfinder_organism.txt"
    audit_file = tmp_path / "015_amrfinder_organism_decision.tsv"

    exit_code = gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "015",
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
        "--sample", "006",
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
        ("006_illumina", lineage("Klebsiella", "Klebsiella pneumoniae")),
        ("006_ont", lineage("Klebsiella", "Klebsiella oxytoca")),
    ])
    audit_file = tmp_path / "audit.tsv"
    gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "006",
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
        ("006", lineage("Salmonella", "Salmonella enterica")),
    ])
    organism_file = tmp_path / "08.mobilome" / "006" / "organism.txt"
    audit_file = tmp_path / "08.mobilome" / "006" / "audit.tsv"
    gao.main([
        "--gtdbtk-summary", summary,
        "--sample", "006",
        "--out-organism", str(organism_file),
        "--out-audit", str(audit_file),
    ])
    assert organism_file.read_text().strip() == "Salmonella"


# ── Regression anchor against a real GTDB-Tk file ────────────────────────────

# A real v2 validation run: sample 006 is Pseudomonas_E, classified from both
# the Illumina and the ONT assembly. The expected answer is "no organism".
# Skipped automatically on any machine that does not have the validation data.
REAL_SUMMARY = ("/media/data/antonielli_dir/BacFlux_v2_validation/"
                "hybrid_screen_batch2/output_dir/03.taxonomy/006/"
                "gtdbtk.bac120.summary.tsv")


def test_real_gtdbtk_summary_gives_no_organism(tmp_path):
    if not os.path.isfile(REAL_SUMMARY):
        import pytest
        pytest.skip("validation data not present on this machine")

    organism, classification, reason = gao.decide_organism(REAL_SUMMARY, "006")
    assert organism == ""
    assert "g__Pseudomonas_E" in classification
    assert "Pseudomonas_E" in reason
    # Both hybrid rows were seen and agreed.
    assert "006_illumina" in reason and "006_ont" in reason
