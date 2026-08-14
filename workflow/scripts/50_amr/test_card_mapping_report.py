#!/usr/bin/env python3
"""Unit tests for workflow/scripts/50_amr/card_mapping_report.py.

No tools, no CARD download: the fixtures below are a handful of real BBMap
covstats lines (copied verbatim from a finished run, trailing space and all) and
a cut-down aro_index.tsv holding only the accessions those lines mention.

The script joins TWO BBMap passes over the same sample — one at a near-exact read
identity, one relaxed — so that a divergent member of a resistance family is
reported as `divergent` instead of vanishing. It then attaches CARD's own
classification, because plenty of CARD entries are efflux subunits, regulators or
porins rather than acquired resistance genes.

What is actually being checked here is the part that can silently go wrong. The
join key is an ARO accession dug out of a free-text defline; the classification
comes from CARD's own controlled vocabulary; and the "divergent" call depends on
two coverage figures being compared the right way round. Each of those is tested
against a case where getting it wrong would produce a plausible-looking but wrong
report rather than a crash.

Run either way:
    pytest workflow/scripts/50_amr/test_card_mapping_report.py
    python workflow/scripts/50_amr/test_card_mapping_report.py
"""

import os
import sys
import tempfile
import unittest

# The module under test sits in this same directory, so both pytest and a direct
# `python test_*.py` already have it on the import path — no sys.path juggling.
import card_mapping_report as cmr  # noqa: E402


# ── Real covstats rows at two identity filters, and a cut-down aro_index ────
# Real covstats rows from a Pseudomonas isolate, plus two invented ones so the
# fixture covers the categories that isolate happens not to carry. Note the
# trailing space after the closing bracket of the organism: BBMap keeps whatever
# whitespace CARD's FASTA defline had, and the parser must survive it.
COVSTATS_STRICT = "\n".join([
    "#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev",
    "gb|AE004091.2|+|2810008-2813197|ARO:3000804|MexF [Pseudomonas aeruginosa PAO1] \t39.57\t3189\t0.65\t100.0000\t3189\t448\t448\t0.67\t44\t15.72",
    "gb|AE004091.2|+|179521-182569|ARO:3003681|TriC [Pseudomonas aeruginosa PAO1] \t18.11\t3048\t0.65\t95.7677\t2919\t195\t193\t0.65\t10\t17.56",
    "gb|AE004091.2|+|4292297-4292998|ARO:3000829|CpxR [Escherichia coli] \t12.00\t702\t0.53\t99.0000\t695\t60\t60\t0.53\t12\t4.00",
    # An acquired beta-lactamase: the category the report exists to keep visible.
    "gb|FJ234049.1|+|0-861|ARO:3001864|OXA-58 [Acinetobacter baumannii] \t25.00\t861\t0.40\t98.0000\t844\t70\t70\t0.40\t25\t6.00",
    # An outer-membrane porin whose LOSS confers resistance. Present here, which
    # means the susceptible state, not resistance.
    "gb|CP000647.1|+|0-1104|ARO:3003480|OmpK36 [Klebsiella pneumoniae] \t30.00\t1104\t0.55\t99.5000\t1098\t95\t95\t0.55\t30\t5.00",
    # Below the coverage threshold at both settings: must never reach the report.
    "gb|XX000000.1|+|0-800|ARO:3009999|GhostGene [Nowhere bacterium] \t2.00\t800\t0.50\t11.0000\t88\t5\t5\t0.50\t2\t1.00",
    "",
])

# Same sample at the relaxed identity filter. Two deliberate differences:
#   * OXA-58's coverage is unchanged (it matched exactly, so it stays `exact`);
#   * a second beta-lactamase, CTX-M-15, only clears the threshold here. That is
#     the `divergent` case — a variant of the family is present, but not this
#     reference allele.
COVSTATS_RELAXED = "\n".join([
    "#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev",
    "gb|AE004091.2|+|2810008-2813197|ARO:3000804|MexF [Pseudomonas aeruginosa PAO1] \t41.00\t3189\t0.65\t100.0000\t3189\t460\t460\t0.67\t45\t15.00",
    "gb|AE004091.2|+|179521-182569|ARO:3003681|TriC [Pseudomonas aeruginosa PAO1] \t19.00\t3048\t0.65\t97.0000\t2956\t200\t200\t0.65\t11\t17.00",
    "gb|AE004091.2|+|4292297-4292998|ARO:3000829|CpxR [Escherichia coli] \t12.50\t702\t0.53\t99.0000\t695\t62\t62\t0.53\t12\t4.00",
    "gb|FJ234049.1|+|0-861|ARO:3001864|OXA-58 [Acinetobacter baumannii] \t25.00\t861\t0.40\t98.0000\t844\t70\t70\t0.40\t25\t6.00",
    "gb|CP000647.1|+|0-1104|ARO:3003480|OmpK36 [Klebsiella pneumoniae] \t30.00\t1104\t0.55\t99.5000\t1098\t95\t95\t0.55\t30\t5.00",
    "gb|AY044436.1|+|0-876|ARO:3001872|CTX-M-15 [Escherichia coli] \t14.00\t876\t0.52\t92.0000\t806\t40\t40\t0.52\t14\t7.00",
    "gb|XX000000.1|+|0-800|ARO:3009999|GhostGene [Nowhere bacterium] \t3.00\t800\t0.50\t19.0000\t152\t8\t8\t0.50\t3\t1.00",
    "",
])

# Cut-down aro_index.tsv. The column NAMES matter (the parser looks them up by
# name, not position); the column ORDER here is deliberately not CARD's, to prove
# that the lookup really is by name.
ARO_INDEX = "\n".join([
    "ARO Accession\tARO Name\tAMR Gene Family\tDrug Class\tResistance Mechanism",
    "ARO:3000804\tMexF\tresistance-nodulation-cell division (RND) antibiotic efflux pump\tfluoroquinolone antibiotic\tantibiotic efflux",
    "ARO:3003681\tTriC\tresistance-nodulation-cell division (RND) antibiotic efflux pump\ttriclosan\tantibiotic efflux",
    "ARO:3000829\tCpxR\tresistance-nodulation-cell division (RND) antibiotic efflux pump\taminoglycoside antibiotic\tantibiotic efflux",
    "ARO:3001864\tOXA-58\tOXA beta-lactamase\tcarbapenem\tantibiotic inactivation",
    "ARO:3003480\tOmpK36\tgeneral bacterial porin with reduced permeability to beta-lactams\tcephalosporin\tresistance by absence",
    "ARO:3001872\tCTX-M-15\tCTX-M beta-lactamase\tcephalosporin\tantibiotic inactivation",
    "ARO:3009999\tGhostGene\tsome family\tsome drug\tantibiotic target alteration",
    "",
])


def write_temp(content, suffix=".tsv"):
    """Write content to a throwaway temp file and return its path."""
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, encoding="utf-8"
    )
    handle.write(content)
    handle.close()
    return handle.name


def rows_by_name(rows):
    """Index a list of report rows by ARO name, for readable assertions."""
    return {row["aro_name"]: row for row in rows}


def build_default_rows():
    """Run the whole join over the fixtures with the workflow's real settings."""
    strict = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
    relaxed = cmr.parse_covstats(write_temp(COVSTATS_RELAXED))
    index = cmr.parse_aro_index(write_temp(ARO_INDEX))
    return cmr.build_rows(strict, relaxed, index, 70.0, "99", "95")


# ── Reading BBMap covstats and CARD's index ─────────────────────────────────
# Both parsers stop the run on a file of the wrong shape rather than returning
# nothing, because an empty covstats or a mis-read index would still produce a
# perfectly formatted report — one saying the genome carries no AMR gene at all.
# The single thing they do tolerate is a truncated FINAL line: losing that one
# row beats losing the report.
class TestParsers(unittest.TestCase):

    def test_covstats_ids_keep_the_full_defline_minus_trailing_space(self):
        # The join key is dug out of this string later, so it must arrive intact.
        # BBMap leaves CARD's trailing space in place; strip() is what makes the
        # strict and relaxed dictionaries share keys at all.
        coverage = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
        key = "gb|AE004091.2|+|2810008-2813197|ARO:3000804|MexF [Pseudomonas aeruginosa PAO1]"
        self.assertIn(key, coverage)
        self.assertEqual(coverage[key], 100.0)

    def test_header_and_blank_lines_are_not_rows(self):
        coverage = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
        self.assertEqual(len(coverage), 6)

    def test_a_truncated_final_row_is_skipped_not_fatal(self):
        # A killed BBMap can leave a half-written last line. Losing that one row
        # is preferable to losing the whole report.
        truncated = COVSTATS_STRICT + "gb|AE004091.2|+|1-2|ARO:30008"
        coverage = cmr.parse_covstats(write_temp(truncated))
        self.assertEqual(len(coverage), 6)

    def test_an_empty_covstats_is_fatal_not_an_empty_report(self):
        # map_amr_db runs without `set -e`, so a BBMap pass that dies still lets
        # the rule exit 0 with an empty covstats behind it. Reading that leniently
        # would produce a tidy report saying the genome carries no AMR gene at
        # all — wrong, plausible and silent. It has to stop here instead.
        with self.assertRaises(ValueError):
            cmr.parse_covstats(write_temp(""))

    def test_a_covstats_with_rows_but_no_header_is_also_fatal(self):
        # Same failure, caught one step later: a half-written file.
        headerless = "\n".join(COVSTATS_STRICT.split("\n")[1:])
        with self.assertRaises(ValueError):
            cmr.parse_covstats(write_temp(headerless))

    def test_a_genome_with_no_covered_gene_is_still_a_valid_table(self):
        # The opposite case, and a completely normal result: BBMap ran fine and
        # every reference came back at low coverage. That must parse.
        coverage = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
        self.assertTrue(coverage)

    def test_aro_index_is_read_by_column_name_not_position(self):
        # The fixture's column order is not CARD's. If the parser were positional
        # it would put the gene family where the drug class belongs.
        index = cmr.parse_aro_index(write_temp(ARO_INDEX))
        self.assertEqual(index["ARO:3001864"]["gene_family"], "OXA beta-lactamase")
        self.assertEqual(index["ARO:3001864"]["drug_class"], "carbapenem")
        self.assertEqual(index["ARO:3001864"]["mechanism"], "antibiotic inactivation")

    def test_a_file_that_is_not_aro_index_fails_loudly(self):
        # Pointing the rule at the wrong CARD file must not yield an empty index
        # and a report where every row says NA.
        with self.assertRaises(ValueError):
            cmr.parse_aro_index(write_temp("col1\tcol2\nfoo\tbar\n"))


# ── Telling a resistance gene from a pump subunit or a porin ────────────────
# The rules are ORDERED, and the order is the whole design. "Resistance by
# absence" is checked first because it inverts the meaning of the row: those
# genes confer resistance when they are LOST, so detecting one by presence is
# evidence of the susceptible state. After that, CARD's stated mechanism outranks
# the naming convention that a trailing R means a regulator. An environmental
# Gram-negative encodes whole RND repertoires, so lumping those in with acquired
# genes badly overstates what the genome carries.
class TestClassify(unittest.TestCase):

    def test_resistance_by_absence_wins_over_every_other_rule(self):
        # OmpK36 is a porin, so it would also match the efflux/structural wording
        # if the order were wrong. Its presence means the SUSCEPTIBLE state, and
        # that has to be the label the reader sees.
        category, notes = cmr.classify(
            "OmpK36",
            "general bacterial porin with reduced permeability to beta-lactams",
            "resistance by absence",
        )
        self.assertEqual(category, "presence_indicates_susceptibility")
        self.assertIn("ABSENCE", notes[0])

    def test_rnd_pump_subunit_is_flagged_as_an_efflux_component(self):
        category, notes = cmr.classify(
            "MexF",
            "resistance-nodulation-cell division (RND) antibiotic efflux pump",
            "antibiotic efflux",
        )
        self.assertEqual(category, "efflux_component")
        self.assertIn("multi-subunit", notes[0])

    def test_a_name_ending_in_R_is_flagged_as_a_regulator(self):
        category, _ = cmr.classify("CpxR", "some efflux family", "antibiotic efflux")
        self.assertEqual(category, "regulator")

    def test_a_determinant_whose_name_ends_in_R_is_not_demoted(self):
        # Found on a real Bacillus from this project's test set: vmlR is an ABC-F
        # ribosomal protection protein and a genuine resistance gene, and the bare
        # -R convention filed it as a regulator. CARD's mechanism is what settles
        # it — target protection, not efflux, so the convention does not apply.
        category, notes = cmr.classify(
            "vmlR",
            "Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal "
            "protection proteins",
            "antibiotic target protection",
        )
        self.assertEqual(category, "resistance_determinant")
        self.assertEqual(notes, [])

    def test_the_regulator_rule_does_not_fire_mid_name(self):
        # Only a trailing R counts. "OXA-58" and "arnA" must not be regulators,
        # and neither must a family name that merely contains an R.
        self.assertEqual(
            cmr.classify("OXA-58", "OXA beta-lactamase", "antibiotic inactivation")[0],
            "resistance_determinant",
        )
        self.assertEqual(
            cmr.classify("arnA", "pmr phosphoethanolamine transferase",
                         "antibiotic target alteration")[0],
            "resistance_determinant",
        )

    def test_efflux_outside_the_known_families_is_still_marked_as_efflux(self):
        # tet(A) is an efflux pump and genuinely mobile. It must not be filed with
        # the chromosomal RND machinery, but it must not be lost either.
        category, notes = cmr.classify(
            "tet(A)", "major facilitator superfamily (MFS) antibiotic efflux pump",
            "antibiotic efflux",
        )
        self.assertEqual(category, "efflux_component")
        category, notes = cmr.classify(
            "vgaC", "ABC-F ATP-binding cassette ribosomal protection protein",
            "antibiotic efflux",
        )
        self.assertEqual(category, "efflux_other")
        self.assertEqual(notes, [])

    def test_an_acquired_beta_lactamase_is_a_plain_determinant(self):
        category, notes = cmr.classify(
            "CTX-M-15", "CTX-M beta-lactamase", "antibiotic inactivation")
        self.assertEqual(category, "resistance_determinant")
        self.assertEqual(notes, [])

    def test_several_semicolon_separated_mechanisms_are_all_seen(self):
        # CARD entries can carry more than one mechanism; a plain equality test
        # against the whole string would miss the second one.
        category, _ = cmr.classify(
            "SomeGene", "some family",
            "antibiotic target alteration;antibiotic efflux")
        self.assertEqual(category, "efflux_other")


# ── Joining the two passes into one row per reference ───────────────────────
# A reference reaching the coverage threshold at EITHER identity setting gets a
# row, carrying both coverage figures. Seen at the strict setting it is `exact`;
# seen only at the relaxed one it is `divergent` — a variant of the family is
# present, but not this reference allele.
class TestBuildRows(unittest.TestCase):

    def test_only_sequences_clearing_the_threshold_are_reported(self):
        # Five from the strict pass plus CTX-M-15, which only the relaxed pass
        # sees. GhostGene clears neither and must not appear.
        rows = rows_by_name(build_default_rows())
        self.assertNotIn("GhostGene", rows)
        self.assertEqual(len(rows), 6)

    def test_a_gene_seen_only_at_the_relaxed_filter_is_called_divergent(self):
        rows = rows_by_name(build_default_rows())
        self.assertIn("CTX-M-15", rows)
        self.assertEqual(rows["CTX-M-15"]["detection"], "divergent")
        self.assertEqual(rows["CTX-M-15"]["covered_percent_id99"], "0.00")
        self.assertEqual(rows["CTX-M-15"]["covered_percent_id95"], "92.00")

    def test_a_gene_seen_at_the_strict_filter_is_called_exact(self):
        rows = rows_by_name(build_default_rows())
        self.assertEqual(rows["OXA-58"]["detection"], "exact")
        self.assertEqual(rows["OXA-58"]["covered_percent_id99"], "98.00")

    def test_the_aro_accession_is_pulled_from_the_defline_by_pattern(self):
        # Positional splitting would break the moment CARD changed the number of
        # pipe-separated fields, and would do so silently.
        rows = rows_by_name(build_default_rows())
        self.assertEqual(rows["MexF"]["aro_accession"], "ARO:3000804")
        self.assertEqual(rows["CTX-M-15"]["aro_accession"], "ARO:3001872")

    def test_the_reference_organism_is_reported_and_is_not_the_sample(self):
        # Kept because a report on a soil isolate that says "MexF" is much easier
        # to read when it also says the reference came from P. aeruginosa.
        rows = rows_by_name(build_default_rows())
        self.assertEqual(rows["MexF"]["reference_organism"], "Pseudomonas aeruginosa PAO1")

    def test_cards_own_classification_is_carried_through_unchanged(self):
        rows = rows_by_name(build_default_rows())
        self.assertEqual(rows["OXA-58"]["drug_class"], "carbapenem")
        self.assertEqual(rows["OXA-58"]["resistance_mechanism"], "antibiotic inactivation")
        self.assertEqual(rows["MexF"]["category"], "efflux_component")
        self.assertEqual(rows["OmpK36"]["category"], "presence_indicates_susceptibility")

    def test_an_accession_missing_from_the_index_still_produces_a_row(self):
        # CARD's FASTA and its index come from the same release, but a user can
        # point the two at different ones. That must degrade to an unclassified
        # row, not to a KeyError halfway through a 56-genome run.
        strict = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
        empty_index = cmr.parse_aro_index(write_temp("ARO Accession\tARO Name\n"))
        rows = rows_by_name(cmr.build_rows(strict, {}, empty_index, 70.0, "99", "95"))
        self.assertIn("MexF", rows)
        self.assertEqual(rows["MexF"]["amr_gene_family"], "NA")
        # With no mechanism to go on it falls through to the neutral label rather
        # than guessing.
        self.assertEqual(rows["MexF"]["category"], "resistance_determinant")

    def test_the_two_coverage_columns_are_named_after_the_filters(self):
        strict = cmr.parse_covstats(write_temp(COVSTATS_STRICT))
        relaxed = cmr.parse_covstats(write_temp(COVSTATS_RELAXED))
        index = cmr.parse_aro_index(write_temp(ARO_INDEX))
        rows = cmr.build_rows(strict, relaxed, index, 70.0, "98", "90")
        self.assertIn("covered_percent_id98", rows[0])
        self.assertIn("covered_percent_id90", rows[0])


# ── End to end, through main() ──────────────────────────────────────────────
# These drive the CLI exactly as the Snakemake rule does, and most of what they
# check is the SORT: the file is read from the top, so acquired determinants have
# to sit above the chromosomal efflux machinery, and a divergent hit must not be
# pushed to the bottom by the empty strict-coverage column it has by definition.
class TestMain(unittest.TestCase):

    def run_main(self, min_covered="70"):
        """Drive the CLI exactly as the Snakemake rule does; return the TSV lines."""
        out = write_temp("")
        argv = sys.argv
        sys.argv = [
            "card_mapping_report.py",
            "--sample", "sampleA",
            "--covstats-strict", write_temp(COVSTATS_STRICT),
            "--covstats-relaxed", write_temp(COVSTATS_RELAXED),
            "--aro-index", write_temp(ARO_INDEX),
            "--strict-id", "99",
            "--relaxed-id", "95",
            "--min-covered", min_covered,
            "--out", out,
        ]
        try:
            cmr.main()
        finally:
            sys.argv = argv
        with open(out, encoding="utf-8") as handle:
            return [line.rstrip("\n") for line in handle if line.strip()]

    def test_the_header_is_the_documented_column_set(self):
        lines = self.run_main()
        self.assertEqual(lines[0].split("\t"), [
            "aro_accession", "aro_name", "detection",
            "covered_percent_id99", "covered_percent_id95", "category",
            "resistance_mechanism", "amr_gene_family", "drug_class",
            "reference_organism", "note",
        ])

    def test_real_determinants_are_sorted_above_efflux_machinery(self):
        # The whole point of the sort: someone reading only the top of the file
        # must see the acquired genes, not the pump subunits.
        lines = self.run_main()
        names = [line.split("\t")[1] for line in lines[1:]]
        self.assertEqual(names[:2], ["OXA-58", "CTX-M-15"])
        self.assertEqual(names[-1], "OmpK36")

    def test_a_divergent_hit_is_not_buried_by_its_empty_strict_column(self):
        # CTX-M-15 has a strict coverage of 0 by definition. Ranking on the strict
        # column alone would drop it below every exact hit in its category, which
        # is the opposite of what the relaxed pass was added for. Here it beats
        # nothing else in its category, so the check is that it stays adjacent to
        # OXA-58 rather than falling to the end of the file.
        lines = self.run_main()
        names = [line.split("\t")[1] for line in lines[1:]]
        self.assertLess(names.index("CTX-M-15"), names.index("MexF"))

    def test_an_isolate_with_no_hits_writes_a_header_only_file(self):
        # A genome carrying nothing is a normal, common result. It must produce a
        # valid empty table, so the rule succeeds and the run continues.
        lines = self.run_main(min_covered="101")
        self.assertEqual(len(lines), 1)


if __name__ == "__main__":
    unittest.main()
