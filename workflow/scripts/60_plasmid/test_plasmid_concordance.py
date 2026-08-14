#!/usr/bin/env python3
"""Unit tests for workflow/scripts/60_plasmid/plasmid_concordance.py (D9).

These tests need no tools and no databases — they feed tiny in-memory fixtures
(Platon TSV, Platon chromosome FASTA, verified_plasmids.txt, geNomad plasmid
summary) through the parsers and the join, and assert the resulting rows and
confidence tiers. This is the one Stage-2b deliverable that is fully verifiable
in a plain Python environment.

What the script decides is how much a plasmid call is worth: Platon and geNomad
reach it by different routes, so where they agree the evidence is genuinely
stronger, and where they clash the row is FLAGGED rather than dropped. The tests
below walk one contig through each of those outcomes.

Run:
    pytest workflow/scripts/60_plasmid/test_plasmid_concordance.py

Running the file directly with python also works but executes only 18 of the 23
tests: the unittest.main() guard sits ABOVE the last class in this file, so it
exits before TestPlatonFailureIsNotASecondOpinion is even defined. pytest ignores
the guard and collects everything, which is why it is the way to run this.
"""

import os
import sys
import tempfile
import unittest

# The module under test sits in this same directory, so both pytest and a direct
# `python test_*.py` already have it on the import path — no sys.path juggling.
import plasmid_concordance as pc  # noqa: E402


def write_temp(content):
    """Write content to a throwaway temp file and return its path."""
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, encoding="utf-8"
    )
    handle.write(content)
    handle.close()
    return handle.name


# ── Four contigs, one per agreement tier ─────────────────────────────────────
# Minimal but realistic fixtures. Four contigs, chosen to exercise every tier:
#   contig_1 -> Platon plasmid + geNomad plasmid            -> both / high
#   contig_2 -> Platon plasmid only (geNomad absent)        -> platon_only / medium
#   contig_3 -> geNomad plasmid only (Platon never called)  -> genomad_only / medium
#   contig_4 -> Platon CHROMOSOME vs geNomad plasmid        -> conflict / low
# contig_5 exists only in the chromosome FASTA and is the negative control: both
# tools called it chromosomal, so it must not appear in the output at all.
PLATON_TSV = (
    "ID\tLength\tCoverage\t# ORFs\tRDS\tCircular\tInc Type(s)\t# Plasmid Hits\n"
    "contig_1\t45000\t12.0\t60\t0.85\tyes\tIncF\t3\n"
    "contig_2\t8000\t10.0\t12\t0.41\tno\t-\t0\n"
)

# geNomad calls contig_1 (agrees) and contig_3 (alone) and contig_4 (clash).
# contig_3's fdr is NA, contig_1's is populated — both must survive intact.
GENOMAD_SUMMARY = (
    "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\n"
    "contig_1\t45000\tno_terminal_repeats\t60\t11\t0.95\t0.012\t2\n"
    "contig_3\t3000\tno_terminal_repeats\t5\t11\t0.77\tNA\t0\n"
    "contig_4\t2000000\tno_terminal_repeats\t1800\t11\t0.60\tNA\t1\n"
)

# Platon called contig_4 a chromosome (so contig_4 is a real conflict, not just
# "not_called"). contig_3 appears nowhere on the Platon side.
CHROMOSOME_FASTA = (
    ">contig_4 length=2000000\n"
    "ACGTACGTACGT\n"
    ">contig_5 length=1500000\n"
    "TTTTGGGGCCCC\n"
)

# The kept v1 BLAST-text check: one hit, one no-hit (only Platon-plasmid contigs).
VERIFIED_PLASMIDS = (
    "sampleA: contig_1 is a plasmid.\n"
    "sampleA: contig_2 was not verified by BLAST search.\n"
)


# ── Parsing what each tool actually writes ───────────────────────────────────
# One shape has to fail loudly: a Platon TSV whose columns changed. A parser that
# quietly returned nothing there would produce a tidy concordance table saying
# the genome carries no plasmid — a wrong answer nobody would question.
# The other odd shapes below are all NORMAL and must be tolerated: an empty
# Platon table, the "no plasmid found" line, and a geNomad summary with no fdr
# column, which is what geNomad writes without score calibration.

class TestParsers(unittest.TestCase):
    def test_parse_platon_tsv_by_name(self):
        path = write_temp(PLATON_TSV)
        result = pc.parse_platon_tsv(path)
        self.assertEqual(result, {"contig_1": "0.85", "contig_2": "0.41"})

    def test_parse_platon_tsv_empty_file(self):
        # An isolate with no plasmids: Platon may write nothing / header-only.
        self.assertEqual(pc.parse_platon_tsv(write_temp("")), {})

    def test_parse_platon_tsv_missing_column_raises(self):
        # A changed format must fail loudly, not silently mis-join.
        bad = write_temp("ID\tLength\tCoverage\ncontig_1\t100\t9\n")
        with self.assertRaises(ValueError):
            pc.parse_platon_tsv(bad)

    def test_parse_fasta_ids_first_token(self):
        path = write_temp(CHROMOSOME_FASTA)
        self.assertEqual(pc.parse_fasta_ids(path), {"contig_4", "contig_5"})

    def test_parse_verified_plasmids_states(self):
        path = write_temp(VERIFIED_PLASMIDS)
        self.assertEqual(
            pc.parse_verified_plasmids(path),
            {"contig_1": "hit", "contig_2": "no_hit"},
        )

    def test_parse_verified_plasmids_no_plasmid_line(self):
        # The "none found" line carries no per-contig state and is skipped.
        path = write_temp("Platon found no plasmid in sample sampleA.\n")
        self.assertEqual(pc.parse_verified_plasmids(path), {})

    def test_parse_genomad_plasmids_with_and_without_fdr(self):
        path = write_temp(GENOMAD_SUMMARY)
        result = pc.parse_genomad_plasmids(path)
        self.assertEqual(result["contig_1"], ("0.95", "0.012"))
        self.assertEqual(result["contig_3"], ("0.77", "NA"))
        self.assertEqual(result["contig_4"], ("0.60", "NA"))

    def test_parse_genomad_plasmids_fdr_column_absent(self):
        # Without score calibration geNomad may omit the fdr column entirely.
        no_fdr = write_temp(
            "seq_name\tlength\tplasmid_score\n"
            "contig_1\t45000\t0.95\n"
        )
        result = pc.parse_genomad_plasmids(no_fdr)
        self.assertEqual(result["contig_1"], ("0.95", "NA"))


# ── Turning two calls into an agreement and a confidence ─────────────────────
# The tier is driven by AGREEMENT, not by either tool's own score, because the
# two reach a plasmid call by different routes: Platon from replicon-distribution
# scores, geNomad from gene content and markers. An unforeseen combination has to
# come back as a sentinel rather than raise, so one odd contig cannot lose the
# whole table.

class TestClassify(unittest.TestCase):
    def test_all_tiers(self):
        self.assertEqual(pc.classify("plasmid", "plasmid"), ("both", "high"))
        self.assertEqual(pc.classify("plasmid", "absent"), ("platon_only", "medium"))
        self.assertEqual(pc.classify("not_called", "plasmid"), ("genomad_only", "medium"))
        self.assertEqual(pc.classify("chromosome", "plasmid"), ("conflict", "low"))

    def test_unexpected_pair_is_sentinel_not_crash(self):
        self.assertEqual(pc.classify("chromosome", "absent"), ("undetermined", "low"))


# ── The joined table: which contigs appear, and what they carry ──────────────
# One row per plasmid CANDIDATE, meaning the union of the two tools' plasmid
# calls — so a contig both tools called chromosomal is simply absent, and a
# conflict is present and labelled rather than resolved.

class TestBuildRows(unittest.TestCase):
    def setUp(self):
        self.platon = pc.parse_platon_tsv(write_temp(PLATON_TSV))
        self.chrom = pc.parse_fasta_ids(write_temp(CHROMOSOME_FASTA))
        self.blast = pc.parse_verified_plasmids(write_temp(VERIFIED_PLASMIDS))
        self.genomad = pc.parse_genomad_plasmids(write_temp(GENOMAD_SUMMARY))
        self.rows = pc.build_rows(
            "sampleA", self.platon, self.chrom, self.blast, self.genomad
        )
        self.by_contig = {row["contig"]: row for row in self.rows}

    def test_row_set_is_union_of_plasmid_calls(self):
        # contig_5 is chromosome-only in both tools -> NOT a plasmid candidate.
        self.assertEqual(
            set(self.by_contig), {"contig_1", "contig_2", "contig_3", "contig_4"}
        )

    def test_rows_sorted_for_determinism(self):
        # Two runs on the same assembly must produce byte-identical tables, or a
        # re-run looks like a changed result when nothing changed.
        self.assertEqual(
            [row["contig"] for row in self.rows],
            ["contig_1", "contig_2", "contig_3", "contig_4"],
        )

    def test_both_high(self):
        row = self.by_contig["contig_1"]
        self.assertEqual(row["platon_call"], "plasmid")
        self.assertEqual(row["platon_rds"], "0.85")
        self.assertEqual(row["platon_blast_hit"], "hit")
        self.assertEqual(row["genomad_call"], "plasmid")
        self.assertEqual(row["genomad_score"], "0.95")
        self.assertEqual(row["genomad_fdr"], "0.012")
        self.assertEqual((row["agreement"], row["confidence"]), ("both", "high"))

    def test_platon_only_medium(self):
        row = self.by_contig["contig_2"]
        self.assertEqual(row["platon_call"], "plasmid")
        self.assertEqual(row["platon_blast_hit"], "no_hit")
        self.assertEqual(row["genomad_call"], "absent")
        self.assertEqual(row["genomad_score"], "NA")
        self.assertEqual((row["agreement"], row["confidence"]), ("platon_only", "medium"))

    def test_genomad_only_medium(self):
        row = self.by_contig["contig_3"]
        self.assertEqual(row["platon_call"], "not_called")
        self.assertEqual(row["platon_rds"], "NA")
        self.assertEqual(row["platon_blast_hit"], "NA")  # check never ran for it
        self.assertEqual(row["genomad_call"], "plasmid")
        self.assertEqual((row["agreement"], row["confidence"]), ("genomad_only", "medium"))

    def test_conflict_low_is_kept_not_dropped(self):
        row = self.by_contig["contig_4"]
        self.assertEqual(row["platon_call"], "chromosome")
        self.assertEqual(row["genomad_call"], "plasmid")
        self.assertEqual((row["agreement"], row["confidence"]), ("conflict", "low"))


# ── Writing the TSV and reading it back ──────────────────────────────────────
# The header must be written even when there is nothing to report: an isolate
# with no plasmid is an ordinary result, and a header-only file says "ran, found
# nothing" where a zero-byte file reads as a step that died — a distinction this
# module cares about elsewhere too (see the last class in this file).
# platon_replicons.py reads this table when geNomad was run, and copes with both
# the absent and the header-only case.

class TestWriteAndReadBack(unittest.TestCase):
    def test_roundtrip_tsv(self):
        rows = pc.build_rows(
            "sampleA",
            pc.parse_platon_tsv(write_temp(PLATON_TSV)),
            pc.parse_fasta_ids(write_temp(CHROMOSOME_FASTA)),
            pc.parse_verified_plasmids(write_temp(VERIFIED_PLASMIDS)),
            pc.parse_genomad_plasmids(write_temp(GENOMAD_SUMMARY)),
        )
        out = write_temp("")
        pc.write_rows(out, rows)
        with open(out, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], "\t".join(pc.OUTPUT_COLUMNS))
        self.assertEqual(len(lines), 1 + 4)  # header + 4 candidate contigs

    def test_empty_table_still_has_header(self):
        out = write_temp("")
        pc.write_rows(out, [])
        with open(out, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines, ["\t".join(pc.OUTPUT_COLUMNS)])


# This guard sits mid-file, so `python <this file>` stops here and never defines
# the class below it — see the note in the module docstring. Under pytest the
# guard is false and everything is collected.
if __name__ == "__main__":
    unittest.main(verbosity=2)


# ── A crashed Platon is not a second opinion ─────────────────────────────────

class TestPlatonFailureIsNotASecondOpinion(unittest.TestCase):
    """A crashed Platon must not be scored as a genuine second opinion.

    The distinction that matters, and that an earlier attempt at this fix got
    wrong: an EMPTY Platon directory is NOT a failure. A completed, exit-0 Platon
    run writes almost nothing when every contig exceeds its 500 kb size filter,
    which is the ordinary outcome for a closed genome. Only the explicit
    "Platon exited with status N" line that the plasmid_search rule writes means
    the tool actually crashed.
    """

    def test_crash_marker_is_detected(self):
        crashed = write_temp("sampleA: Platon exited with status 1; see the log.\n")
        self.assertTrue(pc.platon_run_failed(crashed))

    def test_a_normal_no_plasmid_run_is_not_a_failure(self):
        # What plasmid_search writes for a genome where Platon found no plasmid.
        normal = write_temp("Platon found no plasmid in sample sampleA.\n")
        self.assertFalse(pc.platon_run_failed(normal))

    def test_an_empty_platon_directory_is_not_a_failure(self):
        # The closed-genome case: nothing written at all, exit 0.
        empty = write_temp("")
        self.assertFalse(pc.platon_run_failed(empty))
        self.assertFalse(pc.platon_run_failed("/nonexistent/verified_plasmids.txt"))

    def test_a_crash_downgrades_genomad_only_to_platon_unavailable(self):
        # Without this, geNomad's calls would be reported as agreement=genomad_only
        # at medium confidence: a two-tool tier that was never assessed.
        rows = pc.build_rows(
            "sampleA", {}, set(), {},
            pc.parse_genomad_plasmids(write_temp(GENOMAD_SUMMARY)),
            platon_assessed=False,
        )
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual(row["platon_call"], "not_assessed")
            self.assertEqual(row["agreement"], "platon_unavailable")
            self.assertEqual(row["confidence"], "low")

    def test_classify_refuses_to_grade_an_unassessed_contig(self):
        self.assertEqual(pc.classify("not_assessed", "plasmid"),
                         ("platon_unavailable", "low"))
        self.assertEqual(pc.classify("not_assessed", "absent"),
                         ("platon_unavailable", "low"))
