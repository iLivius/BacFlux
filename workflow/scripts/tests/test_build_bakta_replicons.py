#!/usr/bin/env python3
"""Unit tests for workflow/scripts/build_bakta_replicons.py.

These need no tools and no databases: tiny in-memory fixtures (a FASTA, a Flye
assembly_info.txt and a dnaapler reorientation summary) go through the parsers
and the join, and the resulting rows are checked.

Two of the assertions are INVARIANTS rather than behaviour checks, because
breaking either one makes Bakta exit with a single unhelpful line:
  * every replicon row has EXACTLY five fields;
  * the replicon file is never empty while the assembly has at least one record.

Run either way:
    pytest workflow/scripts/tests/test_build_bakta_replicons.py
    python workflow/scripts/tests/test_build_bakta_replicons.py
"""

import os
import sys
import tempfile
import unittest

# Import the module under test regardless of where the runner starts from.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import build_bakta_replicons as br  # noqa: E402


def write_temp(content, suffix=".txt"):
    """Write content to a throwaway temp file and return its path."""
    handle = tempfile.NamedTemporaryFile(
        mode="w", suffix=suffix, delete=False, encoding="utf-8"
    )
    handle.write(content)
    handle.close()
    return handle.name


# ── Fixtures ─────────────────────────────────────────────────────────────────
# Five contigs, each exercising one decision path:
#   contig_1  circular + strong dnaA        -> chromosome / circular
#   contig_2  circular + weak repA          -> contig / circular
#   contig_3  linear   + strong dnaA        -> chromosome / linear
#   contig_4  circular + terL (phage)       -> contig / circular
#   contig_5  in the FASTA only             -> contig / linear + no_flye_row
#
# The FASTA header of contig_1 deliberately carries a description, to prove the
# join keys on the first token only.
FASTA = (
    ">contig_1 length=5000000 circular=true\n"
    "ACGTACGTACGT\n"
    ">contig_2\n"
    "ACGTACGTACGT\n"
    ">contig_3\n"
    "ACGTACGTACGT\n"
    ">contig_4\n"
    "ACGTACGTACGT\n"
    ">contig_5\n"
    "ACGTACGTACGT\n"
)

# Real Flye 2.9 header, verified against an assembly_info.txt on disk.
FLYE_INFO = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "contig_1\t5000000\t96\tY\tN\t1\t*\t1\n"
    "contig_2\t45000\t120\tY\tN\t1\t*\t2\n"
    "contig_3\t300000\t80\tN\tN\t1\t*\t3\n"
    "contig_4\t40000\t60\tY\tN\t1\t*\t4\n"
)

DNAAPLER_HEADER = (
    "Contig\tGene_Reoriented\tStart\tStrand\tTop_Hit\tTop_Hit_Length\t"
    "Covered_Length\tCoverage\tIdentical_AAs\tIdentity_Percentage\t"
    "Overlapping_Contig_End\n"
)

# The numbers are the ones observed on this project's real isolates:
# a genuine DnaA covers ~100% of the reference; the marginal repA/terL hits do not.
DNAAPLER_SUMMARY = DNAAPLER_HEADER + (
    "contig_1 length=5000000\tdnaA\t1024681\treverse\tsp~P21173~DNAA_MICLU\t515\t518\t100.58\t354\t68.34\tFalse\n"
    "contig_2\trepA\t120\tforward\tUniRef90_Q4K3U4\t300\t43\t14.35\t22\t51.04\tFalse\n"
    "contig_3\tdnaA\t900\tforward\tsp~A5CLT3~DNAA_CLAM3\t477\t483\t101.26\t383\t79.30\tFalse\n"
    "contig_4\tterL\t77\tforward\tphrog_2\t600\t108\t17.98\t62\t57.53\tFalse\n"
)


def build(fasta=FASTA, flye=FLYE_INFO, dnaapler=DNAAPLER_SUMMARY,
          min_coverage=80.0, min_identity=40.0):
    """Run the whole parse+join chain on the given fixtures -> {contig: row}."""
    contig_ids = br.read_contig_ids(write_temp(fasta, ".fasta"))
    topology = br.parse_flye_info(write_temp(flye))
    markers = br.parse_dnaapler_summary(write_temp(dnaapler, ".tsv"))
    rows = br.build_rows("sampleA", contig_ids, topology, markers,
                         min_coverage, min_identity)
    return {row["contig"]: row for row in rows}, rows


class TestTypeAndTopology(unittest.TestCase):

    def test_circular_plus_strong_dnaa_is_chromosome(self):
        rows, _ = build()
        self.assertEqual(rows["contig_1"]["type"], "chromosome")
        self.assertEqual(rows["contig_1"]["topology"], "circular")
        self.assertIn("strong_dnaA", rows["contig_1"]["reason"])

    def test_circular_plus_weak_repa_falls_back_to_contig(self):
        # 14.35% coverage is far below the 80% floor, so the plasmid call is not
        # made — but the circularity, which is a measurement, is still reported.
        rows, _ = build()
        self.assertEqual(rows["contig_2"]["type"], "contig")
        self.assertEqual(rows["contig_2"]["topology"], "circular")
        self.assertIn("marker_below_threshold", rows["contig_2"]["reason"])

    def test_linear_plus_strong_dnaa_is_chromosome_linear(self):
        rows, _ = build()
        self.assertEqual(rows["contig_3"]["type"], "chromosome")
        self.assertEqual(rows["contig_3"]["topology"], "linear")

    def test_terl_is_never_typed(self):
        # Bakta has no phage replicon type, so terL stays neutral by design.
        rows, _ = build()
        self.assertEqual(rows["contig_4"]["type"], "contig")
        self.assertIn("marker_not_typed(terL)", rows["contig_4"]["reason"])

    def test_contig_absent_from_both_tables(self):
        rows, _ = build()
        self.assertEqual(rows["contig_5"]["type"], "contig")
        self.assertEqual(rows["contig_5"]["topology"], "linear")
        self.assertIn("no_flye_row", rows["contig_5"]["reason"])
        self.assertIn("no_dnaapler_row", rows["contig_5"]["reason"])

    def test_join_uses_first_token_of_the_dnaapler_contig_column(self):
        # dnaapler stores the full FASTA description; contig_1's row above is
        # "contig_1 length=5000000". If the first token were not used, contig_1
        # would come back untyped.
        rows, _ = build()
        self.assertEqual(rows["contig_1"]["marker"], "dnaA")


class TestDnaaplerStatusStrings(unittest.TestCase):
    """Each status string dnaapler can write into EVERY column must be safe."""

    def _row_for_status(self, status):
        summary = DNAAPLER_HEADER + "\t".join(["contig_1"] + [status] * 10) + "\n"
        rows, _ = build(dnaapler=summary)
        return rows["contig_1"]

    def test_contig_ignored(self):
        row = self._row_for_status("Contig_ignored")
        self.assertEqual(row["type"], "contig")
        self.assertIn("no_marker(Contig_ignored)", row["reason"])

    def test_no_mmseqs2_hits(self):
        row = self._row_for_status("No_MMseqs2_hits")
        self.assertEqual(row["type"], "contig")
        self.assertIn("no_marker(No_MMseqs2_hits)", row["reason"])

    def test_contig_already_reoriented(self):
        # The documented blind spot: the marker identity is lost, so no type.
        row = self._row_for_status("Contig_already_reoriented")
        self.assertEqual(row["type"], "contig")
        self.assertIn("no_marker(Contig_already_reoriented)", row["reason"])

    def test_autocomplete_method(self):
        row = self._row_for_status("autocomplete_method_mystery")
        self.assertEqual(row["type"], "contig")
        self.assertIn("no_marker(autocomplete_method_mystery)", row["reason"])


class TestGracefulDegradation(unittest.TestCase):

    def test_zero_overlap_join_is_fatal(self):
        # Simulates the silent failure this script exists to catch: a later step
        # renamed the contigs, so neither table joins.
        #
        # This case USED to be only a warning, and the row-building half below still
        # shows why that was not enough: with nothing joined, build_rows produces a
        # perfectly well-formed all-linear, all-"contig" table that Bakta would
        # accept without complaint, annotating exactly as if no replicon table had
        # been passed at all. Nothing downstream can tell that apart from a genome
        # that really is all linear contigs. So report_join_health was deliberately
        # changed to STOP the run on zero overlap; this test pins that down, because
        # a regression here would be invisible in the output.
        renamed_flye = FLYE_INFO.replace("contig_", "polished_")
        renamed_dnaapler = DNAAPLER_SUMMARY.replace("contig_", "polished_")
        rows, all_rows = build(flye=renamed_flye, dnaapler=renamed_dnaapler)
        for row in all_rows:
            self.assertEqual(row["type"], "contig")
            self.assertEqual(row["topology"], "linear")

        contig_ids = [row["contig"] for row in all_rows]
        topology = br.parse_flye_info(write_temp(renamed_flye))
        markers = br.parse_dnaapler_summary(write_temp(renamed_dnaapler, ".tsv"))
        with self.assertRaises(SystemExit) as caught:
            br.report_join_health(contig_ids, topology, markers)
        # The message must name the cause (renamed contigs), not just fail.
        self.assertIn("NONE of their IDs match", str(caught.exception))

    def test_partial_overlap_join_is_only_a_warning(self):
        # The counterpart to the test above, and the reason zero-overlap has to be
        # judged separately: dnaapler only reports the contigs it could reorient, so
        # a PARTIAL join is a normal biological outcome, not a broken pipeline. It
        # must stay non-fatal, or every ordinary genome would fail the run.
        contig_ids = br.read_contig_ids(write_temp(FASTA, ".fasta"))
        topology = br.parse_flye_info(write_temp(FLYE_INFO))
        partial_markers = br.parse_dnaapler_summary(
            write_temp(DNAAPLER_SUMMARY, ".tsv")
        )
        partial_markers.pop(next(iter(partial_markers)))  # drop one reoriented contig
        flye_matches, dnaapler_matches = br.report_join_health(
            contig_ids, topology, partial_markers
        )
        self.assertGreater(flye_matches, 0)
        self.assertGreater(dnaapler_matches, 0)

    def test_missing_input_files_do_not_crash(self):
        contig_ids = br.read_contig_ids(write_temp(FASTA, ".fasta"))
        topology = br.parse_flye_info("/nonexistent/assembly_info.txt")
        markers = br.parse_dnaapler_summary("/nonexistent/summary.tsv")
        rows = br.build_rows("sampleA", contig_ids, topology, markers, 80.0, 40.0)
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(row["type"] == "contig" for row in rows))
        self.assertTrue(all(row["topology"] == "linear" for row in rows))


class TestFileInvariants(unittest.TestCase):

    def test_every_replicon_row_has_exactly_five_fields(self):
        # Bakta unpacks each row into five variables inside a bare try/except;
        # any other count exits with "ERROR: wrong replicon table file format!".
        _, rows = build()
        out_path = write_temp("", ".tsv")
        br.write_replicons(out_path, rows)
        with open(out_path, encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle if line.strip()]
        self.assertEqual(len(lines), 5)
        for line in lines:
            self.assertEqual(len(line.split("\t")), 5, msg=line)

    def test_replicon_file_is_not_empty_when_the_assembly_has_records(self):
        # An empty file makes csv.Sniffer raise inside Bakta, which then hard-exits.
        _, rows = build()
        out_path = write_temp("", ".tsv")
        br.write_replicons(out_path, rows)
        self.assertGreater(os.path.getsize(out_path), 0)

    def test_no_header_line_in_the_replicon_file(self):
        _, rows = build()
        out_path = write_temp("", ".tsv")
        br.write_replicons(out_path, rows)
        with open(out_path, encoding="utf-8") as handle:
            first_line = handle.readline()
        self.assertTrue(first_line.startswith("contig_1\t"))

    def test_audit_file_has_a_header_and_one_row_per_contig(self):
        _, rows = build()
        out_path = write_temp("", ".tsv")
        br.write_audit(out_path, rows)
        with open(out_path, encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle if line.strip()]
        self.assertEqual(lines[0].split("\t"), br.AUDIT_COLUMNS)
        self.assertEqual(len(lines), 1 + 5)


if __name__ == "__main__":
    unittest.main()
