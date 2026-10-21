#!/usr/bin/env python3
"""Unit tests for workflow/scripts/10_decontam/select_contigs_by_taxonomy.py.

The script under test decides which assembly contigs belong to the isolate and
which are contamination, so a mistake here removes genome or keeps a
contaminant before anything else runs. These tests cover the parts that pick a
genus:

  * the auto vote, on base pairs (auto_vote: bp, the default) and on contig
    counts (auto_vote: contigs);
  * the genus written to the kept-genus file, which Bakta receives as the
    organism name;
  * the composition report's line shape.

No tools and no databases: a made-up BlobTools table and FASTA are written to a
temporary folder and the script is run on them, the way rule select_contigs
runs it.

Run either way:
    pytest workflow/scripts/10_decontam/test_select_contigs_by_taxonomy.py
    python workflow/scripts/10_decontam/test_select_contigs_by_taxonomy.py
"""

import os
import subprocess
import sys
import tempfile
import unittest

# Import the module under test regardless of where the runner starts from.
import select_contigs_by_taxonomy as sel  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "select_contigs_by_taxonomy.py")


# ── A made-up assembly: one long chromosome, three short plasmids ────────────
# The shape of a closed long-read genome. BLAST assigns the chromosome to
# GenusK and the three plasmids to GenusE, so the two votes disagree:
#   by bp:      GenusK 5,000 bp  vs  GenusE 3 x 100 = 300 bp   -> GenusK
#   by contigs: GenusK 1 contig  vs  GenusE 3 contigs          -> GenusE
CONTIGS = {
    "chrom_1": ("GenusK", 5000),
    "plasmid_1": ("GenusE", 100),
    "plasmid_2": ("GenusE", 100),
    "plasmid_3": ("GenusE", 100),
}


def blobtools_row(contig, genus):
    """One line of a `blobtools view --rank all` table: the contig name in
    column 1 and the genus in column 22. The other columns are filler, because
    the selector reads only those two."""
    fields = ["0"] * 24
    fields[0] = contig
    fields[21] = genus
    return "\t".join(fields)


def write_inputs(folder, contigs):
    """Write the BlobTools table and the FASTA for a {name: (genus, length)}
    dict, and return their two paths."""
    table = os.path.join(folder, "bestscore.table.txt")
    fasta = os.path.join(folder, "contigs.fasta")
    with open(table, "w", encoding="utf-8") as handle:
        handle.write("# name\tlength\t...\n")
        for name, (genus, _) in contigs.items():
            handle.write(blobtools_row(name, genus) + "\n")
    with open(fasta, "w", encoding="utf-8") as handle:
        for name, (_, length) in contigs.items():
            handle.write(f">{name} length={length}\n{'A' * length}\n")
    return table, fasta


def run_selector(folder, contigs, extra_args):
    """Run the script as rule select_contigs does and return the output paths
    plus the completed process, so a test can read whichever file it checks."""
    table, fasta = write_inputs(folder, contigs)
    out = {
        "list": os.path.join(folder, "contigs.list"),
        "fasta": os.path.join(folder, "kept.fasta"),
        "composition": os.path.join(folder, "composition.txt"),
        "decisions": os.path.join(folder, "decisions.tsv"),
        "kept_genus": os.path.join(folder, "kept_genus.txt"),
    }
    command = [
        sys.executable, SCRIPT,
        "--bestscore", table,
        "--contigs", fasta,
        "--sample", "SAMPLE1",
        "--output-list", out["list"],
        "--output-fasta", out["fasta"],
        "--composition", out["composition"],
        "--decisions", out["decisions"],
    ] + extra_args
    out["process"] = subprocess.run(command, capture_output=True, text=True)
    return out


def read_text(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


# ── The vote itself ──────────────────────────────────────────────────────────

class TestChooseAutoGenus(unittest.TestCase):

    def test_largest_vote_wins(self):
        # The same function serves both votes; only the numbers differ.
        self.assertEqual(sel.choose_auto_genus({"GenusK": 5000, "GenusE": 300}, True), "GenusK")
        self.assertEqual(sel.choose_auto_genus({"GenusK": 1, "GenusE": 3}, True), "GenusE")

    def test_exact_tie_goes_to_the_alphabetically_first_genus(self):
        # Ties are broken by name only so that the same input always gives the
        # same answer.
        self.assertEqual(sel.choose_auto_genus({"Zeta": 10, "Alpha": 10}, True), "Alpha")

    def test_no_hit_cannot_win_when_discarded(self):
        votes = {"no-hit": 9000, "GenusK": 10}
        self.assertEqual(sel.choose_auto_genus(votes, discard_no_hit=True), "GenusK")
        self.assertEqual(sel.choose_auto_genus(votes, discard_no_hit=False), "no-hit")

    def test_nothing_to_vote_for(self):
        self.assertIsNone(sel.choose_auto_genus({}, True))
        self.assertIsNone(sel.choose_auto_genus({"no-hit": 5}, True))


class TestChooseLabelGenus(unittest.TestCase):

    def test_no_hit_never_names_the_genome(self):
        self.assertEqual(sel.choose_label_genus({"no-hit": 9000, "GenusK": 10}), "GenusK")

    def test_only_no_hit_kept_gives_no_name(self):
        self.assertIsNone(sel.choose_label_genus({"no-hit": 9000}))


# ── The whole script, the way the rule runs it ───────────────────────────────

class TestSelectorRun(unittest.TestCase):

    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def test_bp_vote_keeps_the_chromosome(self):
        out = run_selector(self.folder, CONTIGS, ["--mode", "auto", "--auto-vote", "bp",
                                                  "--kept-genus", os.path.join(self.folder, "kept_genus.txt")])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertEqual(read_text(out["list"]).split(), ["chrom_1"])
        self.assertEqual(read_text(out["kept_genus"]), "GenusK\n")

    def test_bp_is_the_default_vote(self):
        out = run_selector(self.folder, CONTIGS, ["--mode", "auto",
                                                  "--kept-genus", os.path.join(self.folder, "kept_genus.txt")])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertEqual(read_text(out["list"]).split(), ["chrom_1"])

    def test_contigs_vote_keeps_the_majority_of_contigs(self):
        out = run_selector(self.folder, CONTIGS, ["--mode", "auto", "--auto-vote", "contigs",
                                                  "--kept-genus", os.path.join(self.folder, "kept_genus.txt")])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertEqual(read_text(out["list"]).split(), ["plasmid_1", "plasmid_2", "plasmid_3"])
        self.assertEqual(read_text(out["kept_genus"]), "GenusE\n")

    def test_composition_lines_carry_bp_and_contig_shares(self):
        out = run_selector(self.folder, CONTIGS, ["--mode", "auto"])
        lines = read_text(out["composition"]).splitlines()
        # Sorted by bp share: 5,000 of 5,300 bp is 0.94; 1 of 4 contigs is 0.25.
        self.assertEqual(lines[0], "GenusK: bp 0.94  contigs 0.25")
        self.assertEqual(lines[1], "GenusE: bp 0.06  contigs 0.75")

    def test_mode_off_names_the_genus_holding_most_bp(self):
        # Mode off keeps everything; the name is the kept genus with most bp.
        out = run_selector(self.folder, CONTIGS, ["--mode", "off",
                                                  "--kept-genus", os.path.join(self.folder, "kept_genus.txt")])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertEqual(len(read_text(out["list"]).split()), 4)
        self.assertEqual(read_text(out["kept_genus"]), "GenusK\n")

    def test_all_no_hit_kept_writes_an_empty_genus_file(self):
        # When BLAST placed nothing and nothing is filtered, there is no genus
        # to give Bakta: the file exists and is empty.
        unplaced = {"contig_1": ("no-hit", 500), "contig_2": ("no-hit", 400)}
        out = run_selector(self.folder, unplaced, ["--mode", "off",
                                                   "--kept-genus", os.path.join(self.folder, "kept_genus.txt")])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertEqual(read_text(out["kept_genus"]), "")

    def test_kept_genus_file_is_optional(self):
        # The long-read final screen runs the script without --kept-genus.
        out = run_selector(self.folder, CONTIGS, ["--mode", "off"])
        self.assertEqual(out["process"].returncode, 0, out["process"].stderr)
        self.assertFalse(os.path.exists(out["kept_genus"]))


if __name__ == "__main__":
    unittest.main()
