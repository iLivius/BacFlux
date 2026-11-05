#!/usr/bin/env python3
"""Start-up test: an unquoted on/off in the config stops the run, never misreads it.

Snakemake reads the config file as YAML 1.1, where an UNQUOTED on, yes or true is
the boolean True and an unquoted off, no or false is the boolean False, not a word.
So `caller: off` reaches BacFlux as False. Left unchecked, several keys would
quietly turn such a value into their default (`phage.caller: off` would run
VirSorter2; `decontamination: off`, written as one value instead of a block,
would run auto), number keys would read it as 1 or 0 (`contig_boundary_bp: off`
would be a 0 bp window), and others would crash with a bare Python error. The
checks in workflow/rules/shared/00_common.smk (section 0, and _stop_on_yaml_bool
or _config_number where each key is read) stop the run instead, with a [BacFlux]
message and exit status 1. A number key also takes its default when left blank,
and stops on text that is not a number.

How the test works: it builds a small illumina config from the shipped template,
config/config.yaml, in a temporary directory that holds one empty stub sample
(SAMPLE1) and stub database paths. Each test writes that config with one value
changed, runs `snakemake -n` (a dry run: nothing is installed or executed) on it,
and checks the exit status and the message. The bad values are written into the
file as bare words, exactly as a user would type them, so the test goes through
the same YAML reading as a real run.

The last tests check the other direction: the untouched template still parses,
and the places where a bare on/off IS meaningful (decontamination `mode: off`,
`long_read_mode: off`, `medaka_model: off`) are still accepted, and a blank
`long_read_mode:` gets its default, off, as a blank `mode:` gets its own.

The opt-in layers (Baktfold, the trait profile) get tests of their own: left off,
as the template ships them, they add no job; switched on with stub databases,
they plan their jobs; and their database settings stop a run only when the layer
is on.

Snakemake is looked for next to the Python running the tests, then on PATH; the
tests are skipped when it is not found. Each dry run takes about a second.
Run with:
    pytest miscellaneous/test_config_yaml_booleans.py
"""

import copy
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO, "workflow", "Snakefile")
TEMPLATE = os.path.join(REPO, "config", "config.yaml")

# A placeholder the config is written with, then swapped for the bare word under
# test. PyYAML would quote a string such as "off" when writing it, which is the
# safe form and not the mistake this test is about, so the word is put into the
# file text afterwards instead.
BARE_WORD_SLOT = "BARE_WORD_SLOT"

# Every number key in the shipped template, by its dotted name. Each one must stop
# the run on a bare on/off. ConfigNumberKeyListTest checks this list against the
# template, and against the keys 00_common.smk and 80_mobilome.smk actually check.
NUMBER_KEYS = (
    "resources.threads", "resources.ram_gb", "resources.runtime_hours",
    "parameters.long_read_qc.min_length", "parameters.long_read_qc.keep_percent",
    "parameters.long_read_qc.length_weight",
    "mobilome.max_composite_span_bp", "mobilome.coverage_profile",
    "mobilome.tncentral.min_identity", "mobilome.tncentral.min_reference_coverage",
    "mobilome.iceberg.min_identity", "mobilome.iceberg.min_overlap_fraction",
    "mobilome.isosdb.min_covered_percent", "mobilome.isosdb.min_copy_number",
    "mobilome.contig_boundary_bp",
    "traits.max_gene_gap",
)

# The trait gene dictionary that ships with BacFlux. The trait-profile tests build
# their stub EcoFoldDB from the entry names it uses (see make_trait_stubs below).
GENE_TABLE = os.path.join(REPO, "workflow", "scripts", "40_annotation", "trait_genes_v1.tsv")


def ecofolddb_names_in_gene_table():
    """The EcoFoldDB entry names the shipped gene dictionary uses, sorted.

    Each name is Gene@Pathway/Activity, as EcoFoldDB's descriptions table builds it;
    a cell of the table's ecofolddb column holds none, one, or several separated by
    ";". Lines starting with "#" (the table's note on what it is) are skipped."""
    names = set()
    header = None
    with open(GENE_TABLE, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            for name in fields[header.index("ecofolddb")].split(";"):
                if name.strip():
                    names.add(name.strip())
    return sorted(names)


def find_snakemake():
    """The snakemake next to this Python (same conda env), or the one on PATH."""
    beside_python = os.path.join(os.path.dirname(sys.executable), "snakemake")
    if os.path.exists(beside_python):
        return beside_python
    return shutil.which("snakemake")


SNAKEMAKE = find_snakemake()


def set_key(config, dotted_name, value):
    """Set config["a"]["b"]["c"] from "a.b.c", creating empty blocks on the way."""
    keys = dotted_name.split(".")
    node = config
    for key in keys[:-1]:
        if not isinstance(node.get(key), dict):
            node[key] = {}
        node = node[key]
    node[keys[-1]] = value


@unittest.skipIf(SNAKEMAKE is None, "snakemake not found next to this Python or on PATH")
class TestConfigYamlBooleans(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # One temporary directory holds everything a dry run touches: the stub
        # inputs, the stub databases, the output path, the config variants, and
        # Snakemake's own .snakemake/ folder. Snakemake runs with this directory as
        # its working directory, so nothing is written into the repository and a
        # real run going on in the repository is never disturbed.
        cls.tmp = tempfile.TemporaryDirectory()
        root = cls.tmp.name

        # One empty stub sample per technology. Sample discovery only globs the
        # file names, so empty files are enough for a dry run.
        illumina_dir = os.path.join(root, "in", "illumina")
        nanopore_dir = os.path.join(root, "in", "nanopore")
        os.makedirs(illumina_dir)
        os.makedirs(nanopore_dir)
        for name in ("SAMPLE1_R1.fastq.gz", "SAMPLE1_R2.fastq.gz"):
            open(os.path.join(illumina_dir, name), "w").close()
        open(os.path.join(nanopore_dir, "SAMPLE1_ont.fastq.gz"), "w").close()

        # Database paths are not opened while the workflow is read, with one
        # exception: the BlobTools rule names two NCBI taxonomy files in blast_db
        # as inputs, so a dry run that builds the whole job list needs them to exist.
        db_dir = os.path.join(root, "db")
        blast_db = os.path.join(db_dir, "blast_db")
        os.makedirs(blast_db)
        for name in ("nodes.dmp", "names.dmp"):
            open(os.path.join(blast_db, name), "w").close()

        # The shipped template, with its empty placeholders filled in. The links
        # are only checked for being set, so placeholder addresses do.
        with open(TEMPLATE, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        config["mode"] = "illumina"
        config["input"]["illumina_dir"] = illumina_dir
        config["input"]["nanopore_dir"] = nanopore_dir
        config["directories"]["output_dir"] = os.path.join(root, "out")
        for key in ("bakta_db", "eggnog_db", "gtdbtk_db", "platon_db"):
            config["directories"][key] = os.path.join(db_dir, key)
        config["directories"]["blast_db"] = blast_db
        config["links"]["phix_link"] = "https://example.org/phix.fna.gz"
        config["links"]["card_link"] = "https://example.org/card.tar.bz2"
        cls.base_config = config
        cls.variant_count = 0

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def dry_run(self, changes, bare_word=None, extra_args=(), targets=()):
        """Write the base config with `changes` applied, dry-run it, return (status, output).

        changes:    {dotted key: value}. A value of BARE_WORD_SLOT is replaced in the
                    written file by `bare_word`, unquoted, e.g. `caller: off`.
        extra_args: more Snakemake options, e.g. ["-p"] to print the shell
                    command of every job the dry run lists.
        targets:    what to build instead of the whole workflow, e.g.
                    ["download_all"]. Written right after `snakemake`, because
                    --configfile takes several values and would read a target
                    written after it as one more config file.
        Returns:    Snakemake's exit status and everything it printed.
        """
        config = copy.deepcopy(self.base_config)
        for dotted_name, value in changes.items():
            set_key(config, dotted_name, value)
        text = yaml.safe_dump(config, sort_keys=False)
        if bare_word is not None:
            text = text.replace(BARE_WORD_SLOT, bare_word)

        type(self).variant_count += 1
        path = os.path.join(self.tmp.name, f"variant_{self.variant_count}.yaml")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)

        result = subprocess.run(
            [SNAKEMAKE, *targets, "-n", "--cores", "1", "--snakefile", SNAKEFILE,
             "--configfile", path, *extra_args],
            cwd=self.tmp.name, capture_output=True, text=True, timeout=600,
        )
        return result.returncode, result.stdout + result.stderr

    def assert_stops(self, changes, bare_word, *expected_fragments):
        """The run must stop with exit status 1 and a [BacFlux] message holding every fragment."""
        status, output = self.dry_run(changes, bare_word)
        self.assertEqual(status, 1, output[-3000:])
        for fragment in expected_fragments:
            self.assertIn(fragment, output, output[-3000:])

    # ── WORD keys: one word from a fixed list, never a boolean ──

    def test_phage_caller_off_stops_and_says_there_is_no_off_switch(self):
        self.assert_stops({"phage.caller": BARE_WORD_SLOT}, "off",
                          "[BacFlux] phage.caller is an unquoted off", "virsorter2 or genomad",
                          "Phage calling has no off switch")

    def test_phage_caller_quoted_off_or_other_word_stops_and_says_there_is_no_off_switch(self):
        # Quoted, "off" stays a word and passes the boolean check, so this reaches
        # the word check instead. It must give the same reminder, so the advice to
        # quote a word never leads to a shorter, less helpful message. "none" is
        # an ordinary word in YAML (only ~, null and an empty value are null).
        for word in ("off", "none"):
            self.assert_stops({"phage.caller": word}, None,
                              "[BacFlux] phage.caller must be virsorter2 or genomad",
                              f"(got: '{word}')", "Phage calling has no off switch")

    def test_phage_caller_on_stops(self):
        self.assert_stops({"phage.caller": BARE_WORD_SLOT}, "yes",
                          "[BacFlux] phage.caller is an unquoted on", "virsorter2 or genomad")

    def test_long_read_mode_on_stops_like_mode_on(self):
        self.assert_stops({"parameters.decontamination.long_read_mode": BARE_WORD_SLOT}, "on",
                          "[BacFlux] parameters.decontamination.long_read_mode is an unquoted on")
        self.assert_stops({"parameters.decontamination.mode": BARE_WORD_SLOT}, "on",
                          "[BacFlux] parameters.decontamination.mode is an unquoted on")

    def test_auto_vote_off_stops(self):
        self.assert_stops({"parameters.decontamination.auto_vote": BARE_WORD_SLOT}, "no",
                          "[BacFlux] parameters.decontamination.auto_vote", "bp or contigs")

    def test_spades_kmers_off_stops(self):
        self.assert_stops({"parameters.spades_kmers": BARE_WORD_SLOT}, "off",
                          "[BacFlux] parameters.spades_kmers is an unquoted off")

    def test_nt_version_off_stops(self):
        self.assert_stops({"parameters.nt_version": BARE_WORD_SLOT}, "off",
                          "[BacFlux] parameters.nt_version is an unquoted off")

    def test_pipeline_mode_off_stops(self):
        self.assert_stops({"mode": BARE_WORD_SLOT}, "off",
                          "[BacFlux] mode is an unquoted off",
                          "illumina, nanopore, hybrid or contigs")

    def test_flye_input_mode_off_stops_in_nanopore_mode(self):
        self.assert_stops({"mode": "nanopore", "parameters.nanopore.flye_input_mode": BARE_WORD_SLOT},
                          "off",
                          "[BacFlux] parameters.nanopore.flye_input_mode is an unquoted off",
                          "auto, nano-raw or nano-hq")

    # ── BLOCKS written as one value ──

    def test_decontamination_block_written_as_off_stops(self):
        self.assert_stops({"parameters.decontamination": BARE_WORD_SLOT}, "off",
                          "[BacFlux] parameters.decontamination is a block",
                          "write mode: off indented beneath it")

    def test_decontamination_block_written_as_quoted_off_stops(self):
        # Quoted, so YAML gives the word "off", not a boolean: still not a block.
        self.assert_stops({"parameters.decontamination": "off"}, None,
                          "[BacFlux] parameters.decontamination is a block",
                          "write mode: off indented beneath it")

    def test_phage_block_written_as_off_stops(self):
        self.assert_stops({"phage": BARE_WORD_SLOT}, "off",
                          "[BacFlux] phage is a block", "Phage calling has no off switch")

    def test_mobilome_block_written_as_on_stops(self):
        self.assert_stops({"mobilome": BARE_WORD_SLOT}, "on",
                          "[BacFlux] mobilome is a block", "run: true or run: false")

    def test_mobilome_layer_block_written_as_on_stops(self):
        self.assert_stops({"mobilome.tncentral": BARE_WORD_SLOT}, "on",
                          "[BacFlux] mobilome.tncentral is a block")

    # ── NUMBER keys: a number, never a boolean ──

    def test_keep_percent_on_stops_and_says_there_is_no_off(self):
        # keep_percent is handed to Filtlong unconverted, so without the check a
        # bare on would reach Filtlong as the text True.
        self.assert_stops({"mode": "nanopore",
                           "parameters.long_read_qc.keep_percent": BARE_WORD_SLOT}, "on",
                          "[BacFlux] parameters.long_read_qc.keep_percent is an unquoted on",
                          "not as a number", "such as 95",
                          "There is no off: 100 rejects no read by percentage.")

    def test_min_length_off_stops(self):
        # A whole-number key: int() would read off as 0.
        self.assert_stops({"mode": "nanopore",
                           "parameters.long_read_qc.min_length": BARE_WORD_SLOT}, "off",
                          "[BacFlux] parameters.long_read_qc.min_length is an unquoted off",
                          "not as a number", "a whole number of bases, such as 1000")

    def test_length_weight_no_stops(self):
        # A decimal key: float() would read no as 0.0. Checked in hybrid, the one
        # mode that hands it to Filtlong.
        self.assert_stops({"mode": "hybrid",
                           "parameters.long_read_qc.length_weight": BARE_WORD_SLOT}, "no",
                          "[BacFlux] parameters.long_read_qc.length_weight is an unquoted off",
                          "not as a number", "a number, such as 1")

    def test_coverage_profile_on_stops_even_with_the_module_off(self):
        # A fraction: float() would read on as 1.0, which is inside (0, 1] and so
        # would pass the range check. The template leaves mobilome.run false, and
        # the run must stop anyway, as it does for a mobilome block written as one
        # value.
        self.assert_stops({"mobilome.coverage_profile": BARE_WORD_SLOT}, "on",
                          "[BacFlux] mobilome.coverage_profile is an unquoted on",
                          "a fraction greater than 0 and at most 1, such as 0.5")

    def test_coverage_profile_off_names_yaml_not_the_range(self):
        # With the module on, off must get the YAML message, not the range
        # check's "it is 0.0", which would not say how the value was misread.
        status, output = self.dry_run({"mobilome.run": True,
                                       "mobilome.coverage_profile": BARE_WORD_SLOT}, "off")
        self.assertEqual(status, 1, output[-3000:])
        self.assertIn("[BacFlux] mobilome.coverage_profile is an unquoted off", output)
        self.assertNotIn("but it is 0.0", output)

    def test_resources_ram_gb_off_stops_with_the_yaml_message(self):
        self.assert_stops({"resources.ram_gb": BARE_WORD_SLOT}, "off",
                          "[BacFlux] resources.ram_gb is an unquoted off",
                          "a whole number, 1 or more")

    def test_every_number_key_stops_on_a_bare_word(self):
        # One dry run per number key in the template, alternating the two
        # spellings. ConfigNumberKeyListTest below checks that this list and the
        # template agree, so a new number key cannot be left out by accident.
        for index, key in enumerate(NUMBER_KEYS):
            word, read_as = (("off", "an unquoted off"), ("on", "an unquoted on"))[index % 2]
            with self.subTest(key=key, word=word):
                self.assert_stops({key: BARE_WORD_SLOT}, word,
                                  f"[BacFlux] {key} is {read_as}", "not as a number")

    # ── What must still parse ──

    def test_template_config_parses(self):
        status, output = self.dry_run({})
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("Phage caller: virsorter2", output)

    def test_bare_off_still_means_off_for_both_decontamination_modes(self):
        # mode and long_read_mode read a bare off as off. long_read_mode is checked
        # in every pipeline mode, though only hybrid acts on it.
        status, output = self.dry_run(
            {"parameters.decontamination.mode": BARE_WORD_SLOT,
             "parameters.decontamination.long_read_mode": BARE_WORD_SLOT},
            "off",
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("Decontamination mode: off", output)

    def test_blank_long_read_mode_gets_its_default_off_like_a_blank_mode(self):
        # A blank long_read_mode (nothing after the colon, which YAML reads as "no
        # value") or an empty "" names no policy, so it gets the default, off, just
        # as a blank decontamination mode gets its default. Checked in hybrid, the
        # one mode that acts on the key. -p prints the shell command of every job
        # in the dry run; the second screen (rule screen_final_contigs) is the only
        # selector call run with --check-only, and its --mode must be off.
        for written in ("", "''"):
            with self.subTest(long_read_mode=written or "(blank)"):
                status, output = self.dry_run(
                    {"mode": "hybrid",
                     "parameters.decontamination.long_read_mode": BARE_WORD_SLOT},
                    written, extra_args=["-p"],
                )
                self.assertEqual(status, 0, output[-3000:])
                screen_commands = [line for line in output.splitlines()
                                   if "--check-only" in line]
                self.assertTrue(screen_commands, output[-3000:])
                for command in screen_commands:
                    self.assertIn("--mode off", command)

    def test_valid_numbers_reach_filtlong_unchanged(self):
        # The checks only stop booleans and text: a number, written as an
        # integer, a decimal or in quotes, is converted exactly as before. -p prints each
        # job's shell command, so the values Filtlong receives can be read off.
        status, output = self.dry_run(
            {"mode": "hybrid",
             "parameters.long_read_qc.min_length": 500,
             "parameters.long_read_qc.keep_percent": 97.5,
             "parameters.long_read_qc.length_weight": 0},
            extra_args=["-p"],
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("--min_length 500", output)
        self.assertIn("--keep_percent 97.5", output)
        self.assertIn("--length_weight 0.0", output)
        status, output = self.dry_run(
            {"mode": "nanopore", "parameters.long_read_qc.keep_percent": "95"},
            extra_args=["-p"],
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("--keep_percent 95 ", output)

    def test_blank_numbers_take_their_defaults(self):
        # A number key with nothing after the colon (YAML reads it as None) takes
        # its default, as a missing key does. Before, min_length and length_weight
        # crashed with a bare Python error and keep_percent reached Filtlong as
        # the text None.
        for key, expected in (("parameters.long_read_qc.min_length", "--min_length 1000"),
                              ("parameters.long_read_qc.keep_percent", "--keep_percent 95 "),
                              ("parameters.long_read_qc.length_weight", "--length_weight 1.0")):
            with self.subTest(key=key):
                status, output = self.dry_run({"mode": "hybrid", key: BARE_WORD_SLOT}, "",
                                              extra_args=["-p"])
                self.assertEqual(status, 0, output[-3000:])
                self.assertIn(expected, output)

    def test_text_in_a_number_key_stops(self):
        # Text that is not a number stops the run with a [BacFlux] message,
        # instead of a bare Python error (float("ninety")) or, for keep_percent,
        # being handed to Filtlong as it is.
        self.assert_stops({"mobilome.tncentral.min_identity": "ninety"}, None,
                          "[BacFlux] mobilome.tncentral.min_identity is 'ninety', "
                          "which is not a number", "a percent identity, such as 90.0")
        self.assert_stops({"mode": "nanopore", "parameters.long_read_qc.keep_percent": "all"},
                          None,
                          "[BacFlux] parameters.long_read_qc.keep_percent is 'all', "
                          "which is not a number", "There is no off")

    def test_medaka_model_off_still_skips_medaka(self):
        # medaka_model is the one word key where a boolean is a real choice.
        status, output = self.dry_run(
            {"mode": "nanopore", "parameters.nanopore.medaka_model": BARE_WORD_SLOT}, "off"
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("Medaka disabled via config", output)

    # ── The opt-in Baktfold layer ──

    def test_baktfold_off_by_default_adds_no_job(self):
        # The template ships baktfold.run: false, so no Baktfold rule is planned.
        status, output = self.dry_run({}, extra_args=["-p"])
        self.assertEqual(status, 0, output[-3000:])
        self.assertNotIn("structure_annotation", output)

    def test_baktfold_on_with_a_local_database_adds_its_jobs(self):
        # With baktfold.run: true and a database the user holds, the plan gains the
        # symlink view of that database and a Baktfold job per sample, on the CPU;
        # baktfold.fast adds --fast. A stub folder holding the one file the
        # start-up check looks for stands in for the real 22 GB database.
        stub = os.path.join(self.tmp.name, "baktfold_db_stub")
        os.makedirs(stub, exist_ok=True)
        open(os.path.join(stub, "swissprot.dbtype"), "w").close()
        status, output = self.dry_run(
            {"baktfold.run": True, "baktfold.fast": True, "directories.baktfold_db": stub},
            extra_args=["-p"],
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("rule structure_annotation_db_local", output)
        self.assertIn("rule structure_annotation:", output)
        self.assertIn("baktfold run", output)
        self.assertIn("--cpu --fast", output)

    def test_baktfold_database_path_is_checked_only_when_the_layer_is_on(self):
        # A wrong directories.baktfold_db stops a run that uses Baktfold, and is
        # ignored by a run that does not.
        missing = os.path.join(self.tmp.name, "no_such_baktfold_db")
        self.assert_stops({"baktfold.run": True, "directories.baktfold_db": missing}, None,
                          "[BacFlux] directories.baktfold_db points at")
        status, output = self.dry_run({"directories.baktfold_db": missing})
        self.assertEqual(status, 0, output[-3000:])

    # ── The opt-in trait profile (EcoFoldDB) ──

    # Every rule the trait layer adds. None may be planned while traits.run is false.
    TRAIT_RULES = ("trait_db", "trait_db_local", "trait_structure_query",
                   "trait_structure_search", "trait_profile", "traits_summary")

    def make_trait_stubs(self, folder="ecofolddb_stub", leave_out=0):
        """Stub databases for a dry run with traits.run: true; returns (ecofolddb, baktfold).

        A dry run opens no database file, with one exception: when
        directories.ecofolddb_db is set, the start-up check reads its
        EcoFoldDB_descriptions.txt and compares the entry names with the gene
        dictionary that ships with BacFlux. So the EcoFoldDB stub holds the file the
        start-up check looks for (EcoFoldDB.dbtype) and a descriptions table in the
        real column layout, written here from the names that dictionary uses, with
        invented values in every other column. leave_out drops that many names
        (the first ones, sorted), as in a copy that names some entries
        differently. The Baktfold stub holds swissprot.dbtype, as in the Baktfold
        tests above.
        """
        ecofolddb = os.path.join(self.tmp.name, folder)
        os.makedirs(ecofolddb, exist_ok=True)
        open(os.path.join(ecofolddb, "EcoFoldDB.dbtype"), "w").close()
        names = ecofolddb_names_in_gene_table()[leave_out:]
        with open(os.path.join(ecofolddb, "EcoFoldDB_descriptions.txt"), "w", encoding="utf-8") as out:
            out.write("ID\tGene\tProtein\tCategory\tSub-category\tPathway/Activity\tEC\tKO\tUniProt\n")
            for number, name in enumerate(names, start=1):
                gene, pathway = name.split("@", 1)
                out.write(f"STUB{number}\t{gene}\tinvented\tinvented\tinvented\t{pathway}\tNA\tNA\tNA\n")
        baktfold = os.path.join(self.tmp.name, "trait_baktfold_db_stub")
        os.makedirs(baktfold, exist_ok=True)
        open(os.path.join(baktfold, "swissprot.dbtype"), "w").close()
        return ecofolddb, baktfold

    def test_traits_off_by_default_adds_no_job(self):
        # The template ships traits.run: false: no trait rule is planned, the
        # Baktfold database is not asked for on the trait profile's account, and
        # the trait heatmap is not among MultiQC's inputs.
        status, output = self.dry_run({}, extra_args=["-p"])
        self.assertEqual(status, 0, output[-3000:])
        for rule_name in self.TRAIT_RULES:
            self.assertNotIn("rule " + rule_name, output)
        self.assertNotIn("structure_annotation", output)
        self.assertNotIn("traits_summary_mqc", output)
        self.assertNotIn("Trait profile: ON", output)

    def test_traits_on_with_local_databases_adds_its_jobs(self):
        # traits.run on, Baktfold off, both databases held locally: the plan gains
        # the two database views and the four per-sample and summary rules, with
        # EcoFoldDB's own settings, and the heatmap joins MultiQC's inputs.
        # Baktfold itself stays off: its database is linked only for the ProstT5
        # model, no Baktfold annotation runs, and none is passed to the calls.
        ecofolddb, baktfold = self.make_trait_stubs()
        status, output = self.dry_run(
            {"traits.run": True, "directories.ecofolddb_db": ecofolddb,
             "directories.baktfold_db": baktfold},
            extra_args=["-p"],
        )
        self.assertEqual(status, 0, output[-3000:])
        for planned in ("rule structure_annotation_db_local:", "rule trait_db_local:",
                        "rule trait_structure_query:", "rule trait_structure_search:",
                        "rule trait_profile:", "rule traits_summary:"):
            self.assertIn(planned, output)
        self.assertNotIn("rule trait_db:", output)
        self.assertNotIn("rule structure_annotation:", output)
        self.assertNotIn("--baktfold-tsv", output)
        self.assertIn("baktfold proteins-predict", output)
        self.assertIn("-s 8", output)
        self.assertIn("--max-protein-length 4000", output)
        self.assertIn("--max-gene-gap 10", output)
        # MultiQC is handed the heatmap as a path relative to the output root
        # (after a space in its -d list), unlike the absolute path in the plan
        # of rule traits_summary.
        self.assertIn(" 09.report/traits_summary_mqc.tsv", output)

    def test_traits_on_with_baktfold_on_passes_baktfold_names(self):
        # With both layers on, Baktfold runs, and its names reach the trait calls.
        ecofolddb, baktfold = self.make_trait_stubs()
        status, output = self.dry_run(
            {"traits.run": True, "baktfold.run": True,
             "directories.ecofolddb_db": ecofolddb, "directories.baktfold_db": baktfold},
            extra_args=["-p"],
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("rule structure_annotation:", output)
        self.assertIn("--baktfold-tsv", output)

    def test_download_all_with_traits_on_fetches_both_databases(self):
        # snakemake download_all must fetch everything the run downloads. With
        # traits.run on and baktfold.run off, that includes Baktfold's database
        # (for the ProstT5 model) as well as EcoFoldDB, and no analysis job.
        status, output = self.dry_run({"traits.run": True}, targets=["download_all"])
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("rule structure_annotation_db:", output)
        self.assertIn("rule trait_db:", output)
        self.assertNotIn("rule trait_structure_query:", output)

    def test_traits_database_paths_are_checked_only_when_the_layer_is_on(self):
        # A wrong directories.ecofolddb_db stops a run that uses the trait profile,
        # and is ignored by a run that does not.
        missing = os.path.join(self.tmp.name, "no_such_ecofolddb")
        self.assert_stops({"traits.run": True, "directories.ecofolddb_db": missing}, None,
                          "[BacFlux] directories.ecofolddb_db points at")
        status, output = self.dry_run({"directories.ecofolddb_db": missing})
        self.assertEqual(status, 0, output[-3000:])
        # The trait profile needs Baktfold's database too, so with traits.run on a
        # wrong directories.baktfold_db stops the run even with Baktfold off.
        ecofolddb, _baktfold = self.make_trait_stubs()
        self.assert_stops({"traits.run": True, "directories.ecofolddb_db": ecofolddb,
                           "directories.baktfold_db": os.path.join(self.tmp.name, "no_such_baktfold_db")},
                          None, "[BacFlux] directories.baktfold_db points at")

    def test_ecofolddb_copy_with_other_entry_names_stops(self):
        # The copy of EcoFoldDB in its authors' GitHub repository names some entries
        # differently from the Zenodo release; with it, trait_profile.py would stop
        # on the missing names only after the 3Di step and the search, about
        # 1 h 50 min per genome. A copy missing one name the gene dictionary uses
        # must stop the run at start-up instead, naming that entry and pointing to
        # the Zenodo archive.
        ecofolddb, baktfold = self.make_trait_stubs(folder="ecofolddb_stub_renamed", leave_out=1)
        left_out = ecofolddb_names_in_gene_table()[0]
        self.assert_stops({"traits.run": True, "directories.ecofolddb_db": ecofolddb,
                           "directories.baktfold_db": baktfold}, None,
                          "[BacFlux] directories.ecofolddb_db is",
                          "lacks 1 of the EcoFoldDB entries the trait rules use", left_out,
                          "Zenodo")

    def test_traits_on_with_no_ecofolddb_source_stops(self):
        # EcoFoldDB has no downloader of its own: with no link and no local copy
        # there is nothing to search, so the run stops at start-up, naming both keys.
        self.assert_stops({"traits.run": True, "links.ecofolddb_link": ""}, None,
                          "[BacFlux] traits.run is true, but neither links.ecofolddb_link "
                          "nor directories.ecofolddb_db is set")

    def test_ecofolddb_link_that_is_not_an_archive_stops_only_when_on(self):
        # Zenodo's API address for the same file ends in /content, not .tar.gz, and
        # would give the downloaded file the wrong name. It stops a run that uses
        # the link, and is ignored by a run that does not.
        api_link = "https://zenodo.org/api/records/16908168/files/EcoFoldDB_v2.0.tar.gz/content"
        self.assert_stops({"traits.run": True, "links.ecofolddb_link": api_link}, None,
                          "[BacFlux] Invalid links.ecofolddb_link", "got 'content'")
        status, output = self.dry_run({"links.ecofolddb_link": api_link})
        self.assertEqual(status, 0, output[-3000:])

    def test_max_gene_gap_below_one_stops_only_when_on(self):
        # trait_profile.py refuses a gap below 1 (next-door genes are 1 apart, so
        # with 0 no two different proteins could be co-localised), but only when
        # rule trait_profile runs, after the 3Di step of about 1 h 50 min per
        # genome. The start-up check must stop both 0 and a negative gap first,
        # at the script's limit. With the layer off the value is not used, and the
        # run goes ahead.
        for gap in (0, -1):
            self.assert_stops({"traits.run": True, "traits.max_gene_gap": gap}, None,
                              "[BacFlux] traits.max_gene_gap must be a whole number of "
                              "genes, 1 or more", f"but it is {gap}.")
            status, output = self.dry_run({"traits.max_gene_gap": gap})
            self.assertEqual(status, 0, output[-3000:])

    def test_traits_block_written_as_on_stops(self):
        self.assert_stops({"traits": BARE_WORD_SLOT}, "on",
                          "[BacFlux] traits is a block", "run: true or run: false")


class ConfigBlockListTest(unittest.TestCase):
    """The block check covers every block the shipped template has.

    00_common.smk keeps its list of blocks (_CONFIG_BLOCKS) by hand. If a new block
    is added to config/config.yaml and not to that list, a user who writes it as a
    single value gets no [BacFlux] stop, only the old silent default or a bare
    Python error. This test reads both and fails when they differ. It needs no
    Snakemake: it compares the template's nested settings with the list's text.
    """

    def test_every_template_block_is_checked(self):
        # Every key in the template whose value is itself a set of settings (a
        # mapping), named by its dotted path, e.g. "parameters.decontamination".
        with open(TEMPLATE) as handle:
            template = yaml.safe_load(handle)
        template_blocks = set()
        pending = [("", template)]
        while pending:
            prefix, mapping = pending.pop()
            for key, value in mapping.items():
                if isinstance(value, dict):
                    dotted = prefix + key
                    template_blocks.add(dotted)
                    pending.append((dotted + ".", value))

        # The names inside the _CONFIG_BLOCKS = ( ... ) tuple in 00_common.smk.
        common = os.path.join(REPO, "workflow", "rules", "shared", "00_common.smk")
        with open(common) as handle:
            text = handle.read()
        start = text.index("_CONFIG_BLOCKS = (")
        end = text.index(")", start)
        listed_blocks = set(re.findall(r'"([A-Za-z0-9_.]+)"', text[start:end]))

        self.assertEqual(template_blocks, listed_blocks)


class ConfigNumberKeyListTest(unittest.TestCase):
    """Every number key in the template is checked against a bare on/off.

    The rule files check each number key by name, one call per key. If a number key
    is added to config/config.yaml without that call, a bare on/off in it is read as
    1 or 0 without a word. This test compares three lists: the template's number
    keys, NUMBER_KEYS above (what the dry-run test loops over), and the keys named
    in the rule files' checks. It needs no Snakemake.
    """

    def test_every_template_number_key_is_checked(self):
        # Every setting in the template whose value is a number. bool is excluded
        # by name because Python counts True and False as numbers too.
        with open(TEMPLATE) as handle:
            template = yaml.safe_load(handle)
        template_numbers = set()
        pending = [("", template)]
        while pending:
            prefix, mapping = pending.pop()
            for key, value in mapping.items():
                if isinstance(value, dict):
                    pending.append((prefix + key + ".", value))
                elif isinstance(value, (int, float)) and not isinstance(value, bool):
                    template_numbers.add(prefix + key)

        # The keys the rule files check: the first argument of every
        # _config_number(...) call, plus resources.<key> for every
        # _whole_number("<key>", ...) call, which runs the same check.
        rules = os.path.join(REPO, "workflow", "rules", "shared")
        text = ""
        for name in ("00_common.smk", "80_mobilome.smk"):
            with open(os.path.join(rules, name)) as handle:
                text += handle.read()
        checked = set(re.findall(r'_config_number\(\s*"([a-z0-9_.]+)"', text))
        checked |= {"resources." + key for key in re.findall(r'_whole_number\("([a-z_]+)"', text)}

        self.assertEqual(template_numbers, set(NUMBER_KEYS))
        self.assertEqual(template_numbers, checked)


if __name__ == "__main__":
    unittest.main()
