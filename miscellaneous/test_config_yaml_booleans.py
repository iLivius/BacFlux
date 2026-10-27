#!/usr/bin/env python3
"""Start-up test: an unquoted on/off in the config stops the run, never misreads it.

Snakemake reads the config file as YAML 1.1, where an UNQUOTED on, yes or true is
the boolean True and an unquoted off, no or false is the boolean False, not a word.
So `caller: off` reaches BacFlux as False. Left unchecked, several keys would
quietly turn such a value into their default (`phage.caller: off` would run
VirSorter2; `decontamination: off`, written as one value instead of a block,
would run auto) and others would crash with a bare Python error. The checks in
workflow/rules/shared/00_common.smk (section 0, and _stop_on_yaml_bool where each
key is read) stop the run instead, with a [BacFlux] message and exit status 1.

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

    def dry_run(self, changes, bare_word=None, extra_args=()):
        """Write the base config with `changes` applied, dry-run it, return (status, output).

        changes:    {dotted key: value}. A value of BARE_WORD_SLOT is replaced in the
                    written file by `bare_word`, unquoted, e.g. `caller: off`.
        extra_args: more Snakemake options, e.g. ["-p"] to print the shell
                    command of every job the dry run lists.
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
            [SNAKEMAKE, "-n", "--cores", "1", "--snakefile", SNAKEFILE, "--configfile", path,
             *extra_args],
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

    def test_medaka_model_off_still_skips_medaka(self):
        # medaka_model is the one word key where a boolean is a real choice.
        status, output = self.dry_run(
            {"mode": "nanopore", "parameters.nanopore.medaka_model": BARE_WORD_SLOT}, "off"
        )
        self.assertEqual(status, 0, output[-3000:])
        self.assertIn("Medaka disabled via config", output)


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


if __name__ == "__main__":
    unittest.main()
