"""Plan the two workbook jobs through the real workflow, without running any tool.

Rules amr_workbook and traits_workbook (shared/90_report.smk) pass different
optional tables to final_reports.py depending on the input mode, the phage
caller, the mobilome module and the trait profile. For each of the 32
combinations, this test asks Snakemake for a dry run that prints the shell
commands, then checks two things: that the command is one final_reports.py
accepts, with exactly the optional tables that combination has, and that the
rule's inputs are the folders and files the producing rules declare, so
Snakemake links each workbook to the jobs that make its inputs.

It needs Snakemake (the BacFlux snakemake environment) and is skipped without it.
Run: pytest miscellaneous/test_final_report_integration.py -q
"""

import copy
import importlib.util
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

import pytest
import yaml

import test_config_yaml_booleans as config_tests


REPO = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO / "workflow/Snakefile"
REPORT_SCRIPT = REPO / "workflow/scripts/90_report/final_reports.py"
SNAKEMAKE = config_tests.find_snakemake()
pytestmark = pytest.mark.skipif(SNAKEMAKE is None, reason="Snakemake is not installed")


@pytest.fixture(scope="module")
def report_parser():
    """final_reports.py's own command-line parser, imported from the script."""
    sys.path.insert(0, str(REPORT_SCRIPT.parent))
    try:
        spec = importlib.util.spec_from_file_location("report_integration_cli", REPORT_SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module.make_parser()


@pytest.fixture(scope="module")
def input_config(tmp_path_factory):
    """Use real filenames and small database stubs; dry runs never open genomes."""
    root = tmp_path_factory.mktemp("final_report_inputs")
    config = yaml.safe_load((REPO / "config/config.yaml").read_text())
    for technology, names in {
        "illumina": ("SAMPLE1_R1.fastq.gz", "SAMPLE1_R2.fastq.gz"),
        "nanopore": ("SAMPLE1_ont.fastq.gz",),
        "contigs": ("SAMPLE1.fasta",),
    }.items():
        folder = root / technology
        folder.mkdir()
        for name in names:
            (folder / name).touch()
        config["input"][technology + "_dir"] = str(folder)

    databases = root / "databases"
    databases.mkdir()
    for key in ("bakta_db", "eggnog_db", "gtdbtk_db", "platon_db", "blast_db"):
        folder = databases / key
        folder.mkdir()
        config["directories"][key] = str(folder)
    for name in ("nodes.dmp", "names.dmp"):
        (databases / "blast_db" / name).touch()
    config["links"]["phix_link"] = "https://example.org/phix.fna.gz"
    config["links"]["card_link"] = "https://example.org/card.tar.bz2"

    # Startup checks require catalogue names, but no structure database content.
    ecofolddb = databases / "ecofolddb"
    ecofolddb.mkdir()
    (ecofolddb / "EcoFoldDB.dbtype").touch()
    catalogue = ["ID\tGene\tProtein\tCategory\tSub-category\tPathway/Activity\tEC\tKO\tUniProt"]
    for index, name in enumerate(config_tests.ecofolddb_names_in_gene_table(), 1):
        gene, pathway = name.split("@", 1)
        catalogue.append(f"STUB{index}\t{gene}\tstub\tstub\tstub\t{pathway}\tNA\tNA\tNA")
    (ecofolddb / "EcoFoldDB_descriptions.txt").write_text("\n".join(catalogue) + "\n")
    baktfold = databases / "baktfold"
    baktfold.mkdir()
    (baktfold / "swissprot.dbtype").touch()
    config["directories"]["ecofolddb_db"] = str(ecofolddb)
    config["directories"]["baktfold_db"] = str(baktfold)
    return config


def rendered_report_commands(output):
    """Find the final_reports.py commands in the dry run's output, as word lists."""
    pattern = r"(?m)^\s*python\b[^\n]*final_reports\.py\s+(amr|traits)\b[\s\S]*?--out[^\n]*"
    commands = {}
    for match in re.finditer(pattern, output):
        words = shlex.split(match.group(0).replace("\\\n", " "))
        # Shell redirection belongs to the job, not to the Python arguments.
        if ">" in words:
            words = words[:words.index(">")]
        commands[match.group(1)] = words[2:]
    return commands


def report_inputs(output, rule):
    """The input paths the dry run lists for one rule."""
    match = re.search(rf"(?:local)?rule {rule}:\n\s+input: ([^\n]*)", output)
    assert match is not None, f"No planned input line for {rule}:\n{output[-5000:]}"
    return {path.strip() for path in match.group(1).split(",")}


@pytest.mark.parametrize("mode", ["contigs", "illumina", "nanopore", "hybrid"])
@pytest.mark.parametrize("caller", ["virsorter2", "genomad"])
@pytest.mark.parametrize("mobilome", [False, True])
@pytest.mark.parametrize("traits", [False, True])
def test_report_dag_and_cli_matrix(tmp_path, input_config, report_parser,
                                  mode, caller, mobilome, traits):
    config = copy.deepcopy(input_config)
    config["mode"] = mode
    config["phage"]["caller"] = caller
    config["mobilome"]["run"] = mobilome
    config["traits"]["run"] = traits
    output_root = tmp_path / "output with spaces"
    config["directories"]["output_dir"] = str(output_root)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    targets = [str(output_root / "10.workbooks/SAMPLE1_amr.xlsx")]
    if traits:
        targets.append(str(output_root / "10.workbooks/SAMPLE1_traits.xlsx"))
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = "/tmp/bacflux_test_cache"
    result = subprocess.run(
        [SNAKEMAKE, *targets, "--dry-run", "--printshellcmds", "--cores", "1",
         "--snakefile", str(SNAKEFILE), "--configfile", str(config_path)],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-10000:]
    commands = rendered_report_commands(output)
    assert set(commands) == ({"amr", "traits"} if traits else {"amr"}), output[-5000:]

    for command, words in commands.items():
        # The command parses, names the sample and caller, and carries each
        # optional table exactly when its step is part of the run.
        args = report_parser.parse_args(words)
        assert args.sample == "SAMPLE1" and args.phage_caller == caller
        assert args.output_root == str(output_root)
        assert bool(args.plasmid_concordance) == (caller == "genomad")
        assert bool(args.replicons) == mobilome
        if command == "amr":
            assert bool(args.ice_table) == mobilome
            assert bool(args.mobility) == mobilome
            assert bool(args.mobility_audit) == mobilome
        assert "{sample}" not in " ".join(words)

        # The rule's inputs are the folders the producing rules declare as
        # directory() outputs; the command reads files inside them.
        inputs = report_inputs(output, command + "_workbook")
        viral_directory = output_root / "07.phages" / (
            "genomad" if caller == "genomad" else "virsorter"
        ) / "SAMPLE1"
        checkv_directory = output_root / "07.phages/checkv/SAMPLE1"
        assert str(viral_directory) in inputs
        assert str(checkv_directory) in inputs
        assert args.viral_table not in inputs
        assert Path(args.checkv_quality).parent == checkv_directory
        assert Path(args.checkv_contamination).parent == checkv_directory
        assert Path(args.viral_table).is_relative_to(viral_directory)
        assert args.contigs in inputs
        if command == "traits":
            assert args.provenance in inputs and args.gene_dictionary in inputs
            assert str(Path(args.catalogue).parent) in inputs
            assert args.catalogue not in inputs
            # rule trait_profile writes the provenance file this command reads.
            assert "--out-provenance " in output
            assert args.provenance in output
