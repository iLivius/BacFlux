#!/usr/bin/env python3
"""Write the tool versions for the Software Versions section of the MultiQC report.

What it is for
--------------
The Software Versions section at the end of the MultiQC report says which
version of each tool made the results it shows. MultiQC finds some versions
itself, but several report files carry none (CheckM's table, for one). This
script fills that in, from the best source there is:
  1. the version the tool wrote into its own output for the samples in the
     report (a log line, a JSON field, an HTML banner), when it wrote one;
  2. otherwise the version pinned in the tool's environment file in
     workflow/envs/, labelled "(workflow pin)". A pin says which version this
     copy of BacFlux installs; it does not prove which version made a result
     computed earlier, which is why it is labelled.
It reads files only: it never runs a tool or builds an environment.

Takes in (from rule multiqc in shared/90_report.smk): the folder of the
environment files, the BacFlux output folder, the sample names, and which
kinds of tools the run used (--short-reads, --long-reads, --ont-trimming,
--traits).
Gives back: 09.report/multiqc_versions.yaml, a MultiQC configuration file with a
software_versions entry, which rule multiqc hands to MultiQC with a second
--config.

Runs in envs/multiqc.yaml, MultiQC's own environment, which is why MultiQC's
own version can be read directly. Tests: test_multiqc_versions.py.
"""

import argparse
import importlib.metadata
import json
import re
from pathlib import Path


# The tools of every run's report, each as (name shown in the report,
# environment file without .yaml, package name inside that environment).
# Read-mode tools and the trait profile's structure tools are added in
# collect_versions() when the run used them.
SHARED_TOOLS = (
    ("QUAST", "quast", "quast"),
    ("CheckM", "checkm", "checkm-genome"),
    ("GTDB-Tk", "gtdbtk", "gtdbtk"),
    ("Bakta", "bakta", "bakta"),
    # The adapter screen of the finished genome (Adapters table, every mode). Its
    # BLAST table records no version, so the pin is shown.
    ("BLAST+", "blast", "blast"),
)

# How each tool writes its version into its own output, as a regular expression
# whose group is the version number. fastp and fastplong are read from their JSON
# reports instead.
VERSION_PATTERNS = {
    "QUAST": r"^Version:\s*([^\s]+)",
    "GTDB-Tk": r"\bGTDB-Tk v([^\s]+)",
    "Bakta": r"\bUTILS - version=([^\s]+)",
    "NanoPlot": r"\bNanoPlot(?:\s+version)?[:\s]+v?([0-9][\w.]+)",
    "QualiMap": r"\bQualimap(?:\s+version)?\s+(?:v\.?)?([0-9][\w.]+)",
}


def configured_version(env_dir, environment, package):
    """Return the version of a package pinned for one environment.

    Takes in: the environment folder (workflow/envs), the environment's name
              (for example "checkm" for checkm.yaml) and the package's name.
    Reads the exact version from the environment file ("- checkm-genome=1.2.5")
    or, when the file does not pin it, from the package file name in the pin
    list beside it (checkm.linux-64.pin.txt), which lists every package build
    the environment is made from.
    Gives back: the version as text, so "1.10" stays "1.10". Stops with an error
    when neither file names the package.
    """
    recipe = Path(env_dir) / f"{environment}.yaml"
    pattern = re.compile(r"^\s*-\s*" + re.escape(package) + r"=([^=\s#]+)(?:=|\s|$)")
    for line in recipe.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return match.group(1)
    pin_list = recipe.with_suffix(".linux-64.pin.txt")
    if pin_list.is_file():
        for line in pin_list.read_text(encoding="utf-8").splitlines():
            # A pin-list line ends in a file name such as
            # foldseek-10.941cd33-h5021889_0.conda: name, version, build.
            filename = line.rsplit("/", 1)[-1].split("#", 1)[0]
            for suffix in (".conda", ".tar.bz2"):
                if filename.startswith(package + "-") and filename.endswith(suffix):
                    return filename[len(package) + 1:-len(suffix)].rsplit("-", 1)[0]
    raise ValueError(f"No exact version for {package} in {recipe} or {pin_list}")


def version_sources(root, sample, tool):
    """The files of one sample in which a tool writes its version.

    Takes in: the BacFlux output folder, a sample name and a tool name.
    Gives back: the paths to look in, all inside that sample's own folders, so a
    sample left over from an earlier run in the same output folder is never
    read. An empty list for a tool that writes no version anywhere (CheckM).
    """
    assembly = root / "02.assembly" / sample / "eval"
    if tool == "QUAST":
        return list((assembly / "quast").rglob("quast.log"))
    if tool == "GTDB-Tk":
        return [root / "03.taxonomy" / sample / "gtdbtk.log"]
    if tool == "Bakta":
        return [root / "04.annotation" / "bakta" / sample / f"{sample}.log"]
    if tool == "fastp":
        return list((root / "01.reads" / sample / "illumina").glob("*.json"))
    if tool == "fastplong":
        return list((root / "01.reads" / sample / "ont" / "adapters").glob("*_fastplong.json"))
    if tool == "NanoPlot":
        return list((root / "01.reads" / sample / "ont").rglob("*NanoPlot*.log"))
    if tool == "QualiMap":
        return [assembly / "qualimap" / "qualimapReport.html"]
    return []


def observed_versions(root, samples, tool):
    """Return every version a tool wrote for the samples in the report.

    Takes in: the output folder, the sample names and a tool name.
    Gives back: the distinct versions, sorted. Several are possible when samples
    were analysed at different times, and all are shown. An empty list when the
    tool wrote none.
    """
    versions = set()
    for sample in samples:
        for path in version_sources(root, sample, tool):
            if not path.is_file():
                continue
            if tool in {"fastp", "fastplong"}:
                data = json.loads(path.read_text(encoding="utf-8"))
                version = data.get("summary", {}).get(f"{tool}_version")
                if version:
                    versions.add(str(version))
                continue
            # A tool writes its version near the start of its log, so the first
            # 128 kB is enough and a long log is never read in full.
            with path.open(encoding="utf-8", errors="replace") as handle:
                text = handle.read(128 * 1024)
            matches = re.findall(VERSION_PATTERNS[tool], text, flags=re.I | re.M)
            versions.update(matches)
    return sorted(versions)


def collect_versions(env_dir, root, samples, short_reads=False, long_reads=False, traits=False,
                     ont_trimming=False):
    """Decide the version(s) to show for every tool of this run's report.

    Takes in: the environment folder, the output folder, the samples, and which
              kinds of tools the run used (ont_trimming: fastplong trimmed the ONT
              reads, which the Adapters table reports).
    Gives back: {tool: [version, ...]}, each from the tool's own output when it
              recorded one, otherwise the pin labelled "(workflow pin)".
    """
    tools = list(SHARED_TOOLS)
    if short_reads:
        tools.append(("fastp", "fastp", "fastp"))
    if long_reads:
        tools.append(("NanoPlot", "nanoplot", "nanoplot"))
    if short_reads or long_reads:
        tools.append(("QualiMap", "qualimap", "qualimap"))
    if ont_trimming:
        tools.append(("fastplong", "fastplong", "fastplong"))
    if traits:
        # The trait profile's structure search: Baktfold's environment brings
        # Foldseek, and Baktfold runs ProstT5 to give each protein its 3Di
        # letters (see docs/analysis/traits.md).
        tools.extend((("Baktfold", "baktfold", "baktfold"), ("Foldseek", "baktfold", "foldseek")))

    versions = {}
    for name, environment, package in tools:
        observed = observed_versions(Path(root), samples, name)
        if observed:
            versions[name] = observed
        else:
            versions[name] = [configured_version(env_dir, environment, package) + " (workflow pin)"]
    # This script runs in MultiQC's own environment, so the installed MultiQC is
    # the one that writes the report.
    try:
        versions["MultiQC"] = [importlib.metadata.version("multiqc")]
    except importlib.metadata.PackageNotFoundError:
        versions["MultiQC"] = [configured_version(env_dir, "multiqc", "multiqc") + " (workflow pin)"]
    return versions


def versions_yaml(versions):
    """Write the versions as MultiQC's software_versions configuration.

    MultiQC groups the Software Versions table by report section, so each tool
    goes under the section that shows its results: NanoPlot under NanoStat (the
    MultiQC module that reads NanoPlot's statistics), BLAST+ and fastplong under
    Adapters, the structure tools under Trait profile, and QualiMap's version under
    the name BamQC, the key MultiQC's
    own QualiMap module uses, so the table gets one column for it, not two.
    Every value is written as a quoted string (JSON strings are valid YAML), so
    a version such as 1.10 is not read back as the number 1.1.
    """
    groups = {}
    for tool, values in versions.items():
        if tool == "NanoPlot":
            group = "NanoStat"
        elif tool in {"BLAST+", "fastplong"}:
            group = "Adapters"
        elif tool in {"Baktfold", "Foldseek"}:
            group = "Trait profile"
        else:
            group = tool
        software = "BamQC" if tool == "QualiMap" else tool
        groups.setdefault(group, {})[software] = values

    lines = ["software_versions:"]
    for group, tools in groups.items():
        lines.append("  " + json.dumps(group) + ":")
        for tool, values in tools.items():
            lines.append("    " + json.dumps(tool) + ": " + json.dumps(values))
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Read the command line from rule multiqc, collect the versions, write the YAML."""
    parser = argparse.ArgumentParser(description="Write the tool versions for MultiQC.")
    parser.add_argument("--env-dir", required=True, help="workflow/envs")
    parser.add_argument("--output-root", required=True, help="The BacFlux output folder.")
    parser.add_argument("--samples", nargs="+", required=True, help="The samples in the report.")
    parser.add_argument("--short-reads", action="store_true", help="The run used Illumina reads.")
    parser.add_argument("--long-reads", action="store_true", help="The run used ONT reads.")
    parser.add_argument("--ont-trimming", action="store_true",
                        help="fastplong trimmed the ONT adapters and barcodes.")
    parser.add_argument("--traits", action="store_true", help="The trait profile ran.")
    parser.add_argument("--out", required=True, help="The YAML file to write.")
    args = parser.parse_args(argv)
    versions = collect_versions(args.env_dir, args.output_root, args.samples,
                                args.short_reads, args.long_reads, args.traits,
                                args.ont_trimming)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(versions_yaml(versions), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
