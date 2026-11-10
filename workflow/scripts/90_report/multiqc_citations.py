#!/usr/bin/env python3
"""Write the References section of the MultiQC report.

What it is for
--------------
The MultiQC report shows results of several tools; a reader who uses them
should cite each tool's paper. This script writes that list as a small HTML
page that MultiQC shows as its own section, titled References. The references
are written into the report itself, so they can be read offline, and only the
tools this run's report shows are listed: the read-QC tools in the read modes,
fastplong when it trimmed the ONT reads, the structure-search methods with the
trait profile.

Takes in (from rule multiqc in shared/90_report.smk): which kinds of tools the
run used (--short-reads, --long-reads, --ont-trimming, --traits).
Gives back: 09.report/bacflux_references_mqc.html. MultiQC reads it as custom
content (a page another program wrote for MultiQC to show) because its name
ends in _mqc.html, and the comment block at its top gives the section's ID and
title.

MultiQC's own modules (QUAST, CheckM, ...) and the trait section also register
each paper's DOI with MultiQC's citation button; this page does not register
them again, because a second registration replaces the tool name MultiQC shows
beside the DOI. Runs in envs/multiqc.yaml. Tests: test_multiqc_citations.py.
"""

import argparse
from dataclasses import dataclass
from html import escape
from pathlib import Path


@dataclass(frozen=True)
class Reference:
    """One paper: the tool it describes, its authors, year, title, journal and
    DOI, and an optional note shown in bold (such as "Preprint")."""
    tool: str
    authors: str
    year: int
    title: str
    journal: str
    doi: str
    note: str = ""


# The papers of the tools every report shows. The details match
# docs/about/citation.md, where the full author lists are; the 2012 QualiMap
# paper below is the one MultiQC's own QualiMap module cites. GTDB is the
# database GTDB-Tk classifies against, and is cited beside the tool.
SHARED_REFERENCES = (
    Reference("MultiQC", "Ewels et al.", 2016,
              "MultiQC: summarize analysis results for multiple tools and samples in a single report",
              "Bioinformatics", "10.1093/bioinformatics/btw354"),
    Reference("QUAST", "Gurevich et al.", 2013,
              "QUAST: quality assessment tool for genome assemblies",
              "Bioinformatics", "10.1093/bioinformatics/btt086"),
    Reference("CheckM", "Parks et al.", 2015,
              "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes",
              "Genome Research", "10.1101/gr.186072.114"),
    Reference("GTDB-Tk", "Chaumeil et al.", 2022,
              "GTDB-Tk v2: memory friendly classification with the genome taxonomy database",
              "Bioinformatics", "10.1093/bioinformatics/btac672"),
    Reference("GTDB", "Parks et al.", 2022,
              "GTDB: an ongoing census of bacterial and archaeal diversity through a "
              "phylogenetically consistent, rank normalized and complete genome-based taxonomy",
              "Nucleic Acids Research", "10.1093/nar/gkab776"),
    Reference("Bakta", "Schwengers et al.", 2021,
              "Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification",
              "Microbial Genomics", "10.1099/mgen.0.000685"),
    # The adapter screen of the finished genome (the Adapters table).
    Reference("BLAST+", "Camacho et al.", 2009,
              "BLAST+: architecture and applications",
              "BMC Bioinformatics", "10.1186/1471-2105-10-421"),
)

# The read-QC tools, listed only in the modes that have reads.
FASTP_REFERENCE = Reference(
    "fastp", "Chen et al.", 2018,
    "fastp: an ultra-fast all-in-one FASTQ preprocessor",
    "Bioinformatics", "10.1093/bioinformatics/bty560",
)
# fastplong trims the ONT adapters and barcodes when parameters.<mode>.ont_adapter is
# not off. It has no paper of its own; its README asks for the two fastp papers, the
# 2018 one above and this one.
FASTPLONG_REFERENCE = Reference(
    "fastplong (cite fastp)", "Chen", 2023,
    "Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp",
    "iMeta", "10.1002/imt2.107",
)
NANOPACK_REFERENCE = Reference(
    "NanoStat / NanoPlot", "De Coster et al.", 2018,
    "NanoPack: visualizing and processing long-read sequencing data",
    "Bioinformatics", "10.1093/bioinformatics/bty149",
)
QUALIMAP_REFERENCES = (
    Reference("QualiMap", "García-Alcalde et al.", 2012,
              "Qualimap: evaluating next-generation sequencing alignment data",
              "Bioinformatics", "10.1093/bioinformatics/bts503"),
    Reference("QualiMap 2", "Okonechnikov et al.", 2016,
              "Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data",
              "Bioinformatics", "10.1093/bioinformatics/btv566"),
)
# The methods of the trait profile (traits.run): EcoFoldDB's structure catalogue,
# and Baktfold, ProstT5 and Foldseek, which predict and search the structures.
TRAIT_REFERENCES = (
    Reference("EcoFoldDB", "Ghaly et al.", 2025,
              "EcoFoldDB: protein structure-guided functional profiling of ecologically relevant microbial traits at the metagenome scale",
              "Environmental Microbiology", "10.1111/1462-2920.70178"),
    Reference("Baktfold", "Bouras et al.", 2026,
              "Baktfold: sensitive protein functional annotation across the microbial tree of life using structural information",
              "bioRxiv", "10.64898/2026.03.31.715528", "Preprint"),
    Reference("ProstT5", "Heinzinger et al.", 2024,
              "Bilingual language model for protein sequence and structure",
              "NAR Genomics and Bioinformatics", "10.1093/nargab/lqae150"),
    Reference("Foldseek", "van Kempen et al.", 2024,
              "Fast and accurate protein structure search with Foldseek",
              "Nature Biotechnology", "10.1038/s41587-023-01773-0"),
)


def references_for_run(short_reads=False, long_reads=False, traits=False, ont_trimming=False):
    """Choose the papers of the tools this run's report shows.

    Gives back two lists: the report's tools, and the trait profile's methods
    (empty without the trait profile), which get their own heading.
    """
    methods = list(SHARED_REFERENCES)
    if short_reads or ont_trimming:
        methods.append(FASTP_REFERENCE)
    if ont_trimming:
        methods.append(FASTPLONG_REFERENCE)
    if long_reads:
        methods.append(NANOPACK_REFERENCE)
    if short_reads or long_reads:
        methods.extend(QUALIMAP_REFERENCES)
    trait_methods = list(TRAIT_REFERENCES) if traits else []
    return methods, trait_methods


def reference_list(references):
    """Write the papers as the lines of an HTML list, one paper per item.

    Every text is escaped (html.escape), so a character such as "<" in a title
    is shown as itself instead of being read as HTML.
    """
    lines = ["<ul>"]
    for reference in references:
        doi = escape(reference.doi, quote=True)
        note = f" <strong>{escape(reference.note)}</strong>." if reference.note else ""
        lines.append(
            f"<li><strong>{escape(reference.tool)}</strong>: "
            f"{escape(reference.authors)} ({reference.year}). "
            f"{escape(reference.title)}. <em>{escape(reference.journal)}</em>.{note} "
            f'<a href="https://doi.org/{doi}" target="_blank" rel="noopener noreferrer">'
            f"{doi}</a>.</li>"
        )
    lines.append("</ul>")
    return lines


def render_references(short_reads=False, long_reads=False, traits=False, ont_trimming=False):
    """Write the whole page: MultiQC's header comment, then the reference lists.

    The header comment (between <!-- and -->) tells MultiQC the section's ID,
    title and kind. It carries no DOI, for the reason given in the module
    docstring.
    """
    methods, trait_methods = references_for_run(short_reads, long_reads, traits, ont_trimming)
    lines = [
        "<!--",
        "id: bacflux_references",
        "section_name: References",
        "plot_type: html",
        "-->",
        "<p>References for the methods represented in this MultiQC report. "
        "The complete BacFlux analysis can require additional tool and database "
        "citations; see the workflow documentation.</p>",
        "<h3>Report methods</h3>",
        *reference_list(methods),
    ]
    if trait_methods:
        lines.extend([
            "<h3>Trait profile methods</h3>",
            *reference_list(trait_methods),
            "<p>The papers behind each trait's rule are listed under "
            "<strong>References and assays</strong> in the Trait profile section, "
            "and in each strain's traits workbook (10.workbooks). They support "
            "the genomic prediction; a laboratory assay is needed to test the "
            "trait itself.</p>",
        ])
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Read the command line from rule multiqc and write the page."""
    parser = argparse.ArgumentParser(description="Write the References section for MultiQC.")
    parser.add_argument("--out", required=True,
                        help="Output custom-content HTML, normally ending in _mqc.html.")
    parser.add_argument("--short-reads", action="store_true")
    parser.add_argument("--long-reads", action="store_true")
    parser.add_argument("--ont-trimming", action="store_true")
    parser.add_argument("--traits", action="store_true")
    args = parser.parse_args(argv)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_references(args.short_reads, args.long_reads, args.traits,
                                      args.ont_trimming),
                    encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
