"""Tests for multiqc_citations.py: the References section of the MultiQC report.

The section must list the papers of exactly the tools the run's report shows,
each with a working DOI link, and register no DOI with MultiQC itself.

Run: pytest workflow/scripts/90_report/test_multiqc_citations.py -q
"""

from html.parser import HTMLParser
import re

import pytest

import multiqc_citations as citations


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.targets.append(dict(attrs)["href"])


@pytest.mark.parametrize("short_reads,long_reads,extra", [
    (False, False, set()),
    (True, False, {"fastp", "QualiMap", "QualiMap 2"}),
    (False, True, {"NanoStat / NanoPlot", "QualiMap", "QualiMap 2"}),
    (True, True, {"fastp", "NanoStat / NanoPlot", "QualiMap", "QualiMap 2"}),
])
def test_only_methods_represented_by_the_mode_are_selected(short_reads, long_reads, extra):
    methods, trait_methods = citations.references_for_run(short_reads, long_reads)
    assert {reference.tool for reference in methods} == {
        "MultiQC", "QUAST", "CheckM", "GTDB-Tk", "GTDB", "Bakta", "BLAST+",
    } | extra
    assert trait_methods == []


@pytest.mark.parametrize("short_reads", [False, True])
def test_ont_trimming_adds_the_two_fastp_papers_once(short_reads):
    # fastplong has no paper of its own; its README asks for both fastp papers. In
    # hybrid mode the 2018 one is already there for fastp itself, and is not repeated.
    methods, _ = citations.references_for_run(short_reads=short_reads, long_reads=True,
                                              ont_trimming=True)
    tools = [reference.tool for reference in methods]
    assert tools.count("fastp") == 1
    assert "fastplong (cite fastp)" in tools
    assert citations.FASTPLONG_REFERENCE.doi == "10.1002/imt2.107"


def test_traits_add_the_four_structural_method_references_and_explain_rule_sources():
    content = citations.render_references(traits=True)
    assert "Trait profile methods" in content
    for reference in citations.TRAIT_REFERENCES:
        assert reference.title in content
        assert reference.authors in content
        assert str(reference.year) in content
    assert "<strong>Preprint</strong>" in content
    assert "References and assays" in content
    assert "laboratory assay" in content
    assert "Trait profile methods" not in citations.render_references()


def test_visible_links_are_unique_without_relabelling_native_tool_citations():
    content = citations.render_references(short_reads=True, long_reads=True, traits=True)
    header = content.split("-->", 1)[0]
    assert "id: bacflux_references" in header
    assert "section_name: References" in header
    assert "plot_type: html" in header
    assert "doi:" not in header
    links = Links()
    links.feed(content)
    dois = [url.removeprefix("https://doi.org/") for url in links.targets]
    # 15: the seven shared references (MultiQC, QUAST, CheckM, GTDB-Tk, GTDB, Bakta,
    # BLAST+), fastp, NanoPack, the two QualiMap papers and the four trait methods.
    assert len(dois) == len(set(dois)) == 15
    assert all(re.fullmatch(r"10\.\d{4,9}/\S+", doi) for doi in dois)


def test_bibliographic_text_is_escaped_instead_of_becoming_markup():
    reference = citations.Reference('<script>unsafe</script>', 'A & B "quoted"', 2026,
                                     "Title <tag>", "Journal & more", "10.1234/test", "<Note>")
    content = "\n".join(citations.reference_list([reference]))
    assert "<script>" not in content
    assert "&lt;script&gt;unsafe&lt;/script&gt;" in content
    assert "A &amp; B &quot;quoted&quot;" in content
    assert "Title &lt;tag&gt;" in content
    assert "&lt;Note&gt;" in content


def test_command_writes_deterministic_utf8_html_and_creates_parent(tmp_path):
    output = tmp_path / "09.report" / "bacflux_references_mqc.html"
    args = ["--out", str(output), "--short-reads", "--traits"]
    assert citations.main(args) == 0
    first = output.read_bytes()
    assert "García-Alcalde" in first.decode("utf-8")
    assert citations.main(args) == 0
    assert output.read_bytes() == first
