"""Tests for final_reports.py: the per-strain AMR and traits workbooks.

The workbooks answer three context questions for every gene (on a plasmid?
inside an ICE or IME? inside a phage region?) with yes, no, uncertain or
not_assessed. Each answer has its own rule, and each rule is tested here state
by state, with the evidence source that drives it: Platon and geNomad for the
plasmid, the mobilome module's ICE table for ICE/IME, VirSorter2 or geNomad
with CheckV's viral and host parts for the phage. The grouping of database hits
into loci, the copying of values (the minus strand, GAPS as text, AMRFinderPlus
STRESS rows) and the paths shown in the workbook are tested too.

Everything is invented and written into pytest's tmp_path with the columns the
tools write; the sample is SAMPLE1. Writing the .xlsx file itself needs
XlsxWriter (envs/report.yaml); that one test is skipped where it is missing.

Run: pytest workflow/scripts/90_report/test_final_reports.py -q
"""

import hashlib
import json
import os
import zipfile

import pytest

import final_reports as fr


# ── Small builders for the input tables ──────────────────────────────────────

def write_tsv(path, header, rows):
    """Write a tab-separated table; rows are dicts, missing cells become NA."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(str(row.get(column, "NA")) for column in header))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


ABRICATE_HEADER = ["#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
                   "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE",
                   "ACCESSION", "PRODUCT", "RESISTANCE"]


def abricate_row(contig, start, end, gene, database, strand="+", gaps="0/0"):
    return {"#FILE": "contigs_final.fasta", "SEQUENCE": contig, "START": start, "END": end,
            "STRAND": strand, "GENE": gene, "COVERAGE": f"1-{end - start + 1}/{end - start + 1}",
            "COVERAGE_MAP": "===============", "GAPS": gaps, "%COVERAGE": "100.00",
            "%IDENTITY": "99.50", "DATABASE": database, "ACCESSION": "AB000001.1:1-900",
            "PRODUCT": f"{gene} product", "RESISTANCE": "Aminoglycoside"}


def hit(contig, start, end, database="CARD", gene="geneX", strand="+", number=1):
    """A hit as read_abricate returns it, for the grouping tests."""
    return {"hit_id": f"{database}:{number:04d}", "database": database, "contig": contig,
            "start": start, "end": end, "strand": strand, "gene": gene,
            "pct_identity": 99.0, "pct_coverage": 100.0, "reference_accession": "AB1",
            "product": "p", "resistance": "r"}


def element(element_type, start, end, contig="contig_1", mge_id=None, name="NA",
            mobility="predicted self-transmissible"):
    """An element as read_ice_elements returns it."""
    return {"mge_id": mge_id or f"{contig}|{element_type}-{start}:{end}", "mge_name": name,
            "element_type": element_type, "mobility": mobility, "contig": contig,
            "start": start, "end": end}


def region(region_id, contig, start, end, viral_parts, host_parts):
    """A virus region with CheckV's parts already in contig coordinates."""
    return {"region_id": region_id, "caller": "VirSorter2", "contig": contig,
            "start": start, "end": end, "viral_parts": viral_parts, "host_parts": host_parts}


# ── Grouping the database hits of one gene into one locus ────────────────────

def test_hits_of_one_gene_become_one_locus_anchored_on_the_longest():
    """Three databases report one gene with slightly different ends: one locus,
    with the coordinates of the longest hit and the databases in a fixed order."""
    hits = [hit("contig_1", 1010, 1800, "ResFinder", "aadA2", number=1),
            hit("contig_1", 1000, 1810, "CARD", "ANT(3'')-Ia", number=1),
            hit("contig_1", 1005, 1805, "ARG-ANNOT", "(AGly)aadA2", number=1)]
    loci = fr.group_hits_into_loci(hits, "SAMPLE1")
    assert len(loci) == 1
    assert (loci[0]["start"], loci[0]["end"]) == (1000, 1810)
    assert loci[0]["locus_id"] == "SAMPLE1|contig_1:1000-1810"
    assert loci[0]["databases"] == "CARD; ResFinder; ARG-ANNOT"
    assert loci[0]["hits"][0]["database"] == "CARD"      # the anchor comes first
    assert all(h["locus_id"] == loci[0]["locus_id"] for h in hits)


def test_neighbouring_genes_stay_apart():
    hits = [hit("contig_1", 1000, 2000, number=1), hit("contig_1", 1990, 3000, number=2)]
    assert len(fr.group_hits_into_loci(hits, "SAMPLE1")) == 2


def test_a_short_hit_cannot_chain_two_genes_together():
    """A short hit lying in the overlap of two genes joins one of them only; the
    two genes stay two loci. Grouping every overlapping pair would chain them."""
    hits = [hit("contig_1", 1000, 2000, number=1), hit("contig_1", 1800, 2800, number=2),
            hit("contig_1", 1850, 1950, "ResFinder", number=1)]
    loci = fr.group_hits_into_loci(hits, "SAMPLE1")
    assert len(loci) == 2
    assert sorted(len(locus["hits"]) for locus in loci) == [1, 2]


def test_hits_on_different_contigs_are_never_grouped():
    hits = [hit("contig_1", 1000, 2000, number=1), hit("contig_2", 1000, 2000, number=2)]
    assert len(fr.group_hits_into_loci(hits, "SAMPLE1")) == 2


def test_abricate_values_are_copied_as_written(tmp_path):
    """The minus strand stays "-", and GAPS stays text, so 6/4 is not a date."""
    path = write_tsv(tmp_path / "card.tsv", ABRICATE_HEADER,
                     [abricate_row("contig_1", 100, 900, "geneX", "card", strand="-", gaps="6/4")])
    hits = fr.read_abricate(path, "CARD", "card", str(tmp_path))
    assert hits[0]["strand"] == "-"
    assert hits[0]["gaps"] == "6/4"
    assert hits[0]["source_file"] == "card.tsv"


def test_line_numbers_count_every_line_of_the_file(tmp_path):
    """Blank lines and leading "#" lines are skipped, but a row keeps the number
    of its line in the file, so an error message points at the right line."""
    path = tmp_path / "rules.tsv"
    path.write_text("# a comment\ntrait_id\ttrait\n\nt1\tTrait one\nt2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 5 has 1 fields"):
        fr.read_tsv(str(path), skip_comment_lines=True)
    path.write_text("# a comment\ntrait_id\ttrait\n\nt1\tTrait one\n", encoding="utf-8")
    _header, rows = fr.read_tsv(str(path), skip_comment_lines=True)
    assert rows[0]["__line__"] == 4


def test_a_table_given_as_the_wrong_database_stops_the_run(tmp_path):
    path = write_tsv(tmp_path / "card.tsv", ABRICATE_HEADER,
                     [abricate_row("contig_1", 100, 900, "geneX", "resfinder")])
    with pytest.raises(ValueError, match="given as the CARD table"):
        fr.read_abricate(path, "CARD", "card", str(tmp_path))


# ── Plasmid: the answer for each replicon call ───────────────────────────────

def replicon(replicon_word, source="platon", mobility="NA", evidence="evidence text"):
    return {"contig_1": {"contig": "contig_1", "replicon": replicon_word,
                         "plasmid_mobility": mobility, "mobility_evidence": evidence,
                         "replicon_call_source": source}}


@pytest.mark.parametrize("replicon_word,source,crashed,expected", [
    ("plasmid", "both", False, "yes"),
    ("plasmid", "platon", False, "yes"),
    ("plasmid", "genomad", False, "yes"),        # geNomad alone; the source says so
    ("chromosome", "platon", False, "no"),
    ("chromosome", "conflict", False, "uncertain"),
    ("unknown", "platon", False, "uncertain"),
    ("unknown", "platon", True, "not_assessed"),  # Platon failed
    ("plasmid", "genomad", True, "yes"),          # geNomad's call stands alone
])
def test_plasmid_answer_for_every_replicon_call(replicon_word, source, crashed, expected):
    answer = fr.plasmid_answer("contig_1", replicon(replicon_word, source), crashed)
    assert answer["plasmid"] == expected
    assert answer["plasmid_call_source"] == source


def test_a_closed_genome_without_platon_files_is_still_assessed(tmp_path):
    """Platon writes nothing when every contig is 500 kb or longer. The contig
    lengths decide: 2 Mb or more is a chromosome (no), between is uncertain."""
    platon_dir = tmp_path / "platon"
    platon_dir.mkdir()
    lengths = {"chromosome": 5000000, "megaplasmid": 1123075}
    calls = fr.read_replicon_calls(str(platon_dir), "contigs_final", lengths, None, None)
    assert fr.plasmid_answer("chromosome", calls, False)["plasmid"] == "no"
    assert fr.plasmid_answer("megaplasmid", calls, False)["plasmid"] == "uncertain"


def test_a_megaplasmid_in_platons_chromosome_fasta_is_not_called_chromosome(tmp_path):
    """With a small plasmid in the genome, Platon writes the skipped 1.1 Mb contig
    into its chromosome FASTA; that is not a call, so the answer is uncertain."""
    platon_dir = tmp_path / "platon"
    platon_dir.mkdir()
    (platon_dir / "contigs_final.chromosome.fasta").write_text(">chromosome\nA\n>mega\nA\n")
    (platon_dir / "contigs_final.plasmid.fasta").write_text(">small\nA\n")
    write_tsv(platon_dir / "contigs_final.tsv",
              ["ID", "# Conjugation", "# Mobilization", "# OriT"],
              [{"ID": "small", "# Conjugation": 0, "# Mobilization": 1, "# OriT": 0}])
    lengths = {"chromosome": 5000000, "mega": 1123075, "small": 60000}
    calls = fr.read_replicon_calls(str(platon_dir), "contigs_final", lengths, None, None)
    assert fr.plasmid_answer("mega", calls, False)["plasmid"] == "uncertain"
    assert fr.plasmid_answer("small", calls, False)["plasmid"] == "yes"
    assert fr.plasmid_answer("small", calls, False)["plasmid_mobility"] == "mobilisable"


def test_the_mobilome_replicon_table_is_read_when_given(tmp_path):
    table = write_tsv(tmp_path / "replicons.tsv",
                      ["contig", "replicon", "replicon_id", "plasmid_mobility",
                       "mobility_evidence", "replicon_call_source"],
                      [{"contig": "contig_1", "replicon": "plasmid", "replicon_id": "contig_1",
                        "plasmid_mobility": "conjugative", "mobility_evidence": "from table",
                        "replicon_call_source": "both"}])
    calls = fr.read_replicon_calls(str(tmp_path), "contigs_final", {"contig_1": 1000}, None, table)
    answer = fr.plasmid_answer("contig_1", calls, False)
    assert (answer["plasmid"], answer["plasmid_evidence"]) == ("yes", "from table")


def test_platon_crash_is_read_from_verified_plasmids(tmp_path):
    (tmp_path / "verified_plasmids.txt").write_text(
        "SAMPLE1: Platon exited with status 1; see the log.\n")
    assert fr.platon_crashed(str(tmp_path)) is True
    (tmp_path / "verified_plasmids.txt").write_text("Platon found no plasmid in sample SAMPLE1.\n")
    assert fr.platon_crashed(str(tmp_path)) is False


# ── Phage: the caller's regions and CheckV's viral and host parts ────────────

VS2_HEADER = ["seqname", "trim_bp_start", "trim_bp_end", "seqname_new",
              "final_max_score", "final_max_score_group"]
CONTAMINATION_HEADER = ["contig_id", "contig_length", "total_genes", "viral_genes",
                        "host_genes", "provirus", "proviral_length", "host_length",
                        "region_types", "region_lengths", "region_coords_bp"]
QUALITY_HEADER = ["contig_id", "contig_length", "provirus", "proviral_length", "gene_count",
                  "viral_genes", "host_genes", "checkv_quality", "miuvig_quality",
                  "completeness", "completeness_method", "contamination", "kmer_freq",
                  "warnings"]


def vs2_regions_with_checkv(tmp_path, contamination_rows=True):
    """A VirSorter2 region of 163,418 bp starting at base 2,102,153 of contig_3,
    which CheckV splits into host (1-3603), viral (3604-48893) and host parts,
    as in a real hybrid Klebsiella genome."""
    viral_table = write_tsv(tmp_path / "final-viral-boundary.tsv", VS2_HEADER, [
        {"seqname": "contig_3", "trim_bp_start": 2102153, "trim_bp_end": 2265570,
         "seqname_new": "contig_3||1_partial", "final_max_score": 0.907,
         "final_max_score_group": "dsDNAphage"}])
    quality = write_tsv(tmp_path / "quality_summary.tsv", QUALITY_HEADER, [
        {"contig_id": "contig_3||1_partial", "provirus": "Yes", "checkv_quality": "High-quality",
         "completeness": "98.1", "contamination": "72.3", "warnings": ""}])
    contamination = tmp_path / "contamination.tsv"
    if contamination_rows:
        write_tsv(contamination, CONTAMINATION_HEADER, [
            {"contig_id": "contig_3||1_partial", "contig_length": 163418, "provirus": "Yes",
             "region_types": "host,viral,host",
             "region_coords_bp": "1-3603,3604-48893,48894-163418"}])
    regions = fr.read_viral_regions(viral_table, "virsorter2", {"contig_3": 5579416})
    fr.attach_checkv(regions, quality, str(contamination))
    return regions


def test_checkv_parts_are_moved_to_contig_coordinates(tmp_path):
    regions = vs2_regions_with_checkv(tmp_path)
    assert regions[0]["viral_parts"] == [(2105756, 2151045)]
    assert regions[0]["host_parts"] == [(2102153, 2105755), (2151046, 2265570)]


def test_a_gene_in_the_viral_part_is_in_a_phage(tmp_path):
    regions = vs2_regions_with_checkv(tmp_path)
    answer, names, reason = fr.phage_answer("contig_3", [(2120000, 2121000)], regions, "VirSorter2")
    assert (answer, names) == ("yes", "contig_3||1_partial")
    assert "2105756-2151045" in reason


def test_a_gene_in_the_host_part_of_a_region_is_uncertain(tmp_path):
    """VirSorter2 put the gene in its region, but CheckV calls that part host DNA:
    the two tools do not agree, so the answer is uncertain, not yes."""
    regions = vs2_regions_with_checkv(tmp_path)
    answer, _names, reason = fr.phage_answer("contig_3", [(2200000, 2201000)], regions,
                                             "VirSorter2")
    assert answer == "uncertain"
    assert "host DNA" in reason


def test_a_region_checkv_did_not_assess_is_uncertain(tmp_path):
    regions = vs2_regions_with_checkv(tmp_path, contamination_rows=False)
    answer, _names, reason = fr.phage_answer("contig_3", [(2120000, 2121000)], regions,
                                             "VirSorter2")
    assert answer == "uncertain"
    assert "did not assess" in reason


def test_a_gene_outside_every_region_is_not_in_a_phage(tmp_path):
    regions = vs2_regions_with_checkv(tmp_path)
    assert fr.phage_answer("contig_3", [(10, 900)], regions, "VirSorter2")[0] == "no"
    assert fr.phage_answer("contig_9", [(2120000, 2121000)], regions, "VirSorter2")[0] == "no"


def test_a_virsorter2_full_region_is_the_whole_contig(tmp_path):
    """VirSorter2 runs with --keep-original-seq, so a contig||full region is the
    whole contig, and CheckV's coordinates for it start at the contig's first
    base, whatever trim_bp_start and trim_bp_end say."""
    viral_table = write_tsv(tmp_path / "final-viral-boundary.tsv", VS2_HEADER, [
        {"seqname": "NODE_15", "trim_bp_start": 1, "trim_bp_end": 167784,
         "seqname_new": "NODE_15||full", "final_max_score": 0.95,
         "final_max_score_group": "dsDNAphage"}])
    quality = write_tsv(tmp_path / "quality_summary.tsv", QUALITY_HEADER, [
        {"contig_id": "NODE_15||full", "provirus": "Yes"}])
    contamination = write_tsv(tmp_path / "contamination.tsv", CONTAMINATION_HEADER, [
        {"contig_id": "NODE_15||full", "contig_length": 172524, "provirus": "Yes",
         "region_types": "host,viral,host",
         "region_coords_bp": "1-100180,100181-162355,162356-172524"}])
    regions = fr.read_viral_regions(viral_table, "virsorter2", {"NODE_15": 172524})
    fr.attach_checkv(regions, quality, contamination)
    assert (regions[0]["start"], regions[0]["end"]) == (1, 172524)
    assert regions[0]["viral_parts"] == [(100181, 162355)]
    assert fr.phage_answer("NODE_15", [(170000, 170900)], regions, "VirSorter2")[0] == "uncertain"


def test_checkv_and_caller_lengths_must_agree(tmp_path):
    """If CheckV's length for a region is not the caller's, the two do not
    describe the same DNA, and the workbook stops instead of guessing."""
    viral_table = write_tsv(tmp_path / "final-viral-boundary.tsv", VS2_HEADER, [
        {"seqname": "contig_3", "trim_bp_start": 1001, "trim_bp_end": 2000,
         "seqname_new": "contig_3||0_partial", "final_max_score": 0.9,
         "final_max_score_group": "dsDNAphage"}])
    quality = write_tsv(tmp_path / "quality_summary.tsv", QUALITY_HEADER, [
        {"contig_id": "contig_3||0_partial", "provirus": "No"}])
    contamination = write_tsv(tmp_path / "contamination.tsv", CONTAMINATION_HEADER, [
        {"contig_id": "contig_3||0_partial", "contig_length": 5000, "provirus": "No",
         "region_types": "NA", "region_coords_bp": "NA"}])
    regions = fr.read_viral_regions(viral_table, "virsorter2", {"contig_3": 9000})
    with pytest.raises(ValueError, match="5000 bp, but the caller's region is 1000 bp"):
        fr.attach_checkv(regions, quality, contamination)


def test_parts_left_by_a_failed_checkv_run_are_not_used(tmp_path):
    """CheckV writes contamination.tsv first; when it fails later, the rule leaves
    a header-only quality_summary.tsv. Without a quality row the region counts as
    not assessed, so a gene in it is uncertain, not yes."""
    regions = vs2_regions_with_checkv(tmp_path)
    write_tsv(tmp_path / "quality_summary.tsv", QUALITY_HEADER, [])
    fr.attach_checkv(regions, str(tmp_path / "quality_summary.tsv"),
                     str(tmp_path / "contamination.tsv"))
    assert regions[0]["viral_parts"] is None
    assert fr.phage_answer("contig_3", [(2120000, 2121000)], regions, "VirSorter2")[0] == \
        "uncertain"


def test_genomad_regions_whole_contig_and_prophage(tmp_path):
    """geNomad names a prophage contig|provirus_START_END with coordinates
    START-END, and a whole contig called a virus with coordinates NA."""
    table = write_tsv(tmp_path / "virus_summary.tsv",
                      ["seq_name", "length", "topology", "coordinates", "virus_score", "taxonomy"],
                      [{"seq_name": "contig_1|provirus_1000_5000", "topology": "Provirus",
                        "coordinates": "1000-5000", "virus_score": 0.97, "taxonomy": "Viruses"},
                       {"seq_name": "contig_2", "topology": "DTR", "coordinates": "NA",
                        "virus_score": 0.99, "taxonomy": "Viruses"}])
    regions = fr.read_viral_regions(table, "genomad", {"contig_1": 90000, "contig_2": 40000})
    assert [(r["contig"], r["start"], r["end"]) for r in regions] == [
        ("contig_1", 1000, 5000), ("contig_2", 1, 40000)]


def test_a_region_without_host_part_is_viral_as_a_whole(tmp_path):
    table = write_tsv(tmp_path / "virus_summary.tsv",
                      ["seq_name", "coordinates", "topology", "virus_score", "taxonomy"],
                      [{"seq_name": "contig_1|provirus_1000_5000", "coordinates": "1000-5000",
                        "topology": "Provirus", "virus_score": 0.97, "taxonomy": "Viruses"}])
    quality = write_tsv(tmp_path / "quality_summary.tsv", QUALITY_HEADER, [
        {"contig_id": "contig_1|provirus_1000_5000", "provirus": "No"}])
    contamination = write_tsv(tmp_path / "contamination.tsv", CONTAMINATION_HEADER, [
        {"contig_id": "contig_1|provirus_1000_5000", "contig_length": 4001, "provirus": "No",
         "region_types": "NA", "region_coords_bp": "NA"}])
    regions = fr.read_viral_regions(table, "genomad", {"contig_1": 90000})
    fr.attach_checkv(regions, quality, contamination)
    assert regions[0]["viral_parts"] == [(1000, 5000)]
    assert fr.phage_answer("contig_1", [(4000, 4500)], regions, "geNomad")[0] == "yes"


def test_a_gene_across_the_origin_is_tested_as_two_pieces():
    """Bakta writes a gene across the origin of a circular contig with start
    greater than end; it covers the contig's end and its start."""
    assert fr.gene_intervals(9900, 100, 10000) == [(9900, 10000), (1, 100)]
    regions = [region("r1", "contig_1", 9950, 10000, [(9950, 10000)], [])]
    assert fr.phage_answer("contig_1", fr.gene_intervals(9900, 100, 10000),
                           regions, "geNomad")[0] == "yes"


# ── ICE and IME: decided by coordinates for every gene ───────────────────────

@pytest.mark.parametrize("element_type", ["ice", "ime", "aice"])
def test_a_gene_inside_an_ice_ime_or_aice_is_yes(element_type):
    elements = [element(element_type, 10000, 60000)]
    answer, names, reason = fr.ice_ime_answer("contig_1", 20000, 21000, elements, True)
    assert answer == "yes"
    assert names == elements[0]["mge_id"]
    assert reason.startswith(f"inside {element_type}")


def test_a_gene_inside_a_conjugative_region_is_uncertain():
    elements = [element("conjugative_region", 10000, 60000,
                        mobility="conjugative plasmid region, not an ICE")]
    answer, _names, reason = fr.ice_ime_answer("contig_1", 20000, 21000, elements, True)
    assert answer == "uncertain"
    assert "conjugative plasmid region, not an ICE" in reason


def test_a_gene_near_an_ice_is_uncertain_with_its_distance():
    """The interval of an ICE is the span of its machinery genes, so a gene just
    outside it may still be inside the element."""
    elements = [element("ice", 10000, 60000)]
    answer, _names, reason = fr.ice_ime_answer("contig_1", 63001, 64000, elements, True)
    assert answer == "uncertain"
    assert "3,000 bp from ice" in reason


def test_a_gene_half_inside_an_ice_is_uncertain():
    """Inside means at least 90 % of the gene, as in colocalise.py; a gene half in
    the interval is near it, at distance 0."""
    elements = [element("ice", 10000, 60000)]
    answer, _names, reason = fr.ice_ime_answer("contig_1", 59500, 60500, elements, True)
    assert answer == "uncertain"
    assert "0 bp from ice" in reason


def test_a_gene_inside_a_genomic_island_is_no_with_the_island_named():
    elements = [element("genomic_island", 10000, 60000, mobility="passive")]
    answer, names, reason = fr.ice_ime_answer("contig_1", 20000, 21000, elements, True)
    assert answer == "no"
    assert names == elements[0]["mge_id"]
    assert "cannot be moved by conjugation on its own" in reason


def test_a_gene_near_a_conjugative_region_only_is_no_with_the_region_named():
    elements = [element("conjugative_region", 10000, 60000)]
    answer, _names, reason = fr.ice_ime_answer("contig_1", 62001, 63000, elements, True)
    assert answer == "no"
    assert "2,000 bp from conjugative_region" in reason


def test_a_gene_far_from_every_element_is_no():
    elements = [element("ice", 10000, 60000), element("ice", 10000, 60000, contig="contig_2")]
    answer, names, reason = fr.ice_ime_answer("contig_1", 80000, 81000, elements, True)
    assert (answer, names) == ("no", "NA")
    assert "no ICE, IME, AICE or conjugative region within 5,000 bp" in reason


def test_without_the_mobilome_module_ice_ime_is_not_assessed():
    answer, _names, reason = fr.ice_ime_answer("contig_1", 20000, 21000, [], False)
    assert answer == "not_assessed"
    assert "mobilome.run: false" in reason


def test_the_rules_are_the_mobilome_modules_own():
    assert fr.INSIDE_FRACTION == fr.colocalise.MIN_AMR_FRACTION_INSIDE_ELEMENT
    assert fr.NEAR_ELEMENT_BP == fr.colocalise.IS_ADJACENT_MAX_BP


# ── The mobilome module's mobility tiers ─────────────────────────────────────

def mobility_row(gene, element_type, start, end, tier, context="none", line=2,
                 method="BLASTP"):
    return {"sample": "SAMPLE1", "contig": "contig_1", "amr_gene": gene,
            "amr_element_type": element_type, "amrfinder_method": method,
            "amr_class": "CLASS", "amr_subclass": "SUB",
            "start": start, "end": end, "mge_context": context, "mobility_tier": str(tier),
            "mobility_tier_label": f"tier {tier}", "confidence": "medium", "__line__": line}


def test_a_stress_gene_gets_its_tier():
    """qacEdelta1 is a STRESS gene for AMRFinderPlus and a resistance gene for
    the databases; it is matched, and "none" is copied as written."""
    locus = {"contig": "contig_1", "start": 1000, "end": 1347}
    columns, matched = fr.mobility_columns(
        locus, [mobility_row("qacEdelta1", "STRESS", 1000, 1347, 1, "is_adjacent")], True)
    assert columns["amrfinderplus_genes"] == "qacEdelta1"
    assert columns["mobility_tier"] == 1
    assert columns["mge_context"] == "is_adjacent"
    assert len(matched) == 1


def test_a_virulence_gene_is_not_matched_to_a_resistance_locus():
    locus = {"contig": "contig_1", "start": 1000, "end": 2000}
    columns, matched = fr.mobility_columns(
        locus, [mobility_row("ybtP", "VIRULENCE", 1000, 2000, 6)], True)
    assert matched == []
    assert columns["mobility_tier"] == "NA"


def test_a_point_mutation_is_not_matched_to_a_database_hit_on_its_gene():
    """A point mutation is reported with its whole gene's coordinates; a CARD hit
    on the same gene (acrB, say) must not take the mutation's name and tier."""
    locus = {"contig": "contig_1", "start": 1000, "end": 4149}
    columns, matched = fr.mobility_columns(
        locus, [mobility_row("acrB_R717Q", "AMR", 1000, 4149, 1, method="POINTX")], True)
    assert matched == []
    assert columns["amrfinderplus_genes"] == "NA"


def test_the_highest_tier_is_shown_when_two_rows_match():
    locus = {"contig": "contig_1", "start": 1000, "end": 2000}
    rows = [mobility_row("geneA", "AMR", 1000, 2000, 2, line=2),
            mobility_row("geneA", "AMR", 1005, 1995, 5, "plasmid", line=3)]
    columns, _matched = fr.mobility_columns(locus, rows, True)
    assert columns["mobility_tier"] == 5
    assert columns["mge_context"] == "none; plasmid"


def test_without_the_mobilome_module_the_tiers_are_not_assessed():
    columns, matched = fr.mobility_columns({"contig": "c", "start": 1, "end": 2}, [], False)
    assert matched == []
    assert set(columns.values()) == {"not_assessed"}


# ── Copying values, paths and cells ──────────────────────────────────────────

def test_plain_numbers_become_numbers_but_names_stay_text():
    assert fr.typed("1234", "amr_start") == 1234
    assert fr.typed("1.5e-30", "evalue") == 1.5e-30
    assert fr.typed("1", "contig") == "1"
    assert fr.typed("-", "amr_strand") == "-"
    assert fr.typed("none", "mge_context") == "none"
    assert fr.typed("3.1.3.8", "ec") == "3.1.3.8"


def test_paths_are_shown_without_the_computers_folders(tmp_path):
    inside = tmp_path / "out" / "05.amr" / "card.tsv"
    assert fr.shown_path(str(inside), str(tmp_path / "out")) == os.path.join("05.amr", "card.tsv")
    assert fr.shown_path("/somewhere/else/trait_rules_v1.tsv", str(tmp_path / "out")) == \
        "trait_rules_v1.tsv"


def test_an_overlong_cell_is_cut_with_a_visible_mark():
    text = fr.cell_text("x" * 40000)
    assert len(text) == fr.EXCEL_CELL_LIMIT
    assert text.endswith(fr.CUT_MARK)
    assert fr.cell_text(None) == "NA"


# ── The whole AMR workbook, from a small output folder ───────────────────────

def write_fasta(path, lengths):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for name, length in lengths.items():
            handle.write(f">{name}\n")
            for start in range(0, length, 1000):
                handle.write("A" * min(1000, length - start) + "\n")


def make_output_folder(tmp_path, mobilome=True):
    """A BacFlux output folder for SAMPLE1: a 300 kb chromosome and a 60 kb
    plasmid (both short enough for Platon), three ABRicate tables, one VirSorter2
    region with CheckV parts, and, with the mobilome, one ICE and two
    AMRFinderPlus rows. Returns the command line of the amr command."""
    out = tmp_path / "out"
    write_fasta(out / "02.assembly/SAMPLE1/contigs_final.fasta",
                {"chromosome": 300000, "plasmid": 60000})
    write_tsv(out / "02.assembly/SAMPLE1/eval/SAMPLE1_qc_genomes.tsv",
              ["bin_id", "role", "technology", "source_path"],
              [{"bin_id": "SAMPLE1", "role": "primary", "technology": "NA",
                "source_path": str(out / "02.assembly/SAMPLE1/contigs_final.fasta")}])
    write_tsv(out / "03.taxonomy/SAMPLE1/classify/gtdbtk.bac120.summary.tsv",
              ["user_genome", "classification", "classification_method"],
              [{"user_genome": "SAMPLE1",
                "classification": "d__Bacteria;p__Pseudomonadota;g__Klebsiella;s__",
                "classification_method": "ani_screen"}])
    abricate = out / "05.amr/abricate/SAMPLE1"
    write_tsv(abricate / "card.tsv", ABRICATE_HEADER, [
        abricate_row("chromosome", 20000, 20800, "chromGene", "card", strand="-"),
        abricate_row("chromosome", 150000, 150900, "phageGene", "card"),
        abricate_row("plasmid", 5000, 5800, "aadA2", "card", gaps="6/4")])
    write_tsv(abricate / "resfinder.tsv", ABRICATE_HEADER, [
        abricate_row("plasmid", 5010, 5790, "aadA2_1", "resfinder")])
    write_tsv(abricate / "argannot.tsv", ABRICATE_HEADER, [])
    write_tsv(out / "05.amr/abricate/abricate_databases.tsv",
              ["DATABASE", "SEQUENCES", "DBTYPE", "DATE", "ABRICATE_VERSION"],
              [{"DATABASE": name, "SEQUENCES": 1000, "DBTYPE": "nucl", "DATE": "2026-Apr-3",
                "ABRICATE_VERSION": "1.4.0"} for name in ("card", "resfinder", "argannot")])

    platon = out / "06.plasmids/SAMPLE1/platon"
    platon.mkdir(parents=True)
    (platon / "contigs_final.chromosome.fasta").write_text(">chromosome\nA\n")
    (platon / "contigs_final.plasmid.fasta").write_text(">plasmid\nA\n")
    write_tsv(platon / "contigs_final.tsv", ["ID", "# Conjugation", "# Mobilization", "# OriT"],
              [{"ID": "plasmid", "# Conjugation": 3, "# Mobilization": 1, "# OriT": 1}])
    (platon / "verified_plasmids.txt").write_text("SAMPLE1: plasmid is a plasmid.\n")

    write_tsv(out / "07.phages/virsorter/SAMPLE1/final-viral-boundary.tsv", VS2_HEADER, [
        {"seqname": "chromosome", "trim_bp_start": 140001, "trim_bp_end": 180000,
         "seqname_new": "chromosome||0_partial", "final_max_score": 0.95,
         "final_max_score_group": "dsDNAphage"}])
    write_tsv(out / "07.phages/checkv/SAMPLE1/quality_summary.tsv", QUALITY_HEADER, [
        {"contig_id": "chromosome||0_partial", "provirus": "Yes",
         "checkv_quality": "Medium-quality", "completeness": "70.0", "contamination": "25.0",
         "warnings": ""}])
    write_tsv(out / "07.phages/checkv/SAMPLE1/contamination.tsv", CONTAMINATION_HEADER, [
        {"contig_id": "chromosome||0_partial", "contig_length": 40000, "provirus": "Yes",
         "region_types": "host,viral", "region_coords_bp": "1-5000,5001-40000"}])

    args = ["amr", "--sample", "SAMPLE1", "--output-root", str(out),
            "--qc-genomes", str(out / "02.assembly/SAMPLE1/eval/SAMPLE1_qc_genomes.tsv"),
            "--taxonomy-dir", str(out / "03.taxonomy/SAMPLE1"),
            "--contigs", str(out / "02.assembly/SAMPLE1/contigs_final.fasta"),
            "--card", str(abricate / "card.tsv"), "--resfinder", str(abricate / "resfinder.tsv"),
            "--argannot", str(abricate / "argannot.tsv"),
            "--abricate-databases", str(out / "05.amr/abricate/abricate_databases.tsv"),
            "--platon-dir", str(platon), "--platon-prefix", "contigs_final",
            "--phage-caller", "virsorter2",
            "--viral-table", str(out / "07.phages/virsorter/SAMPLE1/final-viral-boundary.tsv"),
            "--checkv-quality", str(out / "07.phages/checkv/SAMPLE1/quality_summary.tsv"),
            "--checkv-contamination", str(out / "07.phages/checkv/SAMPLE1/contamination.tsv"),
            "--out", str(out / "10.workbooks/SAMPLE1_amr.xlsx")]
    if mobilome:
        mobilome_dir = out / "08.mobilome/SAMPLE1"
        write_tsv(mobilome_dir / "SAMPLE1_ice_candidates.tsv",
                  ["sample", "mge_id", "mge_name", "element_type", "mobility", "contig",
                   "start", "end", "boundary_method"],
                  [{"sample": "SAMPLE1", "mge_id": "chromosome|ice-10000:40000",
                    "mge_name": "NA", "element_type": "ice",
                    "mobility": "predicted self-transmissible", "contig": "chromosome",
                    "start": 10000, "end": 40000, "boundary_method": "none"}])
        write_tsv(mobilome_dir / "SAMPLE1_amr_mobility.tsv",
                  ["sample", "contig", "amr_gene", "amr_element_type", "amr_class",
                   "amr_subclass", "amr_start", "amr_end", "amrfinder_method",
                   "mge_context", "mobility_tier", "mobility_tier_label", "confidence"],
                  [{"sample": "SAMPLE1", "contig": "plasmid", "amr_gene": "aadA2",
                    "amr_element_type": "AMR", "amr_class": "AMINOGLYCOSIDE",
                    "amr_subclass": "STREPTOMYCIN", "amr_start": 5000, "amr_end": 5800,
                    "amrfinder_method": "EXACTP",
                    "mge_context": "plasmid", "mobility_tier": 6,
                    "mobility_tier_label": "conjugative plasmid", "confidence": "high"},
                   {"sample": "SAMPLE1", "contig": "chromosome", "amr_gene": "gyrA_S83I",
                    "amr_element_type": "AMR", "amr_class": "QUINOLONE",
                    "amr_subclass": "QUINOLONE", "amr_start": 250000, "amr_end": 252000,
                    "amrfinder_method": "POINTP",
                    "mge_context": "none", "mobility_tier": 1,
                    "mobility_tier_label": "intrinsic candidate", "confidence": "high"}])
        write_tsv(mobilome_dir / "SAMPLE1_amr_mobility_audit.tsv",
                  ["sample", "contig", "amr_gene", "amr_start", "amr_end", "decision",
                   "reason", "detail"],
                  [{"sample": "SAMPLE1", "contig": "chromosome", "amr_gene": "gyrA_S83I",
                    "amr_start": 250000, "amr_end": 252000, "decision": "tier_assigned",
                    "reason": "no_mobile_element", "detail": "none within 5000 bp"}])
        args += ["--ice-table", str(mobilome_dir / "SAMPLE1_ice_candidates.tsv"),
                 "--mobility", str(mobilome_dir / "SAMPLE1_amr_mobility.tsv"),
                 "--mobility-audit", str(mobilome_dir / "SAMPLE1_amr_mobility_audit.tsv")]
    return args


def sheet_rows(sheets, name):
    """The rows of one sheet as dicts keyed by its columns."""
    sheet = next(s for s in sheets if s["name"] == name)
    return [dict(zip(sheet["columns"], row)) for row in sheet["rows"]]


def test_the_amr_workbook_answers_every_question_for_every_locus(tmp_path):
    args = fr.make_parser().parse_args(make_output_folder(tmp_path))
    sheets, counts = fr.build_amr_workbook(args)
    assert counts == {"database hits": 4, "loci": 3}
    loci = {row["genes"].split("; ")[0]: row for row in sheet_rows(sheets, "AMR_loci")}
    assert set(loci) == {"chromGene", "phageGene", "aadA2"}

    chrom = loci["chromGene"]
    assert (chrom["plasmid"], chrom["ice_ime"], chrom["phage"]) == ("no", "yes", "no")
    assert chrom["strand"] == "-"
    assert chrom["mobility_tier_label"].startswith("no AMRFinderPlus")

    # In the VirSorter2 region (140,001-180,000) and in CheckV's viral part.
    phage = loci["phageGene"]
    assert phage["phage"] == "yes"
    assert phage["ice_ime"] == "no"

    plasmid = loci["aadA2"]
    assert (plasmid["plasmid"], plasmid["plasmid_mobility"]) == ("yes", "conjugative")
    assert plasmid["databases"] == "CARD; ResFinder"
    assert plasmid["mobility_tier"] == 6

    # ResFinder, the database of acquired genes, found only aadA2.
    assert (plasmid["in_resfinder"], chrom["in_resfinder"]) == ("yes", "no")
    summary = {row["item"]: row for row in sheet_rows(sheets, "Summary")}
    assert summary["Found by ResFinder"]["value"] == 1
    # The database version and dates are reported with the results.
    assert summary["Databases"]["value"] == "ABRicate 1.4.0"
    assert "CARD, dated 2026-Apr-3, 1000 sequences" in summary["Databases"]["note"]

    hits = sheet_rows(sheets, "Database_hits")
    assert {h["gaps"] for h in hits} == {"0/0", "6/4"}
    assert "-" in {h["strand"] for h in hits}

    # The point mutation has no database locus, and is kept without one.
    mobility = sheet_rows(sheets, "Mobility_evidence")
    gyrase = next(row for row in mobility if row["amr_gene"] == "gyrA_S83I")
    assert gyrase["locus_ids"] == "NA"
    assert gyrase["amr_start"] == 250000           # a number, not text

    # No cell shows the folders of the computer the analysis ran on.
    for sheet in sheets:
        for row in sheet["rows"]:
            assert not any(str(tmp_path) in str(value) for value in row), sheet["name"]


def test_without_the_mobilome_module_the_workbook_says_not_assessed(tmp_path):
    args = fr.make_parser().parse_args(make_output_folder(tmp_path, mobilome=False))
    sheets, _counts = fr.build_amr_workbook(args)
    for row in sheet_rows(sheets, "AMR_loci"):
        assert row["ice_ime"] == "not_assessed"
        assert row["mobility_tier"] == "not_assessed"
    assert sheet_rows(sheets, "Mobility_evidence")[0]["note"].startswith("The mobilome module")


def test_the_three_mobilome_tables_go_together(tmp_path):
    words = make_output_folder(tmp_path)
    words = words[:words.index("--mobility")]
    with pytest.raises(ValueError, match="go together"):
        fr.main(words)


def test_every_column_has_a_meaning_in_the_data_dictionary(tmp_path):
    args = fr.make_parser().parse_args(make_output_folder(tmp_path))
    sheets, _counts = fr.build_amr_workbook(args)
    explained = {(row["sheet"], row["column"]) for row in sheet_rows(sheets, "Data_dictionary")}
    for sheet in sheets:
        if sheet["name"] in ("Data_dictionary", "Mobility_evidence", "Mobility_audit"):
            continue
        for column in sheet["columns"]:
            assert (sheet["name"], column) in explained


def test_the_workbook_file_has_every_sheet(tmp_path):
    pytest.importorskip("xlsxwriter")
    words = make_output_folder(tmp_path)
    (tmp_path / "out/10.workbooks").mkdir()
    assert fr.main(words) == 0
    with zipfile.ZipFile(words[words.index("--out") + 1]) as workbook:
        index = workbook.read("xl/workbook.xml").decode("utf-8")
    for name in ("Summary", "AMR_loci", "Database_hits", "ICE_IME_elements", "Viral_regions",
                 "Mobility_evidence", "Mobility_audit", "Run_QC", "Data_dictionary"):
        assert f'name="{name}"' in index


# ── The traits workbook: the catalogue must be the one the calls used ────────

def test_the_traits_workbook_stops_when_the_catalogue_changed(tmp_path):
    rules = tmp_path / "trait_rules_v1.tsv"
    rules.write_text("trait_id\ttrait\n", encoding="utf-8")
    profile = tmp_path / "SAMPLE1_trait_profile.tsv"
    profile.write_text("sample\ttrait_id\n", encoding="utf-8")
    digest = hashlib.sha256(rules.read_bytes()).hexdigest()
    record = tmp_path / "SAMPLE1_trait_provenance.json"
    record.write_text(json.dumps({
        "schema_version": "1.0", "sample": "SAMPLE1",
        "inputs": {"rules": {"file": "trait_rules_v1.tsv", "sha256": digest}},
        "outputs": {"profile": {"file": profile.name,
                                "sha256": hashlib.sha256(profile.read_bytes()).hexdigest()}},
    }), encoding="utf-8")
    inputs = {"rules": str(rules)}
    outputs = {"profile": str(profile)}
    assert fr.verify_trait_provenance(str(record), "SAMPLE1", inputs, outputs)["sample"] == "SAMPLE1"
    rules.write_text("trait_id\ttrait\nnew\tNew trait\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not the file the trait calls were made with"):
        fr.verify_trait_provenance(str(record), "SAMPLE1", inputs, outputs)
    with pytest.raises(ValueError, match="written for sample"):
        fr.verify_trait_provenance(str(record), "SAMPLE2", inputs, outputs)
