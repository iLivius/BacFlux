"""Unit tests for trait_profile.py, the trait call of the EcoFoldDB trait layer.

trait_profile.py turns up to three kinds of evidence per protein (an EcoFoldDB
structure hit, the Bakta name, the Baktfold name) into one decision per trait:
called, partial or not called. A wrong call is a wrong biological claim about a
strain ("this isolate makes a siderophore"), so most tests pin down when a gene
must NOT count and when a trait must NOT be called.

Everything is pure text handling. The Bakta and Baktfold tables, the Foldseek hit
tables, the catalogue, the antiSMASH output and both trait tables are invented and
written into pytest's tmp_path in the real column layouts: no EcoFoldDB rows and
no real genomes. Gene names such as abcA and catalogue IDs such as CAT_ABCA are
made up. The tests at the end check the two shipped tables themselves.

Run: pytest workflow/scripts/40_annotation/test_trait_profile.py -q
"""

import csv
import gzip
import json
import os

import pytest

import trait_profile as tp


# ── Fixtures: invented inputs in the real layouts ─────────────────────────────

# Bakta's own header lines, as Bakta writes them above its {sample}.tsv table.
BAKTA_HEADER = [
    "# Annotated with Bakta",
    "# Software: v1.12.1",
    "# Database: v6.0, full",
    "# DOI: 10.1099/mgen.0.000685",
    "# URL: github.com/oschwengers/bakta",
    "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs",
]

# The real column layout of EcoFoldDB_descriptions.txt.
CATALOGUE_HEADER = ["ID", "Gene", "Protein", "Category", "Sub-category",
                    "Pathway/Activity", "EC", "KO", "UniProt"]

# An invented catalogue: (ID, Gene, Pathway/Activity). CAT_ABCA2 is a second
# structure with the same key (Gene@Pathway/Activity) as CAT_ABCA; xyzA is a
# catalogue gene that no trait uses.
CATALOGUE = [
    ("CAT_ABCA", "abcA", "invented route"),
    ("CAT_ABCA2", "abcA", "invented route"),
    ("CAT_ABCB", "abcB", "invented route"),
    ("CAT_ABCC", "abcC", "invented route"),
    ("CAT_XYZA", "xyzA", "unused route"),
]

# Three more invented catalogue genes, for the few tests that need more than three.
CATALOGUE_SIX = CATALOGUE + [
    ("CAT_ABCD", "abcD", "invented route"),
    ("CAT_ABCE", "abcE", "invented route"),
    ("CAT_ABCF", "abcF", "invented route"),
]


def protein(locus, product="hypothetical protein", gene="", contig="contig_1",
            start=None, strand="+", dbxrefs="", feature="cds"):
    """One Bakta table row: contig, type (cds or sorf), start, stop, strand, locus
    tag, gene symbol, product, DbXrefs. Without a start, the gene sits at 1000
    times its locus number (SAMPLE1_00003 at 3000), so the locus order is the gene
    order; every gene is 900 bases long."""
    if start is None:
        start = int(locus.split("_")[-1]) * 1000
    return [contig, feature, str(start), str(start + 899), strand, locus, gene,
            product, dbxrefs]


def fillers(first, last, contig="contig_1"):
    """Unnamed proteins SAMPLE1_<first>..SAMPLE1_<last>, to space genes apart."""
    return [protein(f"SAMPLE1_{number:05d}", contig=contig)
            for number in range(first, last + 1)]


def gene_entry(gene, **changes):
    """One row of the gene dictionary (trait_genes_v1.tsv): recognised by its own
    catalogue entry and by the word itself in the product (case ignored);
    s_alone no unless changed."""
    entry = {
        "gene": gene,
        "protein": f"invented protein {gene}",
        "ecofolddb": f"{gene}@invented route",
        "ko": "",
        "ec": "",
        "bakta_gene": "",
        "symbol_needs_product": "no",
        "product_include": rf"\b{gene}\b",
        "product_exclude": "",
        "s_alone": "no",
        "notes": "Invented entry for a test.",
    }
    entry.update(changes)
    return entry


def rule_entry(trait_id, required, **changes):
    """One row of the rules table, with invented text in the descriptive columns."""
    rule = {
        "trait_id": trait_id,
        "trait": f"Invented trait {trait_id}",
        "application": "biostimulant",
        "required": required,
        "supporting": "",
        "colocalised": "",
        "key_genes": "",
        "caveat": "Genes present, activity not shown.",
        "assay": "An invented assay.",
        "reference": "Invented et al. 2020, doi:10.0000/invented",
    }
    rule.update(changes)
    return rule


def write_table(path, columns, rows):
    """Write a trait table: a '#' provenance line, the header, the rows."""
    lines = ["# Invented test table in the real column layout."]
    lines.append("\t".join(columns))
    lines += ["\t".join(row[column] for column in columns) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def write_hits(path, hits):
    """Write a Foldseek table: (query, target, evalue, qcov, tcov) per line,
    gzip-compressed when the name ends in .gz."""
    text = "".join("\t".join(hit) + "\n" for hit in hits)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return str(path)


def write_inputs(tmp_path, proteins, genes, rules, hits=(), catalogue=CATALOGUE,
                 lengths=None, baktfold=None, antismash=None, hits_name="hits.tsv"):
    """Write every input file into tmp_path; return the arguments for main().

    The .faa gives every protein 300 amino acids unless lengths says otherwise,
    and the catalogue rows get filler in the columns the script does not read.
    The Baktfold and antiSMASH files are written only when given, so a test can
    run with or without them.
    """
    bakta = tmp_path / "SAMPLE1.tsv"
    bakta.write_text("\n".join(BAKTA_HEADER + ["\t".join(row) for row in proteins])
                     + "\n", encoding="utf-8")

    lengths = lengths or {}
    faa_lines = []
    for row in proteins:
        faa_lines.append(f">{row[5]} {row[7]}")
        faa_lines.append("M" * lengths.get(row[5], 300))
    faa = tmp_path / "SAMPLE1.faa"
    faa.write_text("\n".join(faa_lines) + "\n", encoding="utf-8")

    catalogue_lines = ["\t".join(CATALOGUE_HEADER)]
    for target, gene, pathway in catalogue:
        catalogue_lines.append("\t".join([target, gene, "invented protein", "Invented",
                                          "Invented", pathway, "", "", ""]))
    catalogue_file = tmp_path / "EcoFoldDB_descriptions.txt"
    catalogue_file.write_text("\n".join(catalogue_lines) + "\n", encoding="utf-8")

    args = [
        "--sample", "SAMPLE1",
        "--bakta-tsv", str(bakta),
        "--bakta-faa", str(faa),
        "--hits", write_hits(tmp_path / hits_name, hits),
        "--catalogue", str(catalogue_file),
        "--rules", write_table(tmp_path / "rules.tsv", tp.RULE_COLUMNS, rules),
        "--genes", write_table(tmp_path / "genes.tsv", tp.GENE_COLUMNS, genes),
        "--out-genes", str(tmp_path / "out" / "SAMPLE1_trait_genes.tsv"),
        "--out-profile", str(tmp_path / "out" / "SAMPLE1_trait_profile.tsv"),
        "--out-audit", str(tmp_path / "out" / "SAMPLE1_trait_profile_audit.tsv"),
    ]
    if baktfold is not None:
        baktfold_file = tmp_path / "SAMPLE1_baktfold.tsv"
        baktfold_file.write_text(
            "\n".join(BAKTA_HEADER + ["\t".join(row) for row in baktfold]) + "\n",
            encoding="utf-8")
        args += ["--baktfold-tsv", str(baktfold_file)]
    if antismash is not None:
        antismash_file = tmp_path / "SAMPLE1.json"
        antismash_file.write_text(json.dumps({"version": "8.0.4",
                                              "records": antismash}), encoding="utf-8")
        args += ["--antismash-json", str(antismash_file)]
    return args


def read_tsv(path):
    """Return (header, rows as dicts) of an output table, split on tabs only, and
    check on the way that every line has as many fields as the header."""
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle]
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    for line in lines:
        assert len(line.split("\t")) == len(header), line
    return header, rows


def run(tmp_path, proteins, genes, rules, extra_args=(), **inputs):
    """Write the inputs, run main(), and return the three outputs.

    Returns a dict: "genes" (gene-table rows), "profile" and "audit" (rows keyed
    by trait_id).
    """
    args = write_inputs(tmp_path, proteins, genes, rules, **inputs) + list(extra_args)
    assert tp.main(args) == 0
    out = tmp_path / "out"
    _, gene_rows = read_tsv(out / "SAMPLE1_trait_genes.tsv")
    _, profile_rows = read_tsv(out / "SAMPLE1_trait_profile.tsv")
    _, audit_rows = read_tsv(out / "SAMPLE1_trait_profile_audit.tsv")
    return {
        "genes": gene_rows,
        "profile": {row["trait_id"]: row for row in profile_rows},
        "audit": {row["trait_id"]: row for row in audit_rows},
    }


def rows_for(result, locus, gene=None):
    """The gene-table rows of one protein (and one gene, when given)."""
    return [row for row in result["genes"]
            if row["locus_tag"] == locus and (gene is None or row["gene"] == gene)]


# Three genes and one trait that needs all three, used by many tests below.
ABC_GENES = [gene_entry("abcA"), gene_entry("abcB"), gene_entry("abcC")]


def abc_rule(**changes):
    """The rule abc_trait, which needs abcA+abcB+abcC, with any column changed."""
    return rule_entry("abc_trait", "abcA+abcB+abcC", **changes)


def named_abc(contig_of_c="contig_1"):
    """Three adjacent proteins that Bakta names AbcA, AbcB and AbcC."""
    return [protein("SAMPLE1_00001", "AbcA protein"),
            protein("SAMPLE1_00002", "AbcB protein"),
            protein("SAMPLE1_00003", "AbcC protein", contig=contig_of_c)]


# ── Decisions: called, partial, not called ───────────────────────────────────
# How the genes counted in a genome become a decision: complete sets,
# co-localisation, half sets, key_genes and the order of alternative sets. Most
# proteins here are named by Bakta (evidence Q), which keeps the inputs short.

def test_complete_colocalised_set_is_called(tmp_path):
    """All three genes named, next to each other: the trait is called, and the
    profile lists each gene with its protein and evidence."""
    result = run(tmp_path, named_abc(), ABC_GENES, [abc_rule(colocalised="all")])

    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "called"
    assert audit["reason"] == "set_complete"
    assert audit["detail"] == ("abcA+abcB+abcC complete; co-localised: "
                               "abcA+abcB+abcC")
    profile = result["profile"]["abc_trait"]
    assert profile["genes"] == ("abcA:SAMPLE1_00001(Q); abcB:SAMPLE1_00002(Q); "
                                "abcC:SAMPLE1_00003(Q)")
    assert profile["missing"] == "NA"
    assert profile["borderline"] == "no"
    assert profile["structure_only"] == "NA"
    assert profile["caveat"] == "Genes present, activity not shown."


def test_complete_set_on_two_contigs_is_partial_when_it_must_be_colocalised(tmp_path):
    """A pathway whose genes must form one cluster is not called when one gene sits
    on another contig: that could be a stray paralogue, so the call is partial."""
    result = run(tmp_path, named_abc(contig_of_c="contig_2"), ABC_GENES,
                 [abc_rule(colocalised="all")])

    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "set_complete_not_colocalised"
    assert audit["detail"].startswith("abcA+abcB+abcC complete but not co-localised")
    profile = result["profile"]["abc_trait"]
    assert profile["missing"] == "NA"
    assert "abcC:SAMPLE1_00003(Q)" in profile["genes"]


def test_without_a_colocalisation_condition_contigs_do_not_matter(tmp_path):
    """Without a colocalised condition, a complete set is called wherever its genes
    lie, here with abcC on another contig."""
    result = run(tmp_path, named_abc(contig_of_c="contig_2"), ABC_GENES, [abc_rule()])
    assert result["audit"]["abc_trait"]["decision"] == "called"
    assert result["audit"]["abc_trait"]["detail"] == "abcA+abcB+abcC complete"


def test_colocalised_chain_counts_each_step_not_the_whole_span(tmp_path):
    """Co-localised means each gene at most max_gene_gap genes from the previous
    one, as in an operon interrupted by other genes: abcA, abcB and abcC eight
    genes apart each are co-localised with the default gap of 10, although abcA
    and abcC are sixteen apart. With a gap of 7 they are not."""
    proteins = ([protein("SAMPLE1_00001", "AbcA protein")] + fillers(2, 8)
                + [protein("SAMPLE1_00009", "AbcB protein")] + fillers(10, 16)
                + [protein("SAMPLE1_00017", "AbcC protein")])
    rules = [abc_rule(colocalised="all")]

    assert run(tmp_path, proteins, ABC_GENES, rules)["audit"]["abc_trait"][
        "decision"] == "called"
    narrow = run(tmp_path, proteins, ABC_GENES, rules,
                 extra_args=["--max-gene-gap", "7"])
    assert narrow["audit"]["abc_trait"]["decision"] == "partial"


def test_short_orfs_count_in_the_gene_distance(tmp_path):
    """The CDS index ranks cds AND sorf rows, the proteins of Bakta's .faa. Ten
    steps apart is co-localised with the default gap; one sORF in between makes it
    eleven, which is not."""
    genes = [gene_entry("abcA"), gene_entry("abcB")]
    rules = [rule_entry("ab_trait", "abcA+abcB", colocalised="all")]
    base = ([protein("SAMPLE1_00001", "AbcA protein")] + fillers(2, 10)
            + [protein("SAMPLE1_00011", "AbcB protein")])

    assert run(tmp_path, base, genes, rules)["audit"]["ab_trait"]["decision"] == "called"

    with_sorf = base + [protein("SAMPLE1_00099", start=5500, feature="sorf")]
    result = run(tmp_path, with_sorf, genes, rules)
    assert result["audit"]["ab_trait"]["decision"] == "partial"
    assert rows_for(result, "SAMPLE1_00011", "abcB")[0]["cds_index"] == "11"


def test_half_of_a_set_is_partial_and_less_is_not(tmp_path):
    """Two of three genes (at least half, rounded up) give a partial call that names
    the missing gene; one of three is not called, and the audit names the closest
    set and what is missing from it."""
    two = run(tmp_path, named_abc()[:2], ABC_GENES, [abc_rule()])
    audit = two["audit"]["abc_trait"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "half_of_set_found"
    assert audit["detail"] == "missing abcC (2 of 3 genes of abcA+abcB+abcC found)"
    assert two["profile"]["abc_trait"]["missing"] == "abcC"

    one = run(tmp_path, named_abc()[:1], ABC_GENES, [abc_rule()])
    audit = one["audit"]["abc_trait"]
    assert audit["decision"] == "not called"
    assert audit["reason"] == "too_few_genes"
    assert audit["detail"] == ("1 of 3 genes found in the closest set abcA+abcB+abcC "
                               "(missing abcB, abcC)")
    assert one["profile"]["abc_trait"]["genes"] == "NA"


def test_one_gene_sets_have_no_half(tmp_path):
    """A trait made of single-gene alternatives has no half of one gene: with no
    pseudogene in play, it is either called or not. (A gene found only as a
    pseudogene makes it partial; see the pseudogene tests below.)"""
    genes = [gene_entry("abcA"), gene_entry("abcB")]
    rules = [rule_entry("single", "abcA | abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "unrelated enzyme")], genes, rules)
    assert result["audit"]["single"]["decision"] == "not called"
    assert result["audit"]["single"]["reason"] == "no_gene_found"

    called = run(tmp_path, [protein("SAMPLE1_00001", "AbcB protein")], genes, rules)
    assert called["audit"]["single"]["decision"] == "called"
    assert called["audit"]["single"]["detail"] == "abcB complete"


def test_key_genes_withhold_a_partial_without_a_key_gene(tmp_path):
    """key_genes names the genes that set a pathway apart from look-alikes, such as
    the enzyme that commits it: abcB and abcC alone, without abcA, are not
    reported as partial."""
    rules = [abc_rule(key_genes="abcA")]
    result = run(tmp_path, named_abc()[1:], ABC_GENES, rules)
    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "not called"
    assert audit["reason"] == "partial_withheld"
    assert audit["detail"] == ("no key gene (abcA) found (missing abcA (2 of 3 genes "
                               "of abcA+abcB+abcC found))")
    assert result["profile"]["abc_trait"]["genes"] == "NA"

    with_key_gene = run(tmp_path, named_abc()[:2], ABC_GENES, rules)
    assert with_key_gene["audit"]["abc_trait"]["decision"] == "partial"


def test_withheld_detail_names_every_key_gene(tmp_path):
    """The detail of a withheld partial names the rule's key genes as
    alternatives: 'abcA or abcB' for two, 'abcA, abcB or abcC' for three. (Rule
    rest is there only to use the gene entries abcE and abcF: every entry must
    belong to a trait.)"""
    six = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC",
                                         "abcD", "abcE", "abcF")]
    two_keys = [rule_entry("r", "abcA+abcB+abcC+abcD", key_genes="abcA | abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00003", "AbcC protein"),
                            protein("SAMPLE1_00004", "AbcD protein")],
                 six, two_keys + [rule_entry("rest", "abcE+abcF")],
                 catalogue=CATALOGUE_SIX)
    assert result["audit"]["r"]["detail"] == (
        "no key gene (abcA or abcB) found (missing abcA, abcB (2 of 4 genes of "
        "abcA+abcB+abcC+abcD found))")

    three_keys = [rule_entry("r", "abcA+abcB+abcC+abcD+abcE+abcF",
                             key_genes="abcA | abcB | abcC")]
    result = run(tmp_path, [protein("SAMPLE1_00004", "AbcD protein"),
                            protein("SAMPLE1_00005", "AbcE protein"),
                            protein("SAMPLE1_00006", "AbcF protein")],
                 six, three_keys, catalogue=CATALOGUE_SIX)
    assert result["audit"]["r"]["detail"] == (
        "no key gene (abcA, abcB or abcC) found (missing abcA, abcB, abcC (3 of 6 "
        "genes of abcA+abcB+abcC+abcD+abcE+abcF found))")


def test_key_genes_accept_any_of_their_genes(tmp_path):
    """key_genes lists alternatives, and any one found is enough: with abcC counted,
    the half-found abcB+abcC of abcA+abcB+abcC is reported as partial, although
    abcA, the other key gene, is absent."""
    rules = [abc_rule(key_genes="abcA | abcC")]
    result = run(tmp_path, named_abc()[1:], ABC_GENES, rules)
    assert result["audit"]["abc_trait"]["decision"] == "partial"


def test_colocalisation_group_of_another_set_is_ignored(tmp_path):
    """A group applies only to a set that holds all its genes. The group abcA+abcB
    belongs to the first alternative; the second alternative, abcA+abcC, has no
    condition and is called with its genes on two contigs."""
    rules = [rule_entry("alt", "abcA+abcB | abcA+abcC", colocalised="abcA+abcB")]
    proteins = [protein("SAMPLE1_00001", "AbcA protein"),
                protein("SAMPLE1_00003", "AbcC protein", contig="contig_2")]
    result = run(tmp_path, proteins, ABC_GENES, rules)
    assert result["audit"]["alt"]["decision"] == "called"
    assert result["audit"]["alt"]["detail"] == "abcA+abcC complete"


def test_first_set_in_table_order_wins(tmp_path):
    """When two alternative sets are both complete, the first in table order is the
    call, so the same genome always gives the same detail."""
    rules = [rule_entry("alt", "abcA+abcB | abcA+abcC")]
    result = run(tmp_path, named_abc(), ABC_GENES, rules)
    assert result["audit"]["alt"]["detail"] == "abcA+abcB complete"


def test_a_copy_on_another_contig_does_not_block_a_cluster(tmp_path):
    """A second abcA on another contig, beside a full abcA-abcB-abcC cluster on
    contig_1: the chain uses the copy inside the cluster, so a stray paralogue
    cannot hide a real operon."""
    proteins = named_abc() + [protein("SAMPLE1_00009", "AbcA protein",
                                      contig="contig_2")]
    result = run(tmp_path, proteins, ABC_GENES, [abc_rule(colocalised="all")])
    assert result["audit"]["abc_trait"]["decision"] == "called"
    assert result["profile"]["abc_trait"]["genes"].startswith(
        "abcA:SAMPLE1_00001(Q),SAMPLE1_00009(Q); ")


def test_contigs_are_treated_as_linear(tmp_path):
    """Bakta's table does not say which contigs are circular, so the first and the
    last gene of a contig are not neighbours. Here abcA is the first gene and
    abcB, abcC the last two: across the origin of a circular contig they would be
    next-door genes, but they read as not co-localised, a partial, never a wrong
    call."""
    proteins = ([protein("SAMPLE1_00001", "AbcA protein")] + fillers(2, 20)
                + [protein("SAMPLE1_00021", "AbcB protein"),
                   protein("SAMPLE1_00022", "AbcC protein")])
    result = run(tmp_path, proteins, ABC_GENES, [abc_rule(colocalised="all")])
    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "set_complete_not_colocalised"


def test_detail_names_the_group_that_is_not_colocalised(tmp_path):
    """Two groups, abcA+abcB side by side and abcC+abcD eighteen genes apart: the
    partial's detail names only the group that failed."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    rules = [rule_entry("r", "abcA+abcB+abcC+abcD",
                        colocalised="abcA+abcB; abcC+abcD")]
    proteins = ([protein("SAMPLE1_00001", "AbcA protein"),
                 protein("SAMPLE1_00002", "AbcB protein"),
                 protein("SAMPLE1_00003", "AbcC protein")] + fillers(4, 20)
                + [protein("SAMPLE1_00021", "AbcD protein")])
    result = run(tmp_path, proteins, genes, rules, catalogue=CATALOGUE_SIX)
    audit = result["audit"]["r"]
    assert audit["reason"] == "set_complete_not_colocalised"
    assert audit["detail"] == ("abcA+abcB+abcC+abcD complete but not co-localised "
                               "(abcC+abcD: not within 10 genes of each other on "
                               "one contig)")


def test_key_gene_may_count_outside_the_partial_set(tmp_path):
    """key_genes asks whether a gene counts anywhere in the genome. Here abcD
    counts only in the second alternative, and that is enough to report the
    partial of the first."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    rules = [rule_entry("r", "abcA+abcB+abcC | abcA+abcD", key_genes="abcD")]
    proteins = [protein("SAMPLE1_00002", "AbcB protein"),
                protein("SAMPLE1_00003", "AbcC protein"),
                protein("SAMPLE1_00004", "AbcD protein")]
    result = run(tmp_path, proteins, genes, rules, catalogue=CATALOGUE_SIX)
    audit = result["audit"]["r"]
    assert audit["decision"] == "partial"
    assert audit["detail"] == "missing abcA (2 of 3 genes of abcA+abcB+abcC found)"


def test_closest_set_tie_goes_to_the_first_set(tmp_path):
    """Not called, one gene found in each of two sets: the audit names the first
    set in table order as the closest."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC",
                                           "abcD", "abcE", "abcF")]
    rules = [rule_entry("r", "abcA+abcB+abcC | abcD+abcE+abcF")]
    proteins = [protein("SAMPLE1_00001", "AbcA protein"),
                protein("SAMPLE1_00004", "AbcD protein")]
    result = run(tmp_path, proteins, genes, rules, catalogue=CATALOGUE_SIX)
    assert result["audit"]["r"]["detail"] == ("1 of 3 genes found in the closest "
                                              "set abcA+abcB+abcC (missing abcB, "
                                              "abcC)")


@pytest.mark.parametrize("gap", ["0", "-1"])
def test_max_gene_gap_below_one_is_refused(tmp_path, capsys, gap):
    """Next-door genes are 1 apart, so a gap of 0 would let no two different
    proteins be co-localised. The option is refused before any file is read."""
    args = write_inputs(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                        [rule_entry("a", "abcA")])
    with pytest.raises(SystemExit):
        tp.main(args + ["--max-gene-gap", gap])
    assert "--max-gene-gap must be 1 or more" in capsys.readouterr().err


def test_max_gene_gap_of_one_joins_only_next_door_genes(tmp_path):
    """The smallest gap, 1: next-door genes are co-localised, and one gene in
    between breaks the chain."""
    rules = [abc_rule(colocalised="all")]
    next_door = run(tmp_path, named_abc(), ABC_GENES, rules,
                    extra_args=["--max-gene-gap", "1"])
    assert next_door["audit"]["abc_trait"]["decision"] == "called"

    one_between = ([protein("SAMPLE1_00001", "AbcA protein"),
                    protein("SAMPLE1_00002", "AbcB protein")] + fillers(3, 3)
                   + [protein("SAMPLE1_00004", "AbcC protein")])
    spaced = run(tmp_path, one_between, ABC_GENES, rules,
                 extra_args=["--max-gene-gap", "1"])
    assert spaced["audit"]["abc_trait"]["decision"] == "partial"


# ── Evidence: names ──────────────────────────────────────────────────────────
# When a Bakta or Baktfold name points to a gene: product pattern, KO, EC
# number and symbol, and which proteins can carry a Baktfold name.

def test_ko_ec_and_symbol_each_recognise_a_gene(tmp_path):
    """Besides a product pattern, a KEGG ortholog, an EC number or the exact gene
    symbol recognises a gene, each read from Bakta's DbXrefs or Gene column; the
    detail names which part matched."""
    genes = [gene_entry("abcA", product_include="", ko="K99001"),
             gene_entry("abcB", product_include="", ec="9.9.9.9"),
             gene_entry("abcC", product_include="", bakta_gene="abcC")]
    proteins = [protein("SAMPLE1_00001", "enzyme", dbxrefs="KEGG:K99001, SO:0001217"),
                protein("SAMPLE1_00002", "enzyme", dbxrefs="EC:9.9.9.9"),
                protein("SAMPLE1_00003", "enzyme", gene="abcC")]
    result = run(tmp_path, proteins, genes, [abc_rule()])
    assert result["audit"]["abc_trait"]["decision"] == "called"
    assert "KO K99001" in rows_for(result, "SAMPLE1_00001", "abcA")[0]["detail"]
    assert "EC 9.9.9.9" in rows_for(result, "SAMPLE1_00002", "abcB")[0]["detail"]
    assert "symbol abcC" in rows_for(result, "SAMPLE1_00003", "abcC")[0]["detail"]


def test_symbol_that_needs_product_counts_only_with_product_or_ko(tmp_path):
    """Some symbols are given by Bakta to unrelated proteins too; for those the
    symbol counts only when the product or the KO agrees."""
    genes = [gene_entry("abcA", bakta_gene="abcA", symbol_needs_product="yes",
                        product_include="real AbcA", ko="K99001")]
    rules = [rule_entry("a_trait", "abcA")]

    alone = run(tmp_path, [protein("SAMPLE1_00001", "other enzyme", gene="abcA")],
                genes, rules)
    assert alone["audit"]["a_trait"]["decision"] == "not called"

    with_ko = run(tmp_path, [protein("SAMPLE1_00001", "other enzyme", gene="abcA",
                                     dbxrefs="KEGG:K99001")], genes, rules)
    assert with_ko["audit"]["a_trait"]["decision"] == "called"


def test_product_pattern_ignores_case(tmp_path):
    """Product patterns ignore case: a product can be written with different capitals
    in different annotation sources."""
    genes = [gene_entry("abcA", product_include="abca synthase")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "ABCA Synthase")], genes,
                 [rule_entry("a_trait", "abcA")])
    assert result["audit"]["a_trait"]["decision"] == "called"


def test_baktfold_name_counts_as_f_for_a_hypothetical_protein(tmp_path):
    """Baktfold names some proteins Bakta could only call "hypothetical protein";
    such a name is evidence F and counts like a Bakta name."""
    proteins = [protein("SAMPLE1_00001")]
    baktfold = [protein("SAMPLE1_00001", "AbcA protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA")],
                 [rule_entry("a_trait", "abcA")], baktfold=baktfold)

    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert row["evidence"] == "F"
    assert row["reason"] == "name_only"
    assert row["baktfold_product"] == "AbcA protein"
    assert row["bakta_product"] == "hypothetical protein"
    assert result["profile"]["a_trait"]["genes"] == "abcA:SAMPLE1_00001(F)"


def test_baktfold_name_of_a_pseudogene_is_shown_but_not_counted(tmp_path):
    """Baktfold renames the proteins Bakta flagged hypothetical, pseudogenes
    included, and writes Bakta's "(pseudo) " prefix too. Both products are read
    without it, so the Baktfold name gives F evidence (here through a pattern
    anchored with '^'), but the protein is a pseudogene: the row is not counted,
    and the one-gene trait is partial."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) hypothetical protein")]
    baktfold = [protein("SAMPLE1_00001", "(pseudo) AbcA protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA", product_include="^AbcA protein$")],
                 [rule_entry("a_trait", "abcA")], baktfold=baktfold)
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["evidence"], row["counted"], row["reason"]) == ("F", "no", "pseudogene")
    assert row["bakta_product"] == "(pseudo) hypothetical protein"
    assert row["baktfold_product"] == "(pseudo) AbcA protein"
    assert row["detail"] == ("Baktfold name matches by product; no structure hit; "
                             "Bakta marks the gene as a pseudogene (disrupted reading "
                             "frame), so it does not count")
    audit = result["audit"]["a_trait"]
    assert (audit["decision"], audit["reason"]) == ("partial",
                                                    "set_complete_but_pseudogene")


def test_baktfold_row_of_a_protein_bakta_named_is_ignored(tmp_path):
    """A protein Bakta named keeps Bakta's name only: a Baktfold row for it gives
    no F evidence and no baktfold_product."""
    proteins = [protein("SAMPLE1_00001", "unrelated enzyme")]
    baktfold = [protein("SAMPLE1_00001", "AbcA protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA")],
                 [rule_entry("a_trait", "abcA")], baktfold=baktfold)
    assert result["audit"]["a_trait"]["decision"] == "not called"
    assert result["genes"] == []


def test_baktfold_hypothetical_protein_is_no_name_whatever_its_prefix(tmp_path, capsys):
    """Baktfold writes Bakta's prefix too, so its product is read without it as
    well. A Baktfold product that is itself "hypothetical protein", prefix aside,
    is no name: the protein shows no baktfold_product, and the log does not count
    it among the Baktfold names. Its gene-table row comes from a KO, which Bakta
    gave to a pseudogene it could only call "hypothetical protein"."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) hypothetical protein",
                        dbxrefs="KEGG:K99001")]
    baktfold = [protein("SAMPLE1_00001", "(pseudo) hypothetical protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA", ko="K99001")],
                 [rule_entry("a", "abcA")], baktfold=baktfold)
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["evidence"], row["reason"]) == ("Q", "pseudogene")
    assert row["baktfold_product"] == "NA"
    assert "Baktfold names: 0." in capsys.readouterr().out


def test_one_protein_can_count_for_two_genes(tmp_path):
    """A fused protein carries two activities; it counts for both genes, since one
    name, Bakta's, names both, and the two are trivially co-localised."""
    genes = [gene_entry("abcA"), gene_entry("abcB")]
    rules = [rule_entry("ab_trait", "abcA+abcB", colocalised="all")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA-AbcB fusion protein")],
                 genes, rules)
    assert result["audit"]["ab_trait"]["decision"] == "called"
    assert result["profile"]["ab_trait"]["genes"] == ("abcA:SAMPLE1_00001(Q); "
                                                      "abcB:SAMPLE1_00001(Q)")


def test_several_proteins_of_one_gene_are_listed_together(tmp_path):
    """Two proteins of one gene are both listed in the genes column, joined by ',', so
    the reader sees every copy behind a call."""
    proteins = [protein("SAMPLE1_00001", "AbcA protein"),
                protein("SAMPLE1_00002", "AbcA protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA")], [rule_entry("a", "abcA")])
    assert result["profile"]["a"]["genes"] == "abcA:SAMPLE1_00001(Q),SAMPLE1_00002(Q)"


def test_named_short_orf_counts_as_a_protein(tmp_path):
    """sORF rows are proteins of Bakta's .faa just as CDS rows are, so a named sORF
    counts for its gene."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein", feature="sorf")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")])
    assert rows_for(result, "SAMPLE1_00001", "abcA")[0]["counted"] == "yes"
    assert result["audit"]["a"]["decision"] == "called"


# ── Evidence: structure ──────────────────────────────────────────────────────
# EcoFoldDB's filter on the best Foldseek hit, the s_alone column, the veto that
# cancels a structure hit, and how the hit table is read and checked.

def kept(query, target, evalue="1.000E-30", qcov="0.950", tcov="0.950"):
    """A Foldseek line; by default a hit that passes EcoFoldDB's filter easily."""
    return (query, target, evalue, qcov, tcov)


def test_structure_alone_counts_only_where_s_alone_is_yes(tmp_path):
    """A structure hit with no agreeing name counts only for genes whose s_alone
    column says yes; for the others it is reported and not counted."""
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]
    rules = [rule_entry("a_trait", "abcA")]

    trusted = run(tmp_path, [protein("SAMPLE1_00001")],
                  [gene_entry("abcA", s_alone="yes")], rules, hits=hits)
    row = rows_for(trusted, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["evidence"], row["reason"]) == (
        "yes", "S", "structure_counts_alone")
    assert trusted["audit"]["a_trait"]["decision"] == "called"
    assert trusted["profile"]["a_trait"]["structure_only"] == "abcA"

    needs_name = run(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                     rules, hits=hits)
    row = rows_for(needs_name, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["decision"], row["reason"]) == (
        "no", "not counted", "structure_only_needs_name")
    assert needs_name["audit"]["a_trait"]["decision"] == "not called"


def test_name_and_structure_agree(tmp_path):
    """Bakta's name and the kept structure hit point to the same gene: evidence S+Q,
    reason name_and_structure_agree. The row shows the hit's target, catalogue gene
    and pathway, and the e-value as Foldseek printed it."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert row["evidence"] == "S+Q"
    assert row["reason"] == "name_and_structure_agree"
    assert row["ecofolddb_target"] == "CAT_ABCA"
    assert row["ecofolddb_gene"] == "abcA"
    assert row["ecofolddb_pathway"] == "invented route"
    assert row["evalue"] == "1.000E-30"            # as Foldseek printed it
    assert row["borderline"] == "no"
    assert result["profile"]["a"]["genes"] == "abcA:SAMPLE1_00001(S+Q)"
    assert result["profile"]["a"]["structure_only"] == "NA"


def test_veto_cancels_structure_evidence_too(tmp_path):
    """product_exclude names proteins that look like the gene but are something
    else. Its veto removes every kind of evidence on that protein, so a structure
    hit cannot bring back a protein the name has ruled out, even where structure
    alone would count."""
    genes = [gene_entry("abcA", s_alone="yes", product_exclude="decoy")]
    rules = [rule_entry("a_trait", "abcA")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    result = run(tmp_path, [protein("SAMPLE1_00001", "Decoy enzyme")], genes, rules,
                 hits=hits)
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["evidence"], row["reason"]) == (
        "no", "S", "vetoed_by_name")
    assert "Bakta product 'Decoy enzyme'" in row["detail"]
    assert result["audit"]["a_trait"]["decision"] == "not called"


def test_baktfold_product_can_veto(tmp_path):
    """A Baktfold product vetoes a gene as a Bakta product does, here cancelling a
    structure hit that would count alone (s_alone yes)."""
    genes = [gene_entry("abcA", s_alone="yes", product_exclude="decoy")]
    result = run(tmp_path, [protein("SAMPLE1_00001")], genes,
                 [rule_entry("a_trait", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")],
                 baktfold=[protein("SAMPLE1_00001", "decoy enzyme")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert row["reason"] == "vetoed_by_name"
    assert "Baktfold product 'decoy enzyme'" in row["detail"]
    assert result["audit"]["a_trait"]["decision"] == "not called"


def test_baktfold_product_vetoes_a_bakta_ko_or_symbol_match(tmp_path):
    """The Baktfold product's veto cancels Bakta's own evidence on that protein
    too: here a KO on one protein and a symbol on the other, both of which Bakta
    gave to proteins it could only call "hypothetical protein"."""
    genes = [gene_entry("abcA", ko="K99001", bakta_gene="abcA", product_exclude="decoy")]
    proteins = [protein("SAMPLE1_00001", dbxrefs="KEGG:K99001"),
                protein("SAMPLE1_00002", gene="abcA")]
    baktfold = [protein("SAMPLE1_00001", "decoy enzyme"),
                protein("SAMPLE1_00002", "decoy enzyme")]
    result = run(tmp_path, proteins, genes, [rule_entry("a", "abcA")], baktfold=baktfold)
    for locus in ("SAMPLE1_00001", "SAMPLE1_00002"):
        row = rows_for(result, locus, "abcA")[0]
        assert (row["evidence"], row["reason"]) == ("Q", "vetoed_by_name"), locus
        assert "Baktfold product 'decoy enzyme'" in row["detail"], locus
    assert result["audit"]["a"]["decision"] == "not called"


def test_baktfold_veto_reads_the_product_without_the_prefix(tmp_path):
    """The Baktfold product is tested by product_exclude without the prefix
    Baktfold copies from Bakta, as the Bakta product is: a pattern anchored with
    '^' vetoes a gene Bakta marks as running off its sequence (3' truncated),
    structure hit included, although that hit alone would count (s_alone yes).
    Without the veto, this decoy would call the trait."""
    genes = [gene_entry("abcA", s_alone="yes", product_exclude="^decoy")]
    proteins = [protein("SAMPLE1_00001", "(3' truncated) hypothetical protein")]
    baktfold = [protein("SAMPLE1_00001", "(3' truncated) Decoy enzyme")]
    result = run(tmp_path, proteins, genes, [rule_entry("a", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")], baktfold=baktfold)
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["evidence"], row["reason"]) == (
        "no", "S", "vetoed_by_name")
    assert row["detail"] == ("product_exclude of abcA matches the Baktfold product "
                             "'Decoy enzyme', so the S evidence does not count; "
                             "Bakta: runs off its sequence at the gene's 3' end "
                             "(3' truncated)")
    assert row["baktfold_product"] == "(3' truncated) Decoy enzyme"
    assert result["audit"]["a"]["decision"] == "not called"


def test_a_veto_with_nothing_to_cancel_writes_no_row(tmp_path):
    """An exclude pattern that matches a protein with no evidence for the gene writes
    no row: a veto row appears only when it cancelled evidence, so unrelated
    proteins do not crowd the gene table."""
    genes = [gene_entry("abcA", product_exclude="decoy")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "Decoy enzyme")], genes,
                 [rule_entry("a", "abcA")])
    assert result["genes"] == []


def test_filter_is_strict_at_its_cut_offs(tmp_path):
    """EcoFoldDB's filter is e-value below 1e-10 and coverage above 0.8, both
    strict: a hit at exactly 1e-10, or at exactly 0.8 query or target coverage,
    does not count."""
    genes = [gene_entry("abcA", s_alone="yes")]
    rules = [rule_entry("a", "abcA")]
    proteins = [protein("SAMPLE1_00001")]

    at_evalue = run(tmp_path, proteins, genes, rules,
                    hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="1.000E-10")])
    assert at_evalue["audit"]["a"]["decision"] == "not called"
    at_coverage = run(tmp_path, proteins, genes, rules,
                      hits=[kept("SAMPLE1_00001", "CAT_ABCA", qcov="0.800")])
    assert at_coverage["audit"]["a"]["decision"] == "not called"
    at_target_coverage = run(tmp_path, proteins, genes, rules,
                             hits=[kept("SAMPLE1_00001", "CAT_ABCA", tcov="0.800")])
    assert at_target_coverage["audit"]["a"]["decision"] == "not called"
    inside = run(tmp_path, proteins, genes, rules,
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="9.900E-11",
                            qcov="0.801", tcov="0.801")])
    assert inside["audit"]["a"]["decision"] == "called"


def test_best_hit_outside_the_catalogue_blocks_structure_evidence(tmp_path):
    """The filter takes the single best hit across the whole database. When that
    is a background structure, not a catalogue entry, the protein has no structure
    evidence, even if a weaker hit to a catalogue entry passes the thresholds."""
    genes = [gene_entry("abcA", s_alone="yes")]
    hits = [kept("SAMPLE1_00001", "BACKGROUND_1", evalue="1.000E-40"),
            kept("SAMPLE1_00001", "CAT_ABCA", evalue="1.000E-30")]
    result = run(tmp_path, [protein("SAMPLE1_00001")], genes, [rule_entry("a", "abcA")],
                 hits=hits)
    assert result["audit"]["a"]["decision"] == "not called"
    assert result["genes"] == []


def test_name_row_shows_a_failed_structure_hit_and_why(tmp_path):
    """A protein counted from its name whose best structure hit failed EcoFoldDB's
    filter: the row still shows the hit, and the detail says which thresholds it
    failed. borderline is NA, since only kept hits are flagged."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="3.000E-08",
                            qcov="0.500")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert row["reason"] == "name_only"
    assert row["counted"] == "yes"
    assert row["ecofolddb_target"] == "CAT_ABCA"
    assert row["evalue"] == "3.000E-08"
    assert row["borderline"] == "NA"
    assert row["detail"].endswith("best structure hit failed the filter: e-value at "
                                  "or above 1e-10 and query coverage at or below 0.8")


def test_name_row_says_when_the_structure_points_to_another_gene(tmp_path):
    """Bakta's name points to abcA and the kept structure hit to abcB. The abcA row
    counts from the name, and its detail says where the structure points; the abcB
    row, structure alone with s_alone no, is not counted. The two genes share no
    rule, and the abcB row does not count anyway, so the two are not compared (see
    the one-gene-per-rule tests below)."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA"), gene_entry("abcB")],
                 [rule_entry("a", "abcA"), rule_entry("b", "abcB")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCB")])
    a_row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert a_row["reason"] == "name_only"
    assert a_row["detail"].endswith("best structure hit abcB@invented route is not "
                                    "in the ecofolddb list of abcA")
    b_row = rows_for(result, "SAMPLE1_00001", "abcB")[0]
    assert b_row["reason"] == "structure_only_needs_name"


def test_best_hit_tie_goes_to_the_first_listed(tmp_path):
    """Two hits with the same e-value: the first one in the table is the best hit,
    as in EcoFoldDB's own filter, so the same table always gives the same answer."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB", s_alone="yes")]
    rules = [rule_entry("a", "abcA"), rule_entry("b", "abcB")]
    a_first = [kept("SAMPLE1_00001", "CAT_ABCA"), kept("SAMPLE1_00001", "CAT_ABCB")]

    result = run(tmp_path, [protein("SAMPLE1_00001")], genes, rules, hits=a_first)
    assert result["audit"]["a"]["decision"] == "called"
    assert result["audit"]["b"]["decision"] == "not called"

    result = run(tmp_path, [protein("SAMPLE1_00001")], genes, rules,
                 hits=list(reversed(a_first)))
    assert result["audit"]["a"]["decision"] == "not called"
    assert result["audit"]["b"]["decision"] == "called"


def test_kept_hit_to_an_unused_catalogue_gene_gets_its_own_row(tmp_path):
    """A kept hit to a catalogue entry that no dictionary gene uses gets a row of its
    own (gene NA, catalogue_gene_unused), so no kept structure hit goes unreported.
    The detail names the entry's key."""
    result = run(tmp_path, [protein("SAMPLE1_00001")], ABC_GENES, [abc_rule()],
                 hits=[kept("SAMPLE1_00001", "CAT_XYZA")])
    row = rows_for(result, "SAMPLE1_00001")[0]
    assert row["gene"] == "NA"
    assert row["reason"] == "catalogue_gene_unused"
    assert row["decision"] == "not counted"
    assert row["ecofolddb_gene"] == "xyzA"
    assert row["detail"] == ("kept structure hit xyzA@unused route, which no gene of "
                             "the gene dictionary uses")


def test_protein_too_long_for_the_structure_search_is_listed(tmp_path):
    """Proteins over the length limit never reach ProstT5, so the absence of a
    structure hit means nothing for them; one row says so, and their names still
    count."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 lengths={"SAMPLE1_00001": 4500})
    rows = rows_for(result, "SAMPLE1_00001")
    assert [row["reason"] for row in rows] == ["name_only", "not_searched_too_long"]
    assert rows[0]["detail"].endswith("not searched for structure: 4500 aa, longer "
                                      "than the 4000 aa limit")
    assert rows[1]["decision"] == "not searched"
    assert rows[1]["gene"] == "NA"
    assert rows[1]["detail"].startswith("4500 aa, longer than the 4000 aa limit")
    assert result["audit"]["a"]["decision"] == "called"


def test_protein_exactly_at_the_length_limit_was_searched(tmp_path):
    """Rule trait_structure_query keeps proteins of at most the limit, so one of
    exactly 4000 aa was searched: no not_searched_too_long row, and its name row
    says it simply had no structure hit."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 lengths={"SAMPLE1_00001": 4000})
    rows = rows_for(result, "SAMPLE1_00001")
    assert [row["reason"] for row in rows] == ["name_only"]
    assert rows[0]["detail"].endswith("; no structure hit")


def test_protein_lengths_are_read_from_a_wrapped_fasta(tmp_path):
    """A .faa may wrap a sequence over several lines; the length is the sum of the
    line lengths, as rule trait_structure_query counts it."""
    path = tmp_path / "SAMPLE1.faa"
    path.write_text(">SAMPLE1_00001 AbcA protein\nMKVLA\nAGTW\nM\n"
                    ">SAMPLE1_00002 hypothetical protein\nMKV\n", encoding="utf-8")
    assert tp.read_protein_lengths(str(path)) == {"SAMPLE1_00001": 10,
                                                  "SAMPLE1_00002": 3}


def test_gzipped_hit_table_is_read(tmp_path):
    """A gzip-compressed hit table, named .gz, is read like a plain one."""
    result = run(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA", s_alone="yes")],
                 [rule_entry("a", "abcA")], hits=[kept("SAMPLE1_00001", "CAT_ABCA")],
                 hits_name="hits.tsv.gz")
    assert result["audit"]["a"]["decision"] == "called"


def test_hits_of_another_genome_stop_the_run(tmp_path):
    """If no query of the hit table is a protein of this genome, the table belongs
    to another sample; carrying on would call every trait from names alone."""
    args = write_inputs(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                        [rule_entry("a", "abcA")],
                        hits=[kept("OTHER_00001", "CAT_ABCA")])
    with pytest.raises(ValueError, match="hits belong to another genome"):
        tp.main(args)


def test_empty_hit_table_warns_that_names_decide_alone(tmp_path, capsys):
    """An empty hit table is accepted, so a hand-run without the structure search
    still works. For a genome with proteins short enough to be searched it most
    likely means the search failed, so a warning says so; the names still count.
    A genome whose proteins are all too long for the search gets no warning."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")])
    assert result["audit"]["a"]["decision"] == "called"
    warning = capsys.readouterr().err
    assert warning.startswith("WARNING: the hit table ")
    assert ("holds no hits, although 1 protein(s) are short enough for the "
            "structure search (at most 4000 aa). Every trait is decided from names "
            "alone.") in warning

    run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")], [gene_entry("abcA")],
        [rule_entry("a", "abcA")], lengths={"SAMPLE1_00001": 4500})
    assert "holds no hits" not in capsys.readouterr().err


def test_evalue_that_is_not_a_finite_number_stops_the_run(tmp_path):
    """Python reads 'nan' as a number, but it would pass the e-value filter
    unnoticed: no comparison with it is ever true."""
    args = write_inputs(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                        [rule_entry("a", "abcA")],
                        hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="nan")])
    with pytest.raises(ValueError, match="is not a finite number"):
        tp.main(args)


@pytest.mark.parametrize("qcov, tcov", [("nan", "0.950"), ("0.950", "inf")])
def test_coverage_that_is_not_a_finite_number_stops_the_run(tmp_path, qcov, tcov):
    """The same for the coverage of a best hit: 'nan' would pass the coverage
    filter unnoticed, since it is never "at or below 0.8", and 'inf' is above
    every cut-off."""
    args = write_inputs(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                        [rule_entry("a", "abcA")],
                        hits=[kept("SAMPLE1_00001", "CAT_ABCA", qcov=qcov, tcov=tcov)])
    with pytest.raises(ValueError, match="is not a finite number"):
        tp.main(args)


# ── Borderline flags ─────────────────────────────────────────────────────────
# A kept structure hit close to a cut-off or to a rival is flagged for a second
# look; the flag never changes a decision.

def borderline_of(tmp_path, hits):
    """Run one protein with a structure hit for abcA; return its borderline cell
    and the whole result."""
    result = run(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA", s_alone="yes")],
                 [rule_entry("a", "abcA")], hits=hits)
    return rows_for(result, "SAMPLE1_00001", "abcA")[0]["borderline"], result


def test_borderline_near_the_evalue_cut_off(tmp_path):
    """A kept hit at e-value 1e-11, within a factor of 10 of the 1e-10 cut-off, is
    flagged, and the profile's borderline column says yes; one at 9.9e-12 is not
    flagged."""
    flag, result = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA",
                                                 evalue="1.000E-11")])
    assert flag == "near e-value cut-off"
    assert result["profile"]["a"]["borderline"] == "yes"
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA",
                                            evalue="9.900E-12")])
    assert flag == "no"


def test_borderline_near_the_coverage_cut_off(tmp_path):
    """A coverage at or below 0.85 is flagged, here the target coverage at 0.850;
    0.851 for both coverages is not."""
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", tcov="0.850")])
    assert flag == "near coverage cut-off"
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", qcov="0.851",
                                            tcov="0.851")])
    assert flag == "no"


def test_borderline_near_tie_with_another_gene(tmp_path):
    """A second-best hit to a different catalogue gene within a factor of 10 means
    the structure barely tells the two apart."""
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", "1.000E-30"),
                                       kept("SAMPLE1_00001", "CAT_ABCB", "5.000E-30")])
    assert flag == "near tie with abcB"


def test_borderline_near_tie_at_exactly_ten_times(tmp_path):
    """Foldseek prints e-values to four digits, so a rival at exactly ten times the
    best e-value happens, and it is a near tie. In floating point 1.000E-23 x 10
    falls just below 1.000E-22, hence the boundary test; just above ten times is
    not a tie."""
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", "1.000E-23"),
                                       kept("SAMPLE1_00001", "CAT_ABCB", "1.000E-22")])
    assert flag == "near tie with abcB"
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", "1.000E-23"),
                                       kept("SAMPLE1_00001", "CAT_ABCB", "1.001E-22")])
    assert flag == "no"


def test_borderline_near_tie_with_a_target_outside_the_catalogue(tmp_path):
    """A background structure almost as close as the best catalogue hit is a near tie
    too: the protein may belong to a family the catalogue does not hold."""
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", "1.000E-30"),
                                       kept("SAMPLE1_00001", "BACKGROUND_1", "2.000E-30")])
    assert flag == "near tie with a target outside the catalogue"


def test_no_tie_with_another_structure_of_the_same_gene(tmp_path):
    """A second structure of the same catalogue gene at the same e-value is no rival:
    both point to abcA. The rival, abcB at 1e-20 against 1e-30, is far behind."""
    flag, _ = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA", "1.000E-30"),
                                       kept("SAMPLE1_00001", "CAT_ABCA2", "1.000E-30"),
                                       kept("SAMPLE1_00001", "CAT_ABCB", "1.000E-20")])
    assert flag == "no"


def test_rival_is_found_whatever_the_order_of_the_hits(tmp_path):
    """The best rival must not depend on where the best hit sits in the table. Here
    the best hit (abcA, 1e-30) comes last, after two weaker hits to abcB (1e-20,
    then 4e-30) with a background hit between them; the rival is the closer abcB
    hit, 4e-30."""
    catalogue = {target: {"gene": gene, "pathway": pathway, "key": gene + "@" + pathway}
                 for target, gene, pathway in CATALOGUE}
    path = write_hits(tmp_path / "hits.tsv", [
        kept("P1", "CAT_ABCB", "1.000E-20"),
        kept("P1", "BACKGROUND_1", "1.000E-10"),
        kept("P1", "CAT_ABCB", "4.000E-30"),
        kept("P1", "CAT_ABCA", "1.000E-30"),
    ])
    best, rival = tp.read_best_hits(path, catalogue)
    assert best["P1"].target == "CAT_ABCA"
    assert rival["P1"].target == "CAT_ABCB"
    assert rival["P1"].evalue_text == "4.000E-30"


def test_borderline_never_changes_a_decision(tmp_path):
    """A hit flagged on two counts still makes the call: the flag asks for a second
    look and changes no decision."""
    flag, result = borderline_of(tmp_path, [kept("SAMPLE1_00001", "CAT_ABCA",
                                                 evalue="5.000E-11", qcov="0.810")])
    assert flag == "near e-value cut-off; near coverage cut-off"
    assert result["audit"]["a"]["decision"] == "called"


def test_borderline_is_reported_on_a_partial_call(tmp_path):
    """A partial call carries the borderline column too: here its one counted
    gene rests on a structure hit near the e-value cut-off."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001")], genes,
                 [rule_entry("r", "abcA+abcB")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="5.000E-11")])
    assert result["audit"]["r"]["decision"] == "partial"
    assert result["profile"]["r"]["borderline"] == "yes"


# ── Bakta's marks: pseudogenes and genes that run off their sequence ─────────
# Bakta writes "(pseudo) " before the product of a pseudogene (aligned to its
# reference protein, the gene is broken by a frameshift, an internal stop codon or
# a lost start codon) and "(5' truncated) ", "(3' truncated) " or "(partial) "
# before that of a gene that runs off an end of its sequence, which with
# BacFlux's settings happens only on a circular replicon. Names are read without
# the prefix; a truncated gene counts as usual; a pseudogene never counts.

def test_pseudogene_is_listed_but_never_counted(tmp_path):
    """A pseudogene's evidence is shown, but it does not count: the row reads not
    counted, reason pseudogene, and the detail says what matched and why it does
    not count. The product keeps Bakta's prefix in the bakta_product column."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) AbcA protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["decision"], row["reason"]) == (
        "no", "not counted", "pseudogene")
    assert row["evidence"] == "S+Q"
    assert row["bakta_product"] == "(pseudo) AbcA protein"
    assert row["detail"] == ("Bakta name matches by product; structure hit "
                             "abcA@invented route; Bakta marks the gene as a "
                             "pseudogene (disrupted reading frame), so it does not "
                             "count")


def test_one_gene_set_whose_gene_is_a_pseudogene_is_partial(tmp_path):
    """A set of one gene has no half, but its gene can be in the genome and broken:
    the trait is partial, and the audit and the missing column name the broken
    gene with its locus tag. Nothing counts, so genes stays NA, and borderline and
    structure_only come from counted proteins only. A later set that is called
    still wins over this partial."""
    genes = [gene_entry("abcA"), gene_entry("abcB")]
    rules = [rule_entry("single", "abcA | abcB")]
    broken = protein("SAMPLE1_00001", "(pseudo) AbcA protein")
    result = run(tmp_path, [broken], genes, rules,
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA", evalue="5.000E-11")])
    audit = result["audit"]["single"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "set_complete_but_pseudogene"
    assert audit["detail"] == "abcA complete except abcA:SAMPLE1_00001(pseudogene)"
    profile = result["profile"]["single"]
    assert profile["genes"] == "NA"
    assert profile["missing"] == "abcA:SAMPLE1_00001(pseudogene)"
    assert profile["borderline"] == "no"
    assert profile["structure_only"] == "NA"

    with_abcb = run(tmp_path, [broken, protein("SAMPLE1_00002", "AbcB protein")],
                    genes, rules)
    assert with_abcb["audit"]["single"]["decision"] == "called"
    assert with_abcb["audit"]["single"]["detail"] == "abcB complete"


def test_set_with_two_broken_genes_is_partial_naming_both(tmp_path):
    """Two genes of a set found only as pseudogenes, the third working: the set is
    partial (set_complete_but_pseudogene), and the audit and the missing column name
    both broken genes with their locus tags."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                protein("SAMPLE1_00002", "(pseudo) AbcB protein"),
                protein("SAMPLE1_00003", "AbcC protein")]
    result = run(tmp_path, proteins, ABC_GENES, [abc_rule(colocalised="all")])
    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "set_complete_but_pseudogene"
    assert audit["detail"] == ("abcA+abcB+abcC complete except "
                               "abcA:SAMPLE1_00001(pseudogene), "
                               "abcB:SAMPLE1_00002(pseudogene)")
    profile = result["profile"]["abc_trait"]
    assert profile["genes"] == "abcC:SAMPLE1_00003(Q)"
    assert profile["missing"] == ("abcA:SAMPLE1_00001(pseudogene), "
                                  "abcB:SAMPLE1_00002(pseudogene)")


def test_a_broken_and_an_absent_gene_follow_the_usual_rules(tmp_path):
    """The pseudogene partial needs every gene that is not counted to be a
    pseudogene. With one gene broken and another absent, the usual rules apply,
    with the broken gene not counted: one of three genes is too few, two of four
    is half. The details still name the broken gene with its locus tag."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                protein("SAMPLE1_00003", "AbcC protein")]
    three = run(tmp_path, proteins, ABC_GENES, [abc_rule()])
    audit = three["audit"]["abc_trait"]
    assert audit["decision"] == "not called"
    assert audit["reason"] == "too_few_genes"
    assert audit["detail"] == ("1 of 3 genes found in the closest set abcA+abcB+abcC "
                               "(missing abcA:SAMPLE1_00001(pseudogene), abcB)")

    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    four = run(tmp_path, proteins + [protein("SAMPLE1_00004", "AbcD protein")], genes,
               [rule_entry("r", "abcA+abcB+abcC+abcD")], catalogue=CATALOGUE_SIX)
    audit = four["audit"]["r"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "half_of_set_found"
    assert audit["detail"] == ("missing abcA:SAMPLE1_00001(pseudogene), abcB (2 of 4 "
                               "genes of abcA+abcB+abcC+abcD found)")


def test_key_gene_found_as_a_pseudogene_is_enough(tmp_path):
    """key_genes asks whether a key gene is in the genome, and a broken copy is: a
    cluster whose committing gene abcA is a pseudogene is reported as partial, not
    withheld, since that broken gene is what the partial shows. This holds for a
    half-found set too."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                protein("SAMPLE1_00002", "AbcB protein"),
                protein("SAMPLE1_00003", "AbcC protein")]
    result = run(tmp_path, proteins, ABC_GENES, [abc_rule(key_genes="abcA")])
    audit = result["audit"]["abc_trait"]
    assert audit["decision"] == "partial"
    assert audit["reason"] == "set_complete_but_pseudogene"
    assert audit["detail"] == ("abcA+abcB+abcC complete except "
                               "abcA:SAMPLE1_00001(pseudogene)")

    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    half = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                          protein("SAMPLE1_00003", "AbcC protein"),
                          protein("SAMPLE1_00004", "AbcD protein")],
               genes, [rule_entry("r", "abcA+abcB+abcC+abcD", key_genes="abcA")],
               catalogue=CATALOGUE_SIX)
    assert half["audit"]["r"]["decision"] == "partial"
    assert half["audit"]["r"]["reason"] == "half_of_set_found"


def test_key_genes_never_withhold_a_set_complete_except_for_pseudogenes(tmp_path):
    """key_genes holds back partials that a related pathway could produce from a
    few shared genes. A set complete except for pseudogenes is not one of them:
    all its genes are in the genome. Here the one-gene route abcA is broken and the
    key gene abcB belongs to the other route; the trait reads partial and names
    the broken gene. A half-found set without the key gene is still withheld, and
    it does not block a later set that is complete except for a pseudogene."""
    genes = [gene_entry("abcA"), gene_entry("abcB"), gene_entry("abcC")]
    rules = [rule_entry("r", "abcA | abcB+abcC", key_genes="abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) AbcA protein")], genes,
                 rules)
    audit = result["audit"]["r"]
    assert (audit["decision"], audit["reason"]) == ("partial",
                                                    "set_complete_but_pseudogene")
    assert audit["detail"] == "abcA complete except abcA:SAMPLE1_00001(pseudogene)"
    assert result["profile"]["r"]["missing"] == "abcA:SAMPLE1_00001(pseudogene)"

    four = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    later = [rule_entry("r", "abcA+abcB+abcC | abcD", key_genes="abcA")]
    without_key = [protein("SAMPLE1_00002", "AbcB protein"),
                   protein("SAMPLE1_00003", "AbcC protein")]
    withheld = run(tmp_path, without_key, four, later, catalogue=CATALOGUE_SIX)
    assert withheld["audit"]["r"]["reason"] == "partial_withheld"

    broken_abcd = protein("SAMPLE1_00004", "(pseudo) AbcD protein")
    broken_later = run(tmp_path, without_key + [broken_abcd], four, later,
                       catalogue=CATALOGUE_SIX)
    audit = broken_later["audit"]["r"]
    assert (audit["decision"], audit["reason"]) == ("partial",
                                                    "set_complete_but_pseudogene")
    assert audit["detail"] == "abcD complete except abcD:SAMPLE1_00004(pseudogene)"


def test_no_gene_found_names_a_gene_found_only_as_a_pseudogene(tmp_path):
    """When no gene of any set counts, a pseudogene of a required gene is the one
    trace of the trait: the audit names it, though the trait is not called and the
    profile lists nothing. A genome with no trace at all keeps the plain detail."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) AbcA protein")],
                 ABC_GENES, [abc_rule()])
    audit = result["audit"]["abc_trait"]
    assert (audit["decision"], audit["reason"]) == ("not called", "no_gene_found")
    assert audit["detail"] == ("no gene of any required set counted; found only as a "
                               "pseudogene: abcA:SAMPLE1_00001(pseudogene)")
    profile = result["profile"]["abc_trait"]
    assert (profile["genes"], profile["missing"]) == ("NA", "NA")

    nothing = run(tmp_path, [protein("SAMPLE1_00001", "unrelated enzyme")], ABC_GENES,
                  [abc_rule()])
    plain = nothing["audit"]["abc_trait"]
    assert plain["detail"] == "no gene of any required set found"


def test_partial_detail_names_a_broken_gene_of_another_set(tmp_path):
    """Sets are tried in table order, so the half-found abcA+abcB is the partial,
    although abcC+abcD lacks only a pseudogene. The detail describes the reported
    set and then names the broken gene of the other one, so it shows in the audit
    and not only in the gene table. The missing column keeps to the reported set."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    proteins = [protein("SAMPLE1_00001", "AbcA protein"),
                protein("SAMPLE1_00003", "AbcC protein"),
                protein("SAMPLE1_00004", "(pseudo) AbcD protein")]
    result = run(tmp_path, proteins, genes, [rule_entry("r", "abcA+abcB | abcC+abcD")],
                 catalogue=CATALOGUE_SIX)
    audit = result["audit"]["r"]
    assert (audit["decision"], audit["reason"]) == ("partial", "half_of_set_found")
    assert audit["detail"] == ("missing abcB (1 of 2 genes of abcA+abcB found); "
                               "outside this set, found only as a pseudogene: "
                               "abcD:SAMPLE1_00004(pseudogene)")
    profile = result["profile"]["r"]
    assert (profile["genes"], profile["missing"]) == ("abcA:SAMPLE1_00001(Q)", "abcB")


def test_not_called_detail_names_a_broken_gene_of_another_set(tmp_path):
    """The same note ends the detail of a trait that is not called: after the
    closest set, and inside a partial withheld for want of a key gene. A gene with
    a working copy elsewhere counts, so it is not named."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC",
                                           "abcD", "abcE", "abcF")]
    rules = [rule_entry("r", "abcA+abcB+abcC | abcD+abcE+abcF")]
    broken_abcd = protein("SAMPLE1_00004", "(pseudo) AbcD protein")
    few = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein"), broken_abcd],
              genes, rules, catalogue=CATALOGUE_SIX)
    audit = few["audit"]["r"]
    assert audit["reason"] == "too_few_genes"
    assert audit["detail"] == ("1 of 3 genes found in the closest set abcA+abcB+abcC "
                               "(missing abcB, abcC); outside this set, found only "
                               "as a pseudogene: abcD:SAMPLE1_00004(pseudogene)")

    with_copy = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein"), broken_abcd,
                               protein("SAMPLE1_00009", "AbcD protein")],
                    genes, rules, catalogue=CATALOGUE_SIX)
    assert with_copy["audit"]["r"]["detail"] == (
        "1 of 3 genes found in the closest set abcA+abcB+abcC (missing abcB, abcC)")

    needs_a = [rule_entry("r", "abcA+abcB+abcC | abcD+abcE+abcF",
                          key_genes="abcA")]
    withheld = run(tmp_path, [protein("SAMPLE1_00002", "AbcB protein"),
                              protein("SAMPLE1_00003", "AbcC protein"), broken_abcd],
                   genes, needs_a, catalogue=CATALOGUE_SIX)
    audit = withheld["audit"]["r"]
    assert audit["reason"] == "partial_withheld"
    assert audit["detail"] == ("no key gene (abcA) found (missing abcA (2 of 3 genes "
                               "of abcA+abcB+abcC found); outside this set, found "
                               "only as a pseudogene: abcD:SAMPLE1_00004(pseudogene))")


def test_a_working_copy_elsewhere_still_gives_the_call(tmp_path):
    """A broken abcA beside a working one: the working copy carries the call, and
    the pseudogene keeps its own row in the gene table."""
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                protein("SAMPLE1_00009", "AbcA protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA")], [rule_entry("a", "abcA")])
    assert result["audit"]["a"]["decision"] == "called"
    assert result["audit"]["a"]["reason"] == "set_complete"
    assert result["profile"]["a"]["genes"] == "abcA:SAMPLE1_00009(Q)"
    assert result["profile"]["a"]["missing"] == "NA"
    assert rows_for(result, "SAMPLE1_00001", "abcA")[0]["reason"] == "pseudogene"


def test_anchored_pattern_finds_a_pseudogene_and_lists_it(tmp_path):
    """A pattern anchored with '^' is tested on the product without Bakta's
    prefix, so it finds a pseudogene as an unanchored pattern or a KO would; the
    pseudogene is then listed and not counted."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) AbcA protein")],
                 [gene_entry("abcA", product_include="^AbcA protein$")],
                 [rule_entry("a", "abcA")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["reason"]) == ("no", "pseudogene")


@pytest.mark.parametrize("prefix, note", [
    ("(5' truncated) ", "Bakta: runs off its sequence at the gene's 5' end "
                        "(5' truncated)"),
    ("(3' truncated) ", "Bakta: runs off its sequence at the gene's 3' end "
                        "(3' truncated)"),
    ("(partial) ", "Bakta: runs off both ends of its sequence (partial)"),
])
def test_gene_running_off_its_sequence_counts_with_a_note(tmp_path, prefix, note):
    """A gene Bakta marks as running off a sequence end counts as usual, found here
    by a pattern anchored with ^, and its detail says at which end of the gene.
    With BacFlux's settings the mark appears only on a circular replicon, where
    the gene is not really cut."""
    result = run(tmp_path, [protein("SAMPLE1_00001", prefix + "AbcA protein")],
                 [gene_entry("abcA", product_include="^AbcA protein$")],
                 [rule_entry("a", "abcA")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (row["counted"], row["reason"]) == ("yes", "name_only")
    assert row["detail"] == "Bakta name matches by product; no structure hit; " + note
    assert row["bakta_product"] == prefix + "AbcA protein"
    assert result["audit"]["a"]["decision"] == "called"


def test_veto_reads_the_product_without_bakta_prefix_and_comes_first(tmp_path):
    """product_exclude is tested on the product without Bakta's prefix too, so an
    anchored exclude pattern vetoes a marked protein. The veto comes before the
    pseudogene rule: a vetoed pseudogene reads vetoed_by_name, with its mark noted.
    The detail quotes the product the pattern was tested on."""
    genes = [gene_entry("abcA", s_alone="yes", product_exclude="^decoy")]
    proteins = [protein("SAMPLE1_00001", "(pseudo) Decoy enzyme"),
                protein("SAMPLE1_00002", "(3' truncated) Decoy enzyme")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA"), kept("SAMPLE1_00002", "CAT_ABCA")]
    result = run(tmp_path, proteins, genes, [rule_entry("a", "abcA")], hits=hits)
    pseudo_row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert pseudo_row["reason"] == "vetoed_by_name"
    assert pseudo_row["detail"] == ("product_exclude of abcA matches the Bakta "
                                    "product 'Decoy enzyme', so the S evidence does "
                                    "not count; Bakta marks the gene as a pseudogene "
                                    "(disrupted reading frame)")
    assert rows_for(result, "SAMPLE1_00002", "abcA")[0]["reason"] == "vetoed_by_name"
    assert result["audit"]["a"]["decision"] == "not called"


def test_pseudogene_that_would_not_count_anyway_keeps_its_reason(tmp_path):
    """A pseudogene with only a structure hit, for a gene that needs a name to
    agree, would not count anyway: it keeps reason structure_only_needs_name with
    its mark noted, and it does not make its set a pseudogene partial."""
    result = run(tmp_path, [protein("SAMPLE1_00001", "(pseudo) hypothetical protein")],
                 [gene_entry("abcA")], [rule_entry("a", "abcA")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")])
    row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert row["reason"] == "structure_only_needs_name"
    assert row["detail"].endswith("; Bakta marks the gene as a pseudogene (disrupted "
                                  "reading frame)")
    assert result["audit"]["a"]["reason"] == "no_gene_found"


def test_missing_column_writes_broken_genes_with_their_loci(tmp_path):
    """The missing column lists the genes of the partial set that are not counted:
    a gene found only as pseudogenes as gene:locus_tag(pseudogene), its proteins
    joined by ',' as in the genes column, and a gene not found as its bare name;
    genes joined by ', '."""
    genes = [gene_entry(gene) for gene in ("abcA", "abcB", "abcC", "abcD")]
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA protein"),
                protein("SAMPLE1_00003", "AbcC protein"),
                protein("SAMPLE1_00004", "AbcD protein"),
                protein("SAMPLE1_00005", "(pseudo) AbcA protein")]
    result = run(tmp_path, proteins, genes, [rule_entry("r", "abcA+abcB+abcC+abcD")],
                 catalogue=CATALOGUE_SIX)
    assert result["audit"]["r"]["reason"] == "half_of_set_found"
    assert result["profile"]["r"]["missing"] == (
        "abcA:SAMPLE1_00001(pseudogene),SAMPLE1_00005(pseudogene), abcB")


# ── One protein, one gene per rule ───────────────────────────────────────────
# A protein counts for two genes of one rule only when one source names both: one
# Bakta name (product, symbol, KO and EC together), one Baktfold name or one
# structure hit. When the two genes come from different sources, it counts once:
# for the structure's gene when a structure hit counts alone for it (s_alone yes),
# otherwise for the gene of Bakta's name, then that of Baktfold's. A pseudogene is
# compared in the same way, and the trait calls take it as the gene that ranks
# first. Several tests run the same protein twice, once keeping both genes (or
# with no clash) and once in a clash, so they show both sides of the rule.

def test_name_and_structure_for_two_genes_of_one_set_count_once(tmp_path):
    """The structure hit points to abcA, which counts from structure alone (s_alone
    yes), and Bakta's name to abcB, another gene of the same set: two genes from
    two sources, as when a protein Bakta names as the beta chain of a two-chain
    enzyme has the structure of the alpha chain. The protein counts once, for the
    structure's gene. The abcB row stays in the gene table, not counted, with
    reason counted_for_other_gene, and the protein does not count for abcB in
    abcB's own trait either. A second protein named AbcB completes the set, and
    the trait is called."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB", ko="K99002")]
    rules = [rule_entry("ab_trait", "abcA+abcB"), rule_entry("b_trait", "abcB")]
    shared = protein("SAMPLE1_00001", "AbcB protein", dbxrefs="KEGG:K99002")
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    alone = run(tmp_path, [shared], genes, rules, hits=hits)
    a_row = rows_for(alone, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["counted"], a_row["evidence"], a_row["reason"]) == (
        "yes", "S", "structure_counts_alone")
    assert a_row["detail"] == ("structure hit abcA@invented route; no name matches, "
                               "and a structure hit counts alone for abcA (s_alone "
                               "yes); the protein counts once, for abcA, not for "
                               "abcB (Q)")
    b_row = rows_for(alone, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["decision"], b_row["reason"]) == (
        "no", "not counted", "counted_for_other_gene")
    assert b_row["evidence"] == "Q"
    assert b_row["detail"] == ("Bakta name matches abcB by product, KO K99002; the "
                               "structure hit is to abcA, which counts from "
                               "structure alone (s_alone yes), so the protein "
                               "counts once, for abcA")
    audit = alone["audit"]["ab_trait"]
    assert (audit["decision"], audit["reason"]) == ("partial", "half_of_set_found")
    assert alone["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S)"
    assert alone["profile"]["ab_trait"]["missing"] == "abcB"
    assert alone["audit"]["b_trait"]["decision"] == "not called"

    with_abcb = run(tmp_path, [shared, protein("SAMPLE1_00002", "AbcB protein")],
                    genes, rules, hits=hits)
    assert with_abcb["audit"]["ab_trait"]["decision"] == "called"
    assert with_abcb["profile"]["ab_trait"]["genes"] == (
        "abcA:SAMPLE1_00001(S); abcB:SAMPLE1_00002(Q)")


def test_fusion_named_by_one_source_counts_for_both_genes(tmp_path):
    """A fused enzyme whose Bakta product names both its genes, as a PqqC/D
    fusion's does, counts for both, with a structure hit to one of them too: one
    source, the Bakta name, names both genes. The same protein named after its
    second part only has its name and its structure pointing to two different
    genes, so it counts once, for the structure's gene."""
    genes = [gene_entry("abcC", s_alone="yes", product_include=r"\bAbcC(/D)?\b"),
             gene_entry("abcD", product_include=r"\bAbcD\b|\bAbcC/D\b")]
    rules = [rule_entry("cd_trait", "abcC+abcD")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCC")]

    fusion = run(tmp_path, [protein("SAMPLE1_00001", "AbcC/D fusion protein")],
                 genes, rules, hits=hits, catalogue=CATALOGUE_SIX)
    assert [(row["gene"], row["evidence"], row["reason"])
            for row in fusion["genes"]] == [("abcC", "S+Q", "name_and_structure_agree"),
                                            ("abcD", "Q", "name_only")]
    assert fusion["audit"]["cd_trait"]["decision"] == "called"
    assert fusion["profile"]["cd_trait"]["genes"] == (
        "abcC:SAMPLE1_00001(S+Q); abcD:SAMPLE1_00001(Q)")

    one_part = run(tmp_path, [protein("SAMPLE1_00001", "AbcD protein")], genes, rules,
                   hits=hits, catalogue=CATALOGUE_SIX)
    assert [(row["gene"], row["evidence"], row["reason"])
            for row in one_part["genes"]] == [("abcC", "S", "structure_counts_alone"),
                                              ("abcD", "Q", "counted_for_other_gene")]
    assert one_part["audit"]["cd_trait"]["decision"] == "partial"


def test_product_and_ko_of_one_bakta_name_count_for_two_genes(tmp_path):
    """Bakta's product points to abcA and its KO to abcB: still one source, the
    Bakta name, so the protein counts for both. When the two pieces come from two
    sources, the KO from Bakta and the product from Baktfold, the protein counts
    once, for the gene of Bakta's name, which ranks before Baktfold's."""
    genes = [gene_entry("abcA"), gene_entry("abcB", product_include="", ko="K99002")]
    rules = [rule_entry("ab_trait", "abcA+abcB")]

    one_name = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein",
                                      dbxrefs="KEGG:K99002")], genes, rules)
    assert one_name["audit"]["ab_trait"]["decision"] == "called"
    assert one_name["profile"]["ab_trait"]["genes"] == (
        "abcA:SAMPLE1_00001(Q); abcB:SAMPLE1_00001(Q)")

    two_names = run(tmp_path, [protein("SAMPLE1_00001", dbxrefs="KEGG:K99002")],
                    genes, rules, baktfold=[protein("SAMPLE1_00001", "AbcA protein")])
    a_row = rows_for(two_names, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["evidence"], a_row["reason"]) == ("F", "counted_for_other_gene")
    assert a_row["detail"] == ("Baktfold name matches abcA by product; Bakta name "
                               "matches abcB by KO K99002, which ranks before a "
                               "Baktfold name, so the protein counts once, for abcB")
    b_row = rows_for(two_names, "SAMPLE1_00001", "abcB")[0]
    assert b_row["detail"] == ("Bakta name matches by KO K99002; no structure hit; "
                               "the protein counts once, for abcB, not for abcA (F)")
    assert two_names["audit"]["ab_trait"]["decision"] == "partial"


def test_structure_hit_that_needs_a_name_takes_no_part(tmp_path):
    """The structure hit points to abcA, which counts from structure only when a
    name agrees (s_alone no), and Bakta's name to abcB. The abcA row does not count
    anyway (structure_only_needs_name), so there is no clash: the protein counts
    for abcB, and its row carries no note. With s_alone yes for abcA, the same
    protein counts once, for abcA."""
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    proteins = [protein("SAMPLE1_00001", "AbcB protein")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    needs_name = run(tmp_path, proteins, [gene_entry("abcA"), gene_entry("abcB")],
                     rules, hits=hits)
    a_row = rows_for(needs_name, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["counted"], a_row["reason"]) == ("no", "structure_only_needs_name")
    b_row = rows_for(needs_name, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["reason"]) == ("yes", "name_only")
    assert b_row["detail"] == ("Bakta name matches by product; best structure hit "
                               "abcA@invented route is not in the ecofolddb list of "
                               "abcB")
    assert needs_name["profile"]["ab_trait"]["genes"] == "abcB:SAMPLE1_00001(Q)"

    alone = run(tmp_path, proteins, [gene_entry("abcA", s_alone="yes"),
                                     gene_entry("abcB")], rules, hits=hits)
    b_row = rows_for(alone, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["reason"]) == ("no", "counted_for_other_gene")
    assert alone["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S)"


def test_structure_and_bakta_name_rank_before_a_baktfold_name(tmp_path):
    """S+Q against F: the structure hit and a KO, which Bakta gave to a protein it
    could only call hypothetical, point to abcA; Baktfold's name points to abcB.
    The protein counts once, for abcA, whatever abcA's s_alone says, but the
    reason differs. With s_alone no the structure hit does not decide, so the
    abcB row says that Bakta's name ranks before Baktfold's; with s_alone yes it
    says that abcA is the structure's gene. This pins down that S+Q keeps the
    protein against F, and the reason given in each case."""
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    proteins = [protein("SAMPLE1_00001", dbxrefs="KEGG:K99001")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]
    baktfold = [protein("SAMPLE1_00001", "AbcB protein")]

    needs_name = run(tmp_path, proteins, [gene_entry("abcA", ko="K99001"),
                                          gene_entry("abcB")],
                     rules, hits=hits, baktfold=baktfold)
    a_row = rows_for(needs_name, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["evidence"], a_row["counted"], a_row["reason"]) == (
        "S+Q", "yes", "name_and_structure_agree")
    assert a_row["detail"] == ("Bakta name matches by KO K99001; structure hit "
                               "abcA@invented route; the protein counts once, for "
                               "abcA, not for abcB (F)")
    b_row = rows_for(needs_name, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["evidence"], b_row["counted"], b_row["reason"]) == (
        "F", "no", "counted_for_other_gene")
    assert b_row["detail"] == ("Baktfold name matches abcB by product; Bakta name "
                               "matches abcA by KO K99001 and the structure hit is "
                               "to abcA; Bakta's name ranks before a Baktfold name, "
                               "so the protein counts once, for abcA")
    assert needs_name["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S+Q)"

    alone = run(tmp_path, proteins, [gene_entry("abcA", ko="K99001", s_alone="yes"),
                                     gene_entry("abcB")],
                rules, hits=hits, baktfold=baktfold)
    b_row = rows_for(alone, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["reason"]) == ("no", "counted_for_other_gene")
    assert b_row["detail"] == ("Baktfold name matches abcB by product; Bakta name "
                               "matches abcA by KO K99001 and the structure hit is "
                               "to abcA, which ranks first as the structure's gene "
                               "(s_alone yes), so the protein counts once, for abcA")
    assert alone["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S+Q)"


def test_structure_alone_ranks_before_both_names(tmp_path):
    """Q+F against S: Bakta's KO and Baktfold's name both point to abcB, the
    structure hit to abcA, for which a structure hit counts alone (s_alone yes).
    Two names together do not outrank the structure: the protein counts once, for
    abcA, as it would against one name, because a structure hit to a gene with
    s_alone yes ranks first, however many names point elsewhere. Counting the
    letters instead would give the protein to abcB."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB", ko="K99002")]
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001", dbxrefs="KEGG:K99002")], genes,
                 rules, hits=[kept("SAMPLE1_00001", "CAT_ABCA")],
                 baktfold=[protein("SAMPLE1_00001", "AbcB protein")])
    a_row = rows_for(result, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["evidence"], a_row["counted"], a_row["reason"]) == (
        "S", "yes", "structure_counts_alone")
    assert a_row["detail"] == ("structure hit abcA@invented route; no name matches, "
                               "and a structure hit counts alone for abcA (s_alone "
                               "yes); the protein counts once, for abcA, not for "
                               "abcB (Q+F)")
    b_row = rows_for(result, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["evidence"], b_row["counted"], b_row["reason"]) == (
        "Q+F", "no", "counted_for_other_gene")
    assert b_row["detail"] == ("Bakta name matches abcB by KO K99002 and Baktfold "
                               "name matches abcB by product; the structure hit is "
                               "to abcA, which counts from structure alone (s_alone "
                               "yes), so the protein counts once, for abcA")
    assert result["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S)"


def test_bakta_name_ranks_before_structure_and_baktfold_name_unless_s_alone(tmp_path):
    """S+F against Q: the structure hit and Baktfold's name point to abcA; a KO,
    which Bakta gave to the same protein while calling it hypothetical, points to
    abcB. With s_alone no for abcA, the structure hit does not decide, and
    Baktfold's name rests on the same predicted structure, so the names decide:
    the protein counts once, for abcB, the gene of Bakta's name, and the abcA row
    says why its structure hit did not decide. With s_alone yes, the structure hit
    decides, and the protein counts once, for abcA. This pins down both sides of
    the s_alone condition on a structure hit that a Baktfold name agrees with."""
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    proteins = [protein("SAMPLE1_00001", dbxrefs="KEGG:K99002")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]
    baktfold = [protein("SAMPLE1_00001", "AbcA protein")]
    abcb = gene_entry("abcB", product_include="", ko="K99002")

    needs_name = run(tmp_path, proteins, [gene_entry("abcA"), abcb], rules,
                     hits=hits, baktfold=baktfold)
    a_row = rows_for(needs_name, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["evidence"], a_row["counted"], a_row["reason"]) == (
        "S+F", "no", "counted_for_other_gene")
    assert a_row["detail"] == ("Baktfold name matches abcA by product and the "
                               "structure hit is to abcA; Bakta name matches abcB by "
                               "KO K99002, which ranks before a Baktfold name; the "
                               "structure hit does not decide, because abcA does "
                               "not count from structure alone (s_alone no), so the "
                               "protein counts once, for abcB")
    b_row = rows_for(needs_name, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["evidence"], b_row["counted"], b_row["reason"]) == (
        "Q", "yes", "name_only")
    assert b_row["detail"].endswith("; the protein counts once, for abcB, not for "
                                    "abcA (S+F)")
    assert needs_name["profile"]["ab_trait"]["genes"] == "abcB:SAMPLE1_00001(Q)"

    alone = run(tmp_path, proteins, [gene_entry("abcA", s_alone="yes"), abcb], rules,
                hits=hits, baktfold=baktfold)
    a_row = rows_for(alone, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["evidence"], a_row["counted"], a_row["reason"]) == (
        "S+F", "yes", "name_and_structure_agree")
    b_row = rows_for(alone, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["evidence"], b_row["counted"], b_row["reason"]) == (
        "Q", "no", "counted_for_other_gene")
    assert b_row["detail"] == ("Bakta name matches abcB by KO K99002; Baktfold name "
                               "matches abcA by product and the structure hit is to "
                               "abcA, which ranks first as the structure's gene "
                               "(s_alone yes), so the protein counts once, for abcA")
    assert alone["profile"]["ab_trait"]["genes"] == "abcA:SAMPLE1_00001(S+F)"


def test_genes_that_share_no_rule_are_not_compared(tmp_path):
    """abcA (structure) and abcB (Bakta name) on one protein, each in a trait of its
    own: nothing is compared, and both traits are called. Once abcB is a
    supporting gene of abcA's trait, the two share a rule: the protein counts once,
    for abcA, and so not for abcB in abcB's own trait either."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB")]
    proteins = [protein("SAMPLE1_00001", "AbcB protein")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    apart = run(tmp_path, proteins, genes,
                [rule_entry("a_trait", "abcA"), rule_entry("b_trait", "abcB")],
                hits=hits)
    assert [(row["gene"], row["reason"]) for row in apart["genes"]] == [
        ("abcA", "structure_counts_alone"), ("abcB", "name_only")]
    assert apart["audit"]["a_trait"]["decision"] == "called"
    assert apart["audit"]["b_trait"]["decision"] == "called"

    sharing = run(tmp_path, proteins, genes,
                  [rule_entry("a_trait", "abcA", supporting="abcB"),
                   rule_entry("b_trait", "abcB")], hits=hits)
    assert [(row["gene"], row["reason"]) for row in sharing["genes"]] == [
        ("abcA", "structure_counts_alone"), ("abcB", "counted_for_other_gene")]
    assert sharing["audit"]["a_trait"]["decision"] == "called"
    assert sharing["profile"]["a_trait"]["supporting"] == "NA"
    assert sharing["audit"]["b_trait"]["decision"] == "not called"


def test_vetoed_rows_take_no_part(tmp_path):
    """Rows that do not count anyway are not compared: a vetoed abcA row leaves
    abcB to count. The same protein, not vetoed, counts once, for abcA. (How a
    pseudogene is compared: test_pseudogene_is_taken_as_one_gene_per_rule.)"""
    genes = [gene_entry("abcA", s_alone="yes", product_exclude="decoy"),
             gene_entry("abcB")]
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    vetoed = run(tmp_path, [protein("SAMPLE1_00001", "AbcB decoy protein")], genes,
                 rules, hits=hits)
    assert [(row["gene"], row["reason"]) for row in vetoed["genes"]] == [
        ("abcA", "vetoed_by_name"), ("abcB", "name_only")]
    assert vetoed["profile"]["ab_trait"]["genes"] == "abcB:SAMPLE1_00001(Q)"

    plain = run(tmp_path, [protein("SAMPLE1_00001", "AbcB protein")], genes, rules,
                hits=hits)
    assert [(row["gene"], row["reason"]) for row in plain["genes"]] == [
        ("abcA", "structure_counts_alone"), ("abcB", "counted_for_other_gene")]


def test_pseudogene_is_taken_as_one_gene_per_rule(tmp_path):
    """A pseudogene is compared as a working protein is. On a broken protein the
    structure hit points to abcA (s_alone yes) and Bakta's name to abcB, two genes
    of one set from different sources: the trait calls take the protein as a
    broken abcA only, the gene that ranks first. Both rows keep reason pseudogene,
    and their details say which gene the protein is taken as. With abcC counted,
    the one broken protein cannot fill two places, so the set is not complete
    except for pseudogenes: 1 of 3 genes is too few. Without abcC, the audit names
    abcA alone as found only as a pseudogene. So a broken protein whose sources
    point to two genes of one rule, and no one source names both, is taken as one
    of them, on the same terms as a working protein counts, and is left out of the
    trait calls for the gene it lost."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB"), gene_entry("abcC")]
    broken = protein("SAMPLE1_00001", "(pseudo) AbcB protein")
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]
    pseudo_note = ("; Bakta marks the gene as a pseudogene (disrupted reading "
                   "frame), so it does not count")

    with_abcc = run(tmp_path, [broken, protein("SAMPLE1_00002", "AbcC protein")],
                    genes, [abc_rule()], hits=hits)
    a_row = rows_for(with_abcc, "SAMPLE1_00001", "abcA")[0]
    assert (a_row["counted"], a_row["decision"], a_row["reason"]) == (
        "no", "not counted", "pseudogene")
    assert a_row["detail"] == ("structure hit abcA@invented route; no name matches, "
                               "and a structure hit counts alone for abcA (s_alone "
                               "yes); the broken protein is taken as abcA, not as "
                               "abcB (Q)" + pseudo_note)
    b_row = rows_for(with_abcc, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["decision"], b_row["reason"]) == (
        "no", "not counted", "pseudogene")
    assert b_row["detail"] == ("Bakta name matches by product; best structure hit "
                               "abcA@invented route is not in the ecofolddb list of "
                               "abcB; the structure hit is to abcA, which counts "
                               "from structure alone (s_alone yes), so the broken "
                               "protein is taken as abcA" + pseudo_note)
    audit = with_abcc["audit"]["abc_trait"]
    assert (audit["decision"], audit["reason"]) == ("not called", "too_few_genes")
    assert audit["detail"] == ("1 of 3 genes found in the closest set abcA+abcB+abcC "
                               "(missing abcA:SAMPLE1_00001(pseudogene), abcB)")

    two_genes = run(tmp_path, [broken], genes[:2],
                    [rule_entry("ab_trait", "abcA+abcB")], hits=hits)
    audit = two_genes["audit"]["ab_trait"]
    assert (audit["decision"], audit["reason"]) == ("not called", "no_gene_found")
    assert audit["detail"] == ("no gene of any required set counted; found only as a "
                               "pseudogene: abcA:SAMPLE1_00001(pseudogene)")


def test_fusion_pseudogene_named_by_one_source_is_taken_as_both_genes(tmp_path):
    """A broken fusion whose Bakta product names both its genes is one source
    naming both, as for a working fusion: with a structure hit to abcA too, the
    trait calls take the protein as abcA and as abcB, and neither row notes the
    other gene. With abcC counted, the set is complete except for pseudogenes,
    and the audit and the missing column name the broken protein for both genes.
    So the "one source names both" test applies to a broken protein too: it is
    taken as just one of two genes of a rule only when no one source names
    both."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB"), gene_entry("abcC")]
    proteins = [protein("SAMPLE1_00001", "(pseudo) AbcA/AbcB fusion protein"),
                protein("SAMPLE1_00002", "AbcC protein")]
    result = run(tmp_path, proteins, genes, [abc_rule()],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")])
    pseudo_note = ("; Bakta marks the gene as a pseudogene (disrupted reading "
                   "frame), so it does not count")
    assert [(row["gene"], row["evidence"], row["reason"], row["detail"])
            for row in rows_for(result, "SAMPLE1_00001")] == [
        ("abcA", "S+Q", "pseudogene",
         "Bakta name matches by product; structure hit abcA@invented route"
         + pseudo_note),
        ("abcB", "Q", "pseudogene",
         "Bakta name matches by product; best structure hit abcA@invented route is "
         "not in the ecofolddb list of abcB" + pseudo_note)]
    audit = result["audit"]["abc_trait"]
    assert (audit["decision"], audit["reason"]) == ("partial",
                                                    "set_complete_but_pseudogene")
    assert audit["detail"] == ("abcA+abcB+abcC complete except "
                               "abcA:SAMPLE1_00001(pseudogene), "
                               "abcB:SAMPLE1_00001(pseudogene)")
    assert result["profile"]["abc_trait"]["missing"] == (
        "abcA:SAMPLE1_00001(pseudogene), abcB:SAMPLE1_00001(pseudogene)")


def test_fusion_named_by_baktfold_counts_for_both_genes(tmp_path):
    """One Baktfold name that names both genes is one source too: a protein Bakta
    could only call hypothetical, which Baktfold names after both parts of a
    fusion, counts for both genes (F), and the trait is called."""
    genes = [gene_entry("abcA"), gene_entry("abcB")]
    rules = [rule_entry("ab_trait", "abcA+abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001")], genes, rules,
                 baktfold=[protein("SAMPLE1_00001", "AbcA/AbcB fusion protein")])
    assert [(row["gene"], row["evidence"], row["counted"], row["reason"])
            for row in result["genes"]] == [("abcA", "F", "yes", "name_only"),
                                            ("abcB", "F", "yes", "name_only")]
    assert result["audit"]["ab_trait"]["decision"] == "called"
    assert result["profile"]["ab_trait"]["genes"] == (
        "abcA:SAMPLE1_00001(F); abcB:SAMPLE1_00001(F)")


def test_structure_hit_to_an_entry_both_genes_list_counts_for_both(tmp_path):
    """One structure hit whose catalogue entry two genes list (as aphA and phoN
    both list aphA@Phosphate ester degradation in the shipped gene dictionary) is
    one source naming both: with s_alone yes for both, the protein counts for both
    (S). A third gene of the set, named by Bakta, comes from a source that names
    neither, so the protein does not count for it. Its row names abcA, the first
    of the two in the ranking (on this tie, the order of the gene dictionary), and
    only the abcA row says the protein does not count for abcC."""
    genes = [gene_entry("abcA", s_alone="yes"),
             gene_entry("abcB", s_alone="yes", ecofolddb="abcA@invented route"),
             gene_entry("abcC")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA")]

    two = run(tmp_path, [protein("SAMPLE1_00001")], genes[:2],
              [rule_entry("ab_trait", "abcA+abcB")], hits=hits)
    assert [(row["gene"], row["evidence"], row["counted"], row["reason"])
            for row in two["genes"]] == [("abcA", "S", "yes", "structure_counts_alone"),
                                         ("abcB", "S", "yes", "structure_counts_alone")]
    assert two["audit"]["ab_trait"]["decision"] == "called"

    three = run(tmp_path, [protein("SAMPLE1_00001", "AbcC protein")], genes,
                [rule_entry("abc_trait", "abcA+abcB+abcC")], hits=hits)
    assert [(row["gene"], row["counted"], row["reason"]) for row in three["genes"]] == [
        ("abcA", "yes", "structure_counts_alone"),
        ("abcB", "yes", "structure_counts_alone"),
        ("abcC", "no", "counted_for_other_gene")]
    assert rows_for(three, "SAMPLE1_00001", "abcA")[0]["detail"].endswith(
        "; the protein counts once, for abcA, not for abcC (Q)")
    assert rows_for(three, "SAMPLE1_00001", "abcB")[0]["detail"] == (
        "structure hit abcA@invented route; no name matches, and a structure hit "
        "counts alone for abcB (s_alone yes)")
    assert rows_for(three, "SAMPLE1_00001", "abcC")[0]["detail"] == (
        "Bakta name matches abcC by product; the structure hit is to abcA, which "
        "counts from structure alone (s_alone yes), so the protein counts once, "
        "for abcA")
    assert three["profile"]["abc_trait"]["genes"] == (
        "abcA:SAMPLE1_00001(S); abcB:SAMPLE1_00001(S)")


def test_a_gene_that_lost_the_protein_takes_no_further_part(tmp_path):
    """A gene that has lost the protein cannot take it from a third gene. The
    structure hit points to abcA (s_alone yes), Bakta's KO to abcB and Baktfold's
    name to abcC; abcA and abcB share trait r1, abcB and abcC share trait r2.
    abcB loses the protein to abcA in r1, and so counts nowhere. abcC shares no
    trait with abcA, so it keeps the protein, and its row carries no note. Each
    trait is partial, with one of its two genes."""
    genes = [gene_entry("abcA", s_alone="yes"),
             gene_entry("abcB", product_include="", ko="K99002"),
             gene_entry("abcC")]
    rules = [rule_entry("r1", "abcA+abcB"), rule_entry("r2", "abcB+abcC")]
    result = run(tmp_path, [protein("SAMPLE1_00001", dbxrefs="KEGG:K99002")], genes,
                 rules, hits=[kept("SAMPLE1_00001", "CAT_ABCA")],
                 baktfold=[protein("SAMPLE1_00001", "AbcC protein")])
    assert [(row["gene"], row["evidence"], row["counted"], row["reason"])
            for row in result["genes"]] == [
        ("abcA", "S", "yes", "structure_counts_alone"),
        ("abcB", "Q", "no", "counted_for_other_gene"),
        ("abcC", "F", "yes", "name_only")]
    assert rows_for(result, "SAMPLE1_00001", "abcC")[0]["detail"] == (
        "Baktfold name matches by product; best structure hit abcA@invented route "
        "is not in the ecofolddb list of abcC")
    assert result["profile"]["r1"]["genes"] == "abcA:SAMPLE1_00001(S)"
    assert result["profile"]["r2"]["genes"] == "abcC:SAMPLE1_00001(F)"
    assert result["audit"]["r1"]["reason"] == "half_of_set_found"
    assert result["audit"]["r2"]["reason"] == "half_of_set_found"


def test_one_gene_can_keep_the_protein_against_two_others(tmp_path):
    """Three genes of one set from three sources: the structure hit points to abcA
    (s_alone yes), Bakta's KO to abcB and Baktfold's name to abcC. The protein
    counts for abcA alone, and the abcA row names both genes it does not count
    for, in the order of the gene dictionary, joined by "or"."""
    genes = [gene_entry("abcA", s_alone="yes"),
             gene_entry("abcB", product_include="", ko="K99002"),
             gene_entry("abcC")]
    result = run(tmp_path, [protein("SAMPLE1_00001", dbxrefs="KEGG:K99002")], genes,
                 [abc_rule()], hits=[kept("SAMPLE1_00001", "CAT_ABCA")],
                 baktfold=[protein("SAMPLE1_00001", "AbcC protein")])
    assert [(row["gene"], row["reason"]) for row in result["genes"]] == [
        ("abcA", "structure_counts_alone"), ("abcB", "counted_for_other_gene"),
        ("abcC", "counted_for_other_gene")]
    assert rows_for(result, "SAMPLE1_00001", "abcA")[0]["detail"] == (
        "structure hit abcA@invented route; no name matches, and a structure hit "
        "counts alone for abcA (s_alone yes); the protein counts once, for abcA, not "
        "for abcB (Q) or abcC (F)")
    assert result["profile"]["abc_trait"]["decision"] == "not called"


def test_marked_protein_in_a_clash_keeps_the_mark_note_on_both_rows(tmp_path):
    """Bakta's mark is noted on every gene row of the protein, clash or not: the
    abcA row that keeps the protein and the abcB row that loses it both end with
    the note of a gene that runs off its sequence."""
    genes = [gene_entry("abcA", s_alone="yes"), gene_entry("abcB")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "(5' truncated) AbcB protein")],
                 genes, [rule_entry("ab_trait", "abcA+abcB")],
                 hits=[kept("SAMPLE1_00001", "CAT_ABCA")])
    note = "; Bakta: runs off its sequence at the gene's 5' end (5' truncated)"
    assert rows_for(result, "SAMPLE1_00001", "abcA")[0]["detail"] == (
        "structure hit abcA@invented route; no name matches, and a structure hit "
        "counts alone for abcA (s_alone yes); the protein counts once, for abcA, not "
        "for abcB (Q)" + note)
    b_row = rows_for(result, "SAMPLE1_00001", "abcB")[0]
    assert (b_row["counted"], b_row["reason"]) == ("no", "counted_for_other_gene")
    assert b_row["detail"] == (
        "Bakta name matches abcB by product; the structure hit is to abcA, which "
        "counts from structure alone (s_alone yes), so the protein counts once, for "
        "abcA" + note)


# ── antiSMASH regions as supporting evidence ─────────────────────────────────
# An antiSMASH region type in a rule's supporting list is reported when a region
# of exactly that type lies over a protein counted for one of the trait's
# required genes.

def area(start, end, products):
    """An antiSMASH area: 0-based start, exclusive end, as antiSMASH writes them."""
    return {"start": start, "end": end, "products": products}


def test_antismash_region_over_a_required_gene_is_reported(tmp_path):
    """An antiSMASH region of the named type on the same contig, over a protein
    counted for a required gene, is reported as supporting evidence."""
    rules = [rule_entry("a", "abcA", supporting="antiSMASH:NRP-metallophore")]
    records = [{"id": "contig_1", "name": "contig_1",
                "areas": [area(500, 2500, ["NRP-metallophore", "NRPS"])]}]
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], rules, antismash=records)
    assert result["profile"]["a"]["supporting"] == ("antiSMASH:NRP-metallophore"
                                                    "(contig_1:501-2500)")


def test_antismash_region_must_have_the_exact_type_and_contig(tmp_path):
    """An NRPS region over the gene is not an NRP-metallophore region, and an
    NRP-metallophore region on another contig is not over the gene: neither is
    reported."""
    rules = [rule_entry("a", "abcA", supporting="antiSMASH:NRP-metallophore")]
    records = [{"id": "contig_1", "name": "contig_1", "areas": [area(0, 5000, ["NRPS"])]},
               {"id": "contig_2", "name": "contig_2",
                "areas": [area(0, 5000, ["NRP-metallophore"])]}]
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")],
                 [gene_entry("abcA")], rules, antismash=records)
    assert result["profile"]["a"]["supporting"] == "NA"


def test_antismash_region_across_the_origin_is_matched(tmp_path):
    """antiSMASH writes a region that runs across the origin of a circular contig
    with a start larger than its end: here 19956..contig end plus 1..1500. A
    protein in either piece is covered; one between the pieces is not."""
    rules = [rule_entry("a", "abcA", supporting="antiSMASH:NRP-metallophore")]
    records = [{"id": "contig_1", "name": "contig_1",
                "areas": [area(19955, 1500, ["NRP-metallophore"])]}]
    reported = "antiSMASH:NRP-metallophore(contig_1:19956-1500)"

    for start in (1000, 20000):        # the piece at the contig start, at its end
        result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein", start=start)],
                     [gene_entry("abcA")], rules, antismash=records)
        assert result["profile"]["a"]["supporting"] == reported, start

    between = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein", start=5000)],
                  [gene_entry("abcA")], rules, antismash=records)
    assert between["profile"]["a"]["supporting"] == "NA"


def test_overlap_of_a_protein_across_the_origin():
    """A protein whose start is larger than its stop runs across the end of a
    circular contig: it covers start..contig end plus 1..stop, so it overlaps an
    area that reaches either piece, or one that runs across the end itself, and
    not an area lying between the two pieces."""
    across = {"start": 19000, "stop": 500}
    assert tp.overlaps(across, 100, 200)            # inside the piece 1..500
    assert tp.overlaps(across, 19500, 19600)        # inside the piece 19000..end
    assert tp.overlaps(across, 18000, 300)          # an area across the end too
    assert not tp.overlaps(across, 5000, 6000)      # between the two pieces


def test_antismash_region_type_is_not_checked_without_antismash(tmp_path):
    """Without antiSMASH output, an antiSMASH supporting type reads "(not checked)",
    so "not looked for" is not mistaken for "not found"; supporting genes are still
    reported."""
    rules = [rule_entry("a", "abcA", supporting="abcB | antiSMASH:NRP-metallophore")]
    proteins = [protein("SAMPLE1_00001", "AbcA protein"),
                protein("SAMPLE1_00002", "AbcB protein")]
    result = run(tmp_path, proteins, [gene_entry("abcA"), gene_entry("abcB")], rules)
    assert result["profile"]["a"]["supporting"] == ("abcB:SAMPLE1_00002(Q); "
                                                    "antiSMASH:NRP-metallophore"
                                                    "(not checked)")


def test_antismash_record_is_matched_by_name_before_id(tmp_path):
    """Bakta writes its GenBank file with Biopython, which takes a contig name with
    one dot followed only by digits for accession.version and drops the leading
    zeros of the version; antiSMASH reads that back as the record id. The record's
    name keeps the contig name whole. The contig name is invented."""
    contig = "NODE_1_length_5000_cov_7.05"
    rules = [rule_entry("a", "abcA", supporting="antiSMASH:NRP-metallophore")]
    records = [{"id": "NODE_1_length_5000_cov_7.5", "name": contig,
                "areas": [area(0, 5000, ["NRP-metallophore"])]}]
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein", contig=contig)],
                 [gene_entry("abcA")], rules, antismash=records)
    assert result["profile"]["a"]["supporting"] == (
        f"antiSMASH:NRP-metallophore({contig}:1-5000)")


def test_supporting_gene_alone_does_not_make_a_call(tmp_path):
    """A supporting gene is reported beside the call but is not part of it: with
    only the supporting gene present, the trait stays not called. (How a
    supporting gene shares a protein with a required gene of its rule is tested
    in test_genes_that_share_no_rule_are_not_compared.)"""
    rules = [rule_entry("x", "abcB+abcC", supporting="abcA")]
    result = run(tmp_path, [protein("SAMPLE1_00001", "AbcA protein")], ABC_GENES,
                 rules)
    assert result["audit"]["x"]["decision"] == "not called"
    assert result["audit"]["x"]["reason"] == "no_gene_found"
    assert result["profile"]["x"]["supporting"] == "abcA:SAMPLE1_00001(Q)"


# ── Command line and output shape ────────────────────────────────────────────
# The headers and vocabulary of the three tables, the order of the gene table,
# and what happens with an empty genome or a Bakta table without its header.

def test_outputs_have_the_documented_headers_and_vocabulary(tmp_path):
    """The three headers are a contract with rule traits_summary and with the
    manual, and every reason must come from the vocabulary in the script."""
    proteins = (named_abc() + [protein("SAMPLE1_00004", "Decoy AbcB"),
                               protein("SAMPLE1_00005"),
                               protein("SAMPLE1_00006", "(pseudo) AbcC protein")])
    genes = [gene_entry("abcA"), gene_entry("abcB", product_exclude="decoy"),
             gene_entry("abcC")]
    hits = [kept("SAMPLE1_00001", "CAT_ABCA"), kept("SAMPLE1_00005", "CAT_XYZA")]
    args = write_inputs(tmp_path, proteins, genes, [abc_rule()], hits=hits,
                        lengths={"SAMPLE1_00003": 5000})
    assert tp.main(args) == 0
    out = tmp_path / "out"

    header, gene_rows = read_tsv(out / "SAMPLE1_trait_genes.tsv")
    assert header == [
        "sample", "locus_tag", "contig", "start", "stop", "strand", "cds_index",
        "gene", "traits", "counted", "evidence",
        "bakta_gene", "bakta_product", "ko", "ec", "baktfold_product",
        "ecofolddb_target", "ecofolddb_gene", "ecofolddb_pathway",
        "evalue", "qcov", "tcov", "borderline",
        "decision", "reason", "detail"]
    gene_reasons = {"name_and_structure_agree", "name_only", "structure_counts_alone",
                    "pseudogene", "counted_for_other_gene", "structure_only_needs_name",
                    "vetoed_by_name", "catalogue_gene_unused", "not_searched_too_long"}
    assert {row["reason"] for row in gene_rows} <= gene_reasons
    assert {row["sample"] for row in gene_rows} == {"SAMPLE1"}
    assert {row["decision"] for row in gene_rows} <= {"counted", "not counted",
                                                     "not searched"}

    header, profile_rows = read_tsv(out / "SAMPLE1_trait_profile.tsv")
    assert header == ["sample", "trait_id", "trait", "application", "decision",
                      "genes", "missing", "supporting", "borderline",
                      "structure_only", "caveat"]

    header, audit_rows = read_tsv(out / "SAMPLE1_trait_profile_audit.tsv")
    assert header == ["sample", "trait_id", "decision", "reason", "detail"]
    assert {row["reason"] for row in audit_rows} <= {
        "set_complete", "set_complete_not_colocalised", "set_complete_but_pseudogene",
        "half_of_set_found", "partial_withheld", "too_few_genes", "no_gene_found"}


def test_gene_table_is_sorted_by_contig_then_position_then_gene(tmp_path):
    """Gene-table rows follow the contigs in the order Bakta lists them, not by name,
    then the position on the contig, then the gene name (a fusion's two rows)."""
    proteins = [protein("SAMPLE1_00009", "AbcB protein", contig="contig_9"),
                protein("SAMPLE1_00002", "AbcC protein"),
                protein("SAMPLE1_00001", "AbcA-AbcB fusion")]
    result = run(tmp_path, proteins, ABC_GENES, [abc_rule()])
    order = [(row["contig"], row["locus_tag"], row["gene"]) for row in result["genes"]]
    assert order == [("contig_9", "SAMPLE1_00009", "abcB"),
                     ("contig_1", "SAMPLE1_00001", "abcA"),
                     ("contig_1", "SAMPLE1_00001", "abcB"),
                     ("contig_1", "SAMPLE1_00002", "abcC")]


def test_empty_genome_still_writes_every_header(tmp_path):
    """A genome with no proteins still gets all three tables, with headers, and every
    trait reads not called, no_gene_found: the rule's outputs always exist."""
    result = run(tmp_path, [], ABC_GENES, [abc_rule()])
    assert result["genes"] == []
    assert result["audit"]["abc_trait"]["reason"] == "no_gene_found"
    assert result["profile"]["abc_trait"]["decision"] == "not called"


def test_bakta_table_without_its_header_stops_the_run(tmp_path):
    """A Bakta table without its '#Sequence Id' header line stops the run: columns are
    found by name, and without the header no column can be told from another."""
    args = write_inputs(tmp_path, [protein("SAMPLE1_00001")], [gene_entry("abcA")],
                        [rule_entry("a", "abcA")])
    bakta = tmp_path / "SAMPLE1.tsv"
    bakta.write_text("contig_1\tcds\t1\t900\t+\tSAMPLE1_00001\t\tx\t\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no '#Sequence Id' header line"):
        tp.main(args)


# ── Checks on the two tables ─────────────────────────────────────────────────
# load_trait_tables, and read_trait_table inside it, stop the run on a slip in a
# trait table, naming the file and, for most slips, the line. Each case below is
# one such slip, except the last, which shows the '#' provenance line is skipped.

def load(tmp_path, genes, rules, catalogue=CATALOGUE):
    """Write two tables and load them as main() does."""
    catalogue_dict = {target: {"gene": gene, "pathway": pathway,
                               "key": gene + "@" + pathway}
                      for target, gene, pathway in catalogue}
    rules_path = write_table(tmp_path / "rules.tsv", tp.RULE_COLUMNS, rules)
    genes_path = write_table(tmp_path / "genes.tsv", tp.GENE_COLUMNS, genes)
    return tp.load_trait_tables(rules_path, genes_path, catalogue_dict)


@pytest.mark.parametrize("genes, rules, message", [
    ([gene_entry("abcA")], [rule_entry("a", "abcA+abcB")],
     "required gene 'abcB' has no entry in the gene dictionary"),
    ([gene_entry("abcA")], [rule_entry("a", "abcA", supporting="abcB")],
     "supporting gene 'abcB' has no entry in the gene dictionary"),
    ([gene_entry("abcA"), gene_entry("abcB")], [rule_entry("a", "abcA")],
     "gene entry 'abcB' is not used by any trait"),
    ([gene_entry("abcA"), gene_entry("abcB")],
     [rule_entry("a", "abcA", supporting="abcB", colocalised="abcA+abcB")],
     "colocalised gene 'abcB' is not in required"),
    (ABC_GENES, [rule_entry("a", "abcA+abcB | abcA+abcC", colocalised="abcB+abcC")],
     r"colocalised group 'abcB\+abcC' is not inside any one required set"),
    ([gene_entry("abcA"), gene_entry("abcB")],
     [rule_entry("a", "abcA", supporting="abcB", key_genes="abcB")],
     r"key gene 'abcB' \(key_genes\) is not in required"),
    ([gene_entry("abcA", product_include="(unclosed")], [rule_entry("a", "abcA")],
     "product_include is not a valid regular expression"),
    ([gene_entry("abcA", product_exclude="[z-a]")], [rule_entry("a", "abcA")],
     "product_exclude is not a valid regular expression"),
    ([gene_entry("abcA", s_alone="true")], [rule_entry("a", "abcA")],
     "s_alone must be yes or no"),
    ([gene_entry("abcA", symbol_needs_product="Yes")], [rule_entry("a", "abcA")],
     "symbol_needs_product must be yes or no"),
    ([gene_entry("abcA", ecofolddb="abcA@another route")], [rule_entry("a", "abcA")],
     "ecofolddb key 'abcA@another route' is not in the EcoFoldDB catalogue"),
    ([gene_entry("abcA")], [rule_entry("a", "abcA+antiSMASH:NRPS")],
     "antiSMASH region types are supporting evidence only"),
    ([gene_entry("abcA"), gene_entry("abcA")], [rule_entry("a", "abcA")],
     "a second entry for gene abcA"),
    ([gene_entry("abcA")], [rule_entry("a", "abcA"), rule_entry("a", "abcA")],
     "trait_id used twice"),
    ([gene_entry("abcA"), gene_entry("abcB")], [rule_entry("a", "abcA+abcA+abcB")],
     r"required set 'abcA\+abcA\+abcB' names gene 'abcA' more than once"),
    ([gene_entry("abcA"), gene_entry("abcB")],
     [rule_entry("a", "abcA", supporting="abcB | abcB")],
     "supporting lists 'abcB' more than once"),
])
def test_table_problems_stop_the_run_naming_the_row(tmp_path, genes, rules, message):
    """A broken table must stop the run: carrying on would give a profile of "not
    called" that looks like a result. The message names the file and the line of
    the row, or, for a gene entry no trait uses, the file and the entry."""
    with pytest.raises(ValueError, match=message) as error:
        load(tmp_path, genes, rules)
    assert " line " in str(error.value) or "gene entry" in str(error.value)


def test_row_with_a_wrong_number_of_fields_stops_the_run(tmp_path):
    """A trait-table row with fewer fields than the header (a lost tab, say) stops the
    run naming the line, instead of shifting the later cells into the wrong columns."""
    path = tmp_path / "genes.tsv"
    path.write_text("\t".join(tp.GENE_COLUMNS) + "\nabcA\tonly two fields\n",
                    encoding="utf-8")
    with pytest.raises(ValueError, match="line 2 has 2 fields"):
        tp.read_trait_table(str(path), tp.GENE_COLUMNS)


def test_table_without_a_needed_column_stops_the_run(tmp_path):
    """A trait table without one of the columns the script reads stops the run."""
    path = tmp_path / "rules.tsv"
    path.write_text("trait_id\trequired\na\tabcA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="header is missing column"):
        tp.read_trait_table(str(path), tp.RULE_COLUMNS)


def test_comment_lines_are_skipped(tmp_path):
    """The '#' provenance line at the top of a trait table is skipped, read neither as
    the header nor as a row."""
    rules, genes = load(tmp_path, [gene_entry("abcA")], [rule_entry("a", "abcA")])
    assert [rule["trait_id"] for rule in rules] == ["a"]
    assert list(genes) == ["abcA"]


# ── The two shipped tables ───────────────────────────────────────────────────
# The tables that ship with BacFlux must pass the same checks, read the same in
# every tool, and say where they come from.

def shipped_catalogue():
    """A stand-in catalogue holding exactly the keys the shipped gene dictionary
    uses, so the shipped tables can be loaded without the EcoFoldDB download. (Whether
    those keys exist in the real catalogue is checked when the layer runs.)"""
    catalogue = {}
    with open(tp.DEFAULT_GENES_TSV, encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if not line.startswith("#")]
    header = lines[0].split("\t")
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        for key in tp.split_list(row["ecofolddb"]):
            gene, pathway = key.split("@", 1)
            catalogue[f"STANDIN_{len(catalogue)}"] = {"gene": gene, "pathway": pathway,
                                                      "key": key}
    return catalogue


def test_shipped_tables_pass_every_check():
    """Every token has one entry and every entry is used, every pattern compiles,
    every yes/no cell is yes or no: the same checks the layer applies at run
    time, which stop it on any failure."""
    rules, genes = tp.load_trait_tables(tp.DEFAULT_RULES_TSV, tp.DEFAULT_GENES_TSV,
                                        shipped_catalogue())
    assert len(rules) == 49
    assert len(genes) == 274


def test_shipped_tables_read_the_same_with_csv_and_with_a_plain_split():
    """No cell may need CSV quoting: a plain split on tabs and Python's csv module
    must see the same fields, so the tables read the same in any tool."""
    for path in (tp.DEFAULT_RULES_TSV, tp.DEFAULT_GENES_TSV):
        with open(path, encoding="utf-8", newline="") as handle:
            text = handle.read()
        assert '"' not in text, path
        lines = text.rstrip("\n").split("\n")
        with_csv = list(csv.reader(lines, delimiter="\t"))
        with_split = [line.split("\t") for line in lines]
        assert with_csv == with_split, path


def test_shipped_tables_start_with_one_provenance_line():
    """Each shipped table opens with exactly one '#' line, which says where the table
    comes from and points to docs/analysis/traits.md, and its file name carries the
    version (_v1)."""
    for path in (tp.DEFAULT_RULES_TSV, tp.DEFAULT_GENES_TSV):
        with open(path, encoding="utf-8") as handle:
            lines = handle.read().split("\n")
        assert lines[0].startswith("# "), path
        assert "docs/analysis/traits.md" in lines[0], path
        assert not lines[1].startswith("#"), path
        assert os.path.basename(path).endswith("_v1.tsv")


def test_shipped_rule_and_gene_booleans_are_yes_or_no():
    """Every s_alone and symbol_needs_product cell of the shipped gene dictionary is
    yes or no, read straight from the file (the rules table has no yes/no column).
    Any other value would stop every run."""
    with open(tp.DEFAULT_GENES_TSV, encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if not line.startswith("#")]
    header = lines[0].split("\t")
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        assert row["s_alone"] in ("yes", "no"), row["gene"]
        assert row["symbol_needs_product"] in ("yes", "no"), row["gene"]
