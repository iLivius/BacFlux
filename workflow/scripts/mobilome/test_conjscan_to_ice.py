"""Unit tests for workflow/scripts/mobilome/conjscan_to_ice.py.

No tools, no databases, no MacSyFinder run: every test builds a small
CONJscan-shaped best_solution.tsv and a small Bakta-shaped GFF3 by hand and
pushes them through the script, so the whole anchor -> cluster -> classify chain
is verifiable in a plain Python environment.

The fixtures follow the real contracts:
  * CONJscan/MacSyFinder 2.1.6 best_solution.tsv - three '#' banner lines, then
    the verbatim 22-column header, then one row per machinery gene, with blank
    lines between systems;
  * Bakta 1.12.0 GFF3 - '##sequence-region' headers, Pyrodigal CDS rows carrying
    ID= and locus_tag= with the same value, percent-encoded product text, and a
    '##FASTA' section at the end.

Two tests go further and use the REAL output of sample 386 (an Arthrobacter
isolate): the saved fixture testdata/conjscan_386_real_best_solution.tsv, and
the five CDS lines quoted verbatim out of that sample's Bakta GFF3. Passing
those means the parser agrees with the tools, not just with itself.

There is also a contract test that feeds this script's output table straight into
colocalise.py's own parser, which is what consumes it in the workflow.

Run:
    python -m pytest workflow/scripts/mobilome/test_conjscan_to_ice.py -q
"""

import csv
import os
import sys

# Import the script under test (and its sibling consumer) regardless of where
# pytest was started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)
import conjscan_to_ice as ci        # noqa: E402
import colocalise as co             # noqa: E402


REAL_CONJSCAN_FIXTURE = os.path.join(
    TEST_DIR, "testdata", "conjscan_386_real_best_solution.tsv"
)


# ── Fixture builders ─────────────────────────────────────────────────────────

# The verbatim 22 columns of MacSyFinder 2.1.6 best_solution.tsv, in order.
CONJSCAN_COLUMNS = [
    "replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id",
    "sys_loci", "locus_num", "sys_wholeness", "sys_score", "sys_occ",
    "hit_gene_ref", "hit_status", "hit_seq_len", "hit_i_eval", "hit_score",
    "hit_profile_cov", "hit_seq_cov", "hit_begin_match", "hit_end_match",
    "counterpart", "used_in",
]

# One realistic machinery row, used as the base for every fixture row: the
# relaxase of a chromosomal MOB system, with a near-full profile alignment.
CONJSCAN_DEFAULTS = {
    "replicon": "S1",
    "hit_id": "S1_00010",
    "gene_name": "T4SS_MOBF",
    "hit_pos": "10",
    "model_fqn": "CONJScan/Chromosome/MOB",
    "sys_id": "S1_MOB_1",
    "sys_loci": "1",
    "locus_num": "1",
    "sys_wholeness": "1.000",
    "sys_score": "1.200",
    "sys_occ": "1",
    "hit_gene_ref": "T4SS_MOBB",
    "hit_status": "mandatory",
    "hit_seq_len": "546",
    "hit_i_eval": "1.9e-43",
    "hit_score": "146.900",
    "hit_profile_cov": "0.991",
    "hit_seq_cov": "0.498",
    "hit_begin_match": "1",
    "hit_end_match": "272",
    "counterpart": "",
    "used_in": "",
}


def conjscan_row(**overrides):
    """Build one tab-separated CONJscan data line in the tool's column order."""
    fields = dict(CONJSCAN_DEFAULTS)
    fields.update({key: str(value) for key, value in overrides.items()})
    return "\t".join(fields[column] for column in CONJSCAN_COLUMNS)


def write_conjscan(path, rows, with_header=True, with_banner=True):
    """Write a complete best_solution.tsv: banner comments, header, data rows."""
    lines = []
    if with_banner:
        lines += [
            "# macsyfinder 2.1.6  (using MacSyLib 1.0.4 )",
            "# models : CONJScan-2.1.0",
            "# Systems found:",
        ]
    if with_header:
        lines.append("\t".join(CONJSCAN_COLUMNS))
    lines += list(rows)
    with open(str(path), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return str(path)


def gff_cds(contig, start, end, strand, locus_tag, product):
    """Build one Bakta-shaped CDS line.

    Bakta writes the same value in ID= and locus_tag=, which is the join key
    CONJscan's hit_id uses, and percent-encodes commas in the product text.
    """
    encoded_product = product.replace(",", "%2C")
    attributes = (
        f"ID={locus_tag};Name={encoded_product};locus_tag={locus_tag};"
        f"product={encoded_product};Dbxref=SO:0001217"
    )
    return f"{contig}\tPyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}"


def write_gff(path, contig_lengths, cds_lines, with_fasta=True):
    """Write a small Bakta-shaped GFF3, optionally with the trailing FASTA block."""
    lines = ["##gff-version 3", "# Annotated with Bakta"]
    for contig, length in contig_lengths.items():
        lines.append(f"##sequence-region {contig} 1 {length}")
    lines += list(cds_lines)
    if with_fasta:
        # Everything after ##FASTA is sequence and must not be parsed as features.
        lines += ["##FASTA", ">contig_1", "ACGTACGTAC"]
    with open(str(path), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return str(path)


def read_tsv(path):
    """Read a written TSV back as (header list, list of row dicts)."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if row]
    header = rows[0]
    return header, [dict(zip(header, row)) for row in rows[1:]]


def run_main(tmp_path, conjscan=None, gff=None, contig_lengths=None, extra=None,
             sample="S1"):
    """Run main() the way the Snakemake rule will and read both outputs back.

    Returns (return_code, table_rows, audit_rows, out_table_path).
    `conjscan` may be None to model the file simply not being there.
    """
    os.makedirs(str(tmp_path), exist_ok=True)
    out_table = str(tmp_path / "ice_candidates.tsv")
    out_audit = str(tmp_path / "ice_decisions.tsv")

    argv = ["--sample", sample, "--bakta-gff", str(gff),
            "--out-table", out_table, "--out-audit", out_audit]
    if conjscan is not None:
        argv += ["--conjscan-tsv", str(conjscan)]
    if contig_lengths is not None:
        argv += ["--contig-lengths", str(contig_lengths)]
    argv += list(extra or [])

    return_code = ci.main(argv)
    _, table_rows = read_tsv(out_table)
    _, audit_rows = read_tsv(out_audit)
    return return_code, table_rows, audit_rows, out_table


def audit_reasons(audit_rows):
    """The set of reason tokens present in an audit file - what most tests assert on."""
    return {row["reason"] for row in audit_rows}


# A standard scene used by most end-to-end tests: one 200 kb contig carrying an
# integrase, a relaxase, a coupling protein and a VirB4, spaced a few kb apart so
# they all fall inside the default 15 kb clustering window. The machinery genes
# alone (relaxase 55000 -> VirB4 65500) already span more than the default 8 kb
# minimum, so the tests that remove the integrase still produce a candidate
# rather than being dropped for length.
SCENE_CONTIGS = {"contig_1": 200000}
SCENE_CDS = [
    gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
    gff_cds("contig_1", 55000, 56600, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
    gff_cds("contig_1", 57000, 58700, "-", "S1_00016",
            "Type IV secretory pathway, VirD4 component, TraG/TraD family ATPase"),
    gff_cds("contig_1", 63000, 65500, "+", "S1_00017", "conjugal transfer protein TraB"),
]

INTEGRASE_HIT = "S1_00010"
RELAXASE_HIT = "S1_00015"
T4CP_HIT = "S1_00016"
VIRB4_HIT = "S1_00017"


# ── The real CONJscan fixture ────────────────────────────────────────────────

def test_real_conjscan_fixture_parses():
    """The saved real output of sample 386 parses into its four machinery hits.

    This is the ground-truth check: two MOB systems on the chromosome, each
    incomplete (sys_wholeness 0.667), each made of a relaxase plus a coupling
    protein. If MacSyFinder ever changes its column names this test is what
    fails, rather than the workflow quietly reporting "no conjugation machinery".
    """
    hits = ci.read_conjscan_hits(REAL_CONJSCAN_FIXTURE)

    assert len(hits) == 4
    assert [hit["hit_id"] for hit in hits] == [
        "386_00040", "386_00044", "386_01146", "386_01151"
    ]
    assert [hit["gene_name"] for hit in hits] == [
        "T4SS_MOBP1", "T4SS_t4cp2", "T4SS_t4cp2", "T4SS_MOBF"
    ]
    assert {hit["sys_id"] for hit in hits} == {"386_MOB_1", "386_MOB_2"}
    assert all(hit["sys_wholeness"] == 0.667 for hit in hits)
    assert all(hit["model_fqn"] == "CONJScan/Chromosome/MOB" for hit in hits)


def test_real_fixture_classifies_the_two_relaxases_and_two_coupling_proteins():
    """The two relaxase families are relaxases, the two t4cp2 hits are coupling
    proteins - and NOT mating-pair components, which is what keeps sample 386 at
    'mobilisable' rather than 'self-transmissible'."""
    hits = ci.read_conjscan_hits(REAL_CONJSCAN_FIXTURE)
    classes = [ci.anchor_class_for_gene_name(hit["gene_name"]) for hit in hits]
    assert classes == [
        ci.ANCHOR_RELAXASE, ci.ANCHOR_T4CP, ci.ANCHOR_T4CP, ci.ANCHOR_RELAXASE
    ]


# Five CDS lines copied verbatim out of sample 386's real Bakta GFF3
# (04.annotation/bakta/386/386.gff3), so the coordinate join is tested against
# the real attribute layout rather than an invented one.
REAL_386_CDS_LINES = [
    "contig_2\tPyrodigal\tCDS\t37327\t38967\t.\t+\t0\t"
    "ID=386_00040;Name=Relaxase/mobilization nuclease family protein;"
    "locus_tag=386_00040;product=Relaxase/mobilization nuclease family protein;"
    "Dbxref=SO:0001217,UniRef:UniRef50_Q8GAN1",
    "contig_2\tPyrodigal\tCDS\t42322\t44091\t.\t-\t0\t"
    "ID=386_00044;Name=Type IV secretory pathway%2C VirD4 component%2C TraG/TraD family ATPase;"
    "locus_tag=386_00044;"
    "product=Type IV secretory pathway%2C VirD4 component%2C TraG/TraD family ATPase;"
    "Dbxref=SO:0001217;gene=virD4",
    "contig_1\tPyrodigal\tCDS\t1166923\t1168035\t.\t-\t0\t"
    "ID=386_01130;Name=Phage integrase family protein;locus_tag=386_01130;"
    "product=Phage integrase family protein;Dbxref=SO:0001217",
    "contig_1\tPyrodigal\tCDS\t1180113\t1181924\t.\t+\t0\t"
    "ID=386_01146;Name=type IV secretory system conjugative DNA transfer family protein;"
    "locus_tag=386_01146;"
    "product=type IV secretory system conjugative DNA transfer family protein;"
    "Dbxref=SO:0001217",
    "contig_1\tPyrodigal\tCDS\t1187809\t1191384\t.\t-\t0\t"
    "ID=386_01151;Name=TrwC relaxase domain-containing protein;locus_tag=386_01151;"
    "product=TrwC relaxase domain-containing protein;Dbxref=SO:0001217;gene=trwC",
]

REAL_386_CONTIGS = {"contig_2": 57976, "contig_1": 4177018}


def test_real_sample_386_end_to_end():
    """The real sample, end to end: one degraded IME, one cluster too short.

    This reproduces the finding recorded in docs/mobilome_wpA_ground_truth.md:

      * contig_1 carries a phage integrase, a coupling protein and a MOBF
        relaxase within the clustering window, and NO mating-pair apparatus, so
        it is an IME - mobilisable with a helper, not self-transmissible;
      * sys_wholeness is 0.667, so the machinery is flagged incomplete and the
        mobility sentence says so;
      * contig_2's relaxase + coupling protein span only 6765 bp, below the 8 kb
        minimum, so it is dropped with its measured span in the audit file.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as work_dir:
        gff = write_gff(os.path.join(work_dir, "386.gff3"),
                        REAL_386_CONTIGS, REAL_386_CDS_LINES)
        out_table = os.path.join(work_dir, "table.tsv")
        out_audit = os.path.join(work_dir, "audit.tsv")

        return_code = ci.main([
            "--sample", "386",
            "--conjscan-tsv", REAL_CONJSCAN_FIXTURE,
            "--bakta-gff", gff,
            "--out-table", out_table,
            "--out-audit", out_audit,
        ])
        assert return_code == 0

        _, rows = read_tsv(out_table)
        _, audit = read_tsv(out_audit)

    assert len(rows) == 1
    element = rows[0]
    assert element["contig"] == "contig_1"
    assert element["mge_class"] == "ime"
    assert element["element_type"] == "ime"
    assert element["start"] == "1166923"
    assert element["end"] == "1191384"
    assert element["has_integrase"] == "TRUE"
    assert element["has_relaxase"] == "TRUE"
    assert element["has_t4cp"] == "TRUE"
    assert element["has_t4ss"] == "FALSE"
    assert element["relaxase_type"] == "MOBF"
    assert element["mobility"] == "mobilisable (needs a helper) - machinery incomplete"
    assert element["machinery_intact"] == "FALSE"
    assert element["sys_wholeness_min"] == "0.667"
    assert element["spans_contigs"] == "FALSE"
    assert element["at_contig_boundary"] == "FALSE"
    assert element["confidence"] == "medium"
    # Phase 3 is not implemented, and the table says so rather than implying the
    # interval is a resolved element boundary.
    assert element["boundary_method"] == "none"
    assert element["attL"] == "NA" and element["attR"] == "NA"

    dropped = [row for row in audit if row["action"] == "dropped"]
    assert len(dropped) == 1
    assert dropped[0]["contig"] == "contig_2"
    assert dropped[0]["reason"] == "cluster_shorter_than_min"
    assert "6765 bp" in dropped[0]["detail"]
    assert "machinery_degraded" in audit_reasons(audit)


# ── Phase 4: the classification table, tested as a pure function ─────────────

def test_classification_table_matches_the_spec():
    """Spec section 8 Phase 4, row by row."""
    assert ci.classify_cluster(True, True, True) == (
        "ice", "predicted self-transmissible")
    assert ci.classify_cluster(True, True, False) == (
        "ime", "mobilisable (needs a helper)")
    assert ci.classify_cluster(True, False, False) == (
        "cime_or_island", "passive")
    assert ci.classify_cluster(False, True, True) == (
        "conjugative_region", "conjugative region, boundaries not established")


def test_classification_never_says_transmissible_without_predicted():
    """Language discipline (spec section 2.6): an ICE is always 'predicted'."""
    _, mobility = ci.classify_cluster(True, True, True)
    assert mobility.startswith("predicted ")
    assert "predicted self-transmissible" in mobility


def test_classification_fallbacks_for_combinations_the_spec_omits():
    """The two combinations the spec's table does not list.

    Integrase + T4SS but no relaxase cannot mobilise anything, so it is grouped
    with the passive island; machinery without an integrase is always reported as
    an unbounded region and never as an ICE.
    """
    assert ci.classify_cluster(True, False, True)[0] == "cime_or_island"
    assert ci.classify_cluster(False, True, False)[0] == "conjugative_region"
    assert ci.classify_cluster(False, False, True)[0] == "conjugative_region"


# ── End to end: one case per class ───────────────────────────────────────────

def test_ice_needs_all_three_anchor_classes(tmp_path):
    """Integrase + relaxase + mating-pair component on one contig -> ICE, high.

    This is the case the whole module exists for: the element is on the
    CHROMOSOME and is still predicted to move itself.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    element = rows[0]
    assert element["mge_class"] == "ice"
    assert element["element_type"] == "ice"          # colocalise.py raises this to tier 6
    assert element["mobility"] == "predicted self-transmissible"
    assert element["anchor_classes"] == "integrase,relaxase,t4cp,t4ss"
    assert element["n_anchor_classes"] == "4"
    assert element["mpf_type"] == "F"
    assert element["mpf_typed_system"] == "TRUE"
    assert element["machinery_intact"] == "TRUE"
    assert element["confidence"] == "high"
    assert element["mge_id"] == "contig_1|ice-50000:65500"
    assert element["length_bp"] == "15501"
    assert element["integrase_products"] == "Phage integrase family protein"
    assert "S1_00010(integrase)" in element["anchor_ids"]
    assert audit_reasons(audit) == set()             # nothing dropped, nothing capped


def test_ime_is_integrase_plus_relaxase_without_mating_pair(tmp_path):
    """Integrase + relaxase + coupling protein, no MPF -> IME, mobilisable.

    The coupling protein must NOT be read as a T4SS: if it were, this would be
    wrongly promoted to a self-transmissible ICE.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBP1"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    element = rows[0]
    assert element["mge_class"] == "ime"
    assert element["element_type"] == "ime"          # colocalise.py raises this to tier 5
    assert element["mobility"] == "mobilisable (needs a helper)"
    assert element["has_t4cp"] == "TRUE"
    assert element["has_t4ss"] == "FALSE"
    assert element["relaxase_type"] == "MOBP1"
    assert element["mpf_type"] == "NA"               # the MOB model has no MPF type
    assert element["mpf_typed_system"] == "FALSE"
    assert element["n_anchor_classes"] == "3"
    assert element["confidence"] == "high"


def test_integrase_only_cluster_is_dropped_not_reported(tmp_path):
    """A lone integrase is not an element.

    Every genome carries several site-specific recombinases. Phase 2 keeps only
    clusters with a relaxase or a mating-pair component, so the integrase on
    contig_2 here is dropped with a stated reason instead of being reported as a
    passive island - which would put a meaningless row in the table for every
    recombinase in the genome. The classifier's island branch is still exercised
    directly in test_classification_table_matches_the_spec.
    """
    contigs = {"contig_1": 200000, "contig_2": 200000}
    cds = SCENE_CDS + [
        gff_cds("contig_2", 10000, 11200, "+", "S1_00500",
                "tyrosine recombinase XerC"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert [row["contig"] for row in rows] == ["contig_1"]
    dropped = [row for row in audit if row["reason"] == "no_conjugation_anchor"]
    assert len(dropped) == 1
    assert dropped[0]["contig"] == "contig_2"
    assert dropped[0]["action"] == "dropped"
    assert "integrase" in dropped[0]["detail"]


def test_machinery_without_an_integrase_is_not_called_an_ice(tmp_path):
    """Relaxase + mating-pair component, no integrase -> conjugative region.

    The spec is explicit: report it, do not call it an ICE. Nothing here says
    the machinery sits in a discrete element with boundaries, so the row must not
    reach colocalise.py as an ICE either.
    """
    cds_without_integrase = [
        line for line in SCENE_CDS if "integrase" not in line
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeT"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeT"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeT"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds_without_integrase)

    return_code, rows, _audit, out_table = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    element = rows[0]
    assert element["mge_class"] == "conjugative_region"
    assert element["element_type"] == "conjugative_region"
    assert element["mobility"] == "conjugative region, boundaries not established"
    assert element["has_integrase"] == "FALSE"
    assert "predicted self-transmissible" not in element["mobility"]

    # And the downstream consumer must not treat it as an ICE either: an
    # element_type it does not recognise is ignored (loudly, in its own audit).
    elements = co.parse_mobile_elements(out_table)
    assert elements[0]["element_type"] is None


def test_transposase_with_an_integrase_domain_is_not_an_integrase_anchor(tmp_path):
    """IS transposases share the rve 'integrase catalytic domain'.

    Counting one as an integrase would promote a plain insertion sequence next to
    a relaxase into an ICE, so a product that also says 'transposase' is
    excluded - and the exclusion is written to the audit file, never silent.
    """
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010",
                "IS3 family transposase, integrase catalytic domain"),
    ] + SCENE_CDS[1:]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert rows[0]["has_integrase"] == "FALSE"
    assert rows[0]["mge_class"] == "conjugative_region"
    assert "integrase_match_is_a_transposase" in audit_reasons(audit)


def test_reca_is_not_mistaken_for_an_integrase():
    """'recombinase RecA' is in every genome and is not a site-specific integrase.

    Only the 'tyrosine recombinase' / 'serine recombinase' phrasings match, so
    the bare word 'recombinase' must not create an anchor.
    """
    features = [
        {"contig": "c1", "start": 1, "end": 100, "strand": "+",
         "feature_id": "g1", "product": "recombinase RecA"},
        {"contig": "c1", "start": 200, "end": 300, "strand": "+",
         "feature_id": "g2", "product": "Recombinase family protein"},
        {"contig": "c1", "start": 400, "end": 500, "strand": "+",
         "feature_id": "g3", "product": "tyrosine recombinase XerD"},
    ]
    anchors, _audit = ci.find_integrase_anchors("S1", features)
    assert [anchor["feature_id"] for anchor in anchors] == ["g3"]


# ── Degradation, length filters, contig honesty ──────────────────────────────

def test_low_wholeness_downgrades_the_mobility_wording(tmp_path):
    """An incomplete system is reported as incomplete, in the sentence itself.

    Decayed conjugative elements are common and are where naive tools overcall,
    so a system below 0.7 wholeness sets machinery_intact=FALSE, appends
    '- machinery incomplete' to the mobility statement and caps confidence.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", sys_wholeness="0.333"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", sys_wholeness="0.333"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", sys_wholeness="0.333"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    element = rows[0]
    assert element["mge_class"] == "ice"
    assert element["mobility"] == "predicted self-transmissible - machinery incomplete"
    assert element["machinery_intact"] == "FALSE"
    assert element["degraded_reason"] == "low_system_wholeness"
    assert element["confidence"] == "medium"
    assert "machinery_degraded" in audit_reasons(audit)
    assert "machinery_not_intact" in audit_reasons(audit)


def test_decayed_model_and_truncated_relaxase_also_count_as_degraded(tmp_path):
    """The other two ways of being broken: CONJscan's dCONJ_* decayed system
    models, and a relaxase whose HMM alignment covers less than 70% of the
    profile (a fragment of a relaxase nicks nothing)."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/dCONJ_typeF",
                     hit_profile_cov="0.21"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/dCONJ_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/dCONJ_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    element = rows[0]
    assert element["machinery_intact"] == "FALSE"
    assert "decayed_system_model" in element["degraded_reason"]
    assert "truncated_core_hit" in element["degraded_reason"]
    assert element["mobility"].endswith(" - machinery incomplete")


def test_cluster_shorter_than_minimum_is_dropped_with_a_reason(tmp_path):
    """A 2.5 kb span of machinery is a relic, not an integrative element."""
    contigs = {"contig_1": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 51500, 52500, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert rows == []                                  # header-only table, no rows
    dropped = [row for row in audit if row["action"] == "dropped"]
    assert len(dropped) == 1
    assert dropped[0]["reason"] == "cluster_shorter_than_min"
    assert "2501 bp" in dropped[0]["detail"]           # the measured span is reported
    assert dropped[0]["start"] == "50000" and dropped[0]["end"] == "52500"


def test_cluster_longer_than_maximum_is_dropped_with_a_reason(tmp_path):
    """A runaway chain of anchors is not one element.

    Single-linkage clustering can chain anchors much further than the window, so
    --max-element-bp is the guard; the drop records the measured span.
    """
    contigs = {"contig_1": 400000}
    cds = [
        gff_cds("contig_1", 10000, 11000, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 20000, 21000, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
        gff_cds("contig_1", 30000, 31000, "+", "S1_00016", "conjugal transfer protein TraB"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF"),
        conjscan_row(hit_id="S1_00016", gene_name="T4SS_virb4"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--max-element-bp", "15000"]
    )

    assert return_code == 0
    assert rows == []
    dropped = [row for row in audit if row["action"] == "dropped"]
    assert len(dropped) == 1
    assert dropped[0]["reason"] == "cluster_longer_than_max"
    assert "21001 bp" in dropped[0]["detail"]


def test_system_spanning_two_contigs_is_capped_at_low_confidence(tmp_path):
    """MacSyFinder sees a draft assembly as one ordered replicon.

    So it can join genes on either side of a contig break into one 'system' that
    does not exist. Anything spanning contigs is capped at LOW regardless of how
    good the rest of the evidence looks - the spec's rule, applied without
    exception.
    """
    contigs = {"contig_1": 200000, "contig_2": 200000}
    cds = SCENE_CDS + [
        gff_cds("contig_2", 80000, 81200, "+", "S1_00600", "Phage integrase family protein"),
        gff_cds("contig_2", 85000, 86600, "+", "S1_00610", "TrwC relaxase domain-containing protein"),
        gff_cds("contig_2", 90000, 92000, "+", "S1_00620", "conjugal transfer protein TraB"),
    ]
    # One sys_id, hits on BOTH contigs - exactly the artefact described above.
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id="S1_00610", gene_name="T4SS_MOBP1",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id="S1_00620", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 2
    for element in rows:
        assert element["spans_contigs"] == "TRUE"
        assert element["confidence"] == "low"         # capped regardless of the rest
    assert "spans_contigs" in audit_reasons(audit)


def test_element_at_a_contig_end_is_flagged_and_capped(tmp_path):
    """Machinery running into a contig end means the element is truncated.

    The rest of it is in sequence the assembler could not place, so the call is
    flagged and capped at low confidence.
    """
    contigs = {"contig_1": 66000}                      # ends 500 bp after the last anchor
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, SCENE_CDS)

    _return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    element = rows[0]
    assert element["contig_length"] == "66000"
    assert element["dist_to_contig_end"] == "500"
    assert element["at_contig_boundary"] == "TRUE"
    assert element["confidence"] == "low"
    assert "at_contig_boundary" in audit_reasons(audit)


def test_unknown_contig_length_caps_confidence_at_medium(tmp_path):
    """No length, no boundary check - and we say so instead of assuming safety."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    # A GFF with no ##sequence-region lines and no contig-length table given.
    gff = write_gff(tmp_path / "S1.gff3", {}, SCENE_CDS)

    _return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    element = rows[0]
    assert element["contig_length"] == "NA"
    assert element["at_contig_boundary"] == "NA"
    assert element["confidence"] == "medium"
    assert "contig_length_unknown" in audit_reasons(audit)


def test_contig_length_table_is_used_when_given(tmp_path):
    """The supplied contig-length table wins over the GFF's own headers, so this
    script measures against exactly the same contigs as the ISEScan leg."""
    lengths_file = tmp_path / "contig_lengths.tsv"
    with open(str(lengths_file), "w", encoding="utf-8") as handle:
        handle.write("contig\tlength\n")               # a header line is tolerated
        handle.write("contig_1\t500000\n")

    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff, contig_lengths=lengths_file
    )

    assert rows[0]["contig_length"] == "500000"        # not the GFF's 200000
    # The IME's last anchor is the coupling protein at 58700.
    assert rows[0]["dist_to_contig_end"] == "441300"


# ── Graceful degradation: the common case is "nothing found" ─────────────────

def test_missing_conjscan_file_writes_an_empty_table_and_exits_zero(tmp_path):
    """CONJscan is opt-in and most isolates have no conjugative system.

    A missing file must produce a complete, empty table plus an audit line
    explaining that the ICE check was not made - 'not looked at' is a different
    statement from 'looked at and found nothing'.
    """
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, out_table = run_main(
        tmp_path, conjscan=str(tmp_path / "does_not_exist.tsv"), gff=gff
    )

    assert return_code == 0
    assert rows == []
    header, _ = read_tsv(out_table)
    assert header == ci.OUTPUT_COLUMNS                 # a well-formed empty table
    assert "conjscan_output_missing" in audit_reasons(audit)


def test_conjscan_argument_omitted_entirely(tmp_path):
    """--conjscan-tsv is optional; leaving it off behaves like a missing file."""
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)
    return_code, rows, audit, _ = run_main(tmp_path, conjscan=None, gff=gff)
    assert return_code == 0
    assert rows == []
    assert "conjscan_output_missing" in audit_reasons(audit)


def test_empty_conjscan_file(tmp_path):
    """A zero-byte file is 'no systems found', not a crash."""
    empty = tmp_path / "best_solution.tsv"
    open(str(empty), "w").close()
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, empty, gff)

    assert return_code == 0
    assert rows == []
    assert "conjscan_found_no_systems" in audit_reasons(audit)


def test_header_only_conjscan_file(tmp_path):
    """Banner + header + no data rows: CONJscan ran and found nothing.

    This is a genuine negative result and the audit says so explicitly.
    """
    header_only = write_conjscan(tmp_path / "best_solution.tsv", [])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, header_only, gff)

    assert return_code == 0
    assert rows == []
    assert "conjscan_found_no_systems" in audit_reasons(audit)


def test_comment_only_conjscan_file(tmp_path):
    """Some runs write only the '#' banner. Same answer: no systems."""
    banner_only = write_conjscan(tmp_path / "best_solution.tsv", [], with_header=False)
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, banner_only, gff)

    assert return_code == 0
    assert rows == []
    assert "conjscan_found_no_systems" in audit_reasons(audit)


def test_missing_bakta_gff_writes_an_empty_table_and_exits_zero(tmp_path):
    """Without the annotation a protein hit has no coordinates, so nothing can be
    clustered - but the rule still produces readable files and exits 0."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
    ])

    return_code, rows, audit, out_table = run_main(
        tmp_path, conjscan, gff=str(tmp_path / "no_such.gff3")
    )

    assert return_code == 0
    assert rows == []
    header, _ = read_tsv(out_table)
    assert header == ci.OUTPUT_COLUMNS
    assert "bakta_gff_missing" in audit_reasons(audit)


def test_hit_id_absent_from_the_annotation_is_audited(tmp_path):
    """A protein ID that is not in the GFF cannot be placed.

    That usually means CONJscan and Bakta were run on different annotations, so
    the hit is skipped with a reason rather than guessed at.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="NOT_IN_GFF_00001", gene_name="T4SS_MOBF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert rows == []
    assert "hit_id_not_in_annotation" in audit_reasons(audit)
    assert "no_machinery_hit_could_be_placed" in audit_reasons(audit)


def test_unrecognised_conjscan_header_stops_the_rule(tmp_path):
    """A tool-version change must NOT be reported as 'no conjugation machinery'.

    Silently returning a wrong negative is worse than failing, so an unparseable
    header is the one case that exits non-zero.
    """
    bad = tmp_path / "best_solution.tsv"
    with open(str(bad), "w", encoding="utf-8") as handle:
        handle.write("# macsyfinder 9.9.9\n")
        handle.write("protein\tsomething_else\n")
        handle.write("S1_00015\tvalue\n")
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    out_table = str(tmp_path / "table.tsv")
    out_audit = str(tmp_path / "audit.tsv")
    return_code = ci.main([
        "--sample", "S1", "--conjscan-tsv", str(bad), "--bakta-gff", gff,
        "--out-table", out_table, "--out-audit", out_audit,
    ])
    assert return_code == 1


def test_directory_instead_of_file_is_resolved(tmp_path):
    """MacSyFinder writes best_solution.tsv into its --out-dir, so a directory is
    accepted in place of the file itself."""
    out_dir = tmp_path / "conjscan_S1"
    os.makedirs(str(out_dir), exist_ok=True)
    write_conjscan(out_dir / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, _audit, _ = run_main(tmp_path, out_dir, gff)

    assert return_code == 0
    assert len(rows) == 1
    assert rows[0]["mge_class"] == "ime"


# ── The contract with the next script in the chain ───────────────────────────

def test_output_is_readable_by_colocalise(tmp_path):
    """The table must be consumable by colocalise.py's --is-table parser as is.

    That parser matches columns BY NAME, so this test is what catches a rename
    that would otherwise make the co-localisation step silently see no elements.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, _rows, _audit, out_table = run_main(tmp_path, conjscan, gff)

    elements = co.parse_mobile_elements(out_table)
    assert len(elements) == 1
    element = elements[0]
    assert element["contig"] == "contig_1"
    assert element["start"] == 50000
    assert element["end"] == 65500
    assert element["element_type"] == "ice"            # -> mobility tier 6
    assert element["id"] == "contig_1|ice-50000:65500"


def test_element_type_values_are_the_ones_colocalise_knows():
    """ICE and IME must map onto colocalise.py's vocabulary; the two non-mobile
    classes must deliberately NOT, so they cannot raise a gene's mobility tier."""
    assert co.ELEMENT_TYPE_SYNONYMS[ci.ELEMENT_TYPE_FOR_CLASS["ice"]] == "ice"
    assert co.ELEMENT_TYPE_SYNONYMS[ci.ELEMENT_TYPE_FOR_CLASS["ime"]] == "ime"
    assert ci.ELEMENT_TYPE_FOR_CLASS["cime_or_island"] not in co.ELEMENT_TYPE_SYNONYMS
    assert ci.ELEMENT_TYPE_FOR_CLASS["conjugative_region"] not in co.ELEMENT_TYPE_SYNONYMS


def test_no_output_column_collides_with_colocalise_field_names():
    """colocalise.py reads ISEScan's `type` as completeness and `family` as the IS
    family. Naming one of our columns that way would be silently misread, so the
    collision is asserted against here rather than left to a code review."""
    for forbidden in ("type", "complete", "completeness", "is_complete",
                      "family", "cluster", "name", "id"):
        assert forbidden not in ci.OUTPUT_COLUMNS


# ── Clustering and parsing details ───────────────────────────────────────────

def test_anchors_further_apart_than_the_window_form_separate_clusters(tmp_path):
    """Two machinery blocks 40 kb apart are two findings, not one element."""
    contigs = {"contig_1": 300000}
    cds = SCENE_CDS + [
        gff_cds("contig_1", 150000, 151200, "+", "S1_00700", "Phage integrase family protein"),
        gff_cds("contig_1", 155000, 156600, "+", "S1_00710", "TrwC relaxase domain-containing protein"),
        gff_cds("contig_1", 160000, 162000, "+", "S1_00720", "conjugal transfer protein TraB"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF", sys_id="S1_MOB_1",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4", sys_id="S1_MOB_1",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id="S1_00710", gene_name="T4SS_MOBP1", sys_id="S1_MOB_2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id="S1_00720", gene_name="T4SS_virb4", sys_id="S1_MOB_2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    assert len(rows) == 2
    assert [row["start"] for row in rows] == ["50000", "150000"]   # sorted along the contig
    assert all(row["mge_class"] == "ice" for row in rows)


def test_window_is_configurable(tmp_path):
    """Shrinking --window-bp splits a cluster; the threshold is a convention, not
    biology, so it has to be adjustable."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    # 3 kb window: the integrase (ends 51200) is 3799 bp from the relaxase and the
    # VirB4 (starts 63000) is 4299 bp from the coupling protein, so the single ICE
    # falls apart into two machinery fragments and neither can be called an ICE.
    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--window-bp", "3000", "--min-element-bp", "1000"]
    )
    assert len(rows) == 2
    assert all(row["has_integrase"] == "FALSE" for row in rows)
    assert all(row["mge_class"] == "conjugative_region" for row in rows)


def test_gff_parsing_reads_coordinates_products_and_lengths(tmp_path):
    """The GFF3 loader: both ID and locus_tag are join keys, percent-encoded
    product text is decoded, sequence-region gives the contig length, and the
    trailing ##FASTA block is not read as features."""
    gff = write_gff(tmp_path / "S1.gff3", {"contig_1": 200000}, [
        gff_cds("contig_1", 100, 400, "-", "S1_00001",
                "Type IV secretory pathway, VirD4 component"),
    ])
    features_by_id, cds_features, lengths = ci.parse_bakta_gff(gff)

    assert lengths == {"contig_1": 200000}
    assert len(cds_features) == 1                       # the FASTA block was skipped
    feature = features_by_id["S1_00001"]
    assert (feature["contig"], feature["start"], feature["end"], feature["strand"]) == (
        "contig_1", 100, 400, "-")
    # %2C decoded back to a comma, so the product reads as a human would write it.
    assert feature["product"] == "Type IV secretory pathway, VirD4 component"


def test_summary_line_names_each_class(tmp_path, capsys):
    """The stdout line has to distinguish an ICE from a passive island - '2
    elements found' would be a misleading thing to read in a log."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    run_main(tmp_path, conjscan, gff)
    printed = capsys.readouterr().out

    assert "Sample S1:" in printed
    assert "1 ICE (predicted self-transmissible)" in printed
    assert "0 IME" in printed
    assert "confidence high 1" in printed
