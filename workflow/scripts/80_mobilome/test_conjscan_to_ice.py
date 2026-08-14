"""Unit tests for workflow/scripts/80_mobilome/conjscan_to_ice.py.

conjscan_to_ice answers the question no plasmid caller can: an ICE sits in the
chromosome and still carries its own conjugation machinery, so a resistance gene
inside one is chromosomal AND predicted transferable. That is the strongest claim
the module makes — mobility tier 6 — and most of what follows is about refusing
to make it on thin evidence.

No tools, no databases, no MacSyFinder run: every test builds a small
CONJscan-shaped best_solution.tsv and a small Bakta-shaped GFF3 by hand and
pushes them through the script, so the whole anchors → clusters → classification
chain is verifiable in a plain Python environment.

The fixtures follow the real contracts:
  * CONJscan/MacSyFinder 2.1.6 best_solution.tsv — three '#' banner lines, then
    the verbatim 22-column header, then one row per machinery gene, with blank
    lines between systems;
  * Bakta 1.12.0 GFF3 — '##sequence-region' headers, Pyrodigal CDS rows carrying
    ID= and locus_tag= with the same value, percent-encoded product text, and a
    '##FASTA' section at the end.

Most tests pass NO genome FASTA, because they are about classification rather
than boundaries. Without sequence the att search cannot run, so those elements
come back with boundary_method='none' and the machinery span as their interval —
that is the expected result there, not a gap. The tests that ARE about Phase 3
plant a repeat pair in a synthetic contig and hand it over with --genome.

Two tests go further and use the REAL output of sample 386 (an Arthrobacter
isolate): the saved fixture testdata/conjscan_386_real_best_solution.tsv, and
the five CDS lines quoted verbatim out of that sample's Bakta GFF3. Passing
those means the parser agrees with the tools, not just with itself.

There is also a contract test that feeds this script's output table straight into
colocalise.py's own parser, which is what consumes it in the workflow.

Run:
    python -m pytest workflow/scripts/80_mobilome/test_conjscan_to_ice.py -q
"""

import csv
import os
import random
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


def gff_trna(contig, start, end, strand="+", name="tRNA-Gly(gcc)"):
    """Build one Bakta-shaped tRNA line.

    Needed because only a tRNA-ANCHORED att pair is allowed to widen an element
    (a de novo repeat is reported but not applied — see
    refine_candidate_boundaries), so any fixture that tests widening has to plant
    a real tRNA for the probe to come from.
    """
    attributes = f"ID={contig}_{start};Name={name};product={name}"
    return f"{contig}\ttRNAscan-SE\ttRNA\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"


def plant_att_pair(sequence, motif, left_start, right_start):
    """Put the same short motif at two 1-based positions in a sequence.

    Used to build att fixtures: the left copy is the 3' end of a planted tRNA
    (so it looks like a reconstituted integration site) and the right copy sits
    in ordinary sequence, which is the one-in-one-out arrangement real
    integration leaves behind.
    """
    for start in (left_start, right_start):
        index = start - 1
        sequence = sequence[:index] + motif + sequence[index + len(motif):]
    return sequence


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
    """The set of reason tokens present in an audit file — what most tests assert
    on, usually as a subset check so an extra unrelated reason does not fail."""
    return {row["reason"] for row in audit_rows}


# A standard scene used by most end-to-end tests: one 200 kb contig carrying an
# integrase, a relaxase, a coupling protein and a VirB4, spaced a few kb apart so
# they all fall inside the default 15 kb clustering window. The machinery genes
# alone (relaxase 55000 → VirB4 65500) already span more than the default 8 kb
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
    proteins — and NOT mating-pair components, which is what keeps sample 386 at
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
        it is an IME — mobilisable with a helper, not self-transmissible;
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
    # No genome FASTA was passed, so the att search never ran — Phase 3 needs
    # sequence. The table says so rather than implying the machinery span is a
    # resolved element boundary.
    assert element["boundary_method"] == "none"
    assert element["attL"] == "NA" and element["attR"] == "NA"

    dropped = [row for row in audit if row["action"] == "dropped"]
    assert len(dropped) == 1
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
    """Integrase + relaxase + mating-pair component on one contig → ICE, high.

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
    # High: every anchor class present, machinery intact, one contig. No genome was
    # given so the element's ENDS are unknown, but by default that is reported in
    # boundary_method rather than folded into the confidence — the two answer
    # different questions. See test_strict_mode_requires_a_trna_boundary_for_high.
    assert element["confidence"] == "high"
    assert element["mge_id"] == "contig_1|ice-50000:65500"
    assert element["length_bp"] == "15501"
    assert element["integrase_products"] == "Phage integrase family protein"
    assert "S1_00010(integrase)" in element["anchor_ids"]
    # Nothing was dropped. Three reasons may appear, and none of them is a
    # rejection:
    #   assembly_contiguity   one row per run on EVERY sample, recording the
    #                         contig count and N50 the calls came off;
    #   no_genome_for_att_search
    #                         these unit tests deliberately pass no genome FASTA,
    #                         so Phase 3 could not look for the element's real
    #                         ends (the real rule always supplies one). By
    #                         default that is reported in boundary_method and
    #                         does NOT lower the confidence — see
    #                         test_strict_mode_requires_a_trna_boundary_for_high;
    #   integrase_attached_beyond_cluster_window
    #                         the integrase in this fixture sits outside the
    #                         machinery clustering window and is attached by the
    #                         wider integrase search, which is the point of that
    #                         step.
    assert audit_reasons(audit) <= {
        "assembly_contiguity",
        "no_genome_for_att_search",
        "integrase_attached_beyond_cluster_window",
    }


def test_ime_is_integrase_plus_relaxase_without_mating_pair(tmp_path):
    """Integrase + relaxase + coupling protein, no MPF → IME, mobilisable.

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


def test_accessory_virb4_in_a_mob_system_does_not_make_an_ice(tmp_path):
    """A lone VirB4 inside a relaxase-only MOB system is not a mating bridge.

    CONJscan's `MOB` model describes a relaxase-only system — DNA that another
    element's machinery can pick up — and it lists VirB4 as an ACCESSORY gene.
    So a MOB system can quite legitimately contain one VirB4 hit while the cell
    has no mating-pair apparatus at all.

    Counting that hit as an MPF would raise this element from tier 5 (mobilisable,
    needs a helper) to tier 6 (predicted self-transmissible) on the strength of a
    single accessory gene, which is the worst overcall this module could make.
    The hit is not hidden — it is still reported in has_t4ss — but the class must
    come from the SYSTEM's own type, and the audit file must say so.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/MOB"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/MOB",
                     hit_gene_ref="T4SS_t4cp1", hit_status="accessory"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/MOB",
                     hit_gene_ref="T4SS_virb4", hit_status="accessory"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    element = rows[0]
    assert element["mge_class"] == "ime"
    assert element["element_type"] == "ime"           # colocalise.py → tier 5, not 6
    assert element["mobility"] == "mobilisable (needs a helper)"
    assert "self-transmissible" not in element["mobility"]
    # The marker is reported, not suppressed; it simply does not carry the call.
    assert element["has_t4ss"] == "TRUE"
    assert element["mpf_typed_system"] == "FALSE"
    assert element["mpf_type"] == "NA"                # the MOB model has no MPF type
    # And the reason it was not counted is written down, per the audit rule.
    assert "mpf_marker_without_typed_system" in audit_reasons(audit)


def test_the_same_virb4_under_a_typed_t4ss_model_does_make_an_ice(tmp_path):
    """The other side of the test above: same genes, same coordinates, but this
    time CONJscan called a typed T4SS system, so the mating-pair apparatus IS
    evidenced by the system itself and the ICE call is earned."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    element = rows[0]
    assert element["mge_class"] == "ice"
    assert element["mobility"] == "predicted self-transmissible"
    assert element["mpf_typed_system"] == "TRUE"
    assert "mpf_marker_without_typed_system" not in audit_reasons(audit)


def test_integrase_only_cluster_is_dropped_not_reported(tmp_path):
    """A lone integrase is not an element.

    Every genome carries several site-specific recombinases. Phase 2 keeps only
    clusters with a relaxase or a mating-pair component, so the integrase on
    contig_2 here is dropped with a stated reason instead of being reported as a
    passive island — which would put a meaningless row in the table for every
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
    # Integrase-only clusters no longer FORM: machinery is clustered on its own
    # and integrases are attached afterwards, so an integrase with no machinery
    # anchors nothing. The decision is still audited, under a reason describing
    # what actually happens now.
    orphan = [row for row in audit
              if row["reason"] == "integrase_without_conjugation_machinery"]
    assert len(orphan) == 1
    # A run-level summary, so contig is NA; the contigs are named in the detail.
    assert "contig_2" in orphan[0]["detail"]
    assert orphan[0]["action"] == "not_applicable"
    assert "integrase" in orphan[0]["detail"]


def test_machinery_without_an_integrase_is_not_called_an_ice(tmp_path):
    """Relaxase + mating-pair component, no integrase → conjugative region.

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

    # And the downstream consumer must not treat it as an ICE either. It DOES
    # recognise the type now — being unrecognised used to mean the element was
    # dropped from every test in colocalise, so an AMR gene sitting inside a
    # predicted conjugative region came out as "intrinsic candidate" at high
    # confidence. Recognised, but never tier-raising, is the correct handling.
    elements = co.parse_mobile_elements(out_table)
    assert elements[0]["element_type"] == "conjugative_region"
    assert elements[0]["element_type"] in co.CONTEXT_ONLY_ELEMENT_TYPES


def test_transposase_with_an_integrase_domain_is_not_an_integrase_anchor(tmp_path):
    """IS transposases share the rve 'integrase catalytic domain'.

    Counting one as an integrase would promote a plain insertion sequence next to
    a relaxase into an ICE, so a product that also says 'transposase' is
    excluded — and the exclusion is written to the audit file, never silent.
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


def test_hit_is_virb4_knows_the_exchangeable_name_but_not_the_f_type_traU():
    """VirB4 can reach us under two names, and only one look-alike must be refused.

    Every CONJscan model lists `T4SS_I_traU` as an exchangeable profile for
    `T4SS_virb4`, so a VirB4 hit can be reported under either name. `T4SS_F_traU`
    is a different gene entirely — an F-type mating-pair component in its own
    right — and must not be mistaken for the ATPase.
    """
    assert ci.hit_is_virb4({"gene_name": "T4SS_virb4", "hit_gene_ref": "T4SS_virb4"})
    # Found through the exchangeable profile; the model's own gene is alongside it.
    assert ci.hit_is_virb4({"gene_name": "T4SS_I_traU", "hit_gene_ref": "T4SS_virb4"})
    # The same hit from an output that carries no hit_gene_ref column at all.
    assert ci.hit_is_virb4({"gene_name": "T4SS_I_traU"})
    assert not ci.hit_is_virb4({"gene_name": "T4SS_F_traU",
                                "hit_gene_ref": "T4SS_F_traU"})
    assert not ci.hit_is_virb4({"gene_name": "T4SS_MOBF", "hit_gene_ref": "T4SS_MOBB"})


def test_real_output_shows_gene_name_differing_from_the_model_gene():
    """Ground truth for why `hit_gene_ref` is read at all.

    In sample 386's real best_solution.tsv the relaxases are reported as
    T4SS_MOBP1 / T4SS_MOBF while the model's own gene is T4SS_MOBB: MacSyFinder
    writes the profile that ACTUALLY matched in gene_name. The same mechanism
    turns a VirB4 into a `T4SS_I_traU` row, which is what the truncation check
    has to survive.
    """
    hits = ci.read_conjscan_hits(REAL_CONJSCAN_FIXTURE)
    assert [hit["hit_gene_ref"] for hit in hits] == [
        "T4SS_MOBB", "T4SS_t4cp1", "T4SS_t4cp1", "T4SS_MOBB"
    ]
    assert hits[0]["gene_name"] == "T4SS_MOBP1"       # not the model's own gene


def test_truncated_virb4_under_its_exchangeable_name_is_still_flagged(tmp_path):
    """A fragment of VirB4 builds no mating bridge, whatever it was called.

    VirB4 is one of the two components the truncation check tests (the other is
    the relaxase), because both have to WORK for transfer to happen. Here the
    ATPase of a type I system was found through its exchangeable `T4SS_I_traU`
    profile and aligns over only 30% of the HMM. The class is still ICE — a typed
    T4SS system was called — but the machinery must be reported as degraded, in
    the mobility sentence itself, not quietly as intact.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeI"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeI"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_I_traU",
                     hit_gene_ref="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeI",
                     hit_profile_cov="0.30"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    element = rows[0]
    assert element["mge_class"] == "ice"
    assert element["machinery_intact"] == "FALSE"
    assert "truncated_core_hit" in element["degraded_reason"]
    assert element["mobility"] == "predicted self-transmissible - machinery incomplete"
    assert "machinery_degraded" in audit_reasons(audit)


def test_a_lone_relaxase_of_a_typed_system_does_not_make_an_ice(tmp_path):
    """A typed system elsewhere on the replicon must not make THIS cluster an ICE.

    MacSyFinder LONER genes may sit anywhere on the replicon, so a relaxase
    declared a loner of a typed T4SS model carries that model's type letter with
    it. Asking only "is the contributing system typed?" therefore let a cluster
    holding exactly ONE integrase and ONE relaxase — has_t4cp FALSE, has_t4ss
    FALSE, the textbook IME signature — be reported as `ice`, "predicted
    self-transmissible".

    That is the worst overcall this script can make: tier 6 is the answer a
    regulator reads. It happened on NC_013929 in the Phase 7 benchmark, where the
    row contradicted itself — missing_components said "coupling protein,
    mating-pair apparatus" beside the tier-6 claim — and the caller's own audit
    had already refused to merge in the real apparatus, 1.2 Mb away.

    The apparatus must be HERE, not merely somewhere on the replicon.
    """
    contigs = {"contig_1": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 55000, 56600, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
    ]
    # The relaxase is the only hit, and it is attributed to a TYPED model — the
    # loner case. No mating-pair gene is in the cluster.
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--min-element-bp", "1000"])

    assert len(rows) == 1
    row = rows[0]
    assert row["has_t4ss"] == "FALSE"
    # The type letter is still REPORTED — we do not hide what CONJscan said...
    assert row["mpf_typed_system"] == "TRUE"
    # ...but it no longer buys a self-transmissibility claim.
    assert row["mge_class"] == "ime"
    assert row["mobility"].startswith("mobilisable")


def test_a_mating_pair_gene_in_the_cluster_still_makes_an_ice(tmp_path):
    """The complement of the test above: real in-cluster apparatus still means ICE.

    Without this, "require the apparatus in the cluster" could be satisfied by
    simply never calling an ICE again.
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
    _return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    assert len(rows) == 1
    assert rows[0]["has_t4ss"] == "TRUE"
    assert rows[0]["mge_class"] == "ice"
    assert rows[0]["mobility"] == "predicted self-transmissible"


def test_ice_architecture_cluster_shorter_than_minimum_is_dropped(tmp_path):
    """A 2.5 kb span of ICE machinery is a relic, not an integrative element.

    ICE machinery is a ~20-gene mating-pair operon, so a couple of kilobases of
    it is a fragment. The IME floor below does NOT apply here, because this
    cluster carries a mating-pair component and so is not IME-architecture.
    """
    contigs = {"contig_1": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 51500, 52500, "+", "S1_00017", "conjugal transfer protein TraB"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00017", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert rows == []                                  # header-only table, no rows
    dropped = [row for row in audit if row["action"] == "dropped"]
    assert len(dropped) == 1
    assert dropped[0]["reason"] == "cluster_shorter_than_min"
    assert "--min-element-bp" in dropped[0]["detail"]  # the ICE floor, not the IME one
    assert "2501 bp" in dropped[0]["detail"]           # the measured span is reported
    assert dropped[0]["start"] == "50000" and dropped[0]["end"] == "52500"


def test_ime_architecture_cluster_survives_the_lower_floor(tmp_path):
    """The same 2.5 kb span, but IME-architecture, is KEPT.

    Measured on the Phase 7 IME pilot: the size floor is applied to the anchor
    cluster's SPAN, and that span scales with the NUMBER of machinery genes. An
    IME carries a relaxase and an integrase — two genes, 1–6 kb — where an ICE
    carries a twenty-gene operon. An 8,000 bp floor therefore selected for ICEs
    by construction: six of the twelve curated IMEs were clustered and
    classified correctly, then dropped for size, Tn4451 at a span of 1,266 bp.

    This is the same scene as the test above with the mating-pair gene swapped
    for a relaxase — which is precisely the difference between the two classes.
    """
    contigs = {"contig_1": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 51500, 52500, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/MOB"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    assert return_code == 0
    assert len(rows) == 1
    assert rows[0]["mge_class"] == "ime"
    assert rows[0]["mobility"].startswith("mobilisable")
    # Still capped below high: an IME has only two anchor classes by definition.
    assert rows[0]["confidence"] != "high"


def test_ime_floor_can_be_raised_back_to_the_ice_floor(tmp_path):
    """--min-ime-element-bp is a knob, and the empirical cut is not load-bearing.

    The default of 2,000 bp was fitted to a handful of Phase 7 observations, not
    taken from a published bound, so a reader who disagrees must be able to turn
    it off. Setting it to the ICE floor restores the old behaviour exactly.
    """
    contigs = {"contig_1": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 51500, 52500, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/MOB"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--min-ime-element-bp", "8000"])

    assert rows == []
    dropped = [row for row in audit if row["action"] == "dropped"]
    assert dropped and dropped[0]["reason"] == "cluster_shorter_than_min"


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
    good the rest of the evidence looks — the spec's rule, applied without
    exception.
    """
    contigs = {"contig_1": 200000, "contig_2": 200000}
    cds = SCENE_CDS + [
        gff_cds("contig_2", 80000, 81200, "+", "S1_00600", "Phage integrase family protein"),
        gff_cds("contig_2", 85000, 86600, "+", "S1_00610", "TrwC relaxase domain-containing protein"),
        gff_cds("contig_2", 90000, 92000, "+", "S1_00620", "conjugal transfer protein TraB"),
    ]
    # One sys_id, hits on BOTH contigs — exactly the artefact described above.
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
    """No length, no boundary check — and we say so instead of assuming safety."""
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
    explaining that the ICE check was not made — 'not looked at' is a different
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
    clustered — but the rule still produces readable files and exits 0."""
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
    assert element["element_type"] == "ice"            # → mobility tier 6
    assert element["id"] == "contig_1|ice-50000:65500"


def test_element_type_values_are_the_ones_colocalise_knows():
    """ICE and IME must map onto colocalise.py's vocabulary; the two non-mobile
    classes must deliberately NOT, so they cannot raise a gene's mobility tier."""
    assert co.ELEMENT_TYPE_SYNONYMS[ci.ELEMENT_TYPE_FOR_CLASS["ice"]] == "ice"
    assert co.ELEMENT_TYPE_SYNONYMS[ci.ELEMENT_TYPE_FOR_CLASS["ime"]] == "ime"
    # The two non-mobile classes ARE recognised by colocalise (so they can set
    # context and cap confidence) but must never map onto "ice" or "ime", which
    # are the only values that raise the tier to 5 or 6.
    for non_mobile_class in ("cime_or_island", "conjugative_region"):
        element_type = ci.ELEMENT_TYPE_FOR_CLASS[non_mobile_class]
        assert co.ELEMENT_TYPE_SYNONYMS[element_type] not in {"ice", "ime"}
        assert co.ELEMENT_TYPE_SYNONYMS[element_type] in co.CONTEXT_ONLY_ELEMENT_TYPES


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


def write_hmmer_extract(directory, profile, hit_ids, profile_coverage=0.8):
    """Write one MacSyFinder {profile}.res_hmm_extract file.

    Real shape, taken from a CONJscan run: four '#' comment lines (the last one
    being the column names) and then one tab-separated row per hit.
    """
    os.makedirs(str(directory), exist_ok=True)
    path = os.path.join(str(directory), f"{profile}.res_hmm_extract")
    with open(path, "w") as handle:
        handle.write(f"# gene: {profile} extract\n")
        handle.write("# profile length= 500\n")
        handle.write("# i_evalue threshold= 0.001\n")
        handle.write("# hit_id\treplicon_name\tposition_hit\thit_sequence_length\t"
                     "gene_name\ti_eval\tscore\tprofile_coverage\t"
                     "sequence_coverage\tbegin\tend\n")
        for position, hit_id in enumerate(hit_ids, start=1):
            handle.write(
                f"{hit_id}\tcontig_1\t{position}\t500\t{profile}\t"
                f"1.0e-100\t300.0\t{profile_coverage}\t0.9\t1\t490\n"
            )
    return path


def test_profile_hits_recover_an_element_when_no_system_was_assembled(tmp_path):
    """Machinery that fails MacSyFinder's quorum is reported, weakly, not lost.

    The case this comes from is ICEVflInd1 on the Phase 7 benchmark: a genuine
    114 kb SXT/R391-family ICE where CONJscan hit twelve mating-pair profiles, a
    coupling protein and VirB4 — but no relaxase. Every conjugative model
    requires a relaxase, so no system was assembled and the module reported
    nothing whatsoever. A real element vanished on one missing component.

    It is now surfaced as the weakest evidence tier, and must say so in three
    places at once: evidence_level, missing_components, and confidence.
    """
    # best_solution.tsv exists but holds only MacSyFinder's "no system" banner.
    conjscan = tmp_path / "best_solution.tsv"
    conjscan.write_text("# macsyfinder 2.1.6\n# No System found\n")

    hmmer_dir = tmp_path / "hmmer_results"
    write_hmmer_extract(hmmer_dir, "T4SS_t4cp2", [T4CP_HIT])
    write_hmmer_extract(hmmer_dir, "T4SS_virb4", [VIRB4_HIT])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--conjscan-hmmer-dir", str(hmmer_dir), "--min-element-bp", "1000"],
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["evidence_level"] == "profile_hits_only"
    # The relaxase is the component that was absent, and naming it is the whole
    # point: "passive" must be readable as "no relaxase found", not as a positive
    # claim about the element.
    assert "relaxase" in row["missing_components"]
    assert row["has_relaxase"] == "FALSE"
    assert row["mobility"] == "passive"
    # Never presented as though a system had been called.
    assert row["confidence"] == "low"
    assert "no_system_using_profile_hits" in audit_reasons(audit)


def test_profile_hit_confidence_cap_survives_the_boundary_pass(tmp_path):
    """The low cap must not be lost when confidence is settled after Phase 3.

    finalise_confidence rebuilds every confidence from scratch once the att
    search has run, so a cap applied earlier and not restated there is silently
    thrown away. That is exactly what happened on the first cut of this feature:
    the fallback element came out 'high'. This pins it.
    """
    conjscan = tmp_path / "best_solution.tsv"
    conjscan.write_text("# No System found\n")
    hmmer_dir = tmp_path / "hmmer_results"
    write_hmmer_extract(hmmer_dir, "T4SS_MOBF", [RELAXASE_HIT])
    write_hmmer_extract(hmmer_dir, "T4SS_t4cp2", [T4CP_HIT])
    write_hmmer_extract(hmmer_dir, "T4SS_virb4", [VIRB4_HIT])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--conjscan-hmmer-dir", str(hmmer_dir), "--min-element-bp", "1000"],
    )
    assert len(rows) == 1
    # All four anchor classes are present here, so without the cap this would be
    # a 'high' call — which is precisely the misleading output being prevented.
    assert rows[0]["n_anchor_classes"] == "4"
    assert rows[0]["evidence_level"] == "profile_hits_only"
    assert rows[0]["confidence"] == "low"


def test_assembled_system_is_not_labelled_as_profile_hits(tmp_path):
    """The normal path is untouched: a real system keeps evidence_level=system."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)
    _return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)
    assert len(rows) == 1
    assert rows[0]["evidence_level"] == "system"
    assert rows[0]["missing_components"] == "none"
    assert rows[0]["confidence"] == "high"


def test_an_icescan_integrase_cannot_upgrade_profile_hits_to_a_system(tmp_path):
    """An attached integrase must not answer "was a conjugation system assembled?".

    The bug this pins. The profile-hit fallback exists for machinery MacSyFinder
    saw but never assembled into a system, and every such row is capped at low
    confidence and labelled evidence_level=profile_hits_only. That judgement used
    to be made by asking whether ANY anchor in the cluster carried a system id —
    and an ICEscan integrase carries one, because ICEscan (unlike CONJScan) does
    have integrase models.

    So the moment ICEscan was switched on, one integrase was enough to make a
    cluster of loose profile hits look system-backed. The row came out
    evidence_level=system at HIGH confidence while its own audit file carried the
    line `no_system_using_profile_hits` saying the machinery had been assembled
    into nothing. The table and the audit contradicted each other.

    The scene below is exactly that: CONJscan found no system, its three
    machinery profiles are recovered from hmmer_results/, and ICEscan supplies a
    Phage_integrase on the same CDS the Bakta product text already calls an
    integrase — which is how the ICEscan system id ends up on the surviving
    anchor (see the corroboration step in main()).
    """
    conjscan = tmp_path / "best_solution.tsv"
    conjscan.write_text("# No System found\n")
    hmmer_dir = tmp_path / "hmmer_results"
    write_hmmer_extract(hmmer_dir, "T4SS_MOBF", [RELAXASE_HIT])
    write_hmmer_extract(hmmer_dir, "T4SS_t4cp2", [T4CP_HIT])
    write_hmmer_extract(hmmer_dir, "T4SS_virb4", [VIRB4_HIT])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id=INTEGRASE_HIT, gene_name="Phage_integrase",
                    model_fqn="ICEscan/Chromosome/IME"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--conjscan-hmmer-dir", str(hmmer_dir),
               "--icescan-tsv", icescan, "--min-element-bp", "1000"],
    )

    assert len(rows) == 1
    row = rows[0]
    # The integrase is real evidence and still counts towards the class and the
    # anchor classes — it is only the "was a SYSTEM assembled?" question it may
    # not answer.
    assert row["has_integrase"] == "TRUE"
    assert row["n_anchor_classes"] == "4"
    # The two things the bug got wrong. Without the fix these read "system" and
    # "high", because nothing else here caps the confidence.
    assert row["evidence_level"] == "profile_hits_only"
    assert row["confidence"] == "low"
    # And the table now agrees with the audit line that was always being written.
    assert "no_system_using_profile_hits" in audit_reasons(audit)


def test_an_icescan_relaxase_system_still_counts_as_an_assembled_system(tmp_path):
    """The fix above narrows the question to machinery — not to CONJscan.

    ICEscan's Gram-positive and IME relaxase families are the reason it is
    unioned in at all, and a relaxase belonging to a real assembled ICEscan
    system IS a system. Fixing the integrase leak must not sweep those up as
    collateral: this row has to keep evidence_level=system.
    """
    conjscan = tmp_path / "best_solution.tsv"
    conjscan.write_text("# No System found\n")
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id=RELAXASE_HIT, gene_name="Relaxase_firmi_MOBL",
                    model_fqn="ICEscan/Chromosome/IME"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--icescan-tsv", icescan, "--min-element-bp", "1000"])

    assert len(rows) == 1
    assert rows[0]["has_relaxase"] == "TRUE"
    assert rows[0]["evidence_level"] == "system"


def test_window_does_not_split_one_conjscan_system(tmp_path):
    """A narrow --window-bp must NOT cut a single CONJscan system into pieces.

    This is the Phase 7 benchmark's clearest finding. Distance clustering used to
    be the only thing deciding what belonged together, so a system whose genes
    were spread slightly wider than the window came out as several "elements".
    On R391 the two halves were 15,010 bp apart against a 15,000 bp window — a
    ten base pair margin — and the reported element lost a third of its length.

    MacSyFinder has already decided these hits form one system, using gene-count
    co-localisation rules rather than base pairs. That decision now wins:
    --window-bp still governs how MACHINERY of DIFFERENT systems is grouped (see
    the next test), but it can no longer split a system apart.
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

    # 3 kb window: VirB4 (starts 63000) is 4299 bp from the coupling protein
    # (ends 58700), so distance clustering alone would break the machinery in two.
    # All three hits carry the same sys_id, so the split is undone.
    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--window-bp", "3000", "--min-element-bp", "1000"]
    )
    assert len(rows) == 1
    row = rows[0]
    # One element with the complete machinery, spanning relaxase to VirB4 — not
    # two fragments, one of which would have been a spurious passive island.
    assert row["mge_class"] == "ice"
    assert row["has_relaxase"] == "TRUE"
    assert row["has_t4cp"] == "TRUE"
    assert row["has_t4ss"] == "TRUE"
    assert row["has_integrase"] == "TRUE"
    # 50000 rather than the relaxase's 55000: the attached integrase at
    # 50000-51200 is part of the element, and the span runs from the first anchor
    # to the last. 65500 is VirB4's end — the half that used to be cut off.
    assert int(row["machinery_start"]) == 50000
    assert int(row["machinery_end"]) == 65500

    # The merge is a filtering decision, so it is audited with its reason.
    merge_reasons = [entry["reason"] for entry in audit]
    assert "clusters_merged_same_conjscan_system" in merge_reasons


def test_system_merge_refused_when_it_would_exceed_max_element(tmp_path):
    """Sharing a sys_id does not license an impossibly large element.

    Caught on the Phase 7 benchmark: on the 10.1 Mb Streptomyces scabiei
    chromosome MacSyFinder assigned hits 2.19 Mb apart to one system. Merging
    them honestly produced a 2,185,966 bp span, which then exceeded
    --max-element-bp and was dropped altogether — so a partial detection became
    a total miss. Refusing the merge keeps the pieces, which is strictly better.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        # Same sys_id (the fixture default), but far away down the contig.
        conjscan_row(hit_id="S1_00090", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    contigs = {"contig_1": 400000}
    cds = list(SCENE_CDS) + [
        gff_cds("contig_1", 350000, 352000, "+", "S1_00090",
                "conjugal transfer protein TraB"),
    ]
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    # A 100 kb ceiling: relaxase at 55000 to the distant VirB4 end at 352000 is
    # ~297 kb, so the merge must be refused.
    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--max-element-bp", "100000", "--min-element-bp", "1000"],
    )
    reasons = [entry["reason"] for entry in audit]
    assert "system_merge_refused_too_long" in reasons
    # The nearby machinery still yields an element rather than nothing at all.
    assert rows
    assert all(int(row["length_bp"]) <= 100000 for row in rows)


def test_window_still_separates_distinct_systems(tmp_path):
    """--window-bp remains meaningful for machinery of DIFFERENT systems.

    The merge above keys on sys_id, so it must not quietly glue together two
    genuinely separate conjugative systems that happen to sit on one contig.
    Here the relaxase belongs to one system and the distant VirB4 to another;
    they stay apart, and only the one with the full complement is called an ICE.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_MOB_1"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_MOB_1"),
        # A second, unrelated system further along the contig.
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeT",
                     sys_id="S1_MOB_2"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _return_code, rows, _audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--window-bp", "3000", "--min-element-bp", "1000"]
    )
    assert len(rows) == 2
    systems = {row["conjscan_systems"] for row in rows}
    assert systems == {"S1_MOB_1", "S1_MOB_2"}


# ── Loners: a shared system id is not always a claim of proximity ────────────

def test_a_loner_hit_does_not_merge_two_distant_clusters(tmp_path):
    """A MacSyFinder LONER must not join two blocks of machinery into one element.

    The case this comes from. CP011419.1 (Streptococcus suis, IME pilot). A MOBT
    relaxase at gene 102 is a loner of system MOB_3, whose only other member is a
    coupling protein 175 genes away at gene 277. The merge rule keyed on the
    system id alone, so those two anchors became one 179,889 bp "IME" — sixteen
    times the curated element — and the genuine 4,959 bp IME inside it was then
    reported a second time.

    Why the system id is not enough here. Merging on a shared system id is
    justified because MacSyFinder has already applied its own co-localisation
    test, counted in genes. A LONER is precisely the gene it exempted from that
    test — the model lets it join from anywhere on the replicon — and MacSyFinder
    says so by writing a NEGATIVE locus_num. So on a loner the system id carries
    no statement about proximity at all.

    The fixture: a relaxase at 55000 and a VirB4 at 350000 in the same system,
    where the relaxase is the loner. They must stay apart.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        # locus_num -1: MacSyFinder admitted this relaxase to the system without
        # requiring it to sit near anything.
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="-1"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="-1"),
        # A real locus member of the same system, 300 kb away.
        conjscan_row(hit_id="S1_00090", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="1"),
    ])
    contigs = {"contig_1": 400000}
    cds = list(SCENE_CDS) + [
        gff_cds("contig_1", 350000, 352000, "+", "S1_00090",
                "conjugal transfer protein TraB"),
    ]
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--min-element-bp", "1000"])

    # No 300 kb blob: nothing reaches from the relaxase to the distant VirB4.
    assert all(int(row["length_bp"]) < 300000 for row in rows)
    assert "system_merge_refused_loner_only_link" in audit_reasons(audit)
    # The refusal names the loner, so a reader can check it in best_solution.tsv.
    refusal = [entry for entry in audit
               if entry["reason"] == "system_merge_refused_loner_only_link"][0]
    assert RELAXASE_HIT in refusal["detail"]
    assert "LONER" in refusal["detail"]


def test_locus_members_of_one_system_are_still_merged(tmp_path):
    """The loner rule must not undo the fix it sits next to.

    Same geometry as the test above — two machinery blocks farther apart than the
    clustering window, one shared system id — but here BOTH sides are genes
    MacSyFinder placed in a real locus (positive locus_num). That is the R391
    case the merge exists for, and it must still produce one element.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="1"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="1"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", locus_num="1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    # 3 kb window: VirB4 (starts 63000) is 4299 bp from the coupling protein
    # (ends 58700), so distance clustering alone would split the machinery.
    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff,
        extra=["--window-bp", "3000", "--min-element-bp", "1000"])

    assert len(rows) == 1
    assert int(rows[0]["machinery_end"]) == 65500          # VirB4's end, not cut off
    assert "clusters_merged_same_conjscan_system" in audit_reasons(audit)
    assert "system_merge_refused_loner_only_link" not in audit_reasons(audit)


def test_an_absent_locus_num_column_still_merges(tmp_path):
    """An older MacSyFinder table without locus_num behaves as it did before.

    The loner rule reads a column we did not use until now. If a future — or
    past — version of the tool does not write it, nothing is assumed about the
    hits and the merge goes ahead, rather than the module silently stopping to
    merge anything.
    """
    columns_without_locus = [c for c in CONJSCAN_COLUMNS if c != "locus_num"]
    rows_out = []
    for hit_id, gene in ((RELAXASE_HIT, "T4SS_MOBF"), (T4CP_HIT, "T4SS_t4cp2"),
                         (VIRB4_HIT, "T4SS_virb4")):
        fields = dict(CONJSCAN_DEFAULTS)
        fields.update({"hit_id": hit_id, "gene_name": gene,
                       "model_fqn": "CONJScan/Chromosome/T4SS_typeF"})
        rows_out.append("\t".join(fields[c] for c in columns_without_locus))
    path = str(tmp_path / "best_solution.tsv")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns_without_locus) + "\n")
        handle.write("\n".join(rows_out) + "\n")

    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)
    _return_code, rows, audit, _ = run_main(
        tmp_path, path, gff,
        extra=["--window-bp", "3000", "--min-element-bp", "1000"])

    assert len(rows) == 1
    assert "clusters_merged_same_conjscan_system" in audit_reasons(audit)


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
    """The stdout line has to distinguish an ICE from a passive island — '2
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


# ── Phase 3: att-site boundaries, and the biology that gates them ────────────
#
# The search only runs for an element that has an integrase, because an att pair
# is the scar that integrase leaves — no integrase, no scar to find. And only a
# tRNA-anchored pair is applied to the coordinates: ~16% of arbitrary chromosomal
# spans carry a de novo repeat by chance, so widening an element onto one would
# turn unrelated chromosomal genes into cargo of something reported as predicted
# self-transmissible.


def test_att_search_is_skipped_when_there_is_no_integrase(tmp_path):
    """An att site is the scar of integrase-mediated recombination, so an element
    with no integrase cannot have one.

    This is the real failure seen on the K. pneumoniae positive control: without
    this gate the search "resolved"
    boundaries for two conjugative_region calls that had no integrase at all —
    one of them on a plasmid, which does not integrate — while the one genuinely
    integrative element got nothing. Any repeat found in that situation is
    something else (an IS end, a duplication, noise), and widening the element to
    it would manufacture a boundary that does not exist.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF", sys_id="S1_T4SS_1"),
    ])
    # SCENE_CDS[0] is the integrase; integrase anchors are found by product regex
    # from the GFF3, not from CONJscan, so it has to be left OUT of the
    # annotation entirely for this element to have none.
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS, SCENE_CDS[1:])
    # A genome carrying a perfectly good direct repeat bracketing the machinery —
    # which must still NOT be reported, because there is no integrase.
    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(99)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = sequence[:39_999] + motif + sequence[40_000 + len(motif) - 1:]
    sequence = sequence[:74_999] + motif + sequence[75_000 + len(motif) - 1:]
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff,
        extra=["--genome", str(genome)])

    assert len(rows) == 1
    assert rows[0]["has_integrase"] == "FALSE"
    assert rows[0]["boundary_method"] == "none"
    assert rows[0]["start"] == rows[0]["machinery_start"]
    assert "att_search_skipped_no_integrase" in audit_reasons(audit)


def test_att_search_widens_an_integrative_element_to_its_real_ends(tmp_path):
    """With an integrase present and a tRNA-anchored att pair bracketing the
    machinery, the element is widened to it — because the cargo between attL and
    attR is what actually travels when the element moves."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    # The tRNA ends at 40024, so its last 25 bp are 40000..40024 — which is where
    # the left att copy is planted. Integration reconstitutes the host tRNA at one
    # end, so exactly one copy lies inside a tRNA and the other does not.
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS,
                    SCENE_CDS + [gff_trna("contig_1", 39_952, 40_024)])

    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(98)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    # Bracket the machinery span (integrase 50000 .. relaxase 56600).
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff,
        extra=["--genome", str(genome)])

    assert len(rows) == 1
    element = rows[0]
    assert element["has_integrase"] == "TRUE"
    assert element["boundary_method"] == "tRNA"
    assert element["attL"] == "40000..40024"
    assert element["attR"] == "70000..70024"
    # start/end are now the ELEMENT, not the machinery — and the machinery span
    # is preserved rather than overwritten.
    assert element["start"] == "40000"
    assert element["end"] == "70024"
    assert element["machinery_start"] == "50000"
    assert int(element["length_bp"]) > (int(element["machinery_end"])
                                        - int(element["machinery_start"]) + 1)
    assert "att_pair_found_tRNA" in audit_reasons(audit)


def test_a_denovo_repeat_is_reported_but_never_moves_the_element(tmp_path):
    """The restraint that keeps fabricated boundaries out of the AMR table.

    Same fixture as above but with NO tRNA, so the bracketing repeat can only be
    found de novo. On a real chromosome ~16% of arbitrary spans yield such a
    repeat by chance, so it is reported for a human to follow up and the element
    interval stays the machinery span — otherwise every gene in the invented
    interval would be called cargo of a self-transmissible element.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS, SCENE_CDS)

    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(98)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff, extra=["--genome", str(genome)])

    element = rows[0]
    # The repeat IS found and IS reported...
    assert element["boundary_method"] == "denovo"
    assert element["attL"] == "40000..40024"
    assert element["attR"] == "70000..70024"
    # ...but the interval is still the machinery span, not the invented element.
    assert element["start"] == element["machinery_start"]
    assert element["end"] == element["machinery_end"]
    assert "denovo_att_reported_not_applied" in audit_reasons(audit)


def test_the_mge_id_is_not_repointed_when_boundaries_move(tmp_path):
    """mge_id is a join key other tables reference, so widening the element must
    not silently change what it points at."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS,
                    SCENE_CDS + [gff_trna("contig_1", 39_952, 40_024)])
    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(97)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, _audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff, extra=["--genome", str(genome)])

    # The id still names the machinery span it was minted from.
    assert rows[0]["mge_id"] == "contig_1|ime-50000:58700"
    assert rows[0]["start"] == "40000"


def test_widening_recomputes_the_contig_distance_flags_and_confidence(tmp_path):
    """Everything derived from the interval must follow it when it moves.

    build_candidates works out dist_to_contig_*, at_contig_boundary and the
    confidence from the MACHINERY span. Phase 3 then widens the element, which
    can carry it to a contig end the machinery never came near. Leaving those
    fields stale is quietly dangerous rather than untidy: at_contig_boundary
    feeds assess_confidence, so a widened element running off the end of its
    contig would keep a high confidence it no longer deserves.

    Here the element widens to 100 bp from the contig end, well inside the
    1000 bp default boundary window, so the flag must flip and the confidence
    must drop.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS,
                    SCENE_CDS + [gff_trna("contig_1", 39_952, 40_024)])

    # A 200 kb contig with the att pair placed so the element ends 100 bp from
    # the far end of the contig.
    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(96)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = plant_att_pair(sequence, motif, 40_000, 199_876)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff,
        extra=["--genome", str(genome), "--att-flank-window-bp", "170000"])

    element = rows[0]
    assert element["boundary_method"] == "tRNA"
    assert element["end"] == "199900"
    # The flags now describe the ELEMENT, not the machinery it grew from.
    assert element["dist_to_contig_end"] == "100"
    assert element["at_contig_boundary"] == "TRUE"
    assert element["confidence"] != "high"
    assert "confidence_settled_after_boundary_search" in audit_reasons(audit)


def test_widening_that_stays_clear_of_the_contig_ends_keeps_its_confidence(tmp_path):
    """The counterpart: rescoring must not fire when nothing actually changed."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS,
                    SCENE_CDS + [gff_trna("contig_1", 39_952, 40_024)])
    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(95)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff, extra=["--genome", str(genome)])

    element = rows[0]
    assert element["boundary_method"] == "tRNA"
    assert element["at_contig_boundary"] == "FALSE"
    assert element["dist_to_contig_start"] == "39999"
    assert "confidence_settled_after_boundary_search" not in audit_reasons(audit)


# ── The strict spec §8 Phase 6 rule, as an opt-in ────────────────────────────

def test_strict_mode_requires_a_trna_boundary_for_high(tmp_path):
    """--require-trna-boundary-for-high applies the spec's literal Phase 6 rule.

    Same evidence as test_ice_needs_all_three_anchor_classes — all four anchor
    classes, intact machinery, one contig — but no genome, so the element's ends
    were never resolved. By default that is reported in boundary_method and the
    call stays high; in strict mode it caps the call at medium.

    The switch exists because the two facts answer different questions ("is this
    an ICE?" vs "where does it stop?"), and on a fragmented short-read assembly
    the second usually cannot be answered at all.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _rc, lenient_rows, _audit, _path = run_main(tmp_path, conjscan, gff)
    assert lenient_rows[0]["confidence"] == "high"
    assert lenient_rows[0]["boundary_method"] == "none"

    _rc, strict_rows, strict_audit, _path = run_main(
        tmp_path, conjscan, gff, extra=["--require-trna-boundary-for-high"])
    assert strict_rows[0]["confidence"] == "medium"
    # The downgrade is explained in the audit rather than left bare.
    settled = [row for row in strict_audit
               if row["reason"] == "confidence_settled_after_boundary_search"]
    assert settled, "the strict downgrade must be audited"
    assert "no att pair was found" in settled[0]["detail"]
    # The element itself is identical either way — only the label changed.
    assert strict_rows[0]["start"] == lenient_rows[0]["start"]
    assert strict_rows[0]["end"] == lenient_rows[0]["end"]


def test_strict_mode_still_allows_high_when_a_trna_boundary_was_found(tmp_path):
    """Strict mode is a requirement, not a blanket downgrade: an element whose
    ends really were fixed by a tRNA-anchored att pair keeps its high call."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS,
                    SCENE_CDS + [gff_trna("contig_1", 39_952, 40_024)])
    genome = tmp_path / "genome.fna"
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    import random as _random
    generator = _random.Random(94)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome.write_text(">contig_1\n" + sequence + "\n")

    _rc, rows, _audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff,
        extra=["--genome", str(genome), "--require-trna-boundary-for-high"])

    assert rows[0]["boundary_method"] == "tRNA"
    assert rows[0]["confidence"] == "high"


def test_an_unknown_contig_length_is_not_read_as_not_at_the_boundary(tmp_path):
    """NA must stay NA through the final confidence pass.

    finalise_confidence re-reads at_contig_boundary from the row it wrote
    earlier. Reading that cell with a plain TRUE/not-TRUE test collapses NA
    ("no contig length was known, so we could not check") into FALSE ("we
    checked and the element is clear of the ends") — an unknown quietly becoming
    a positive claim, which is the one direction this module must not drift in.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    # No sequence-region lines → no contig length is known.
    gff = write_gff(tmp_path / "S1.gff3", {}, SCENE_CDS)

    _rc, rows, _audit, _path = run_main(tmp_path, conjscan, gff)

    assert rows[0]["at_contig_boundary"] == "NA"
    assert rows[0]["confidence"] == "medium"


# ── The ICEscan union ────────────────────────────────────────────────────────
#
# ICEscan is a FORK of CONJScan 2.0.1 by the same Pasteur authors, one minor
# version behind the CONJScan 2.1.0 we run. It ADDS an IME model, an AICE model
# and 21 profiles; it REMOVES MOB.xml, the decayed dCONJ models and the whole
# Plasmids set, and its T4SS quorum is stricter — so swapping to it LOSES
# elements. We union the two instead, and take from ICEscan only its integrase
# anchors, its Gram-positive/IME relaxase families and its IME/AICE model
# classes. These tests pin the parts of that bargain that are easy to break.


def icescan_row(**overrides):
    """One ICEscan best_solution.tsv line.

    Identical column layout to CONJscan's — both are MacSyFinder — so the same
    row builder is reused and only the model namespace differs.
    """
    fields = {"model_fqn": "ICEscan/Chromosome/IME"}
    fields.update(overrides)
    return conjscan_row(**fields)


def test_the_four_trusted_icescan_integrase_profiles_are_integrases():
    """These four are genuine element integrases and must anchor an element."""
    for profile in ("Phage_integrase", "Recombinase", "UPF0236", "PB001819"):
        assert ci.anchor_class_for_gene_name(profile) == ci.ANCHOR_INTEGRASE


def test_the_four_untrusted_integrase_profiles_anchor_nothing():
    """FIX 1-3, and a deliberate divergence from ICEscan's own IME.xml.

    Upstream lists all four as exchangeables of Phage_integrase. None of them is
    an element integrase:

      TIGR02249  IntI1, the class-1 INTEGRON integrase. Measured on CP042858.1:
                 an att-bounded 32,103 bp tier-6 ICE became a 103,303 bp
                 UNBOUNDED one at unchanged 'high' confidence, anchored on IntI1.
      TIGR02224  XerC   } the chromosomal dif-site recombinases every bacterium
      TIGR02225  XerD   } carries. A XerC attached 43,639 bp from a 945 bp
                 relaxase cluster produced an element six times its true size.
      rve        the DDE catalytic domain shared by IS transposases — 66 of its
                 71 hits on the benchmark are Bakta-annotated transposases.

    They must be neither an integrase NOR — via the fall-through default —
    a mating-pair component, which is why None is the required answer.
    """
    for profile in ("TIGR02249", "TIGR02224", "TIGR02225", "rve"):
        assert ci.anchor_class_for_gene_name(profile) is None


def test_icescan_relaxase_families_are_recognised_as_relaxases():
    """The Gram-positive and IME relaxase families are the reason ICEscan sees
    IMEs that CONJScan 2.1.0 cannot. Missing one would silently drop the anchor
    into the mating-pair bucket and could promote an IME to an ICE."""
    for profile in ("Relaxase_firmi_MOBL", "Relaxase_firmi_Rep_2",
                    "Relaxase_firmi_Viral_Rep_A", "Relaxase_firmi_Viral_Rep_B1",
                    "Relaxase_firmi_Viral_Rep_B2", "Relaxase_PHA_IME_A1",
                    "Relaxase_PHA_IME_B", "Relaxase_profile_MOBT", "T4SS_MOBL"):
        assert ci.anchor_class_for_gene_name(profile) == ci.ANCHOR_RELAXASE


def test_aice_machinery_is_never_conjugation_machinery():
    """An AICE translocates double-stranded DNA through a septal pore; it has no
    relaxase and no mating bridge. Classing any of these as T4SS — which the
    fall-through default would have done — would manufacture a mating-pair
    apparatus and promote elements to 'predicted self-transmissible'."""
    for profile in ("FtsK_SpoIIIE", "Prim-Pol", "RepSAv2", "DUF3631"):
        assert ci.anchor_class_for_gene_name(profile) == ci.ANCHOR_AICE


def test_conjscan_profile_names_are_unchanged_by_the_icescan_vocabulary():
    """The new name rules must not have moved any CONJscan profile between
    classes — that would change every existing call."""
    assert ci.anchor_class_for_gene_name("T4SS_MOBF") == ci.ANCHOR_RELAXASE
    assert ci.anchor_class_for_gene_name("T4SS_t4cp2") == ci.ANCHOR_T4CP
    assert ci.anchor_class_for_gene_name("T4SS_tcpA") == ci.ANCHOR_T4CP
    assert ci.anchor_class_for_gene_name("T4SS_virb4") == ci.ANCHOR_T4SS
    assert ci.anchor_class_for_gene_name("T4SS_F_traU") == ci.ANCHOR_T4SS


def test_without_icescan_the_result_is_byte_for_byte_what_it_always_was(tmp_path):
    """THE CONTROL. A user who has not downloaded the ICEscan models — the
    default — must get exactly today's answer.

    The models are CC BY-NC-SA and fetched at runtime, so most runs will not have
    them. Running the same input with and without --icescan-tsv pointed at
    nothing must produce identical tables."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _rc, without_flag, _audit, path_a = run_main(tmp_path / "a", conjscan, gff)
    # The flag present but pointing at a file that does not exist: the module
    # must degrade rather than fail, and must not change its answer.
    _rc, with_absent, _audit_b, path_b = run_main(
        tmp_path / "b", conjscan, gff,
        extra=["--icescan-tsv", str(tmp_path / "not_there.tsv")])

    assert without_flag == with_absent
    with open(path_a) as a, open(path_b) as b:
        assert a.read() == b.read()


def test_an_icescan_integron_integrase_cannot_anchor_an_element(tmp_path):
    """FIX 1 end to end, as the CP042858.1 regression would have arrived.

    A TIGR02249 hit sits 45 kb from the machinery, well inside the 50 kb
    integrase window. Trusting it would attach it, stretch the element to cover
    it, and report an ICE. It must anchor nothing, leaving the cluster with no
    integrase — so the honest "machinery, but no element boundaries" call stands.

    NOTE the CDS is deliberately annotated "hypothetical protein" so that the
    ONLY thing that could make it an integrase is the ICEscan profile, which is
    what this test is about. Bakta usually annotates IntI1 as "class 1 integron
    integrase IntI1", and INTEGRASE_PRODUCT_PATTERN matches that on purpose (see
    its comment) — so the product-text path has its own, separate exposure to
    integron integrases, which the FIX-4 tie-break rather than this rule is what
    keeps in check.
    """
    cds = list(SCENE_CDS)
    # Replace the real integrase with a CDS only ICEscan could call an integrase.
    cds = [line for line in cds if "S1_00010" not in line]
    cds.append(gff_cds("contig_1", 95000, 96000, "+", "S1_00090",
                       "hypothetical protein"))
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id="S1_00090", gene_name="TIGR02249",
                    hit_gene_ref="Phage_integrase", sys_id="S1_IME_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff,
                                       extra=["--icescan-tsv", icescan])

    assert len(rows) == 1
    assert rows[0]["has_integrase"] == "FALSE"
    assert rows[0]["element_type"] != "ice"
    # The element must not have been stretched to reach the IntI1 at 95–96 kb.
    assert int(rows[0]["end"]) < 90000
    assert "untrusted_integrase_profile" in audit_reasons(audit)


def test_an_icescan_rve_hit_cannot_anchor_an_element(tmp_path):
    """FIX 3. `rve` is the transposase catalytic domain; find_integrase_anchors
    already refuses these on the product-text side, and admitting them as an HMM
    hit would let an insertion sequence next to a relaxase become an ICE."""
    # Annotated as a transposase, which is what these really are: 66 of rve's 71
    # hits on the benchmark are Bakta-annotated transposases. The product-text
    # path already refuses it (TRANSPOSASE_PRODUCT_PATTERN); this pins that the
    # HMM path refuses it too, so the two cannot let it in by different doors.
    cds = [line for line in SCENE_CDS if "S1_00010" not in line]
    cds.append(gff_cds("contig_1", 52000, 53000, "+", "S1_00012",
                       "IS3 family transposase"))
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id="S1_00012", gene_name="rve",
                    hit_gene_ref="Phage_integrase", sys_id="S1_IME_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, _audit, _path = run_main(tmp_path, conjscan, gff,
                                        extra=["--icescan-tsv", icescan])

    assert rows[0]["has_integrase"] == "FALSE"
    assert rows[0]["element_type"] != "ice"


def test_an_icescan_relaxase_can_make_an_ime_conjscan_would_have_missed(tmp_path):
    """The measured gain. A Gram-positive relaxase family CONJScan 2.1.0 does
    not model, plus a product-text integrase, is exactly the IME architecture —
    and its machinery is two genes, so it needs the lower IME size floor."""
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "tyrosine recombinase XerC"),
        gff_cds("contig_1", 52000, 53500, "+", "S1_00011",
                "MobV family relaxase"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id="S1_00011", gene_name="Relaxase_firmi_MOBL",
                    sys_id="S1_IME_1", hit_gene_ref="T4SS_MOBV"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, _audit, _path = run_main(tmp_path, conjscan, gff,
                                        extra=["--icescan-tsv", icescan])

    assert len(rows) == 1
    assert rows[0]["mge_class"] == "ime"
    assert rows[0]["has_relaxase"] == "TRUE"
    assert rows[0]["has_integrase"] == "TRUE"
    assert "icescan" in rows[0]["evidence_sources"]


def test_an_icescan_ime_model_cannot_promote_an_island(tmp_path):
    """ICEscan's IME quorum is two genes, both declared loners, so the model can
    fire from hits anywhere on the replicon. It must never turn a cluster with no
    relaxase into an IME — that would be a tier-5 mobility claim on no evidence.
    """
    # The cluster has a mating-pair gene and an integrase but NO relaxase, so our
    # own rules call it a passive island. ICEscan's IME model fires over the same
    # integrase — and must not be allowed to change the answer.
    cds = list(SCENE_CDS)
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_T4SS_typeF_1"),
    ])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id=INTEGRASE_HIT, gene_name="Phage_integrase",
                    sys_id="S1_IME_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff,
                                       extra=["--icescan-tsv", icescan])

    assert len(rows) == 1
    assert rows[0]["mge_class"] != "ime"
    assert rows[0]["has_relaxase"] == "FALSE"
    assert "icescan_ime_model_not_followed" in audit_reasons(audit)


def test_icescan_system_ids_cannot_collide_with_conjscan_ones():
    """Both tools define a model called T4SS_typeF, and MacSyFinder builds a
    system id as {replicon}_{model}_{n}. Run over one genome they therefore emit
    DIFFERENT systems under IDENTICAL ids. Merging them would pool their contigs
    into a false spans_contigs flag and let two clusters be joined that neither
    tool ever said belonged together."""
    assert ci.namespaced_sys_id("X_T4SS_typeF_3", ci.SOURCE_CONJSCAN) == "X_T4SS_typeF_3"
    assert ci.namespaced_sys_id("X_T4SS_typeF_3", ci.SOURCE_ICESCAN) != "X_T4SS_typeF_3"


# ── AICE: a third class, and it is not on the conjugation ladder ─────────────

def aice_scene(tmp_path):
    """A minimal AICE: integrase, FtsK/SpoIIIE translocase and a Rep protein,
    called as ICEscan's AICE model. No relaxase and no mating-pair gene — which
    is what an AICE is, not what is wrong with it."""
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010",
                "site-specific recombinase"),
        gff_cds("contig_1", 52000, 54400, "+", "S1_00011",
                "FtsK/SpoIIIE family DNA translocase"),
        gff_cds("contig_1", 54600, 55900, "+", "S1_00012",
                "replication initiator protein"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id="S1_00011", gene_name="FtsK_SpoIIIE",
                    model_fqn="ICEscan/Chromosome/AICE", sys_id="S1_AICE_1"),
        icescan_row(hit_id="S1_00012", gene_name="RepSAv2",
                    model_fqn="ICEscan/Chromosome/AICE", sys_id="S1_AICE_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)
    return conjscan, icescan, gff


def test_an_aice_is_reported_as_its_own_class(tmp_path):
    """It must not be folded into `ice`: an AICE has no relaxase and no mating
    bridge, so calling it an ICE would assert self-transmissibility it cannot
    have."""
    conjscan, icescan, gff = aice_scene(tmp_path)
    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff,
                                       extra=["--icescan-tsv", icescan])

    assert len(rows) == 1
    assert rows[0]["mge_class"] == "aice"
    assert rows[0]["element_type"] == "aice"
    assert rows[0]["icescan_model_class"] == "AICE"
    assert "class_taken_from_icescan_aice_model" in audit_reasons(audit)


def test_an_aice_never_claims_a_conjugation_tier(tmp_path):
    """THE HEADLINE RULE. Tiers 5 and 6 are both conjugation — "mobilisable by a
    helper" and "self-transmissible" — and an AICE does neither: it moves as
    double-stranded DNA between hyphal compartments by FtsK/SpoIIIE
    translocation. Either tier would be a false claim, so it gets none, and the
    reason is spelled out rather than left blank."""
    conjscan, icescan, gff = aice_scene(tmp_path)
    _rc, rows, _audit, _path = run_main(tmp_path, conjscan, gff,
                                        extra=["--icescan-tsv", icescan])
    row = rows[0]

    assert row["mobility_tier"] == "NA"
    assert row["mobility_tier_reason"].strip()          # never blank
    assert "conjugation" in row["mobility_tier_reason"]
    # The mobility sentence names the real mechanism instead of hedging.
    assert "FtsK/SpoIIIE" in row["mobility"]
    assert "not on the conjugation mobility ladder" in row["mobility"]
    assert "self-transmissible" not in row["mobility"]
    assert "mobilisable" not in row["mobility"]
    # colocalise.py is the enforcement point: recognised, but context only, so
    # it can raise no gene's tier.
    assert co.ELEMENT_TYPE_SYNONYMS["aice"] == "aice"
    assert "aice" in co.CONTEXT_ONLY_ELEMENT_TYPES
    assert co.ELEMENT_TYPE_SYNONYMS["aice"] not in {"ice", "ime"}


def test_an_aice_does_not_read_as_a_degraded_ice(tmp_path):
    """`missing_components` must not list the relaxase and mating bridge as
    missing: they are absent by definition, not broken."""
    conjscan, icescan, gff = aice_scene(tmp_path)
    _rc, rows, _audit, _path = run_main(tmp_path, conjscan, gff,
                                        extra=["--icescan-tsv", icescan])

    assert rows[0]["missing_components"] == ci.AICE_MISSING_COMPONENTS
    assert "none expected" in rows[0]["missing_components"]
    assert rows[0]["machinery_intact"] == "TRUE"
    assert ci.DEGRADED_MOBILITY_SUFFIX not in rows[0]["mobility"]


def test_loose_ftsk_hits_do_not_manufacture_an_aice(tmp_path):
    """FtsK/SpoIIIE is a core chromosome-partitioning ATPase present in
    essentially every bacterium. Without ICEscan's assembled AICE model behind
    it, it must seed nothing at all — otherwise every genome we ever run grows an
    'AICE'."""
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "integrase"),
        gff_cds("contig_1", 52000, 54400, "+", "S1_00011",
                "DNA translocase FtsK"),
    ]
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [])
    # A profile hit with no system behind it: sys_id empty, no AICE model.
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id="S1_00011", gene_name="FtsK_SpoIIIE",
                    model_fqn="", sys_id=""),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff,
                                       extra=["--icescan-tsv", icescan])

    assert rows == []
    assert "no_conjugation_anchor" in audit_reasons(audit)


# ── FIX 4: which integrase belongs to the element ────────────────────────────

def test_the_tie_break_prefers_the_integrase_that_yields_an_att_pair(tmp_path):
    """FIX 4, and the reason the CP042858.1 regression was possible.

    Two integrases are in the running for one cluster. The old rule was
    "closest, first-seen wins a tie", so when both sat inside the machinery span
    — gap 0 for each, which is the normal case once ICEscan's hits join the pool
    — the winner was decided by list order, i.e. by SOURCE rather than by any
    evidence.

    The geometry here makes distance and att support disagree on purpose:

      40000..40024   attL — the last 25 bp of the tRNA ending at 40024
      41000..42200   the REAL integrase, inside the att-bounded interval,
                     23,800 bp from the machinery
      66000..68500   the conjugation machinery
      70000..70024   attR — the second copy of that same 25-mer
      70500..71500   the DECOY, only 2,000 bp from the machinery, but OUTSIDE
                     the att interval

    Attaching the decoy pushes the element's right edge past attR, so attR is
    swallowed by the span instead of flanking it and no pair can be found.
    Attaching the real integrase leaves both copies in the flanks, where the
    scar of integration actually lies. The decoy is more than ten times closer,
    so if distance still ranked first it would win — and the element would come
    out unbounded, which is precisely the CP042858.1 failure.

    Both candidates come from the Bakta product text, so this isolates key 1
    (att support) from key 3 (source) and tests it against key 2 (distance)
    alone.
    """
    import random as _random
    generator = _random.Random(1729)
    sequence = "".join(generator.choice("ACGT") for _ in range(200_000))
    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    sequence = plant_att_pair(sequence, motif, 40_000, 70_000)
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig_1\n" + sequence + "\n")

    cds = [
        # The REAL integrase: far from the machinery, but inside the element.
        gff_cds("contig_1", 41_000, 42_200, "+", "S1_00005",
                "phage integrase family protein"),
        # The machinery, sitting near the element's right-hand edge.
        gff_cds("contig_1", 66_000, 67_600, "+", "S1_00015",
                "TrwC relaxase domain-containing protein"),
        gff_cds("contig_1", 67_800, 68_200, "-", "S1_00016",
                "Type IV secretory pathway, VirD4 component"),
        gff_cds("contig_1", 68_300, 68_500, "+", "S1_00017",
                "conjugal transfer protein TraB"),
        # The DECOY: much closer, but outside the att-bounded interval.
        gff_cds("contig_1", 70_500, 71_500, "+", "S1_00020",
                "tyrosine recombinase XerD"),
    ]
    trna = gff_trna("contig_1", 39_952, 40_024)
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_T4SS_typeF_1"),
        conjscan_row(hit_id="S1_00016", gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_T4SS_typeF_1"),
        conjscan_row(hit_id="S1_00017", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF",
                     sys_id="S1_T4SS_typeF_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds + [trna],
                    with_fasta=False)

    _rc, rows, _audit, _path = run_main(
        tmp_path, conjscan, gff, extra=["--genome", str(genome)])

    assert len(rows) == 1
    element = rows[0]
    # The att-supporting integrase was chosen over the ten-times-closer decoy.
    assert "phage integrase" in element["integrase_products"].lower()
    assert "XerD" not in element["integrase_products"]
    # ...and because it was, the element has real boundaries instead of none.
    assert element["boundary_method"] == "tRNA"
    assert element["attL"] == "40000..40024"
    assert element["attR"] == "70000..70024"


def test_at_equal_distance_an_hmm_hit_does_not_displace_the_annotation(tmp_path):
    """The last clause of FIX 4. When two integrase candidates are the same
    distance away and neither yields an att pair, the Bakta product text wins:
    it names the protein in words a reader can check, and it is the source that
    found the integrase for 12 of the 30 curated pilot elements."""
    assert (ci.INTEGRASE_SOURCE_RANK[ci.SOURCE_BAKTA_PRODUCT]
            < ci.INTEGRASE_SOURCE_RANK[ci.SOURCE_ICESCAN])

    # The same CDS found by both sources is counted once, keeping the readable
    # product text as its label.
    cds = list(SCENE_CDS)
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    icescan = write_conjscan(tmp_path / "icescan.tsv", [
        icescan_row(hit_id=INTEGRASE_HIT, gene_name="Phage_integrase",
                    sys_id="S1_IME_1"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, cds)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff,
                                       extra=["--icescan-tsv", icescan])

    assert len(rows) == 1
    # Counted once, not twice.
    assert rows[0]["anchor_ids"].count(INTEGRASE_HIT) == 1
    # Labelled with the Bakta product text, not the profile name.
    assert "Phage integrase family protein" in rows[0]["integrase_products"]
    assert "integrase_corroborated_by_icescan" in audit_reasons(audit)


# ── One locus, one row: the nested double-report ─────────────────────────────
#
# Two calls of the same class where one sits inside the other describe the same
# neighbourhood twice, and a reader has no way to tell which line to believe.
# resolve_nested_calls keeps one of them — on evidence, never simply the smaller
# one — and writes the other to the audit. The end-to-end test below reproduces
# the CP011419.1 shape that prompted this; the pure-function tests after it pin
# each rung of the decision ladder separately, including the cases where the
# LARGER call is the one that must survive.

def nesting_row(mge_id, contig, start, end, mge_class="ime",
                boundary_method="none", machinery_gap_bp=0, anchor_ids="x(relaxase)"):
    """A minimal element row, just the columns resolve_nested_calls reads.

    Building these by hand rather than through main() keeps each rule of the
    ladder testable on its own; the end-to-end test covers the wiring.
    """
    return {
        "mge_id": mge_id,
        "contig": contig,
        "start": str(start),
        "end": str(end),
        "length_bp": str(end - start + 1),
        "mge_class": mge_class,
        "boundary_method": boundary_method,
        "machinery_gap_bp": str(machinery_gap_bp),
        "n_anchors": "3",
        "anchor_ids": anchor_ids,
    }


def test_one_locus_is_not_reported_as_two_nested_imes(tmp_path):
    """The CP011419.1 defect, end to end: a blob and the honest call inside it.

    The measured case. On CP011419.1 the caller emitted a 179,889 bp "IME"
    spanning sixteen times the curated element, and — inside it — the honest
    4,959 bp IME that sits 63 bp from the curated start. Both were reported, so
    one locus appeared twice and the benchmark still scored the element as a
    16.19x swallow.

    THE CAUSE, established from MacSyFinder's own output rather than guessed:
    the blob was built by merging through a LONER. CP011419_1_00138 (T4SS_MOBT)
    is listed under two systems with locus_num = -1, and it is the sole entry of
    best_solution_loners.tsv. A loner is precisely the gene a model admits from
    anywhere on the replicon WITHOUT the co-localisation test, and MacSyFinder
    signals that with a negative locus_num. So merge_clusters_sharing_a_system's
    justification — "MacSyFinder already decided these genes form one system" —
    is false for a loner, and merging on it invented a 179,889 bp interval.

    The fixture reproduces that shape: a compact, genuine element, and a distant
    loner of the same system 240 kb away. The loner must not drag the two
    together, so only one call comes out and no nesting resolution is needed.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        # The LONER, 240 kb from everything else. locus_num = -1 is
        # MacSyFinder telling us it was admitted without the co-localisation
        # test, so it must not link two clusters into one element.
        conjscan_row(hit_id="S1_00100", gene_name="T4SS_MOBF",
                     sys_id="S1_MOB_2", locus_num="-1"),
        # The compact element in between: its own system, three genes together.
        conjscan_row(hit_id="S1_00200", gene_name="T4SS_MOBP1",
                     sys_id="S1_MOB_2", locus_num="1"),
        conjscan_row(hit_id="S1_00201", gene_name="T4SS_t4cp2",
                     sys_id="S1_MOB_2", locus_num="1"),
    ])
    contigs = {"contig_1": 400000}
    cds = [
        gff_cds("contig_1", 60000, 61600, "+", "S1_00100",
                "TrwC relaxase domain-containing protein"),
        gff_cds("contig_1", 199000, 200200, "+", "S1_00199",
                "Phage integrase family protein"),
        gff_cds("contig_1", 200400, 202000, "+", "S1_00200",
                "MobA/MobL family protein"),
        gff_cds("contig_1", 202400, 204100, "-", "S1_00201",
                "Type IV secretory pathway, VirD4 component, TraG/TraD family ATPase"),
        gff_cds("contig_1", 300000, 301700, "-", "S1_00400",
                "Type IV secretory pathway, VirD4 component, TraG/TraD family ATPase"),
    ]
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--min-element-bp", "1000"])

    # The compact element is reported at its own size, NOT stretched to 240 kb.
    ime_rows = [row for row in rows if row["mge_class"] == "ime"]
    assert len(ime_rows) == 1
    kept = ime_rows[0]
    assert int(kept["start"]) == 199000 and int(kept["end"]) == 204100

    # The orphaned loner is still REPORTED — separately, and as the weakest class
    # its evidence supports. A relaxase with no integrase beside it is a
    # conjugative region, not an element with boundaries. This is the same shape
    # seen on the real CP011419.1, where the loner's own compact locus came out
    # as its own small call rather than being folded into the element 150 kb away.
    assert all(int(row["length_bp"]) < 100000 for row in rows)
    # The blob is never BUILT, so this is a refusal to merge rather than a
    # suppression after the fact — and the audit says which gene caused it.
    assert "system_merge_refused_loner_only_link" in audit_reasons(audit)


def test_two_separate_elements_on_one_contig_are_both_reported(tmp_path):
    """The guard against over-suppression: nesting is not the same as neighbouring.

    Two genuinely separate IMEs on one contig, neither inside the other, must
    both survive. Without this the fix above would quietly become "report at most
    one element per contig", which is the opposite failure.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00100", gene_name="T4SS_MOBF",
                     sys_id="S1_MOB_1", locus_num="1"),
        conjscan_row(hit_id="S1_00101", gene_name="T4SS_t4cp2",
                     sys_id="S1_MOB_1", locus_num="1"),
        conjscan_row(hit_id="S1_00300", gene_name="T4SS_MOBP1",
                     sys_id="S1_MOB_2", locus_num="1"),
        conjscan_row(hit_id="S1_00301", gene_name="T4SS_t4cp2",
                     sys_id="S1_MOB_2", locus_num="1"),
    ])
    contigs = {"contig_1": 400000}
    cds = [
        gff_cds("contig_1", 60000, 61200, "+", "S1_00099",
                "Phage integrase family protein"),
        gff_cds("contig_1", 61500, 63100, "+", "S1_00100",
                "TrwC relaxase domain-containing protein"),
        gff_cds("contig_1", 63400, 65100, "-", "S1_00101",
                "Type IV secretory pathway, VirD4 component, TraG/TraD family ATPase"),
        gff_cds("contig_1", 300000, 301200, "+", "S1_00299",
                "tyrosine-type recombinase/integrase"),
        gff_cds("contig_1", 301500, 303100, "+", "S1_00300",
                "MobA/MobL family protein"),
        gff_cds("contig_1", 303400, 305100, "-", "S1_00301",
                "Type IV secretory pathway, VirD4 component, TraG/TraD family ATPase"),
    ]
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _return_code, rows, audit, _ = run_main(
        tmp_path, conjscan, gff, extra=["--min-element-bp", "1000"])

    assert len(rows) == 2
    assert {row["mge_class"] for row in rows} == {"ime"}
    assert [int(row["start"]) for row in rows] == [60000, 300000]
    assert "nested_call_of_same_class_suppressed" not in audit_reasons(audit)


def test_a_nested_call_with_an_att_boundary_beats_the_larger_one():
    """Key 1: a tRNA-anchored att pair decides it, whichever call is bigger.

    An att site is the scar left where the element recombined into the
    chromosome, so a call whose ends came from one knows where the element
    starts and stops; a call whose ends are only the span of its machinery does
    not. Here the SMALLER call carries the att pair and wins.
    """
    outer = nesting_row("c1|ime-1000:90000", "c1", 1000, 90000, machinery_gap_bp=0)
    inner = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000,
                        boundary_method="tRNA", machinery_gap_bp=0)
    kept, audit = ci.resolve_nested_calls("S1", [outer, inner])

    assert [row["mge_id"] for row in kept] == ["c1|ime-40000:50000"]
    assert audit[0]["reason"] == "nested_call_of_same_class_suppressed"
    assert "tRNA-anchored att pair" in audit[0]["detail"]


def test_the_larger_call_survives_when_it_is_the_one_with_the_att_boundary():
    """Key 1 again, the other way round — because "keep the smaller one" would be
    wrong. A 100 kb ICE genuinely contains smaller blocks of machinery, and when
    the LARGE call is the one with the att evidence it is the element."""
    outer = nesting_row("c1|ice-1000:90000", "c1", 1000, 90000, mge_class="ice",
                        boundary_method="tRNA", machinery_gap_bp=40000)
    inner = nesting_row("c1|ice-40000:50000", "c1", 40000, 50000, mge_class="ice",
                        machinery_gap_bp=0)
    kept, _audit = ci.resolve_nested_calls("S1", [outer, inner])

    assert [row["mge_id"] for row in kept] == ["c1|ice-1000:90000"]


def test_machinery_coherence_is_reported_but_never_decides():
    """The removed key 2: a wide machinery gap must NOT lose a nesting contest.

    An earlier version preferred the call whose machinery "sits together as one
    operon", on the reasoning that conjugation genes form an operon. Measured on
    the benchmark that premise is false for exactly the elements we care about:
    20 of 37 ice calls (54%) have an anchor-free hole wider than the 15 kb
    window, among them R391 (28,354 bp), SPI-7 (42,039) and Tn4371 (15,120) —
    the spec's own positive controls. Large ICEs carry cargo BETWEEN their
    machinery genes. The rule deleted a 193 kb ICE in favour of a 7 kb element
    inside it, so it was removed; machinery_gap_bp is still reported for a reader
    to judge by eye.
    """
    outer = nesting_row("c1|ime-1000:90000", "c1", 1000, 90000, machinery_gap_bp=70000)
    inner = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000, machinery_gap_bp=300)
    kept, audit = ci.resolve_nested_calls("S1", [outer, inner])

    # Key 3 decides instead: with no att evidence, the wider interval is kept.
    assert [row["mge_id"] for row in kept] == ["c1|ime-1000:90000"]
    assert "no evidence separates them" in audit[0]["detail"]


def test_nested_calls_with_nothing_to_separate_them_keep_the_outer_one():
    """Key 3: with no evidence either way, report the wider interval.

    It already contains every base and every anchor the inner call had, so the
    inner one is cargo of it rather than a second finding. This is the EBI
    Mobilome Annotation Pipeline's convention, adopted here as a design decision
    (their code is CC BY-NC-SA and is never copied — see spec §11).
    """
    outer = nesting_row("c1|ime-1000:90000", "c1", 1000, 90000, machinery_gap_bp=200)
    inner = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000, machinery_gap_bp=100)
    kept, audit = ci.resolve_nested_calls("S1", [outer, inner])

    assert [row["mge_id"] for row in kept] == ["c1|ime-1000:90000"]
    assert "no evidence separates them" in audit[0]["detail"]


def test_an_ime_nested_inside_an_ice_is_still_reported():
    """Two DIFFERENT classes nested are two different elements, and both stand.

    An IME sitting inside an ICE is real cargo — and the more mobile of the two
    findings, since it can be picked up by a helper independently. Suppressing it
    would lose the answer a reader most needs. This pass only removes a duplicate
    description of ONE locus, which is what a same-class nest is.
    """
    ice = nesting_row("c1|ice-1000:90000", "c1", 1000, 90000, mge_class="ice",
                      machinery_gap_bp=40000)
    ime = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000, mge_class="ime",
                      machinery_gap_bp=100)
    kept, audit = ci.resolve_nested_calls("S1", [ice, ime])

    assert len(kept) == 2
    assert audit == []


def test_overlapping_calls_that_do_not_nest_are_both_kept():
    """Partial overlap is not containment. Two calls that merely share some bases
    are two findings with a shared neighbourhood, and both are reported — the
    rule is deliberately narrow."""
    left = nesting_row("c1|ime-1000:50000", "c1", 1000, 50000)
    right = nesting_row("c1|ime-40000:90000", "c1", 40000, 90000)
    kept, audit = ci.resolve_nested_calls("S1", [left, right])

    assert len(kept) == 2
    assert audit == []


def test_nested_calls_on_different_contigs_are_never_compared():
    """Two contigs are two pieces of DNA. Coordinates on one say nothing about
    the other, so a call at 40000-50000 on contig_2 is not 'inside' anything on
    contig_1."""
    outer = nesting_row("c1|ime-1000:90000", "contig_1", 1000, 90000,
                        machinery_gap_bp=70000)
    inner = nesting_row("c2|ime-40000:50000", "contig_2", 40000, 50000,
                        machinery_gap_bp=100)
    kept, audit = ci.resolve_nested_calls("S1", [outer, inner])

    assert len(kept) == 2
    assert audit == []


def test_a_chain_of_three_nested_calls_collapses_to_one():
    """A ⊃ B ⊃ C is still one locus described three times.

    Resolved one pair at a time, so the ladder never has to reason about three
    calls at once. Here only the innermost has an att pair, so it is what stands.
    """
    a = nesting_row("c1|ime-1000:90000", "c1", 1000, 90000, machinery_gap_bp=70000)
    b = nesting_row("c1|ime-30000:60000", "c1", 30000, 60000, machinery_gap_bp=20000)
    c = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000,
                    boundary_method="tRNA", machinery_gap_bp=100)
    kept, audit = ci.resolve_nested_calls("S1", [a, b, c])

    assert [row["mge_id"] for row in kept] == ["c1|ime-40000:50000"]
    assert len(audit) == 2


def test_machinery_gap_bp_measures_the_widest_hole_in_the_machinery(tmp_path):
    """The column the nesting rule reads: the biggest anchor-free stretch inside
    a call. Small for a real operon however long the element is; large only when
    the interval was stitched from distant blocks."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)
    _return_code, rows, _audit, _ = run_main(tmp_path, conjscan, gff)

    assert len(rows) == 1
    # SCENE_CDS anchors: 50000-51200, 55000-56600, 57000-58700, 63000-65500.
    # The widest hole is between the coupling protein and VirB4: 63000-58700-1.
    # Every gap here is small — this is one operon, which is the point.
    assert rows[0]["machinery_gap_bp"] == "4299"


def test_two_calls_on_the_same_interval_collapse_to_the_better_evidenced_one():
    """Identical coordinates are one locus reported twice, whatever the classes.

    This appears when two machinery clusters resolve onto the SAME att pair: both
    are widened to the element the repeats define, and two rows come out with the
    same start and end. Measured on ICEEc2 (GU725392) once the att search was
    fixed — an `ime` row and an `ice` row, both 27-92263. The benchmark's scorer
    (phase7_benchmark/score.py, which lives with the benchmark rather than in
    this repo) broke the tie towards the `ime`: tier 5 for an element correctly
    identified as a tier-6 ICE.

    The same-class nesting rule cannot catch this: the classes differ, and the
    intervals nest in neither direction because they are equal.
    """
    weak = nesting_row("c1|ime-27:92263", "c1", 27, 92263, machinery_gap_bp=100)
    weak["mge_class"] = "ime"
    weak["n_anchor_classes"] = "2"
    strong = nesting_row("c1|ice-27:92263", "c1", 27, 92263, machinery_gap_bp=100)
    strong["mge_class"] = "ice"
    strong["n_anchor_classes"] = "4"

    kept, audit = ci.resolve_nested_calls("S1", [weak, strong])

    assert len(kept) == 1
    assert kept[0]["mge_class"] == "ice"          # the better-evidenced call
    assert "identical_interval_reported_twice" in [row["reason"] for row in audit]


def test_an_ime_genuinely_inside_an_ice_is_not_collapsed():
    """Cross-class nesting stays two findings — only EQUAL intervals collapse.

    An IME sitting inside an ICE is two real elements, and the IME is the more
    mobile finding. The identical-interval rule must not become a back door that
    suppresses it.
    """
    ice = nesting_row("c1|ice-1000:90000", "c1", 1000, 90000, machinery_gap_bp=100)
    ice["mge_class"] = "ice"
    ime = nesting_row("c1|ime-40000:50000", "c1", 40000, 50000, machinery_gap_bp=100)
    ime["mge_class"] = "ime"

    kept, _audit = ci.resolve_nested_calls("S1", [ice, ime])

    assert len(kept) == 2


# ── Fragmented (draft) assemblies ────────────────────────────────────────────
#
# Everything below was written after the module was measured on drafts for the
# first time. Until then every validation had been on CLOSED genomes, where
# spans_contigs was TRUE on 0 of 63 calls — so the guards that exist for
# fragmentation had never been exercised by a test or by a benchmark. The
# fragmented-assembly validation cut 40 benchmark genomes to three contiguities
# (~150 kb, ~50 kb and ~20 kb N50), re-annotated all 120 assemblies and re-ran
# the whole chain; the numbers quoted in these tests come from it.
# Full write-up: docs/mobilome_draft_assemblies.md.


def test_assembly_contiguity_reports_n50_and_contig_count():
    """N50 as a reader expects it: half the assembly is in contigs this long or
    longer. Checked on a hand-computable case so the arithmetic is visible.

    100 + 50 + 30 + 20 = 200 kb total, half is 100 kb; the longest contig alone
    reaches it, so N50 is 100 kb — not the median contig length (40 kb), which
    is the usual way of getting this wrong.
    """
    stats = ci.assembly_contiguity({
        "c1": 100000, "c2": 50000, "c3": 30000, "c4": 20000,
    })

    assert stats["n_contigs"] == 4
    assert stats["total_bp"] == 200000
    assert stats["longest_bp"] == 100000
    assert stats["n50_bp"] == 100000


def test_assembly_contiguity_of_a_closed_genome_is_the_genome():
    """A single-contig genome: N50, longest and total are all the same number.
    This is the input every earlier validation used, and it must not be a
    special case in the code."""
    stats = ci.assembly_contiguity({"chromosome": 5400000})

    assert stats == {
        "n_contigs": 1,
        "total_bp": 5400000,
        "longest_bp": 5400000,
        "n50_bp": 5400000,
    }


def test_assembly_contiguity_is_none_when_there_are_no_contigs():
    """No lengths means no claim — the caller stays silent rather than writing a
    row of zeros that reads like a measurement."""
    assert ci.assembly_contiguity({}) is None


def test_contiguity_verdict_changes_at_the_two_validated_levels():
    """The three sentences map onto the three arms the validation actually ran,
    and each one has to say something a reader can act on."""
    trusted = ci.contiguity_verdict(300000)
    class_only = ci.contiguity_verdict(60000)
    poor = ci.contiguity_verdict(20000)

    assert trusted != class_only != poor
    # The lower divider sits BETWEEN the two arms it separates: detection held
    # at the 50 kb arm and broke at the 20 kb one, so an assembly a little under
    # 50 kb must still get the "class is reliable" message, not the bottom one.
    assert ci.contiguity_verdict(49000) == class_only
    # At good contiguity the message is reassuring but still says lengths are a
    # floor; below 150 kb it must say the extent is not reliable; at the bottom
    # it must say detection itself suffers.
    assert "floor" in trusted
    assert "EXTENT is not" in class_only
    assert "DETECTION" in poor


def test_every_run_writes_one_assembly_contiguity_audit_row(tmp_path):
    """The audit's first job is to say what the calls came off.

    A reader opening the element table cannot tell a closed genome from a
    300-contig draft, and every caveat in this module depends on that
    difference — so the contig count and N50 are recorded once per sample,
    unconditionally, for closed genomes too.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _rc, _rows, audit, _path = run_main(tmp_path, conjscan, gff)

    qc = [row for row in audit if row["action"] == "assembly_qc"]
    assert len(qc) == 1
    assert qc[0]["reason"] == "assembly_contiguity"
    assert "1 contig(s)" in qc[0]["detail"]
    assert "N50 200,000 bp" in qc[0]["detail"]


def test_the_contiguity_line_is_written_even_when_nothing_is_found(tmp_path):
    """'No elements' off a closed genome and 'no elements' off a shattered draft
    are different results. The audit has to distinguish them, so the QC row must
    survive the early exit that a missing CONJscan file takes."""
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan=None, gff=gff)

    assert rows == []
    assert "assembly_contiguity" in audit_reasons(audit)


def test_the_summary_line_carries_the_contig_count_and_n50(tmp_path, capsys):
    """The log line is what a user reads without opening a file, so '4 ICEs' and
    'off a 300-contig draft' should not be a file apart."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    run_main(tmp_path, conjscan, gff)
    printed = capsys.readouterr().out

    assert "Assembly: 1 contig(s), N50 200,000 bp." in printed
    assert "Sample S1:" in printed


def test_an_ice_resting_on_another_contigs_apparatus_says_so(tmp_path):
    """The draft-assembly exemption is allowed, but it must be visible.

    On a fragmented assembly the tra operon routinely lands on a different
    contig from the relaxase, so a cluster with no mating-pair gene of its own
    is still called an ICE when the typed system's other hits are elsewhere —
    without that exemption real ICEs are demoted the moment an assembly breaks.

    The cost is a row that contradicts itself to anyone who does not know the
    rule: mobility says "predicted self-transmissible" while missing_components
    says "mating-pair apparatus". mpf_from_other_contig is the column that
    explains it. Measured on the fragmented benchmark: 9 of 191 draft calls,
    every one already at low confidence, and 0 of 68 closed-genome calls.
    """
    contigs = {"contig_1": 200000, "contig_2": 200000}
    cds = [
        gff_cds("contig_1", 50000, 51200, "+", "S1_00010", "Phage integrase family protein"),
        gff_cds("contig_1", 55000, 56600, "+", "S1_00015", "TrwC relaxase domain-containing protein"),
        # The mating bridge, on the OTHER contig — the assembly broke between them.
        gff_cds("contig_2", 90000, 92000, "+", "S1_00620", "conjugal transfer protein TraB"),
    ]
    # One typed system, hits split across the break.
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id="S1_00015", gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id="S1_00620", gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", contigs, cds)

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan, gff, extra=["--min-element-bp", "1000"])

    element = [row for row in rows if row["contig"] == "contig_1"][0]
    # The exemption still does its job: this is an ICE, not a demoted IME...
    assert element["mge_class"] == "ice"
    assert element["has_t4ss"] == "FALSE"          # no mating-pair gene here
    # ...and it now says out loud where that apparatus actually was.
    assert element["mpf_from_other_contig"] == "TRUE"
    assert element["spans_contigs"] == "TRUE"
    assert element["confidence"] == "low"          # unchanged: spans_contigs caps it
    assert "mpf_apparatus_on_another_contig" in audit_reasons(audit)


def test_an_ice_with_its_own_apparatus_is_not_flagged(tmp_path):
    """The complement: when the mating-pair genes are in the cluster, the new
    column must be FALSE. A flag that is TRUE on every ICE tells a reader
    nothing."""
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=RELAXASE_HIT, gene_name="T4SS_MOBF",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
        conjscan_row(hit_id=VIRB4_HIT, gene_name="T4SS_virb4",
                     model_fqn="CONJScan/Chromosome/T4SS_typeF"),
    ])
    gff = write_gff(tmp_path / "S1.gff3", SCENE_CONTIGS, SCENE_CDS)

    _rc, rows, audit, _path = run_main(tmp_path, conjscan, gff)

    assert rows[0]["mge_class"] == "ice"
    assert rows[0]["mpf_from_other_contig"] == "FALSE"
    assert "mpf_apparatus_on_another_contig" not in audit_reasons(audit)


def test_a_denovo_repeat_family_is_rejected_using_the_whole_assembly(tmp_path):
    """A repeat with copies on OTHER contigs is not an att site.

    att_search already refuses a de novo repeat with more than two copies — attL
    and attR and nothing else. But it counts copies on the contig it was handed,
    and on a closed genome the contig IS the assembly, so nobody noticed the two
    questions were different. On a draft they are not: a family with 30 copies
    genome-wide can show only two on one contig.

    Measured on the fragmented benchmark: all 35 de novo repeats reported on
    drafts had exactly 2 copies on their own contig, but 10 had 3–30 across the
    assembly — one of them an 89 bp repeat with 30 copies attached to a
    high-confidence call. On closed genomes all 29 de novo repeats have 2
    assembly-wide copies, so this guard is a no-op there.

    The fixture: the same motif twice on the element's contig (the pair the
    search finds) and twice more on a second contig — the copies a fragmented
    assembly hides from a per-contig count.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=INTEGRASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS, SCENE_CDS)

    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    generator = random.Random(98)
    element_contig = "".join(generator.choice("ACGT") for _ in range(200_000))
    element_contig = plant_att_pair(element_contig, motif, 40_000, 70_000)
    # A second contig carrying two more copies of the same motif. Nothing on it
    # is a candidate; it exists only to make the repeat a family.
    other_contig = "".join(generator.choice("ACGT") for _ in range(60_000))
    other_contig = plant_att_pair(other_contig, motif, 10_000, 30_000)

    genome = tmp_path / "genome.fna"
    genome.write_text(
        ">contig_1\n" + element_contig + "\n>contig_2\n" + other_contig + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff, extra=["--genome", str(genome)])

    element = rows[0]
    # The boundary CLAIM is withdrawn...
    assert element["boundary_method"] == "none"
    assert element["attL"] == "NA"
    assert element["attR"] == "NA"
    assert element["att_sequence"] == "NA"
    assert "denovo_att_is_a_repeat_family" in audit_reasons(audit)
    # ...and nothing else moves. A de novo repeat never widened an element, so
    # rejecting one cannot change a coordinate.
    assert element["start"] == element["machinery_start"]
    assert element["end"] == element["machinery_end"]


def test_a_denovo_repeat_unique_in_the_assembly_is_still_reported(tmp_path):
    """The complement, and the proof this guard is not just switching Phase 3 off.

    Same fixture with the extra copies removed: two copies in the whole
    assembly, which is what a real attL/attR pair looks like. The repeat is
    still reported as a lead for a human to follow — and still does not move the
    element, which was always the rule for de novo boundaries.
    """
    conjscan = write_conjscan(tmp_path / "best_solution.tsv", [
        conjscan_row(hit_id=INTEGRASE_HIT, gene_name="T4SS_MOBF"),
        conjscan_row(hit_id=T4CP_HIT, gene_name="T4SS_t4cp2"),
    ])
    gff = write_gff(tmp_path / "sample.gff3", SCENE_CONTIGS, SCENE_CDS)

    motif = "GGCTCGAACCCAGGACCTCTTGCAT"
    generator = random.Random(98)
    element_contig = "".join(generator.choice("ACGT") for _ in range(200_000))
    element_contig = plant_att_pair(element_contig, motif, 40_000, 70_000)
    other_contig = "".join(generator.choice("ACGT") for _ in range(60_000))

    genome = tmp_path / "genome.fna"
    genome.write_text(
        ">contig_1\n" + element_contig + "\n>contig_2\n" + other_contig + "\n")

    _rc, rows, audit, _path = run_main(
        tmp_path, conjscan=conjscan, gff=gff, extra=["--genome", str(genome)])

    element = rows[0]
    assert element["boundary_method"] == "denovo"
    assert element["attL"] == "40000..40024"
    assert "denovo_att_is_a_repeat_family" not in audit_reasons(audit)
    assert "denovo_att_reported_not_applied" in audit_reasons(audit)
    assert element["start"] == element["machinery_start"]
