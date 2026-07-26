"""Unit tests for colocalise.py — the AMR x mobile-element co-localisation step.

Everything here runs on SYNTHETIC tables written into a temp directory: no
ISEScan, no AMRFinderPlus, no Platon, no database, no network. The point is to
pin down the mobility-ladder rules (spec §2.5) and the graceful-degradation
behaviour so they cannot drift silently.

Each test builds a tiny genome-on-paper — a contig of a known length with an AMR
gene at known coordinates and mobile elements placed around it — and asserts the
tier, the context, the measured distance and the confidence that come back.

Run: /home/antoniellil/miniconda3/envs/snakemake/bin/python -m pytest \
        workflow/scripts/mobilome/test_colocalise.py -q
"""

import csv
import os

import pytest

import colocalise as co


# ---------------------------------------------------------------------------
# Fixtures and helpers: build input files, run the CLI, read the outputs back.
# ---------------------------------------------------------------------------

# The AMRFinderPlus 4.2.7 header, verbatim (22 columns). Written into every
# synthetic AMR table so the tests exercise the real parser contract, including
# the 4.x column names that differ from the 3.x names in the older docs.
AMRFINDER_HEADER = [
    "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
    "Element name", "Scope", "Type", "Subtype", "Class", "Subclass", "Method",
    "Target length", "Reference sequence length", "% Coverage of reference",
    "% Identity to reference", "Alignment length", "Closest reference accession",
    "Closest reference name", "HMM accession", "HMM description",
]


def amr_row(contig="contig_1", start=10000, stop=11000, strand="+",
            symbol="blaTEM-116", method="ALLELEX", amr_class="BETA-LACTAM"):
    """One AMRFinderPlus row with sensible defaults; override what a test needs."""
    return [
        "s1_00001", contig, str(start), str(stop), strand, symbol,
        "broad-spectrum class A beta-lactamase", "core", "AMR", "AMR",
        amr_class, amr_class, method, "286", "286", "100.00", "100.00", "286",
        "AAA1.1", "beta-lactamase TEM-116", "NA", "NA",
    ]


def write_tsv(path, header, rows):
    """Write a header + rows as a tab-separated file."""
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    return str(path)


def read_tsv(path):
    """Read a tab-separated file back into a list of dicts (one per data row)."""
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# Column order used by every synthetic mobile-element table below. It matches
# what the module's loader is expected to hand over: one row per called element.
# Column names here must be the ones isescan_to_table.py ACTUALLY writes. This
# header used to say "family", which no real BacFlux table has ever contained -
# the real one is "is_family" - and because colocalise matches columns by name and
# silently ignores a name it cannot find, every test below passed against a
# contract that did not exist while is_family came out NA on real data and the
# IS26/IS6 tests never fired. If you add a column here, copy the name from a real
# {sample}_is_elements.tsv rather than inventing one.
IS_HEADER = ["contig", "start", "end", "strand", "is_family", "cluster",
             "is_complete", "is_id", "element_type", "name"]


def is_row(contig="contig_1", start=8000, end=9000, strand="+", family="IS3",
           cluster="IS3_1", complete="c", is_id=None, element_type="",
           name=""):
    """One mobile-element row. Defaults describe a complete IS3 on contig_1."""
    if is_id is None:
        is_id = f"{contig}_IS_{start}"
    return [contig, str(start), str(end), strand, family, cluster, complete,
            is_id, element_type, name]


REPLICON_HEADER = ["contig", "replicon", "replicon_id", "plasmid_mobility",
                   "mobility_evidence"]
LENGTH_HEADER = ["contig", "length"]

# The ICE/IME table is a SEPARATE element source, passed as a second --is-table
# exactly as rule amr_mge_colocalisation does. Column names copied from a real
# {sample}_ice_candidates.tsv. The last four are the quality judgements
# conjscan_to_ice.py makes about the element and that the mobility row must
# inherit - they do not exist in the IS table.
ICE_HEADER = ["contig", "start", "end", "strand", "mge_id", "mge_name",
              "element_type", "mobility", "machinery_intact", "degraded_reason",
              "spans_contigs", "confidence"]


def ice_row(contig="contig_1", start=5000, end=25000, strand=".",
            mge_id="ICE_1", mge_name="ICEEc2-like", element_type="ice",
            mobility="self-transmissible", machinery_intact="TRUE",
            degraded_reason="NA", spans_contigs="FALSE", confidence="high"):
    """One ICE/IME candidate row; defaults describe a clean, intact ICE."""
    return [contig, str(start), str(end), strand, mge_id, mge_name,
            element_type, mobility, machinery_intact, degraded_reason,
            spans_contigs, confidence]


def run_colocalise(tmp_path, amr_rows, is_rows=None, replicon_rows=None,
                   length_rows=(("contig_1", 50000),), max_span=20000,
                   is_table="write", replicon_table="write", ice_rows=None):
    """Write the inputs, run the CLI end to end, return (report_rows, audit_rows).

    `is_table` / `replicon_table` accept the string "absent" to test the path
    where the file was never created at all (ISEScan writes nothing when it finds
    no IS), as opposed to being written but empty.
    """
    amr_path = write_tsv(tmp_path / "amr.tsv", AMRFINDER_HEADER, amr_rows)

    if is_table == "absent":
        is_path = str(tmp_path / "missing_is.tsv")
    else:
        is_path = write_tsv(tmp_path / "is.tsv", IS_HEADER, is_rows or [])

    if replicon_table == "absent":
        replicon_path = str(tmp_path / "missing_replicons.tsv")
    else:
        replicon_path = write_tsv(tmp_path / "replicons.tsv", REPLICON_HEADER,
                                  replicon_rows or [])

    length_path = write_tsv(tmp_path / "lengths.tsv", LENGTH_HEADER,
                            [[name, str(length)] for name, length in length_rows])

    out_table = str(tmp_path / "mobility.tsv")
    out_audit = str(tmp_path / "audit.tsv")

    argv = [
        "--sample", "sampleA",
        "--amrfinder", amr_path,
        "--is-table", is_path,
    ]
    # Second element source, same as the real rule: --is-table is repeatable.
    if ice_rows is not None:
        argv += ["--is-table", write_tsv(tmp_path / "ice.tsv", ICE_HEADER, ice_rows)]
    argv += [
        "--replicons", replicon_path,
        "--contig-lengths", length_path,
        "--max-composite-span", str(max_span),
        "--out-table", out_table,
        "--out-audit", out_audit,
    ]
    co.main(argv)
    return read_tsv(out_table), read_tsv(out_audit)


def audit_reasons(audit_rows, gene=None):
    """Collect the reason slugs from the audit file, optionally for one gene."""
    return {row["reason"] for row in audit_rows
            if gene is None or row["amr_gene"] == gene}


# ---------------------------------------------------------------------------
# The pure helper functions. Cheap to test, and every ladder rule rests on them.
# ---------------------------------------------------------------------------


def test_overlap_and_gap_arithmetic_is_1_based_inclusive():
    # 100-200 and 150-300 share 150..200 = 51 bases.
    assert co.overlap_bp(100, 200, 150, 300) == 51
    # Non-overlapping intervals share nothing.
    assert co.overlap_bp(100, 200, 201, 300) == 0
    # A feature ending at 100 and the next starting at 151: bases 101..150 = 50.
    assert co.gap_bp(1, 100, 151, 300) == 50
    # Directly adjacent, and overlapping, both read as zero distance.
    assert co.gap_bp(1, 100, 101, 300) == 0
    assert co.gap_bp(1, 200, 150, 300) == 0


def test_family_matching_uses_family_then_cluster():
    is3_a = {"family": "IS3", "cluster": "IS3_1", "name": ""}
    is3_b = {"family": "is3", "cluster": "IS3_2", "name": ""}
    is21 = {"family": "IS21", "cluster": "IS21_1", "name": ""}
    assert co.families_match(is3_a, is3_b) == (True, "family")
    assert co.families_match(is3_a, is21) == (False, "family")

    # ISEScan's "new" means "could not be classified": two such elements are NOT
    # automatically the same element, so the finer cluster decides.
    new_same = {"family": "new", "cluster": "IS_cl_9", "name": ""}
    new_other = {"family": "new", "cluster": "IS_cl_4", "name": ""}
    assert co.families_match(new_same, dict(new_same)) == (True, "cluster")
    assert co.families_match(new_same, new_other) == (False, "cluster")

    # Nothing informative at all -> we must not claim the copies are the same.
    blank = {"family": "", "cluster": "", "name": ""}
    assert co.families_match(blank, blank) == (False, "unknown")


def test_is6_family_is_exempt_from_the_orientation_rule():
    assert co.is_orientation_exempt({"family": "IS6", "name": ""}) is True
    assert co.is_orientation_exempt({"family": "is26", "name": ""}) is True
    # Family column blank but the naming cascade called it IS26.
    assert co.is_orientation_exempt({"family": "new", "name": "IS26"}) is True
    assert co.is_orientation_exempt({"family": "IS3", "name": ""}) is False


def test_upstream_depends_on_the_genes_reading_direction():
    plus_gene = {"start": 1000, "end": 2000, "strand": "+"}
    minus_gene = {"start": 1000, "end": 2000, "strand": "-"}
    left_element = {"start": 500, "end": 900}
    right_element = {"start": 2100, "end": 2500}
    # A + strand gene is read left to right, so its upstream side is the LOWER
    # coordinates; for a - strand gene it is the other way round.
    assert co.is_upstream_of_gene(left_element, plus_gene) is True
    assert co.is_upstream_of_gene(right_element, plus_gene) is False
    assert co.is_upstream_of_gene(right_element, minus_gene) is True
    assert co.is_upstream_of_gene(left_element, minus_gene) is False
    # Unknown gene strand: we cannot say which side is upstream, so we do not.
    assert co.is_upstream_of_gene(left_element, {"start": 1000, "end": 2000,
                                                 "strand": "NA"}) is False


def test_confidence_can_only_be_capped_downwards():
    assert co.apply_confidence_caps("high", []) == "high"
    assert co.apply_confidence_caps("high", [("medium", "r", "d")]) == "medium"
    assert co.apply_confidence_caps("high", [("medium", "r", "d"),
                                             ("low", "r2", "d2")]) == "low"
    # A "cap" that is looser than the current level must not raise it.
    assert co.apply_confidence_caps("low", [("high", "r", "d")]) == "low"


# ---------------------------------------------------------------------------
# Ladder tier 1 — chromosomal, nothing nearby.
# ---------------------------------------------------------------------------


def test_tier1_chromosomal_no_context(tmp_path):
    """AmpC-like gene in the middle of a chromosomal contig, no IS anywhere."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="ampC", amr_class="BETA-LACTAM")],
        is_rows=[is_row(contig="contig_1", start=40000, end=41000)],  # 29 kb away
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert len(report) == 1
    row = report[0]
    assert row["mobility_tier"] == "1"
    assert row["mobility_tier_label"] == "intrinsic_candidate"
    assert row["mge_context"] == "none"
    assert row["replicon"] == "chromosome"
    assert row["confidence"] == "high"
    assert row["is_inside_amr_cds"] == "no"
    # The no-call must be explained, and the measured distance must be in it.
    assert "nearest_mobile_element_too_far" in audit_reasons(audit, "ampC")


def test_tier1_reports_the_measured_distance_to_a_nearby_but_unqualified_is(tmp_path):
    """An IS 999 bp DOWNSTREAM: real context, but not the tier-2 mechanism."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[is_row(start=12000, end=13000, strand="+")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "1"          # downstream IS is not mobilisation
    assert row["mge_context"] == "is_adjacent"  # but it IS reported as context
    assert row["distance_bp"] == "999"          # the number, next to the tier
    assert row["confidence"] == "medium"
    assert "nearby_is_not_upstream_or_not_oriented" in audit_reasons(audit)


# ---------------------------------------------------------------------------
# Ladder tier 2 — IS upstream and pointing at the gene (expression, not mobility).
# ---------------------------------------------------------------------------


def test_tier2_upstream_oriented_is_gives_expression_modulation(tmp_path):
    """ISEcp1-like: IS 199 bp upstream, same strand -> hybrid promoter candidate."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaCTX-M-15", start=10000, stop=11000, strand="+")],
        is_rows=[is_row(start=8500, end=9800, strand="+", family="ISEcp1_like")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "2"
    assert row["mobility_tier_label"] == "expression_modulation_not_mobilisation"
    assert row["mge_context"] == "is_adjacent"
    assert row["distance_bp"] == "199"
    assert row["orientation"] == "same"
    # The mechanism is inferred from coordinates only, so never "high".
    assert row["confidence"] == "medium"
    assert "promoter_inferred_from_position_only" in audit_reasons(audit)


def test_tier2_needs_the_is_on_the_genes_own_strand(tmp_path):
    """Same position, opposite strand: the IS cannot read into the gene."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[is_row(start=8500, end=9800, strand="-")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "1"
    assert row["mge_context"] == "is_adjacent"   # context, not mechanism
    assert row["orientation"] == "opposite"


def test_tier2_respects_the_promoter_distance_window(tmp_path):
    """Upstream and correctly oriented, but 1499 bp away — too far to assume it."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[is_row(start=7500, end=8500, strand="+")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "1"
    assert row["distance_bp"] == "1499"
    assert "upstream_oriented_is_beyond_promoter_window" in audit_reasons(audit)


def test_a_distant_upstream_is_does_not_clutter_the_audit(tmp_path):
    """An IS 29 kb away is not a "near miss" — that line would be pure noise."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=30000, stop=31000, strand="+")],
        is_rows=[is_row(start=500, end=1000, strand="+")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["mobility_tier"] == "1"
    reasons = audit_reasons(audit)
    assert "upstream_oriented_is_beyond_promoter_window" not in reasons
    assert "nearest_mobile_element_too_far" in reasons


def test_tier2_on_a_minus_strand_gene_looks_at_the_other_side(tmp_path):
    """For a - strand gene, 'upstream' means the HIGHER coordinates."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="-")],
        is_rows=[is_row(start=11200, end=12000, strand="-")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "2"
    assert row["distance_bp"] == "199"


# ---------------------------------------------------------------------------
# Ladder tier 3 — composite transposon, including the IS26/IS6 exception.
# ---------------------------------------------------------------------------


def test_tier3_composite_same_family_same_orientation(tmp_path):
    """Two complete IS3 copies in the same orientation, gene between them."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="aph(3')-Ia", start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", is_id="ISleft"),
            is_row(start=12000, end=13000, strand="+", family="IS3", is_id="ISright"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "3"
    assert row["mobility_tier_label"] == "composite_mobilisable_within_cell"
    assert row["mge_context"] == "composite"
    assert row["mge_id"] == "ISleft+ISright"
    assert row["is_family"] == "IS3"
    assert row["same_orientation"] == "yes"
    assert row["flanking_span_bp"] == "5001"     # 8000..13000 inclusive
    assert row["distance_bp"] == "999"           # nearest of the two flanks
    assert row["n_flanking_is"] == "2"
    assert row["confidence"] == "high"


def test_composite_rejected_when_the_two_copies_are_different_families(tmp_path):
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", cluster="IS3_1"),
            is_row(start=12000, end=13000, strand="+", family="IS21", cluster="IS21_1"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["mobility_tier"] == "1"
    assert "flanking_is_different_family" in audit_reasons(audit)


def test_composite_rejected_when_the_span_is_too_long(tmp_path):
    """Same element both sides, but 30 kb apart — beyond the configured span."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=1000, end=2000, strand="+", family="IS3"),
            is_row(start=29000, end=30000, strand="+", family="IS3"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        max_span=20000,
    )
    assert report[0]["mobility_tier"] == "1"
    assert "flanking_is_pair_span_too_long" in audit_reasons(audit)
    # ...and raising the limit turns exactly the same data into a composite call.
    report_relaxed, _ = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=1000, end=2000, strand="+", family="IS3"),
            is_row(start=29000, end=30000, strand="+", family="IS3"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        max_span=40000,
    )
    assert report_relaxed[0]["mobility_tier"] == "3"
    assert report_relaxed[0]["flanking_span_bp"] == "29001"


def test_is26_direct_orientation_is_still_called_a_composite(tmp_path):
    """IS26 flanks its cargo in DIRECT orientation — the clinically vital case."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaSHV-12", start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=9000, end=9800, strand="+", family="IS6",
                   cluster="IS6_1", is_id="IS26_a", name="IS26"),
            is_row(start=11200, end=12000, strand="+", family="IS6",
                   cluster="IS6_1", is_id="IS26_b", name="IS26"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "3"
    assert row["mge_context"] == "composite"
    assert row["is_family"] == "IS6"
    assert row["same_orientation"] == "yes"


def test_is26_is_exempt_from_the_orientation_rule_in_either_arrangement(tmp_path):
    """The exception itself: an INVERTED IS26 pair must still call composite.

    The generic rule asks for both copies in the same orientation. IS6/IS26 is
    exempted from that test altogether, so the call does not depend on which way
    round the two copies happen to sit — this is the single most important AMR
    architecture and must not be lost to an orientation technicality.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaSHV-12", start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=9000, end=9800, strand="+", family="IS6", is_id="IS26_a"),
            is_row(start=11200, end=12000, strand="-", family="IS6", is_id="IS26_b"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "3"
    assert row["mge_context"] == "composite"
    assert row["same_orientation"] == "no"      # reported honestly, not hidden
    assert "is26_orientation_exemption_applied" in audit_reasons(audit)


def test_non_is6_inverted_pair_is_rejected_and_the_reason_is_audited(tmp_path):
    """Control for the exception: the same geometry in family IS3 is turned down.

    This is where the convention bites. Some genuine composite transposons (Tn10,
    Tn5) do carry their flanking IS in INVERTED orientation, so this rejection is
    a deliberately conservative convention, not a biological law — which is
    exactly why the rejected pair is written to the audit file with its measured
    coordinates instead of being dropped.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", is_id="ISleft"),
            is_row(start=12000, end=13000, strand="-", family="IS3", is_id="ISright"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["mobility_tier"] == "1"
    reasons = audit_reasons(audit)
    assert "flanking_is_pair_inverted_orientation" in reasons
    detail = [row["detail"] for row in audit
              if row["reason"] == "flanking_is_pair_inverted_orientation"][0]
    assert "ISleft" in detail and "ISright" in detail


def test_partial_flanking_is_lowers_composite_confidence(tmp_path):
    """A partial (fragmentary) IS call still counts, but only at medium."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", complete="c"),
            is_row(start=12000, end=13000, strand="+", family="IS3", complete="p"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["mobility_tier"] == "3"
    assert report[0]["confidence"] == "medium"
    assert "flanking_is_partial" in audit_reasons(audit)


def test_the_tightest_flanking_pair_wins(tmp_path):
    """With nested candidate pairs, the innermost boundary is the call."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=2000, end=3000, strand="+", family="IS3", is_id="far_left"),
            is_row(start=9000, end=9500, strand="+", family="IS3", is_id="near_left"),
            is_row(start=11500, end=12000, strand="+", family="IS3", is_id="near_right"),
            is_row(start=18000, end=19000, strand="+", family="IS3", is_id="far_right"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mge_id"] == "near_left+near_right"
    assert row["flanking_span_bp"] == "3001"    # 9000..12000
    assert row["n_flanking_is"] == "4"          # all four are reported as context


# ---------------------------------------------------------------------------
# Ladder tier 4 — a curated named element containing the gene, which overrides
# the pattern-based composite call.
# ---------------------------------------------------------------------------


def test_tier4_named_unit_transposon_overrides_the_composite_pattern(tmp_path):
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaTEM-1", start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", is_id="ISleft"),
            is_row(start=12000, end=13000, strand="+", family="IS3", is_id="ISright"),
            is_row(start=7500, end=14000, strand="+", family="", cluster="",
                   is_id="Tn3_hit", element_type="unit_transposon", name="Tn3"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "4"
    assert row["mge_context"] == "unit_transposon"
    assert row["mge_name"] == "Tn3"
    assert row["mge_id"] == "Tn3_hit"
    assert row["confidence"] == "high"
    # The pattern-based call it displaced is recorded rather than discarded.
    assert "named_element_overrides_composite_pattern" in audit_reasons(audit)


def test_tier4_integron_cassette(tmp_path):
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="aadA2", start=10000, stop=10800, strand="+")],
        is_rows=[is_row(start=9000, end=15000, family="", cluster="",
                        is_id="In2_hit", element_type="integron", name="In2")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "4"
    assert row["mge_context"] == "integron"
    assert row["mge_name"] == "In2"


def test_unrecognised_element_type_is_ignored_and_audited(tmp_path):
    """An element type we do not understand must not quietly drive a tier."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        is_rows=[is_row(start=9000, end=15000, is_id="odd",
                        element_type="something_new")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["mobility_tier"] == "1"
    assert report[0]["mge_context"] == "none"
    assert "unknown_element_type" in audit_reasons(audit)


# ---------------------------------------------------------------------------
# Ladder tiers 5 and 6 — replicon-driven, plus ICE/IME.
# ---------------------------------------------------------------------------


def test_tier5_on_a_mobilisable_plasmid(tmp_path):
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="contig_7", symbol="tetA", start=3000, stop=4000)],
        is_rows=[],
        replicon_rows=[["contig_7", "plasmid", "pA", "mobilisable", "MOBQ"]],
        length_rows=(("contig_7", 40000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "5"
    assert row["mobility_tier_label"] == "mobilisable_needs_helper"
    assert row["mge_context"] == "plasmid"
    assert row["replicon"] == "plasmid"
    assert row["mge_id"] == "pA"
    assert row["plasmid_mobility"] == "mobilisable"
    assert row["mobility_evidence"] == "MOBQ"


def test_tier5_untyped_plasmid_is_capped_at_medium(tmp_path):
    """A plasmid with no CONJscan typing is ACQUIRED but not shown transferable."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="contig_7", start=3000, stop=4000)],
        is_rows=[],
        replicon_rows=[["contig_7", "plasmid", "pA", "NA", "NA"]],
        length_rows=(("contig_7", 40000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "5"
    assert row["plasmid_mobility"] == "unknown"
    assert row["confidence"] == "medium"
    assert "plasmid_mobility_untyped" in audit_reasons(audit)


def test_tier6_conjugative_plasmid_is_predicted_self_transmissible(tmp_path):
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="contig_7", symbol="blaCTX-M-15",
                          start=3000, stop=4000)],
        is_rows=[],
        replicon_rows=[["contig_7", "plasmid", "pB", "conjugative",
                        "MOBF;T4SS_typeF"]],
        length_rows=(("contig_7", 90000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"
    # Language discipline: "predicted", never a bare "transmissible" (spec §2.6).
    assert row["mobility_tier_label"] == "predicted_self_transmissible"
    assert row["mge_context"] == "plasmid"
    assert row["plasmid_mobility"] == "conjugative"
    # NOT high. This tier rests on Platon's conjugation HMM hit COUNT, which does
    # not establish that a complete mating apparatus is present - and on the
    # chromosome the same evidence is refused an ICE call and marked "machinery
    # incomplete". The tier stands; the certainty is capped, and audited.
    assert row["confidence"] == "medium"
    assert "plasmid_conjugation_from_hit_counts_only" in audit_reasons(_audit)


def test_tier6_ice_on_a_chromosome_beats_the_naive_chromosomal_assumption(tmp_path):
    """An ICE carries its own conjugation machinery: chromosomal is not safe."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="tetM", start=10000, stop=11000)],
        is_rows=[is_row(start=5000, end=25000, family="", cluster="",
                        is_id="ICE_1", element_type="ice", name="ICEEc2-like")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"
    assert row["mge_context"] == "ice"
    assert row["mge_name"] == "ICEEc2-like"
    assert row["replicon"] == "chromosome"      # chromosomal AND tier 6


def test_ime_is_tier5_not_tier6(tmp_path):
    """An IME has a relaxase but no conjugation machinery — it needs a helper."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        is_rows=[is_row(start=5000, end=25000, family="", cluster="",
                        is_id="IME_1", element_type="ime", name="IME-1")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "5"
    assert row["mge_context"] == "ime"


def test_ice_wins_over_a_conjugative_plasmid_when_both_apply(tmp_path):
    """Same tier, but the interval-level call is the more specific answer."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="contig_7", start=10000, stop=11000)],
        is_rows=[is_row(contig="contig_7", start=5000, end=25000, family="",
                        cluster="", is_id="ICE_1", element_type="ice")],
        replicon_rows=[["contig_7", "plasmid", "pB", "conjugative", "MOBF"]],
        length_rows=(("contig_7", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"
    assert row["mge_context"] == "ice"


# ---------------------------------------------------------------------------
# The inactivation special case — reported, never counted as mobilisation.
# ---------------------------------------------------------------------------


def test_is_inside_the_amr_cds_is_flagged_as_inactivation_not_mobilisation(tmp_path):
    """An IS that landed inside the gene usually breaks it: report it separately."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaTEM", start=10000, stop=11000, strand="+")],
        # 10200-10900 = 701 bp of the gene's 1001 bp -> above the 0.5 threshold.
        is_rows=[is_row(start=10200, end=10900, strand="+", family="IS1",
                        is_id="IS_in_gene")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "yes"
    assert row["is_amr_overlap_bp"] == "701"
    # Crucially: this did NOT raise the mobility tier or set a mobilisation context.
    assert row["mobility_tier"] == "1"
    assert row["mge_context"] == "none"
    assert "is_inside_amr_cds_likely_inactivation" in audit_reasons(audit)


def test_small_overlap_is_reported_but_not_called_a_disruption(tmp_path):
    """Below the 0.5 overlap threshold: show the number, do not make the call."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        # 10900-11500 overlaps only 101 bp of the 1001 bp gene.
        is_rows=[is_row(start=10900, end=11500, strand="+", family="IS1")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "no"
    assert row["is_amr_overlap_bp"] == "101"


def test_a_disrupted_gene_can_still_sit_in_a_mobile_block(tmp_path):
    """The two questions are independent: 'is it broken' and 'can it move'.

    A resistance gene knocked out by an IS insertion can still be flanked by a
    composite transposon, so the inactivation flag must not suppress the tier,
    and the tier must not hide the inactivation.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaTEM", start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=10200, end=10900, strand="+", family="IS1",
                   cluster="IS1_1", is_id="IS_in_gene"),
            is_row(start=8000, end=9000, strand="+", family="IS3", is_id="ISleft"),
            is_row(start=12000, end=13000, strand="+", family="IS3", is_id="ISright"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "yes"
    assert row["mobility_tier"] == "3"
    assert row["mge_context"] == "composite"
    assert "is_inside_amr_cds_likely_inactivation" in audit_reasons(audit)


def test_stress_and_virulence_rows_are_kept_and_labelled(tmp_path):
    """Mercury resistance (STRESS) inside a transposon is classic Tn21 cargo.

    Those rows are reported alongside AMR rows rather than dropped, with the
    AMRFinderPlus Type in its own column so a reader can filter to AMR only.
    """
    stress = amr_row(symbol="merA", start=10000, stop=11000)
    stress[8] = "STRESS"          # the "Type" column
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="ampC"), stress],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    types = {row["amr_gene"]: row["amr_element_type"] for row in report}
    assert types == {"ampC": "AMR", "merA": "STRESS"}


# ---------------------------------------------------------------------------
# The short-read honesty signals: contig ends and partial hits.
# ---------------------------------------------------------------------------


def test_gene_near_a_contig_end_is_flagged_and_capped_at_low(tmp_path):
    """The gene is 50 bp from the end: a missing flank may just be a missing contig."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="ampC", start=900, stop=1150, strand="+")],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 1200),),
    )
    row = report[0]
    assert row["dist_to_contig_end"] == "50"
    assert row["is_at_contig_boundary"] == "yes"
    assert row["confidence"] == "low"
    assert "amr_gene_near_contig_end" in audit_reasons(audit)


def test_partial_at_contig_end_method_caps_confidence_at_low(tmp_path):
    """AMRFinderPlus's PARTIAL_CONTIG_ENDX: the gene runs off the assembly."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="blaTEM", start=10000, stop=11000,
                          method="PARTIAL_CONTIG_ENDX")],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["amr_partial_at_contig_end"] == "yes"
    assert row["confidence"] == "low"
    assert "amr_hit_partial_at_contig_end" in audit_reasons(audit)


def test_a_bare_partial_method_is_not_the_contig_end_flag(tmp_path):
    """PARTIALX means truncated but INTERNAL — a pseudogene, not fragmentation."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(method="PARTIALX")],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report[0]["amr_partial_at_contig_end"] == "no"


def test_an_element_touching_a_contig_end_caps_the_call_at_low(tmp_path):
    """The composite is called, but one flank sits at the edge of the assembly."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="+")],
        is_rows=[
            is_row(start=8000, end=9000, strand="+", family="IS3", is_id="ISleft"),
            is_row(start=12000, end=12400, strand="+", family="IS3", is_id="ISright"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 12500),),   # ISright ends 100 bp from the end
    )
    row = report[0]
    assert row["mobility_tier"] == "3"
    assert row["confidence"] == "low"
    assert "context_element_at_contig_end" in audit_reasons(audit)


def test_an_element_reported_on_two_contigs_is_capped_at_low(tmp_path):
    """Anything spanning contigs is capped at low regardless of everything else."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        is_rows=[
            is_row(contig="contig_1", start=5000, end=25000, family="", cluster="",
                   is_id="ICE_1", element_type="ice"),
            is_row(contig="contig_2", start=1, end=4000, family="", cluster="",
                   is_id="ICE_1", element_type="ice"),
        ],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000), ("contig_2", 4000)),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"
    assert row["spans_contigs"] == "yes"
    assert row["confidence"] == "low"
    assert "mge_spans_contigs" in audit_reasons(audit)


def test_unknown_contig_length_caps_confidence_and_is_audited(tmp_path):
    """Without a length we cannot check the main short-read caveat at all."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="contig_99", start=10000, stop=11000)],
        is_rows=[],
        replicon_rows=[["contig_99", "chromosome", "contig_99", "NA", "NA"]],
        length_rows=(("contig_1", 50000),),   # contig_99 is not in the table
    )
    row = report[0]
    assert row["dist_to_contig_end"] == "NA"
    assert row["is_at_contig_boundary"] == "NA"
    assert row["confidence"] == "medium"
    assert "contig_length_unknown" in audit_reasons(audit)


# ---------------------------------------------------------------------------
# Graceful degradation: every input except the AMR table may be missing.
# ---------------------------------------------------------------------------


def test_empty_is_table_still_produces_rows_and_says_why(tmp_path):
    """ISEScan found no IS (it then writes no files at all): not an error."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="ampC"), amr_row(symbol="mexE", start=20000,
                                                  stop=21000)],
        is_rows=[],   # header only
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert len(report) == 2
    for row in report:
        assert row["mobility_tier"] == "1"
        assert row["mge_context"] == "none"
        # "No context" here means "not looked at", so it cannot be high confidence.
        assert row["confidence"] == "medium"
    reasons = audit_reasons(audit)
    assert "is_table_absent_or_empty" in reasons
    assert "no_is_calls_available" in reasons


def test_absent_is_file_behaves_like_an_empty_one(tmp_path):
    """The path does not exist at all — the ISEScan no-output case."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row()],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        is_table="absent",
    )
    assert report[0]["mobility_tier"] == "1"
    assert report[0]["confidence"] == "medium"
    assert "is_table_absent_or_empty" in audit_reasons(audit)


def test_empty_replicon_table_leaves_the_replicon_unknown(tmp_path):
    """No Platon call: we must NOT default the gene to 'chromosomal, intrinsic'."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="ampC")],
        is_rows=[],
        replicon_rows=[],   # header only
    )
    row = report[0]
    assert row["replicon"] == "unknown"
    assert row["plasmid_mobility"] == "NA"
    assert row["mobility_tier"] == "1"
    assert row["confidence"] == "medium"     # never "high" without a replicon call
    reasons = audit_reasons(audit)
    assert "replicon_table_absent_or_empty" in reasons
    assert "replicon_call_unavailable" in reasons


def test_absent_replicon_file_behaves_like_an_empty_one(tmp_path):
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row()],
        is_rows=[],
        replicon_table="absent",
    )
    assert report[0]["replicon"] == "unknown"
    assert "replicon_table_absent_or_empty" in audit_reasons(audit)


def test_every_no_context_gene_has_at_least_one_audit_reason(tmp_path):
    """The BacFlux audit convention, checked directly."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[
            amr_row(symbol="ampC", start=10000, stop=11000),
            amr_row(symbol="mexE", start=30000, stop=31000),
            amr_row(symbol="emhB", start=40000, stop=41000),
        ],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    no_context_genes = {row["amr_gene"] for row in report
                        if row["mge_context"] == "none"}
    assert no_context_genes == {"ampC", "mexE", "emhB"}
    explained = {row["amr_gene"] for row in audit
                 if row["decision"] == "no_mge_context"}
    assert no_context_genes <= explained


def test_outputs_always_have_a_header_even_with_no_amr_genes(tmp_path):
    """An isolate with no AMR hits gets empty tables, not missing files."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    assert report == []
    # The audit still carries the "which inputs were empty" lines.
    assert any(row["decision"] == "input_missing" for row in audit)
    with open(tmp_path / "mobility.tsv", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header == co.OUTPUT_COLUMNS


# ---------------------------------------------------------------------------
# Parser robustness.
# ---------------------------------------------------------------------------


def test_amr_row_without_coordinates_is_not_assessable(tmp_path):
    """Protein-only AMRFinderPlus output has no contig: say so, do not guess."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(contig="NA", start=0, stop=0)],
        is_rows=[],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["mobility_tier"] == "NA"
    assert row["mobility_tier_label"] == "not_assessable"
    assert row["confidence"] == "low"
    assert "no_contig_coordinates" in audit_reasons(audit)


def test_legacy_amrfinderplus_3x_header_is_rejected_loudly(tmp_path):
    """A 3.x file would silently parse as 'no AMR genes' — stop instead."""
    legacy_header = list(AMRFINDER_HEADER)
    legacy_header[legacy_header.index("Element symbol")] = "Gene symbol"
    path = write_tsv(tmp_path / "legacy.tsv", legacy_header, [amr_row()])
    with pytest.raises(ValueError, match="3.x"):
        co.parse_amrfinder(path)


def test_missing_amrfinder_file_is_a_hard_error(tmp_path):
    """The AMR table defines the rows: without it there is nothing to report."""
    with pytest.raises(ValueError, match="required"):
        co.parse_amrfinder(str(tmp_path / "does_not_exist.tsv"))


def test_mobile_element_rows_with_unusable_coordinates_are_skipped(tmp_path):
    """One bad line must not take the whole sample down."""
    path = write_tsv(tmp_path / "is.tsv", IS_HEADER, [
        is_row(start=8000, end=9000),
        ["contig_1", "not_a_number", "9000", "+", "IS3", "IS3_1", "c", "bad", "", ""],
    ])
    elements = co.parse_mobile_elements(path)
    assert len(elements) == 1
    assert elements[0]["start"] == 8000


def test_isescan_native_column_names_are_accepted(tmp_path):
    """seqID/isBegin/isEnd and the c/p `type` column, straight from ISEScan."""
    path = write_tsv(tmp_path / "isescan.tsv",
                     ["seqID", "family", "cluster", "isBegin", "isEnd", "strand", "type"],
                     [["contig_1", "IS3", "IS3_1", "8000", "9000", "+", "c"],
                      ["contig_1", "IS6", "IS6_1", "12000", "13000", "-", "p"]])
    elements = co.parse_mobile_elements(path)
    assert [element["start"] for element in elements] == [8000, 12000]
    assert [element["complete"] for element in elements] == ["complete", "partial"]
    # No element_type column -> everything is an insertion sequence.
    assert {element["element_type"] for element in elements} == {"insertion_sequence"}


def test_contig_lengths_accept_a_headerless_samtools_fai(tmp_path):
    fai_path = tmp_path / "contigs.fasta.fai"
    fai_path.write_text("contig_1\t50000\t12\t60\t61\ncontig_2\t4000\t51000\t60\t61\n",
                        encoding="utf-8")
    lengths = co.parse_contig_lengths(str(fai_path))
    assert lengths == {"contig_1": 50000, "contig_2": 4000}


def test_reversed_element_coordinates_are_repaired(tmp_path):
    path = write_tsv(tmp_path / "is.tsv", IS_HEADER,
                     [is_row(start=9000, end=8000)])
    element = co.parse_mobile_elements(path)[0]
    assert (element["start"], element["end"]) == (8000, 9000)


# ---------------------------------------------------------------------------
# Regression anchor against the real AMRFinderPlus output for sample 006, if the
# shared test-data file is present. Sample 006 is a Pseudomonas_E with five
# chromosomal, intrinsic-looking AMR genes (efflux pumps + AmpC) and nothing
# acquired — so the honest answer is five tier-1 rows.
# ---------------------------------------------------------------------------

REAL_AMRFINDER = os.path.join(os.path.dirname(__file__), "testdata",
                              "amrfinderplus_006_real.tsv")


@pytest.mark.skipif(not os.path.exists(REAL_AMRFINDER),
                    reason="shared testdata/amrfinderplus_006_real.tsv not present")
def test_real_sample_006_is_five_intrinsic_candidates(tmp_path):
    out_table = str(tmp_path / "mobility.tsv")
    out_audit = str(tmp_path / "audit.tsv")
    lengths = write_tsv(tmp_path / "lengths.tsv", LENGTH_HEADER,
                        [["contig_1", "5000000"]])
    replicons = write_tsv(tmp_path / "replicons.tsv", REPLICON_HEADER,
                          [["contig_1", "chromosome", "contig_1", "NA", "NA"]])
    co.main([
        "--sample", "006",
        "--amrfinder", REAL_AMRFINDER,
        "--replicons", replicons,
        "--contig-lengths", lengths,
        "--out-table", out_table,
        "--out-audit", out_audit,
    ])
    report = read_tsv(out_table)
    assert len(report) == 5
    assert {row["amr_gene"] for row in report} == {"emhC", "emhB", "mexE",
                                                   "aac(6')", "ampC"}
    assert {row["mobility_tier"] for row in report} == {"1"}
    assert {row["replicon"] for row in report} == {"chromosome"}
    # No IS table was supplied at all, so none of these may claim high confidence.
    assert {row["confidence"] for row in report} == {"medium"}


# ---------------------------------------------------------------------------
# The ICE/IME element's OWN verdict must reach the AMR row.
#
# conjscan_to_ice.py already decides whether a candidate's machinery is intact,
# whether it straddles a contig break, and what confidence it deserves. Those
# judgements used to be computed and then dropped here, so an AMR gene sitting in
# a half-broken, contig-spanning element was still reported as tier 6 at high
# confidence — the strongest claim the module can make, resting on an element the
# module itself had flagged as doubtful. These tests pin that shut.
# ---------------------------------------------------------------------------


def test_degraded_ice_machinery_caps_the_genes_confidence(tmp_path):
    """A decayed ICE must not yield a high-confidence self-transmissible call."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="tetM", start=10000, stop=11000)],
        ice_rows=[ice_row(machinery_intact="FALSE",
                          degraded_reason="low_system_wholeness",
                          mobility="self-transmissible - machinery incomplete",
                          confidence="medium")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"          # the tier itself is still right
    assert row["confidence"] != "high"          # but the certainty is not
    assert "mge_machinery_not_intact" in audit_reasons(audit)


def test_ice_spanning_contigs_is_capped_at_low_and_says_so(tmp_path):
    """Spec §8 phase 6: spanning contigs caps at low, whatever else is true.

    The generic cross-contig detector cannot see this case — an ICE id embeds its
    contig and each element is one row — so the element's own flag is the only
    evidence there is.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="tetM", start=10000, stop=11000)],
        ice_rows=[ice_row(spans_contigs="TRUE", confidence="low")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["confidence"] == "low"
    assert row["spans_contigs"] == "yes"        # must not claim "no"
    assert "mge_spans_contigs" in audit_reasons(audit)


def test_a_clean_intact_ice_still_earns_high_confidence(tmp_path):
    """The counterpart: the caps must not fire on a good element."""
    report, _audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="tetM", start=10000, stop=11000)],
        ice_rows=[ice_row()],                   # intact, single contig, high
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "6"
    assert row["mobility_tier_label"] == "predicted_self_transmissible"
    assert row["confidence"] == "high"
    assert row["spans_contigs"] == "no"


def test_a_non_mobilisable_plasmid_is_not_labelled_mobilisable(tmp_path):
    """Tier 5 is 'on a plasmid'; its generic label must not contradict the typing.

    The gene IS acquired (that is the regulatory question), so the tier stays 5,
    but calling a plasmid we just typed non-mobilisable "mobilisable_needs_helper"
    contradicts our own evidence in the column most people read.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        replicon_rows=[["contig_1", "plasmid", "p1", "non-mobilisable",
                        "no mobility genes detected by Platon"]],
        length_rows=(("contig_1", 50000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "5"
    assert row["mobility_tier_label"] != "mobilisable_needs_helper"
    assert row["mobility_tier_label"] == "on_plasmid_typed_non_mobilisable"
    assert "plasmid_typed_non_mobilisable" in audit_reasons(audit)


def test_an_is_partly_overlapping_the_gene_is_not_reported_as_clean_context(tmp_path):
    """Below the disruption threshold, an overlapping IS used to vanish entirely.

    It is not a disruption call, and the nearest-IS search skips anything that
    overlaps, so the row read "intrinsic candidate, nothing nearby" at high
    confidence while an IS was physically sitting on the gene.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],          # 1001 bp gene
        is_rows=[is_row(start=9800, end=10200)],              # ~20% overlap
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "no"     # correctly NOT a disruption call
    assert int(row["is_amr_overlap_bp"]) > 0    # the overlap is still measured
    assert row["confidence"] != "high"          # but the context is not clean
    assert "is_partially_overlaps_amr_gene" in audit_reasons(audit)


# ---------------------------------------------------------------------------
# "Nothing found" must mean "looked and found nothing", never "did not look".
# ---------------------------------------------------------------------------


def test_an_ice_row_does_not_stand_in_for_missing_is_calls(tmp_path):
    """One ICE row used to make the module believe ISEScan had run.

    ISEScan legitimately writes nothing when a genome has no IS, so the module
    caps confidence when no IS calls exist. That check counted the POOLED rows of
    both element tables, so a single unrelated ICE row silently satisfied it and
    the same genome jumped from medium to high confidence on identical IS
    evidence (none).
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        is_rows=[],                                   # ISEScan found no IS
        ice_rows=[ice_row(start=400000, end=420000)],  # unrelated, far away
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 500000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "1"
    assert row["confidence"] != "high"
    assert "no_is_calls_available" in audit_reasons(audit)


def test_a_gene_inside_a_conjugative_region_is_not_an_intrinsic_candidate(tmp_path):
    """conjugative_region / genomic_island were unrecognised and vanished.

    They deliberately do not raise the tier (spec §8 phase 4 — a conjugative
    region with no integrase has no established boundaries). But being unmapped
    made them invisible to every test, so the gene came back "intrinsic candidate,
    high confidence" with a contradictory "unrecognised element_type" audit line.
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        ice_rows=[ice_row(start=5000, end=25000, mge_id="REGION_1",
                          element_type="conjugative_region",
                          mobility="conjugative region, boundaries not established")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "1"              # correctly NOT raised
    assert row["mge_context"] == "conjugative_region"
    assert row["mge_id"] == "REGION_1"
    assert row["confidence"] != "high"
    assert "inside_context_only_element" in audit_reasons(audit)


def test_a_gene_just_outside_an_ime_is_reported_with_its_distance(tmp_path):
    """The nearest-element scan was insertion-sequence-only.

    The interval an ICE/IME carries is the span of its machinery genes, not a
    resolved boundary, so a gene sitting just outside it may still be carried by
    the element. It used to read "intrinsic candidate, nothing nearby".
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        ice_rows=[ice_row(start=13000, end=33000, mge_id="IME_9",
                          element_type="ime")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["mobility_tier"] == "1"
    assert row["mge_id"] == "IME_9"
    assert int(row["distance_bp"]) > 0              # the measured gap is reported
    assert row["confidence"] != "high"
    assert "near_but_outside_element_machinery_span" in audit_reasons(audit)


def test_the_no_context_audit_line_only_claims_what_was_checked(tmp_path):
    """It said "no mobile element on this contig" having checked only IS calls."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000)],
        is_rows=[is_row(start=300000, end=301000)],   # far away, so no context
        ice_rows=[ice_row(start=400000, end=420000, mge_id="ICE_FAR")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 500000),),
    )
    reasons = audit_reasons(audit)
    assert "no_mobile_element_on_this_contig" not in reasons


# ---------------------------------------------------------------------------
# The IS-into-gene architecture, as it actually appears in the data.
# ---------------------------------------------------------------------------


def test_an_is_abutting_a_partial_hit_is_called_inactivation_not_expression(tmp_path):
    """The canonical architecture the >=50% overlap test could never catch.

    When an IS lands in a resistance gene, AMRFinderPlus reports only the
    surviving fragment, so the gene's coordinates STOP where the IS starts and the
    overlap is ~0. Taken from the reviewer's real-world example: a PARTIALX
    aac(3)-IIa fragment with an IS26 flush against its upstream side. This used to
    come out as tier 2 "expression modulation" — the opposite conclusion (gene
    working harder) from the true one (gene broken).
    """
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="aac(3)-IIa", start=10000, stop=10450,
                          strand="-", method="PARTIALX")],
        is_rows=[is_row(start=10451, end=11670, strand="-", family="IS6",
                        cluster="IS6_1", is_id="IS26_a")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "yes"
    assert row["mobility_tier"] != "2"          # NOT expression modulation
    assert "is_abuts_partial_amr_hit_likely_inactivation" in audit_reasons(audit)


def test_a_partial_hit_truncated_by_the_contig_end_is_not_called_inactivation(tmp_path):
    """The guard: AMRFinderPlus says when the ASSEMBLY, not an IS, did the cutting."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(symbol="aac(3)-IIa", start=10000, stop=10450,
                          strand="-", method="PARTIAL_CONTIG_END")],
        is_rows=[is_row(start=10451, end=11670, strand="-", family="IS6")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "no"
    assert "is_abuts_partial_amr_hit_likely_inactivation" not in audit_reasons(audit)


def test_a_complete_gene_next_to_an_is_is_still_ordinary_context(tmp_path):
    """The other guard: a full-coverage hit beside an IS is not an inactivation."""
    report, audit = run_colocalise(
        tmp_path,
        amr_rows=[amr_row(start=10000, stop=11000, strand="-", method="EXACTX")],
        is_rows=[is_row(start=11001, end=12200, strand="-", family="IS6")],
        replicon_rows=[["contig_1", "chromosome", "contig_1", "NA", "NA"]],
        length_rows=(("contig_1", 200000),),
    )
    row = report[0]
    assert row["is_inside_amr_cds"] == "no"
    assert "is_abuts_partial_amr_hit_likely_inactivation" not in audit_reasons(audit)
