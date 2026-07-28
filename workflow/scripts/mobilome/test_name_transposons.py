"""Tests for name_transposons.py — the curated naming layer that makes tier 4 real.

WHAT IS AT STAKE HERE
    This script hands colocalise.py intervals and says "a named transposon lives
    here". Every AMR gene inside one of those intervals is then reported as cargo
    of that transposon, at tier 4. So an interval that is too WIDE is not a
    cosmetic problem: it silently converts unrelated chromosomal genes into
    "mobilisable, named architecture".

    That is not hypothetical. The first run of this script on the KPNIH1 positive
    control produced a single "Tn7246" spanning 942,502 bp — 129x the length of
    the 7,325 bp transposon itself — because two separate copies of it, a
    megabase apart, were merged into one element. Most of the tests below exist
    because of that.
"""

import name_transposons as nt


# ── Fixtures ─────────────────────────────────────────────────────────────────

def hsp(contig="contig_1", subject="Tn4401b-JX560992", pident="99.9",
        length=5000, qstart=10_000, qend=15_000, sstart=1, send=5000,
        bitscore="9000", slen=5000, qlen=500_000):
    """One BLAST tabular row, as a dict in BLAST_COLUMNS order."""
    return {
        "qseqid": contig, "sseqid": subject, "pident": str(pident),
        "length": str(length), "qstart": str(qstart), "qend": str(qend),
        "sstart": str(sstart), "send": str(send),
        "evalue": "0.0", "bitscore": str(bitscore),
        "slen": str(slen), "qlen": str(qlen),
    }


def build(hits, sample="S1", min_identity=90.0, min_coverage=0.80):
    return nt.build_elements(sample, hits, min_identity, min_coverage)


def reasons(audit_rows):
    return {row["reason"] for row in audit_rows}


# ── Defline parsing and element classification ───────────────────────────────

def test_element_name_splits_on_the_last_hyphen():
    """TnCentral writes <NAME>-<ACCESSION>, and names contain hyphens themselves,
    so only the final hyphen separates the accession."""
    assert nt.parse_element_name("Tn4401b-JX560992") == ("Tn4401b", "JX560992")
    assert nt.parse_element_name("IS1133_Tn10_IS903B-CP000602.1") == \
        ("IS1133_Tn10_IS903B", "CP000602.1")
    # Real TnCentral deflines include entries with an empty accession.
    assert nt.parse_element_name("Tn7246-") == ("Tn7246", "")
    # No hyphen at all: keep the whole thing as the name rather than lose the hit.
    assert nt.parse_element_name("Tn3000") == ("Tn3000", "")


def test_element_kinds_are_told_apart():
    assert nt.classify_element("Tn4401b") == "unit_transposon"
    assert nt.classify_element("In0") == "integron"
    assert nt.classify_element("IS26") == "insertion_sequence"
    assert nt.classify_element("weird_thing") == "unknown"


def test_a_plain_is_is_skipped_with_a_reason():
    """An IS carries only what it needs to move, so it cannot hold a passenger AMR
    gene — and ISEScan already inventories IS elements on these contigs. Naming it
    here would list one element in two tables."""
    elements, audit = build([hsp(subject="IS26-X00011")])
    assert elements == []
    assert "tncentral_hit_is_a_plain_is" in reasons(audit)


# ── THE bug: separate copies must not be merged ──────────────────────────────

def test_two_copies_of_one_transposon_stay_two_elements():
    """The KPNIH1 Tn7246 failure, reduced to its essentials.

    Two full-length copies of a 5 kb transposon, a megabase apart on the same
    contig. Taking min(qstart)/max(qend) across all HSPs would report ONE element
    spanning the whole megabase, making every gene between the two copies look
    like its cargo.
    """
    hits = [
        hsp(qstart=10_000, qend=15_000),           # copy 1
        hsp(qstart=1_010_000, qend=1_015_000),     # copy 2, a megabase away
    ]
    elements, _audit = build(hits)

    assert len(elements) == 2
    spans = sorted(int(e["end"]) - int(e["start"]) + 1 for e in elements)
    assert spans == [5001, 5001]
    assert all(e["mge_name"] == "Tn4401b" for e in elements)


def test_pieces_of_one_copy_broken_by_indels_are_still_merged():
    """The counterpart, so the split above cannot pass by simply never merging.

    A real transposon is often broken into several HSPs by internal indels. Judged
    separately each piece covers too little of the reference and would be thrown
    out by the coverage threshold — so the elements most worth naming, the big
    mosaic ones, would be exactly the ones systematically missed.
    """
    hits = [
        hsp(qstart=10_000, qend=12_000, sstart=1, send=2_000, length=2_000),
        hsp(qstart=12_100, qend=15_000, sstart=2_050, send=5_000, length=2_950),
    ]
    elements, _audit = build(hits)

    assert len(elements) == 1
    assert elements[0]["start"] == "10000"
    assert elements[0]["end"] == "15000"
    # Coverage is summed across the pieces, not taken from the best one.
    assert float(elements[0]["subject_coverage"]) > 0.95


def test_an_implausibly_wide_span_is_dropped_even_if_clustering_let_it_through():
    """Belt and braces behind the copy clustering.

    A copy of a 5 kb transposon occupies about 5 kb. If an interval reaches many
    times that, HSPs were joined that should not have been, and passing it on
    would make everything inside it look like cargo. The guard is independent of
    the clustering so that a bug in one does not silently disable the other.
    """
    # One continuous HSP claiming a 60 kb query span against a 5 kb reference.
    hits = [hsp(qstart=10_000, qend=70_000, length=5_000, sstart=1, send=5_000)]
    elements, audit = build(hits)

    assert elements == []
    assert "element_span_implausible_for_reference" in reasons(audit)


# ── The naming thresholds ────────────────────────────────────────────────────

def test_a_weak_identity_hit_does_not_earn_the_name():
    """Below the identity threshold the sequence may well be a RELATIVE of the
    element, but calling it by that name would claim more than the data shows."""
    elements, audit = build([hsp(pident="82.0")])
    assert elements == []
    assert "identity_below_naming_threshold" in reasons(audit)


def test_a_fragment_of_a_transposon_is_not_that_transposon():
    """Coverage is measured against the REFERENCE: the question is whether the
    whole known element is present, not how much of our contig it covers."""
    hits = [hsp(length=1_500, sstart=1, send=1_500, qstart=10_000, qend=11_500)]
    elements, audit = build(hits)
    assert elements == []
    assert "reference_coverage_below_threshold" in reasons(audit)


def test_coverage_is_against_the_reference_not_the_contig():
    """A 5 kb transposon fully present in a 500 kb contig covers 1% of the contig
    and 100% of itself, and it is unambiguously there."""
    elements, _audit = build([hsp(qlen=500_000)])
    assert len(elements) == 1
    assert float(elements[0]["subject_coverage"]) == 1.0


# ── Nested references ────────────────────────────────────────────────────────

def test_nested_reference_hits_collapse_to_the_best_one():
    """TnCentral entries nest on purpose — a large transposon contains smaller
    ones — so a single real element matches several references at once. Reporting
    them all would multiply one element into many in the AMR table."""
    hits = [
        hsp(subject="Tn4401b-JX560992", qstart=10_000, qend=15_000, bitscore="9000"),
        hsp(subject="Tn2-AY123456", qstart=11_000, qend=14_000, bitscore="5000",
            length=3_000, sstart=1, send=3_000, slen=3_000),
    ]
    elements, audit = build(hits)

    assert len(elements) == 1
    assert elements[0]["mge_name"] == "Tn4401b"       # the better-scoring one
    assert "nested_or_overlapping_tncentral_hit" in reasons(audit)


def test_hits_on_different_contigs_never_collapse():
    """Overlap is only meaningful within one contig; two contigs can carry copies
    at coincidentally similar coordinates."""
    hits = [hsp(contig="contig_1"), hsp(contig="contig_2")]
    elements, _audit = build(hits)
    assert len(elements) == 2


# ── The output contract with colocalise.py ───────────────────────────────────

def test_rows_carry_the_columns_colocalise_joins_on():
    """colocalise.py reads contig/start/end plus element_type and mge_name. If the
    element_type is not one it recognises, tier 4 is never awarded and the whole
    layer is silently pointless."""
    elements, _audit = build([hsp()])
    row = elements[0]
    for column in ("contig", "start", "end", "element_type", "mge_id", "mge_name"):
        assert row[column] not in ("", "NA"), column
    # These two spellings are what colocalise.py's ELEMENT_TYPE_SYNONYMS maps to
    # the tier-4 branch.
    assert row["element_type"] in ("unit_transposon", "integron")


def test_an_integron_is_typed_as_an_integron():
    elements, _audit = build([hsp(subject="In104-AY463797")])
    assert elements[0]["element_type"] == "integron"
    assert elements[0]["mge_name"] == "In104"


# ── Graceful degradation ─────────────────────────────────────────────────────

def test_no_hits_is_a_normal_result_not_an_error():
    """Most genomes carry no characterised transposon."""
    elements, audit = build([])
    assert elements == []
    assert "no_tncentral_hits" in reasons(audit)


def test_a_missing_blast_file_is_not_fatal():
    assert nt.read_blast_hits("/nonexistent/path.tsv") == []
    assert nt.read_blast_hits("") == []
