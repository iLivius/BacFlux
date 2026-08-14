"""Tests for name_ice_elements.py — putting curated ICEberg names on ICE calls.

This layer labels elements; it does not decide them. conjscan_to_ice.py works out
what is an ICE from machinery evidence, and nothing here may create, move, drop
or re-type an element. Several tests below exist purely to hold that line — in
particular that a `conjugative_region`, which has no integrase and is explicitly
NOT an ICE, can never be handed an ICE name.

The BLAST rows are written out by hand below, so none of this needs an ICEberg
download or a blastn run. The naming layer itself is off unless a user configures
mobilome.iceberg.urls or mobilome.iceberg.dir; ICEberg publishes no licence, so
BacFlux ships the URL and never the data.

Run: pytest workflow/scripts/80_mobilome/test_name_ice_elements.py -q
"""

import name_ice_elements as ni


# ── Building an ICEberg BLAST row and an ICE candidate row ───────────────────

def iceberg_hit(contig="contig_1", name="ICEKp1", accession="CP012345.1",
                pident="99.5", length=50_000, qstart=100_000, qend=150_000,
                sstart=1, send=50_000, bitscore="90000", slen=50_000):
    """One BLAST row against an ICEberg-shaped defline.

    q* fields are our contig, s* fields the curated reference, and slen its full
    length — the three the naming decision turns on.
    """
    subject = f"ICEberg|1174|{name}|GenBank|{accession}|{sstart}..{send}"
    return {
        "qseqid": contig, "sseqid": subject, "pident": str(pident),
        "length": str(length), "qstart": str(qstart), "qend": str(qend),
        "sstart": str(sstart), "send": str(send),
        "evalue": "0.0", "bitscore": str(bitscore),
        "slen": str(slen), "qlen": "5000000",
    }


def ice_row(mge_id="contig_1|ice-100000:150000", contig="contig_1",
            start=100_000, end=150_000, element_type="ice"):
    return {
        "sample": "S1", "mge_id": mge_id, "contig": contig,
        "start": str(start), "end": str(end),
        "element_type": element_type, "mge_name": "NA",
    }


def reasons(audit_rows):
    return {row["reason"] for row in audit_rows}


# ── Pulling the element name out of an ICEberg defline ───────────────────────

def test_the_element_name_is_the_third_pipe_field():
    """ICEberg deflines are pipe-delimited and regular. Names contain hyphens
    themselves, which is why this splits on the pipe."""
    defline = ("ICEberg|1174|ICEKpnATCCBAA-2146-1|GenBank|CP006659.2|"
               "4603840..4661887 Klebsiella pneumoniae ...")
    assert ni.parse_iceberg_name(defline) == ("ICEKpnATCCBAA-2146-1", "CP006659.2")


def test_an_unexpected_defline_keeps_the_hit():
    """A malformed defline should cost us the accession, not the element."""
    name, accession = ni.parse_iceberg_name("something_odd")
    assert name == "something_odd"
    assert accession == ""


# ── When a curated name is earned, and when it is not ────────────────────────
# Two independent measures have to agree before a bare name is given: how much of
# OUR element the hit covers (a curated element clipping our edge is not the same
# element) and how much of the REFERENCE is present inside our interval (less than
# most of it, and the element is only "-like"). Both thresholds are the naming
# cascade's convention, not a biological boundary.

def test_an_overlapping_curated_element_supplies_the_name():
    rows, audit = ni.name_elements("S1", [ice_row()], [iceberg_hit()])
    assert rows[0]["mge_name"] == "ICEKp1"
    assert rows[0]["iceberg_accession"] == "CP012345.1"
    assert "iceberg_match" in reasons(audit)


def test_a_partial_match_is_named_like_rather_than_named():
    """Below 80% of the REFERENCE present, the element is clearly related to the
    curated one but is not the whole of it, so the bare name would overclaim."""
    hit = iceberg_hit(length=30_000, slen=100_000)   # 30% of the reference
    rows, _audit = ni.name_elements("S1", [ice_row()], [hit])
    assert rows[0]["mge_name"] == "ICEKp1-like"


def test_a_weak_identity_hit_confers_no_name():
    rows, audit = ni.name_elements("S1", [ice_row()], [iceberg_hit(pident="62.0")])
    assert rows[0]["mge_name"] == "NA"
    assert "no_curated_name_found" in reasons(audit)


def test_a_hit_that_barely_touches_the_element_confers_no_name():
    """Overlap is measured as a fraction of OUR candidate: a curated element
    clipping the edge of ours is not the same element."""
    hit = iceberg_hit(qstart=148_000, qend=152_000, length=4_000)
    rows, audit = ni.name_elements("S1", [ice_row()], [hit])
    assert rows[0]["mge_name"] == "NA"
    assert "no_curated_name_found" in reasons(audit)


def test_hits_on_another_contig_are_not_used():
    rows, _audit = ni.name_elements(
        "S1", [ice_row(contig="contig_1")], [iceberg_hit(contig="contig_2")])
    assert rows[0]["mge_name"] == "NA"


# ── Holding the line: this layer only labels ─────────────────────────────────

def test_a_conjugative_region_never_receives_an_ice_name():
    """The single most important restraint here.

    A conjugative region has a relaxase and a mating apparatus but NO integrase,
    so it is deliberately not called an ICE (spec §8 phase 4: "report it, do not
    call it an ICE"). Hanging an ICE name on it would quietly undo the
    classification the module went to some trouble to make.
    """
    row = ice_row(element_type="conjugative_region")
    rows, audit = ni.name_elements("S1", [row], [iceberg_hit()])
    assert rows[0]["mge_name"] == "NA"
    assert rows[0]["element_type"] == "conjugative_region"
    assert "element_type_not_nameable" in reasons(audit)


def test_naming_never_changes_coordinates_or_type():
    """Even on a perfect match, the element's own extent and class are untouched —
    the curated record is used for its NAME, not to redraw our call."""
    row = ice_row(start=100_000, end=150_000)
    # A curated element much larger than ours, which is what really happens on
    # the K. pneumoniae positive control.
    hit = iceberg_hit(sstart=1, send=200_000, slen=200_000, length=50_000)
    rows, _audit = ni.name_elements("S1", [row], [hit])
    assert rows[0]["start"] == "100000"
    assert rows[0]["end"] == "150000"
    assert rows[0]["element_type"] == "ice"


def test_a_much_larger_reference_is_reported_as_such():
    """The useful side effect of naming: it shows how far our boundaries fall
    short. On the K. pneumoniae positive control (ATCC BAA-2146, CP006659.2) our
    ICE call now spans 54,943 bp against ICEberg's 58,048 bp for the same element,
    ICEKpnATCCBAA-2146-1 — we recover 0.946 of it and stop 3,138 bp inside its far
    end. The 40 kb below is the wider gap that same element showed before the att
    search was reworked (2026-07-31), when the interval was only the machinery
    span; see docs/mobilome_worked_example.md."""
    row = ice_row(start=100_000, end=140_000)          # 40,001 bp
    hit = iceberg_hit(qstart=100_000, qend=140_000, length=40_000,
                      sstart=1, send=58_000, slen=58_000)
    _rows, audit = ni.name_elements("S1", [row], [hit])
    matched = [a for a in audit if a["reason"] == "iceberg_match"]
    assert matched
    assert "floor" in matched[0]["detail"]


# ── One real element matching many ICEberg entries ───────────────────────────

def test_near_identical_alternatives_are_counted():
    """ICEs of one species are near-identical across strains, so a single real
    element routinely matches many ICEberg entries at ~100%. Reporting one name
    without saying that would imply a precision the data does not have."""
    hits = [
        iceberg_hit(name="ICEKp1", bitscore="90000"),
        iceberg_hit(name="ICEKp2", bitscore="89950"),
        iceberg_hit(name="ICEKp3", bitscore="89900"),
    ]
    rows, audit = ni.name_elements("S1", [ice_row()], hits)
    assert rows[0]["mge_name"] == "ICEKp1"            # best bitscore wins
    assert int(rows[0]["iceberg_alternatives"]) == 2
    matched = [a for a in audit if a["reason"] == "iceberg_match"]
    assert "near-identical group" in matched[0]["detail"]


def test_several_hsps_of_one_reference_are_not_counted_as_alternatives():
    """Alternatives are counted per distinct element NAME, so one reference
    matching in several pieces does not look like several candidates."""
    hits = [
        iceberg_hit(name="ICEKp1", bitscore="90000", qstart=100_000, qend=120_000),
        iceberg_hit(name="ICEKp1", bitscore="89990", qstart=125_000, qend=150_000),
    ]
    rows, _audit = ni.name_elements("S1", [ice_row()], hits)
    assert int(rows[0]["iceberg_alternatives"]) == 0


# ── Degrading when ICEberg has nothing to say ────────────────────────────────
# Most environmental isolates carry no catalogued ICE, and the database is off
# altogether unless a URL is configured, so "no name" is the ordinary outcome.
# The element and its mobility tier must survive it untouched.

def test_no_blast_hits_leaves_everything_unnamed_but_intact():
    rows, audit = ni.name_elements("S1", [ice_row()], [])
    assert rows[0]["mge_name"] == "NA"
    assert rows[0]["element_type"] == "ice"
    assert "no_iceberg_hits" in reasons(audit)


def test_a_missing_blast_file_is_not_fatal():
    assert ni.read_blast_hits("/nonexistent.tsv") == []


def test_like_is_judged_on_the_part_inside_the_element():
    """The whole genome is blasted, so a curated element can match well beyond our
    interval. Judging "-like" on the full HSP answers "how much of the curated
    element is anywhere on this contig?" when the question is "how much of it is
    in the thing we are naming?" — and the first reading hands a bare, confident
    name to an element we have only partly found."""
    row = ice_row(start=100_000, end=110_000)          # our call: 10 kb
    # The curated element matches 50 kb of contig, only 10 kb of it inside us.
    hit = iceberg_hit(qstart=100_000, qend=150_000, length=50_000,
                      sstart=1, send=50_000, slen=50_000)
    rows, _audit = ni.name_elements("S1", [row], [hit])
    assert rows[0]["mge_name"] == "ICEKp1-like"
    assert float(rows[0]["iceberg_reference_coverage"]) < 0.3
