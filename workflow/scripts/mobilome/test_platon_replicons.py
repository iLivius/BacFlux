"""Tests for platon_replicons.py — the per-contig replicon call.

WHY THIS FILE MATTERS
    This script decides, for every contig, "chromosome or plasmid?", and that
    single answer drives tiers 5 and 6 of the mobility ladder. Getting it wrong
    in the CHROMOSOME direction is the worst error the module can make: every AMR
    gene on a missed plasmid is then reported as a tier 1 "intrinsic candidate",
    which is precisely the claim a reader would use to conclude the resistance is
    a species trait and not transferable.

    The tests below concentrate on the four ways Platon and geNomad can combine,
    because that is where that error was hiding.
"""

import os

import platon_replicons as pr


# ── Fixtures ─────────────────────────────────────────────────────────────────

def write_platon_dir(tmp_path, prefix="contigs_final",
                     plasmids=(), chromosomes=(), table_rows=()):
    """Build a directory shaped like Platon's output.

    Platon writes three things this script reads: a per-contig TSV, a FASTA of
    the contigs it called chromosomal, and a FASTA of the ones it called plasmid.
    Any of them can legitimately be missing or header-only — a genome with no
    plasmid produces no plasmid FASTA at all — so the fixtures allow that.
    """
    directory = tmp_path / "platon"
    directory.mkdir(exist_ok=True)

    if chromosomes:
        (directory / f"{prefix}.chromosome.fasta").write_text(
            "".join(f">{name}\nACGT\n" for name in chromosomes))
    if plasmids:
        (directory / f"{prefix}.plasmid.fasta").write_text(
            "".join(f">{name}\nACGT\n" for name in plasmids))

    header = ["ID", "Length", "# Conjugation", "# Mobilization", "# OriT", "Inc Type(s)"]
    lines = ["\t".join(header)]
    for row in table_rows:
        lines.append("\t".join(str(row.get(column, "0")) for column in header))
    (directory / f"{prefix}.tsv").write_text("\n".join(lines) + "\n")
    return str(directory)


def write_concordance(tmp_path, rows):
    """Build a {sample}_plasmid_concordance.tsv as rule plasmid_concordance does.

    Its rows are the UNION of the two tools' plasmid calls, so a contig appears
    only when at least one of them thought it was a plasmid.
    """
    header = ["sample", "contig", "platon_call", "genomad_call",
              "genomad_score", "agreement", "confidence"]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(str(row.get(column, "NA")) for column in header))
    path = tmp_path / "concordance.tsv"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def rows_by_contig(rows):
    return {row["contig"]: row for row in rows}


# ── The four Platon x geNomad outcomes ───────────────────────────────────────

def test_both_tools_agree_the_contig_is_a_plasmid(tmp_path):
    """The easy case: Platon called it, geNomad agrees. Nothing changes except
    that the row now records that two tools backed the call."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p1"], chromosomes=["chr1"],
        table_rows=[{"ID": "p1", "# Mobilization": "1"}])
    concordance = write_concordance(tmp_path, [
        {"contig": "p1", "platon_call": "plasmid", "genomad_call": "plasmid",
         "agreement": "both"},
    ])

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", ["p1", "chr1"],
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["p1"]["replicon"] == "plasmid"
    assert rows["p1"]["replicon_call_source"] == "both"
    assert rows["p1"]["plasmid_mobility"] == "mobilisable"


def test_a_plasmid_platon_never_classified_is_rescued_by_genomad(tmp_path):
    """THE case this wiring exists for.

    Platon does not mention the contig at all. Before geNomad was consulted it
    fell through to 'unknown', and every AMR gene on it was reported as
    chromosomal and intrinsic. geNomad is the only opinion available, so using it
    is not overriding anything.
    """
    platon_dir = write_platon_dir(tmp_path, chromosomes=["chr1"])
    concordance = write_concordance(tmp_path, [
        {"contig": "p_missed", "platon_call": "not_called",
         "genomad_call": "plasmid", "genomad_score": "0.98",
         "agreement": "genomad_only"},
    ])

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", ["chr1", "p_missed"],
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["p_missed"]["replicon"] == "plasmid"
    assert rows["p_missed"]["replicon_call_source"] == "genomad"
    # Mobility is deliberately NOT guessed - Platon's conjugation counts are what
    # type it and they do not exist here. CONJscan resolves tier 5 vs 6 instead.
    assert rows["p_missed"]["plasmid_mobility"] == "unknown"
    assert "geNomad" in rows["p_missed"]["mobility_evidence"]


def test_without_genomad_the_same_contig_stays_unknown(tmp_path):
    """The counterpart, proving the rescue really is geNomad's doing and that the
    script is unchanged when geNomad was not run."""
    platon_dir = write_platon_dir(tmp_path, chromosomes=["chr1"])

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", ["chr1", "p_missed"]))

    assert rows["p_missed"]["replicon"] == "unknown"
    assert rows["p_missed"]["replicon_call_source"] == "platon"


def test_a_straight_disagreement_is_flagged_not_silently_resolved(tmp_path):
    """Platon says chromosome, geNomad says plasmid.

    The call is NOT flipped - Platon is the default caller and made an active
    call - but the conflict has to be visible, because if geNomad is right then
    calling the gene intrinsic is exactly wrong.
    """
    platon_dir = write_platon_dir(tmp_path, chromosomes=["c_disputed"])
    concordance = write_concordance(tmp_path, [
        {"contig": "c_disputed", "platon_call": "chromosome",
         "genomad_call": "plasmid", "genomad_score": "0.91",
         "agreement": "conflict"},
    ])

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", ["c_disputed"],
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["c_disputed"]["replicon"] == "chromosome"
    assert rows["c_disputed"]["replicon_call_source"] == "conflict"
    assert "CONFLICT" in rows["c_disputed"]["mobility_evidence"]


def test_platon_only_plasmid_keeps_its_own_call(tmp_path):
    """geNomad ran but did not call this contig a plasmid. Platon's positive call
    stands; the row simply records that only one tool backed it."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p1"], chromosomes=["chr1"],
        table_rows=[{"ID": "p1", "# Conjugation": "4"}])
    concordance = write_concordance(tmp_path, [
        {"contig": "p1", "platon_call": "plasmid", "genomad_call": "absent",
         "agreement": "platon_only"},
    ])

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", ["p1", "chr1"],
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["p1"]["replicon"] == "plasmid"
    assert rows["p1"]["replicon_call_source"] == "platon"
    assert rows["p1"]["plasmid_mobility"] == "conjugative"


# ── Graceful degradation ─────────────────────────────────────────────────────

def test_a_missing_concordance_file_is_a_silent_no_op(tmp_path):
    """geNomad is opt-in, so its absence is normal, not an error."""
    assert pr.read_genomad_concordance("") == {}
    assert pr.read_genomad_concordance(str(tmp_path / "nope.tsv")) == {}


def test_a_header_only_concordance_means_no_plasmid_candidates(tmp_path):
    """Both tools ran and neither found a plasmid - a perfectly normal result for
    a genome that has none."""
    concordance = write_concordance(tmp_path, [])
    assert pr.read_genomad_concordance(concordance) == {}
