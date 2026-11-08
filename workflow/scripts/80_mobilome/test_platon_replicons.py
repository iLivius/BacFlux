"""Tests for platon_replicons.py — the per-contig "chromosome or plasmid?" call.

That single answer drives tiers 5 and 6 of the mobility ladder, and getting it
wrong in the CHROMOSOME direction is the worst error the module can make: every
AMR gene on a missed plasmid is then reported as a tier 1 "intrinsic candidate",
which is precisely the claim a reader would use to conclude the resistance is a
species trait and not transferable.

The tests concentrate on the four ways Platon and geNomad can combine, because
that is where that error was hiding, and on what happens when geNomad was never
run — it is opt-in, so that is the ordinary case, not an error.

No tool and no database is needed: the fixtures below build a directory shaped
like Platon's output in pytest's tmp_path.

Run: pytest workflow/scripts/80_mobilome/test_platon_replicons.py -q
"""

import os

import platon_replicons as pr


# ── Building a Platon output directory and a concordance table ───────────────

def write_platon_dir(tmp_path, prefix="contigs_final",
                     plasmids=(), chromosomes=(), table_rows=()):
    """Build a directory shaped like Platon's output.

    Platon writes three things platon_replicons.py reads: a per-contig TSV, a
    FASTA of the contigs it called chromosomal, and a FASTA of the ones it called
    plasmid. Any of them can legitimately be missing or header-only — a genome
    with no plasmid produces no plasmid FASTA at all — so the fixtures allow that.
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
    only when at least one of them thought it was a plasmid. A contig both tools
    called chromosomal is therefore ABSENT, which is why "not in this file" can
    never be read as "geNomad disagreed".
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


# ── The four ways Platon and geNomad can combine ─────────────────────────────
# Platon stays the default caller and geNomad is only ever consulted, never
# allowed to overrule an active Platon call. The `replicon_call_source` column
# records which of the two decided, so a reader can always see how much backing a
# plasmid call has.

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
    # Mobility is deliberately NOT guessed — the conjugation/mobilization gene
    # counts that type a plasmid come from Platon's own table, and Platon never
    # saw this contig. Whether it is tier 5 or tier 6 is left to CONJscan, which
    # scans the whole proteome independently.
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

    The call is NOT flipped — Platon is the default caller and made an active
    call — but the conflict has to be visible, because if geNomad is right then
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


# ── Degrading gracefully when geNomad never ran ──────────────────────────────

def test_a_missing_concordance_file_is_a_silent_no_op(tmp_path):
    """geNomad is opt-in, so its absence is normal, not an error."""
    assert pr.read_genomad_concordance("") == {}
    assert pr.read_genomad_concordance(str(tmp_path / "nope.tsv")) == {}


def test_a_header_only_concordance_means_no_plasmid_candidates(tmp_path):
    """Both tools ran and neither found a plasmid — a perfectly normal result for
    a genome that has none."""
    concordance = write_concordance(tmp_path, [])
    assert pr.read_genomad_concordance(concordance) == {}


# ── Platon's size filter: a long contig in the chromosome FASTA is not a call ──
# Platon skips every contig of 500 kb or more. When at least one other contig of the
# genome passes the filter, Platon still writes the skipped ones into its
# chromosome FASTA, because that file holds every input contig it did not call
# plasmid. Reading that file as a chromosome call would turn a 0.5-2 Mb
# megaplasmid into "chromosome", and every resistance gene on it into a tier 1
# intrinsic candidate. The contig lengths from the assembly (--contigs) are what
# lets the script tell the two apart.

MEGAPLASMID_BP = 1123075      # the size of a real megaplasmid (R. jostii pRHL1)
CHROMOSOME_BP = 5000000       # a typical bacterial chromosome
SMALL_PLASMID_BP = 80000      # short enough for Platon to classify


def test_a_megaplasmid_in_the_chromosome_fasta_stays_unknown(tmp_path):
    """The 1.1 Mb contig sits in Platon's chromosome FASTA only because Platon
    skipped it. Its length is in the band where a megaplasmid cannot be ruled
    out, so the honest answer is 'unknown', not 'chromosome'."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p_small"], chromosomes=["chr1", "mega"],
        table_rows=[{"ID": "p_small", "# Mobilization": "1"}])
    lengths = {"chr1": CHROMOSOME_BP, "mega": MEGAPLASMID_BP,
               "p_small": SMALL_PLASMID_BP}

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths))

    assert rows["mega"]["replicon"] == "unknown"
    assert "megaplasmid cannot be ruled out" in rows["mega"]["mobility_evidence"]
    # The small plasmid Platon did classify keeps Platon's call.
    assert rows["p_small"]["replicon"] == "plasmid"


def test_a_long_contig_in_the_chromosome_fasta_is_an_inferred_chromosome(tmp_path):
    """Over the 2 Mb ceiling the contig is called chromosome, but the evidence
    text says the call is inferred from length, not made by Platon."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p_small"], chromosomes=["chr1"])
    lengths = {"chr1": CHROMOSOME_BP, "p_small": SMALL_PLASMID_BP}

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths))

    assert rows["chr1"]["replicon"] == "chromosome"
    assert rows["chr1"]["mobility_evidence"].startswith("inferred chromosomal")


def test_genomad_decides_a_megaplasmid_platon_skipped(tmp_path):
    """Platon never looked at the 1.1 Mb contig, so a geNomad plasmid call is
    the only opinion there is: plasmid, source 'genomad'. It is not a conflict,
    because there is no Platon call to conflict with."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p_small"], chromosomes=["chr1", "mega"])
    concordance = write_concordance(tmp_path, [
        {"contig": "mega", "platon_call": "not_called",
         "genomad_call": "plasmid", "genomad_score": "0.97",
         "agreement": "genomad_only"},
    ])
    lengths = {"chr1": CHROMOSOME_BP, "mega": MEGAPLASMID_BP,
               "p_small": SMALL_PLASMID_BP}

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths,
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["mega"]["replicon"] == "plasmid"
    assert rows["mega"]["replicon_call_source"] == "genomad"
    assert "CONFLICT" not in rows["mega"]["mobility_evidence"]


def test_a_short_contig_in_the_chromosome_fasta_is_still_platons_call(tmp_path):
    """For a contig of 400 kb, under the 500 kb limit, the chromosome FASTA is
    read as Platon's call, with or without known lengths."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p_small"], chromosomes=["chr_fragment"])
    lengths = {"chr_fragment": 400000, "p_small": SMALL_PLASMID_BP}

    with_lengths = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths))
    without_lengths = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths)))

    for rows in (with_lengths, without_lengths):
        assert rows["chr_fragment"]["replicon"] == "chromosome"
        assert rows["chr_fragment"]["mobility_evidence"].startswith("chromosomal contig")


def test_a_contig_under_platons_minimum_is_not_a_chromosome_call(tmp_path):
    """Platon also skips every contig under 1,000 bp and writes it into the
    chromosome FASTA. Such a contig is not classified: unknown, unless geNomad
    calls it a plasmid."""
    platon_dir = write_platon_dir(
        tmp_path, plasmids=["p_small"], chromosomes=["chr1", "short", "short_plasmid"])
    concordance = write_concordance(tmp_path, [
        {"contig": "short_plasmid", "platon_call": "not_called",
         "genomad_call": "plasmid", "genomad_score": "0.95",
         "agreement": "genomad_only"},
    ])
    lengths = {"chr1": 400000, "short": 963, "short_plasmid": 980,
               "p_small": SMALL_PLASMID_BP}

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths,
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["short"]["replicon"] == "unknown"
    assert "1000 bp minimum" in rows["short"]["mobility_evidence"]
    assert rows["short_plasmid"]["replicon"] == "plasmid"
    assert rows["short_plasmid"]["replicon_call_source"] == "genomad"
    assert rows["chr1"]["replicon"] == "chromosome"


def test_a_genomad_plasmid_call_on_a_chromid_sized_contig_is_a_conflict(tmp_path):
    """A contig of 2 Mb or more that only geNomad calls a plasmid is most likely
    a secondary chromosome (a chromid): chromosome, flagged conflict, so its
    genes are neither called mobile nor intrinsic at full confidence. A 1.5 Mb
    megaplasmid with the same geNomad call stays a plasmid."""
    platon_dir = write_platon_dir(tmp_path)          # closed genome: no Platon files
    concordance = write_concordance(tmp_path, [
        {"contig": "chromid", "platon_call": "not_called", "genomad_call": "plasmid",
         "genomad_score": "0.94", "agreement": "genomad_only"},
        {"contig": "megaplasmid", "platon_call": "not_called", "genomad_call": "plasmid",
         "genomad_score": "0.97", "agreement": "genomad_only"},
    ])
    lengths = {"chr1": 4895836, "chromid": 3363523, "megaplasmid": 1471779}

    rows = rows_by_contig(pr.build_rows(
        platon_dir, "contigs_final", sorted(lengths), contig_lengths=lengths,
        genomad_calls=pr.read_genomad_concordance(concordance)))

    assert rows["chromid"]["replicon"] == "chromosome"
    assert rows["chromid"]["replicon_call_source"] == "conflict"
    assert "chromid" in rows["chromid"]["mobility_evidence"]
    assert rows["megaplasmid"]["replicon"] == "plasmid"
    assert rows["megaplasmid"]["replicon_call_source"] == "genomad"
    assert rows["chr1"]["mobility_evidence"].startswith("inferred chromosomal")
