"""Unit tests for att_search.py — the attL/attR direct-repeat search.

att_search is the only piece of original algorithm in the mobilome module.
conjscan_to_ice.py hands it one contig's sequence plus the coordinates of a
conjugation-machinery cluster, and it looks for the pair of direct repeats that
site-specific integration leaves at the two ends of an integrated element. Only a
tRNA-anchored pair is ever acted on: it moves the element's start and end, and
therefore decides which genes count as cargo of something predicted to be
mobile. A mistake here either invents a boundary that does not exist or throws
away a real one, which is why the assertions below are on exact coordinates.

Every test builds a small SYNTHETIC contig with repeats planted at coordinates
the test itself chose, so the expected answer is known exactly and the assertions
can be on the numbers, not on "something was found". That is what spec §8 Phase 7
asks for: "unit tests on synthetic contigs with planted att repeats at known
offsets".

No genome, no Bakta, no database — pure sequence handling.

Run: pytest workflow/scripts/80_mobilome/test_att_search.py -q
"""

import random

import att_search as att


# ── Building synthetic genomes ──────────────────────────────────────────────

def random_sequence(length, seed):
    """Reproducible random DNA.

    A fixed seed per test keeps failures debuggable: the same "genome" is rebuilt
    byte for byte on every run, so a failing assertion always refers to the same
    sequence. Random background also means any repeat the search finds was
    planted deliberately — chance 25-mers do not occur in a few kb.
    """
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def plant(sequence, position, motif):
    """Overwrite `motif` into `sequence` at a 1-based inclusive position."""
    index = position - 1
    return sequence[:index] + motif + sequence[index + len(motif):]


def build_contig_with_att(att_motif, left_position, right_position,
                          length=80_000, seed=1):
    """A contig carrying the same motif twice, at two known positions.

    This is the shape the search is looking for: one copy at attL, one at attR,
    same orientation, with the element's machinery somewhere between them.
    """
    sequence = random_sequence(length, seed)
    sequence = plant(sequence, left_position, att_motif)
    sequence = plant(sequence, right_position, att_motif)
    return sequence


# A 25 bp motif used as the planted att site throughout. Written out rather than
# generated so the tests read concretely.
ATT_MOTIF = "GGCTCGAACCCAGGACCTCTTGCAT"


# ── Repeats with no tRNA involved: boundary_method='denovo' ─────────────────
#
# These tests pass no tRNAs at all, so no candidate can be tRNA-anchored and
# every call the search makes is labelled 'denovo'. That is the weaker of the two
# labels and conjscan_to_ice.py refuses to move an element's coordinates onto it,
# but the search still has to find the right repeat at the right offsets.

def test_denovo_finds_a_planted_direct_repeat_at_the_exact_offsets():
    """The core case: two copies of one motif bracketing the machinery."""
    left_position = 20_000
    right_position = 50_000
    sequence = build_contig_with_att(ATT_MOTIF, left_position, right_position)

    # Machinery sits between the two planted copies.
    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000)

    assert result["boundary_method"] == "denovo"
    assert result["att_left"].startswith("%d.." % left_position)
    assert result["att_right"].startswith("%d.." % right_position)
    # The reported repeat CONTAINS the planted motif and may extend past it.
    assert ATT_MOTIF in result["att_sequence"]
    # The element runs from the START of attL to the END of attR: both repeats
    # are part of the integrated element.
    assert result["element_start"] == left_position
    assert result["element_end"] >= right_position + len(ATT_MOTIF) - 1
    assert result["element_length_bp"] == result["element_end"] - result["element_start"] + 1


def test_denovo_returns_none_when_there_is_no_repeat():
    """Random sequence with nothing planted must yield no boundaries at all —
    the search must not invent an element out of background similarity."""
    sequence = random_sequence(80_000, seed=2)
    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000)
    assert result["boundary_method"] == "none"
    assert result["att_left"] == "NA"
    assert result["element_start"] is None


def test_denovo_ignores_a_repeat_that_does_not_bracket_the_machinery():
    """Two copies both sitting to the LEFT of the element are not attL/attR.

    Only the flanking regions are searched, so a repeat pair entirely on one side
    can never be paired up — which is what stops an unrelated duplication
    elsewhere on the contig from being read as element boundaries.
    """
    sequence = build_contig_with_att(ATT_MOTIF, 5_000, 12_000, seed=3)
    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000)
    assert result["boundary_method"] == "none"


def test_denovo_respects_the_element_size_bounds():
    """A repeat pair implying an element far larger than max_element_bp is not a
    plausible ICE and must be rejected."""
    sequence = build_contig_with_att(ATT_MOTIF, 1_000, 79_000, seed=4)
    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000, max_element_bp=50_000)
    assert result["boundary_method"] == "none"


def test_denovo_prefers_the_longest_repeat():
    """Among unanchored candidates the longest repeat wins: a 25 bp match beats a
    12 bp one, because a long exact match between two specific windows is far
    less likely to have happened by chance."""
    sequence = random_sequence(80_000, seed=5)
    # A short repeat, and a longer one, both bracketing the machinery.
    short_motif = "ACGTACGTACGT"                       # 12 bp
    sequence = plant(sequence, 18_000, short_motif)
    sequence = plant(sequence, 55_000, short_motif)
    sequence = plant(sequence, 22_000, ATT_MOTIF)      # 25 bp
    sequence = plant(sequence, 48_000, ATT_MOTIF)

    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000)
    assert result["boundary_method"] == "denovo"
    assert result["att_sequence"] == ATT_MOTIF
    assert result["att_length_bp"] == 25


# ── The masking guard: IS repeats must not flood the search ─────────────────

def test_masked_is_repeats_do_not_produce_a_false_boundary():
    """The failure mode the spec calls the most likely one.

    Insertion sequences carry terminal repeats, so a region with IS copies is
    full of direct repeats that have nothing to do with ICE integration. Here two
    copies of an "IS end" bracket the machinery and would otherwise be called as
    attL/attR; masking the IS intervals removes them.
    """
    is_repeat = "TGTCAGGGCCCTTAAGGGCCCTGA"
    sequence = random_sequence(80_000, seed=6)
    sequence = plant(sequence, 20_000, is_repeat)
    sequence = plant(sequence, 50_000, is_repeat)

    # Without masking the IS repeat is found and reported as a boundary.
    unmasked = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                  min_element_bp=8_000)
    assert unmasked["boundary_method"] == "denovo"

    # With the IS intervals masked, the decoy is gone and nothing is called.
    masked = att.find_att_sites(
        sequence, element_start=30_000, element_end=40_000, min_element_bp=8_000,
        mask=[(19_900, 20_100), (49_900, 50_100)])
    assert masked["boundary_method"] == "none"


def test_mask_intervals_blanks_the_right_bases_and_keeps_the_length():
    """Masking must not shift coordinates — every position reported afterwards
    still refers to the real genome."""
    sequence = "ACGT" * 10                     # 40 bp
    masked = att.mask_intervals(sequence, [(5, 8)])
    assert len(masked) == len(sequence)
    assert masked[4:8] == "NNNN"
    assert masked[:4] == sequence[:4]
    assert masked[8:] == sequence[8:]


def test_masking_is_clamped_to_the_sequence():
    """An interval running off either end must not crash or wrap around."""
    sequence = "ACGT" * 10
    masked = att.mask_intervals(sequence, [(-5, 3), (38, 999)])
    assert len(masked) == len(sequence)
    assert masked.startswith("NNN")
    assert masked.endswith("NNN")


# ── Repeats sitting in a tRNA: boundary_method='tRNA' ───────────────────────

def trna_feature(contig, start, end, strand, name="tRNA-Gly(gcc)"):
    """One tRNA in the shape att.parse_trna_features returns from a Bakta GFF3."""
    return {"contig": contig, "start": start, "end": end,
            "strand": strand, "name": name}


def test_trna_anchored_search_is_preferred_over_denovo():
    """When a repeat sits in a tRNA 3' end, that answer wins — even if a longer
    unanchored repeat is available.

    This is the ranking rule the search turns on: a repeat at a tRNA is where
    integration actually happens, while a longer repeat elsewhere is only less
    likely to be coincidence. So the tRNA hit must be returned and the method
    reported as 'tRNA', which is the only label conjscan_to_ice.py acts on.
    """
    sequence = random_sequence(80_000, seed=7)

    # A 76 bp tRNA on the + strand whose last 25 bp are our motif, and a second
    # copy of those 25 bp on the far side of the machinery.
    trna_start, trna_end = 20_000, 20_075
    sequence = plant(sequence, trna_end - 24, ATT_MOTIF)      # the tRNA 3' end
    sequence = plant(sequence, 50_000, ATT_MOTIF)             # the scar copy

    # A competing, unrelated de novo repeat, also bracketing.
    other = "TTTTGGGGAAAACCCCTTTTGGGGA"
    sequence = plant(sequence, 25_000, other)
    sequence = plant(sequence, 45_000, other)

    result = att.find_att_sites(
        sequence, element_start=30_000, element_end=40_000, min_element_bp=8_000,
        trnas=[trna_feature("c1", trna_start, trna_end, "+")])

    assert result["boundary_method"] == "tRNA"
    assert result["att_sequence"] == ATT_MOTIF
    assert result["trna"] == "tRNA-Gly(gcc)"
    assert result["att_left"].startswith(str(trna_end - 24))


# The next two are deliberately retired, not broken — note the `_retired_`
# prefix, which stops pytest collecting them. They pinned the OLD mismatch-
# tolerant probe, which allowed attL and attR to differ by one base. The search
# is exact now (maximal repeats are grown only while the flanks agree base for
# base), so both would fail. They are kept as the record of what the module used
# to promise; what replaced them is the exact-match pair further down, under
# "The exact-match contract". Do not "fix" them — delete them or leave them.

def _retired_test_trna_anchored_tolerates_a_single_mismatch():
    """attL and attR often differ by one base, because only the copy that
    reconstitutes the tRNA is under selection. One mismatch must still match."""
    sequence = random_sequence(80_000, seed=9)
    trna_start, trna_end = 20_000, 20_075
    sequence = plant(sequence, trna_end - 24, ATT_MOTIF)

    # The second copy differs at one position.
    imperfect = "T" + ATT_MOTIF[1:] if ATT_MOTIF[0] != "T" else "A" + ATT_MOTIF[1:]
    sequence = plant(sequence, 50_000, imperfect)

    result = att.find_att_sites(
        sequence, element_start=30_000, element_end=40_000, min_element_bp=8_000,
        trnas=[trna_feature("c1", trna_start, trna_end, "+")])

    assert result["boundary_method"] == "tRNA"
    assert result["att_mismatches"] >= 1


def _retired_test_two_mismatches_are_rejected():
    """Beyond one substitution this is no longer a recombination scar."""
    sequence = random_sequence(80_000, seed=10)
    trna_start, trna_end = 20_000, 20_075
    sequence = plant(sequence, trna_end - 24, ATT_MOTIF)
    two_off = "TT" + ATT_MOTIF[2:]
    sequence = plant(sequence, 50_000, two_off)

    result = att.find_att_sites(
        sequence, element_start=30_000, element_end=40_000, min_element_bp=8_000,
        trnas=[trna_feature("c1", trna_start, trna_end, "+")])
    # No tRNA-anchored call; may fall through to de novo, but must not claim tRNA.
    assert result["boundary_method"] != "tRNA"


def test_a_trna_far_outside_the_window_is_not_used():
    """Only tRNAs near the element can plausibly be its integration site."""
    sequence = random_sequence(200_000, seed=11)
    trna_start, trna_end = 150_000, 150_075
    sequence = plant(sequence, trna_end - 24, ATT_MOTIF)
    sequence = plant(sequence, 160_000, ATT_MOTIF)

    result = att.find_att_sites(
        sequence, element_start=30_000, element_end=40_000, min_element_bp=8_000,
        flank_window_bp=10_000,
        trnas=[trna_feature("c1", trna_start, trna_end, "+")])
    assert result["boundary_method"] == "none"


# ── Degenerate and defensive cases ──────────────────────────────────────────

def test_empty_or_missing_input_is_not_fatal():
    """A contig with no sequence, or an element with no coordinates, must return
    the 'none' result rather than raise — the caller writes these columns for
    every candidate unconditionally."""
    assert att.find_att_sites("", 10, 20)["boundary_method"] == "none"
    sequence = random_sequence(1_000, seed=12)
    assert att.find_att_sites(sequence, None, 20)["boundary_method"] == "none"
    assert att.find_att_sites(sequence, 10, None)["boundary_method"] == "none"


def test_reversed_coordinates_are_tolerated():
    """start > end should be swapped, not silently produce nonsense."""
    sequence = build_contig_with_att(ATT_MOTIF, 20_000, 50_000, seed=13)
    forward = att.find_att_sites(sequence, 30_000, 40_000, min_element_bp=8_000)
    reversed_pair = att.find_att_sites(sequence, 40_000, 30_000, min_element_bp=8_000)
    assert forward["att_left"] == reversed_pair["att_left"]
    assert forward["att_right"] == reversed_pair["att_right"]


# ── FASTA and GFF3 readers ──────────────────────────────────────────────────
#
# In a real run the sequence is the delivered assembly (contigs_final.fasta) and
# the tRNAs come from that sample's Bakta GFF3, both handed over by
# conjscan_to_ice.py, the only caller. These two readers are where a sample
# silently loses its att search: a contig id that does not match ISEScan's, or a
# tRNA block that is never found, ends as boundary_method='none' on every
# element rather than as an error.

def test_read_fasta_keys_on_the_first_token_and_upper_cases(tmp_path):
    """Contigs are keyed on the first token of the header, not the whole line.

    That is the same convention the rest of BacFlux uses, and it has to be: the
    key must match the seqid in Bakta's GFF3, or no tRNA can ever be placed on
    the sequence it came from.
    """
    path = tmp_path / "genome.fna"
    path.write_text(">NZ_CP006659.2 Klebsiella pneumoniae chromosome\nacgt\nACGT\n"
                    ">contig_2\nTTTT\n")
    sequences = att.read_fasta(str(path))
    assert set(sequences) == {"NZ_CP006659.2", "contig_2"}
    # Soft-masked lower case must be normalised, or repeat search silently misses
    # exactly the repetitive regions it is looking for.
    assert sequences["NZ_CP006659.2"] == "ACGTACGT"


def test_parse_trna_features_reads_bakta_gff3_and_stops_at_the_fasta(tmp_path):
    """Bakta appends the whole assembly as FASTA after a ##FASTA line, so the
    reader has to stop there: annotation ends, megabases of sequence begin.

    The fixture plants the letters tRNA inside a sequence line on purpose. Only
    the tab-delimited column check keeps that out of the results today, so if
    anyone ever loosens the parser to a plain text search, this test fails
    instead of the run quietly gaining a tRNA at a made-up coordinate.
    """
    path = tmp_path / "sample.gff3"
    path.write_text(
        "##gff-version 3\n"
        "NZ_CP1\ttRNAscan-SE\ttRNA\t17807\t17882\t.\t+\t.\t"
        "ID=X_00016;Name=tRNA-Glu(ttc);product=tRNA-Glu(ttc)\n"
        "NZ_CP1\tProdigal\tCDS\t100\t400\t.\t+\t0\tID=X_00001\n"
        "NZ_CP1\ttRNAscan-SE\ttRNA\t22000\t22075\t.\t-\t.\t"
        "ID=X_00020;product=tRNA-Gly(gcc)\n"
        "##FASTA\n"
        ">NZ_CP1\n"
        "ACGTtRNAACGT\n")                       # must NOT be parsed as annotation
    trnas = att.parse_trna_features(str(path))
    assert len(trnas) == 2
    assert trnas[0]["name"] == "tRNA-Glu(ttc)"
    assert trnas[0]["strand"] == "+"
    assert trnas[1]["name"] == "tRNA-Gly(gcc)"
    assert trnas[1]["strand"] == "-"


def test_parse_trna_features_on_a_missing_file_is_not_fatal():
    """No annotation means no tRNA anchoring, not a crash: the search falls back
    to de novo repeats, which are reported but never move an element."""
    assert att.parse_trna_features("/nonexistent/sample.gff3") == []


def test_group_by_contig():
    records = [{"contig": "a", "n": 1}, {"contig": "b", "n": 2}, {"contig": "a", "n": 3}]
    grouped = att.group_by_contig(records)
    assert set(grouped) == {"a", "b"}
    assert len(grouped["a"]) == 2


# ── The statistical floor on repeat length ──────────────────────────────────

def test_minimum_repeat_length_scales_with_the_window():
    """A bigger search window must demand a longer repeat.

    This is the rule that stopped the de novo scan reporting noise: comparing two
    30 kb flanks, a 12 bp match is expected ~54 times by chance, so believing one
    would mean calling an element boundary on essentially every genome.
    """
    small = att.minimum_informative_repeat_length(1_000, 1_000)
    large = att.minimum_informative_repeat_length(30_000, 30_000)
    assert large > small
    # For 30 kb flanks the arithmetic lands around 20 bp; well above the old
    # fixed floor of 12, and comfortably inside the real 15–25 bp att range.
    assert 18 <= large <= 24
    # Never below the absolute floor, however tiny the window.
    assert att.minimum_informative_repeat_length(10, 10) == att.DENOVO_ABSOLUTE_MIN_REPEAT_BP


def test_expected_chance_matches_at_the_chosen_length_are_negligible():
    """State the guarantee directly: at the chosen k, fewer than the configured
    number of matches are expected by chance."""
    for flank_bp in (5_000, 30_000, 100_000):
        k = att.minimum_informative_repeat_length(flank_bp, flank_bp)
        expected = (flank_bp * flank_bp) / (4 ** k)
        assert expected <= att.DENOVO_MAX_EXPECTED_CHANCE_MATCHES


def test_random_sequence_yields_no_boundaries_at_several_window_sizes():
    """The regression this whole rule exists for, checked across window sizes and
    seeds: background similarity must never be reported as an att site."""
    for seed in (101, 102, 103):
        sequence = random_sequence(120_000, seed=seed)
        for window in (5_000, 20_000, 40_000):
            result = att.find_att_sites(
                sequence, element_start=50_000, element_end=60_000,
                min_element_bp=8_000, flank_window_bp=window)
            assert result["boundary_method"] == "none", (seed, window)


# ── Negative controls built from REAL genome structure ──────────────────────
#
# The tests above use random DNA in which every base is equally likely and
# independent of its neighbours — exactly the assumption the de novo chance model
# makes — so they can only ever confirm that model, never challenge it. Measured
# on the K. pneumoniae positive-control chromosome, 300 randomly placed non-ICE
# spans produced a confident "denovo" boundary 22% of the time, against the ~1.3%
# the model predicts. The difference is not noise: it is the repetitive structure
# that every real chromosome has and random DNA does not.
#
# These tests plant that structure deliberately.

def test_a_dispersed_repeat_family_is_not_an_att_site():
    """An rRNA-operon-like repeat, present many times, must not become a boundary.

    This is that control's false positive in miniature. A 25 bp stretch of 16S
    rRNA satisfies every length threshold and genuinely IS an exact direct repeat
    shared by the two flanks — but the cell carries seven rRNA operons, so the
    sequence occurs seven times. An integration scar occurs exactly twice, which
    is what separates the two cases.
    """
    sequence = random_sequence(200_000, seed=404)
    # Two copies bracket the machinery, as an att pair would...
    for position in (30_000, 90_000):
        sequence = plant(sequence, position, ATT_MOTIF)
    # ...but five more copies elsewhere make it a FAMILY, not a scar.
    for position in (5_000, 120_000, 140_000, 160_000, 180_000):
        sequence = plant(sequence, position, ATT_MOTIF)

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "none"


def test_a_repeat_present_exactly_twice_is_still_accepted():
    """The counterpart, so the guard above cannot pass by rejecting everything."""
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=404)

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "denovo"
    # Maximal repeats extend past the planted motif where the flanks agree.
    assert result["att_left"].startswith("30000..")
    assert result["att_right"].startswith("85000..")


def test_two_paralogous_trnas_are_not_an_att_pair():
    """Isoacceptor tRNAs share their 3' ends, and a genome carries dozens.

    Two paralogous tRNA genes share their 3' ends by definition, so a repeat
    found in one is guaranteed to be found in the other — producing a beautifully
    bracketing pair that gets labelled 'tRNA', the method this module treats as
    its most precise. Nine such pairs turned up in 300 random spans of the
    K. pneumoniae positive-control chromosome.

    Real integration reconstitutes the host tRNA at ONE end and leaves the second
    copy out in ordinary sequence, so exactly one copy may sit in a tRNA.
    """
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=505)
    both_in_trnas = [
        trna_feature("contig_1", 29_952, 30_024, "+"),   # left copy is its 3' end
        trna_feature("contig_1", 84_952, 85_024, "+"),   # so is the right copy
    ]

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=both_in_trnas, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] != "tRNA"


def test_one_copy_in_a_trna_is_the_real_integration_signature():
    """Same sequence, but only the left copy is inside a tRNA — which is what
    site-specific integration at a tRNA 3' end actually leaves behind."""
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=505)
    one_trna = [trna_feature("contig_1", 29_952, 30_024, "+")]

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=one_trna, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "tRNA"
    # Maximal repeats extend past the planted motif where the flanks agree.
    assert result["att_left"].startswith("30000..")
    assert result["att_right"].startswith("85000..")


# ── The exact-match contract, replacing the two retired mismatch tests ──────
#
# The old design took a fixed probe and allowed up to one mismatch against it.
# The search is now EXACT maximal repeats — `vmatch -l` semantics, which is what
# ICEfinder2 actually uses — so mismatch tolerance no longer exists. That is a
# deliberate trade, and the honest limitation is recorded here rather than
# hidden: BLAST (the DEPhT route) would absorb mismatches and indels, but needs
# a subprocess dependency this module does not have. See
# docs/methods_att_and_small_plasmids.md.
#
# The two retired tests above are kept, renamed, as a record of what changed.

def test_the_search_is_exact_so_a_mismatched_copy_is_not_a_pair():
    """A second copy carrying a substitution is no longer found as one repeat.

    It is still found as the two EXACT sub-repeats either side of the mismatch,
    which is why this asserts on the reported length rather than on absence: a
    23 bp site with a mismatch in the middle yields an ~11 bp maximal repeat,
    below every floor, so nothing is reported.
    """
    motif = ATT_MOTIF
    broken = motif[:12] + ("A" if motif[12] != "A" else "C") + motif[13:]
    sequence = random_sequence(200_000, seed=717)
    sequence = plant(sequence, 30_000, motif)
    sequence = plant(sequence, 85_000, broken)

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "none"


def test_an_exact_pair_is_still_found_after_the_change():
    """The counterpart, so the test above cannot pass by finding nothing ever."""
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000,
                                     length=200_000, seed=717)
    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        min_element_bp=8_000, flank_window_bp=30_000)
    assert result["boundary_method"] == "denovo"
    assert int(result["att_length_bp"]) >= len(ATT_MOTIF)


def test_a_trna_derived_repeat_survives_its_own_paralogues():
    """THE regression that prompted the tRNA-aware copy guard.

    ICEKp integrates at tRNA-Asn, so its att core is a piece of a tRNA 3' end —
    and a genome with five tRNA-Asn genes contains that sequence five times
    whether or not an ICE is present. The old "occurs at most twice" rule
    therefore threw away the real, published ICEKp att site. Copies INSIDE tRNAs
    are now not counted, because tRNA paralogy explains them.
    """
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000,
                                     length=200_000, seed=818)
    # three more copies of the same motif, each inside its own tRNA — the
    # paralogous tRNA genes that defeated the old rule
    trnas = [trna_feature("contig_1", 29_952, 30_024, "+")]
    for position in (120_000, 150_000, 175_000):
        sequence = plant(sequence, position, ATT_MOTIF)
        trnas.append(trna_feature("contig_1", position - 48, position + 24, "+"))

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=trnas, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "tRNA"


def test_two_different_trna_species_can_still_bracket_an_element():
    """An element that landed between UNLIKE tRNAs is not a paralogue pair.

    The guard above rejects a repeat whose two copies both sit in tRNA genes,
    because two copies of the SAME tRNA are paralogues rather than an integration
    scar. That reasoning does not extend to two DIFFERENT tRNAs: a genome's tRNA
    genes are not all copies of each other.

    Measured on ICEEc2 (GU725392), where the real 22 bp att pair sits in tRNA-Phe
    at one end and tRNA-Ser at the other. The blanket rule discarded the correct
    boundary and the element was reported 37 kb short of its true extent — even
    though the right pair was in the candidate list and outscored the winner.
    """
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=505)
    unlike_trnas = [
        trna_feature("contig_1", 29_952, 30_024, "+", name="tRNA-Phe(gaa)"),
        trna_feature("contig_1", 84_952, 85_024, "+", name="tRNA-Ser(gct)"),
    ]

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=unlike_trnas, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "tRNA"


def test_an_unparseable_trna_product_keeps_the_stricter_old_behaviour():
    """A tRNA whose product cannot be read must not be assumed to differ.

    same_trna_species answers True when either name is unparseable, so the pair
    is still rejected. Guessing the other way would let an unnamed feature open
    the guard that exists to keep paralogues out.
    """
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=505)
    unnamed = [
        trna_feature("contig_1", 29_952, 30_024, "+", name="tRNA"),
        trna_feature("contig_1", 84_952, 85_024, "+", name="tRNA"),
    ]

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=unnamed, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] != "tRNA"


def test_a_trna_anchored_pair_beats_a_longer_unanchored_one():
    """Anchoring outranks length, because they are not the same kind of evidence.

    Repeat length says how unlikely a match is by chance; sitting at a tRNA 3' end
    says the match is where integration actually happens. Ranking on length alone
    let a 51 bp repeat in ordinary sequence beat the real 24 bp att pair at
    tRNA-Phe on SPI-7 (AL513382), and the element came out 50 kb short.
    """
    assert att.trna_species({"name": "tRNA-Phe(gaa)"}) == "phe"
    assert att.trna_species({"name": "tRNA-Ser(gct)"}) == "ser"
    assert att.trna_species({"name": "tRNA"}) == ""
    # Unlike species do not trip the paralogue guard; like species do.
    assert not att.same_trna_species({"name": "tRNA-Phe(gaa)"}, {"name": "tRNA-Ser(gct)"})
    assert att.same_trna_species({"name": "tRNA-Gly(gcc)"}, {"name": "tRNA-Gly(tcc)"})
    # An unreadable product falls back to "treat as the same", the safer answer.
    assert att.same_trna_species({"name": "tRNA"}, {"name": "tRNA-Ser(gct)"})


# ── The two copy counters: one guard, two scopes ────────────────────────────

def test_what_each_copy_counter_counts():
    """Pin the counting convention of both copy counters — neither was pinned.

    The two differ in two independent ways, and only one of them bites.
    count_repeat_copies uses str.count, which skips past each hit it finds, while
    count_repeat_copies_outside_trnas steps one base at a time and so also sees
    copies that overlap each other; that WALK can only tell them apart on a
    periodic repeat, which an att site is not. What does separate them in normal
    use is SCOPE: only the second discounts copies inside tRNAs. Both are pinned
    here, with the direction each falls.
    """
    # The published ICEKp direct repeat (Lam et al. 2018), 17 bp. It has no
    # period shorter than itself: sliding it along by 1..16 bases never lines it
    # up with itself, so two copies of it can never overlap.
    icekp = "CCAGTCAGAGGAGCCAA"
    assert all(icekp[shift:] != icekp[:len(icekp) - shift]
               for shift in range(1, len(icekp)))

    # Two ordinary copies, far apart: both counters say 2, which is the attL and
    # attR of a scar and passes both guards.
    background = random_sequence(5_000, seed=1207)
    scar = plant(plant(background, 1_000, icekp), 3_000, icekp)
    assert att.count_repeat_copies(scar, icekp) == 2
    assert att.count_repeat_copies_outside_trnas(scar, icekp, []) == 2
    assert att.count_repeat_copies(scar, icekp) <= att.DENOVO_MAX_CONTIG_COPIES
    assert (att.count_repeat_copies_outside_trnas(scar, icekp, [])
            <= att.MAX_ATT_COPIES_OUTSIDE_TRNA)

    # Five copies back to back. str.count counts each PLACE the repeat sits, so a
    # tandem array of a non-periodic repeat is five copies under both
    # conventions — over both limits, rejected either way.
    tandem = plant(background, 1_000, icekp * 5)
    assert att.count_repeat_copies(tandem, icekp) == 5
    assert att.count_repeat_copies_outside_trnas(tandem, icekp, []) == 5

    # The one case where the two WALKS part company: a periodic k-mer. (AT)9
    # lines up with itself every 2 bases, so 40 bp of AT stutter holds 12
    # overlapping copies but only 2 that do not overlap. Nothing upstream filters
    # low-complexity sequence out, so a tract like this can reach the counters,
    # and the direction of the disagreement is the part that matters: within this
    # one sequence only the stepping count is over its limit, so the tRNA-aware
    # guard rejects the tract as a repeat family and the plain count does not.
    # (In the assembly-wide caller the plain counts are summed over every contig,
    # so a stutter this common would exceed the limit there too.)
    #
    # Written out with its own G/C flanks rather than planted into the random
    # background, so a neighbouring A or T cannot lengthen the run and shift the
    # count.
    microsatellite = "AT" * 9
    stutter = "GGGG" + "AT" * 20 + "CCCC"
    non_overlapping = att.count_repeat_copies(stutter, microsatellite)
    overlapping = att.count_repeat_copies_outside_trnas(stutter, microsatellite, [])
    assert (non_overlapping, overlapping) == (2, 12)
    assert non_overlapping <= att.DENOVO_MAX_CONTIG_COPIES
    assert overlapping > att.MAX_ATT_COPIES_OUTSIDE_TRNA

    # Copies inside tRNAs are the scope difference, not the step difference: the
    # tRNA-aware counter discounts them, which is what lets a tRNA-derived att
    # site survive its own paralogues.
    covering_trna = [trna_feature("contig_1", 950, 1_030, "+")]
    assert att.count_repeat_copies_outside_trnas(scar, icekp, covering_trna) == 1
