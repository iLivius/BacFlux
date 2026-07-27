"""Unit tests for att_search.py — the attL/attR direct-repeat search.

Every test builds a small SYNTHETIC contig with repeats planted at coordinates
the test itself chose, so the expected answer is known exactly and the assertions
can be on the numbers, not on "something was found". That is what spec §8 Phase 7
asks for: "unit tests on synthetic contigs with planted att repeats at known
offsets".

No genome, no Bakta, no database - pure sequence handling.

Run: pytest workflow/scripts/mobilome/test_att_search.py -q
"""

import random

import att_search as att


# ── Building synthetic genomes ──────────────────────────────────────────────

def random_sequence(length, seed):
    """Reproducible random DNA.

    A fixed seed per test keeps failures debuggable: the same "genome" is rebuilt
    byte for byte on every run, so a failing assertion always refers to the same
    sequence. Random background also means any repeat the search finds was
    planted deliberately - chance 25-mers do not occur in a few kb.
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


# ── Mode B: de novo direct repeats ──────────────────────────────────────────

def test_denovo_finds_a_planted_direct_repeat_at_the_exact_offsets():
    """The core case: two copies of one motif bracketing the machinery."""
    left_position = 20_000
    right_position = 50_000
    sequence = build_contig_with_att(ATT_MOTIF, left_position, right_position)

    # Machinery sits between the two planted copies.
    result = att.find_att_sites(sequence, element_start=30_000, element_end=40_000,
                                min_element_bp=8_000)

    assert result["boundary_method"] == "denovo"
    assert result["att_left"] == "%d..%d" % (left_position,
                                             left_position + len(ATT_MOTIF) - 1)
    assert result["att_right"] == "%d..%d" % (right_position,
                                              right_position + len(ATT_MOTIF) - 1)
    assert result["att_sequence"] == ATT_MOTIF
    # The element runs from the START of attL to the END of attR: both repeats
    # are part of the integrated element.
    assert result["element_start"] == left_position
    assert result["element_end"] == right_position + len(ATT_MOTIF) - 1
    assert result["element_length_bp"] == result["element_end"] - result["element_start"] + 1


def test_denovo_returns_none_when_there_is_no_repeat():
    """Random sequence with nothing planted must yield no boundaries at all -
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
    can never be paired up - which is what stops an unrelated duplication
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
    """k counts down from 25, so a 25 bp match wins over a 12 bp one. A long
    exact match between two specific windows is far less likely by chance."""
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
    """Masking must not shift coordinates - every position reported afterwards
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


# ── Mode A: tRNA-anchored ───────────────────────────────────────────────────

def trna_feature(contig, start, end, strand, name="tRNA-Gly(gcc)"):
    return {"contig": contig, "start": start, "end": end,
            "strand": strand, "name": name}


def test_trna_anchored_search_is_preferred_over_denovo():
    """When a tRNA 3' end brackets the element, that answer wins.

    Mode A tested a prediction made in advance from the biology (ICEs integrate
    at tRNA 3' ends); Mode B merely found the best available repeat. So even
    though a de novo repeat is also present here, the tRNA hit must be returned
    and the method reported as 'tRNA'.
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


def test_trna_probe_is_strand_aware():
    """A minus-strand tRNA's 3' end is at its LOWER coordinate, and the probe
    must be reverse-complemented to be read in the gene's own direction.

    Getting this wrong would search for a sequence that is simply not present,
    so the test checks the probe itself rather than only the end result.
    """
    sequence = random_sequence(2_000, seed=8)
    motif = "AAAACCCCGGGGTTTTAAAACCCCG"          # 25 bp
    # Plant the motif at the START of the gene; for a '-' strand gene that IS
    # the 3' end, and the probe should come back as its reverse complement.
    sequence = plant(sequence, 1_000, motif)

    minus_strand = trna_feature("c1", 1_000, 1_075, "-")
    probe = att.trna_three_prime_probe(sequence, minus_strand)
    assert probe == att.reverse_complement(motif)

    # The same gene on the + strand takes its probe from the far end instead.
    plus_strand = trna_feature("c1", 1_000, 1_075, "+")
    plus_probe = att.trna_three_prime_probe(sequence, plus_strand)
    assert plus_probe == sequence[1_075 - 25:1_075]


def test_trna_anchored_tolerates_a_single_mismatch():
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


def test_two_mismatches_are_rejected():
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
    the 'none' result rather than raise - the caller writes these columns for
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


def test_count_mismatches_gives_up_once_over_budget():
    """The comparison short-circuits, and N never counts as a match."""
    assert att.count_mismatches("ACGT", "ACGT", 1) == 0
    assert att.count_mismatches("ACGT", "ACGA", 1) == 1
    assert att.count_mismatches("ACGT", "TGCA", 1) > 1
    # N is ambiguous, so it can never satisfy a match even against itself.
    assert att.count_mismatches("ANGT", "ANGT", 1) >= 1


def test_a_probe_drawn_from_masked_sequence_is_refused():
    """If the probe itself contains N there is nothing to prove, so no
    occurrences may be reported."""
    sequence = "ACGT" * 100
    assert att.find_probe_occurrences(sequence, "ACNT", 1, len(sequence)) == []


# ── FASTA and GFF3 readers ──────────────────────────────────────────────────

def test_read_fasta_keys_on_the_first_token_and_upper_cases(tmp_path):
    path = tmp_path / "genome.fna"
    path.write_text(">NZ_CP006659.2 Klebsiella pneumoniae chromosome\nacgt\nACGT\n"
                    ">contig_2\nTTTT\n")
    sequences = att.read_fasta(str(path))
    assert set(sequences) == {"NZ_CP006659.2", "contig_2"}
    # Soft-masked lower case must be normalised, or repeat search silently misses
    # exactly the repetitive regions it is looking for.
    assert sequences["NZ_CP006659.2"] == "ACGTACGT"


def test_parse_trna_features_reads_bakta_gff3_and_stops_at_the_fasta(tmp_path):
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
    # fixed floor of 12, and comfortably inside the real 15-25 bp att range.
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
# The tests above use i.i.d. uniform random DNA, which is exactly the assumption
# the de novo chance model makes - so they can only ever confirm that model, never
# challenge it. Measured on the KPNIH1 chromosome, 300 randomly placed non-ICE
# spans produced a confident "denovo" boundary 22% of the time, against the ~1.3%
# the model predicts. The difference is not noise: it is the repetitive structure
# that every real chromosome has and random DNA does not.
#
# These tests plant that structure deliberately.

def test_a_dispersed_repeat_family_is_not_an_att_site():
    """An rRNA-operon-like repeat, present many times, must not become a boundary.

    This is the KPNIH1 false positive in miniature. A 25 bp stretch of 16S rRNA
    satisfies every length threshold and genuinely IS an exact direct repeat
    shared by the two flanks - but the cell carries seven rRNA operons, so the
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
    assert result["att_left"] == "30000..30024"
    assert result["att_right"] == "85000..85024"


def test_two_paralogous_trnas_are_not_an_att_pair():
    """Isoacceptor tRNAs share their 3' ends, and a genome carries dozens.

    Mode A builds its probe FROM a tRNA 3' end, so a second tRNA of the same
    species is guaranteed to match it - producing a beautifully bracketing pair
    that is labelled 'tRNA', the method this module treats as its most precise.
    Nine such pairs turned up in 300 random spans of the KPNIH1 chromosome.

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
    """Same sequence, but only the left copy is inside a tRNA - which is what
    site-specific integration at a tRNA 3' end actually leaves behind."""
    sequence = build_contig_with_att(ATT_MOTIF, 30_000, 85_000, length=200_000, seed=505)
    one_trna = [trna_feature("contig_1", 29_952, 30_024, "+")]

    result = att.find_att_sites(
        sequence, element_start=50_000, element_end=60_000,
        trnas=one_trna, min_element_bp=8_000, flank_window_bp=30_000)

    assert result["boundary_method"] == "tRNA"
    assert result["att_left"] == "30000..30024"
    assert result["att_right"] == "85000..85024"
