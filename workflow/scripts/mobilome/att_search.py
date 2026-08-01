#!/usr/bin/env python3
"""Find the attL/attR direct repeats that mark the real ends of an ICE.

WHAT THIS IS FOR
    CONJscan tells us WHERE the conjugation machinery is, but machinery genes are
    not the edges of the element. An ICE is typically much larger than the cluster
    of tra/relaxase genes that identify it, and everything between the true edges
    - including any AMR gene sitting in the cargo - travels with it when it moves.
    Reporting the machinery span as if it were the element would therefore
    UNDERSTATE what is mobile, which for an AMR report is the dangerous direction.

    This module recovers the real edges from the scar integration leaves behind.

THE BIOLOGY, IN ONE PARAGRAPH
    An ICE lives integrated in the host chromosome. It got there by site-specific
    recombination between a short site on the element (attP) and a matching site
    on the chromosome (attB), catalysed by its own integrase. Because the two
    sites are nearly identical, recombination DUPLICATES them: one copy ends up at
    each end of the integrated element, called attL (left) and attR (right). They
    are DIRECT repeats - same sequence, same orientation - typically 15-25 bp.
    Finding a matching pair that brackets the machinery is therefore direct
    physical evidence of where the element starts and stops.

    The favourite landing site is the 3' end of a tRNA gene. That is not an
    accident: tRNA genes are short, highly conserved, and present in every
    genome, so an element that targets one can integrate into essentially any
    relative without disrupting anything - the recombination reconstitutes the
    tRNA on one side and leaves the second copy as the scar on the other. So a
    repeat pair with ONE copy sitting inside a tRNA is not merely a repeat that
    happens to be there: it is the exact arrangement the biology predicts, and
    the search below treats that as the strongest evidence it can get.

ONE SEARCH, NOT TWO MODES
    An earlier version of this module had two modes: a tRNA-anchored search that
    cut a fixed 25 bp probe from a tRNA's 3' end and looked for a second copy,
    tried first, with a blind repeat scan as the fallback. That design was
    retired because a fixed-length probe cannot work - the ICEKp direct repeat is
    17 bp, and a 25 bp probe cannot match a 17 bp repeat under any mismatch
    budget. The measurements and the reference tools that settled it are written
    out in full above MIN_ATT_REPEAT_BP; read that before changing anything here.

    What runs now is ONE variable-length search, search_maximal_repeat: take
    every exact repeat shared by the left and right flanks, extend each as far as
    the two sides agree, keep the pairs that bracket the machinery, and rank what
    survives - tRNA-anchored first, then longest. The label on the answer says
    which kind of evidence it rests on:

        boundary_method='tRNA'    one copy sits inside an annotated tRNA gene
        boundary_method='denovo'  neither copy does - a repeat with no biological
                                  story attached, so weaker evidence
        boundary_method='none'    no credible pair at all

    Only 'tRNA' is ACTED ON downstream: conjscan_to_ice.py reports a 'denovo'
    pair but refuses to move the element's coordinates onto it. 'none' leaves the
    machinery span standing, which is the honest answer for a fragmented assembly
    whose contig ends before the element does.

THE FAILURE MODE THIS GUARDS AGAINST
    Insertion sequences are themselves flanked by repeats - terminal inverted
    repeats, plus the short direct repeat of target DNA they duplicate on
    transposition. A genome region containing several IS copies is therefore FULL
    of direct repeats that have nothing whatsoever to do with ICE integration. Run
    the de novo scan over raw sequence and those repeats swamp every real signal.
    So the flanks are MASKED against the ISEScan calls first (see mask_intervals);
    the spec calls this out as the single most likely way to get this wrong.

WHERE IT RUNS
    Imported by conjscan_to_ice.py, which already has the candidate machinery
    spans, the Bakta GFF3 (for tRNAs) and the genome. Pure text/sequence handling,
    standard library only - a plain dict of k-mers is far faster than anything
    these window sizes need, which is why no suffix-array tool (vmatch and
    friends) has to be installed. Unit-testable with no tool or database present.

COORDINATES
    Everything crossing this module's boundary is 1-based inclusive, matching
    GFF3 and every other BacFlux table. Python slicing is 0-based half-open, so
    conversions happen at the edges of the helpers below and nowhere else.
"""

import argparse
import csv
import math
import os
import re
import sys


# ── Tunables ────────────────────────────────────────────────────────────────

# How far beyond the machinery span to look for the element's real ends.
#
# 50 kb rather than the original 30 kb, measured: on SPI-7 (AL513382) the correct
# attR sits in tRNA-Phe 5,800 bp OUTSIDE a 30 kb window, so the element was
# reported 50 kb short for want of somewhere to look. Widening to 50 kb recovers
# it. Going further does not help - 80 kb, 120 kb and 200 kb windows change no
# element's answer on the benchmark - so this is the point where the curve flattens
# rather than an arbitrary larger number.
#
# The cost of a wider window is more candidate repeats to rank, which is why the
# ranking rule in search_maximal_repeat matters more at 50 kb than it did at 30.
DEFAULT_FLANK_WINDOW_BP = 50_000

# NOTE ON THE `DENOVO_` PREFIX below. It is historical: these three constants
# were written for the blind "Mode B" scan back when there were two modes. There
# is one search now and they apply to all of it. The names are kept because
# DENOVO_MAX_CONTIG_COPIES is imported by conjscan_to_ice.py, and renaming one
# but not the others would read worse than leaving all three alone.

# Absolute floor, whatever the arithmetic below says. Real att sites are ~15-25 bp
# and nothing shorter than this is worth reporting even in a tiny window.
DENOVO_ABSOLUTE_MIN_REPEAT_BP = 12

# How many chance matches we are willing to tolerate when choosing the shortest
# acceptable repeat length. See minimum_informative_repeat_length: comparing two
# 50 kb flanks (the shipped window), a 12 bp "repeat" is expected about 150 times
# purely by chance, so a fixed 12 bp floor would report noise as an element
# boundary on essentially every genome. This was caught by the unit tests below,
# which expected random sequence to yield no boundaries and got a confident call
# instead.
DENOVO_MAX_EXPECTED_CHANCE_MATCHES = 0.05

# Sanity bounds on the element the boundaries imply. Defaults mirror the spec's
# ice.min_element_bp / ice.max_element_bp; the caller passes the configured ones.
DEFAULT_MIN_ELEMENT_BP = 8_000
DEFAULT_MAX_ELEMENT_BP = 500_000

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ── Sequence and annotation input ───────────────────────────────────────────

def read_fasta(path):
    """Read a FASTA into {sequence_id: sequence}, upper-cased.

    Input:  any nucleotide FASTA - here the genome Bakta annotated, so the IDs
            match the GFF3's seqids and every coordinate lines up.
    Output: dict of id -> sequence string. The id is the first whitespace-
            separated token of the header, the same convention the rest of
            BacFlux uses, so 'NZ_CP006659.2 Klebsiella...' keys on
            'NZ_CP006659.2'.

    Upper-casing matters: some assemblers soft-mask repeats in lower case, and a
    repeat search that treated 'acgt' and 'ACGT' as different would silently miss
    exactly the repetitive regions this module is looking at.
    """
    sequences = {}
    sequence_id = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_id is not None:
                    sequences[sequence_id] = "".join(chunks).upper()
                sequence_id = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                chunks.append(line)
    if sequence_id is not None:
        sequences[sequence_id] = "".join(chunks).upper()
    return sequences


def parse_trna_features(gff3_path):
    """Pull every tRNA gene out of a Bakta GFF3.

    Input:  {sample}.gff3 as Bakta writes it. Bakta runs tRNAscan-SE and records
            the results as `tRNA` features with a `product=tRNA-Glu(ttc)` style
            attribute.
    Output: list of {contig, start, end, strand, name} dicts, 1-based inclusive,
            in file order.

    Bakta appends its own FASTA after a `##FASTA` line; parsing stops there so
    sequence lines are never mistaken for annotation. Anything malformed is
    skipped rather than raising - a missing tRNA costs precision on one element,
    while a crash costs the whole sample.
    """
    trnas = []
    if not os.path.isfile(gff3_path):
        return trnas

    with open(gff3_path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "tRNA":
                continue
            start = to_int(fields[3])
            end = to_int(fields[4])
            if start is None or end is None:
                continue

            # The human-readable tRNA name, for the report. Bakta writes both
            # Name= and product=; either is fine, product= is the more reliable.
            name = ""
            for attribute in fields[8].split(";"):
                key, _, value = attribute.partition("=")
                if key.strip() in ("product", "Name") and value:
                    name = value.strip()
                    if key.strip() == "product":
                        break

            trnas.append({
                "contig": fields[0],
                "start": min(start, end),
                "end": max(start, end),
                "strand": fields[6],
                "name": name or "tRNA",
            })
    return trnas


def to_int(value):
    """Parse an integer, or return None when the field is not one."""
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def reverse_complement(sequence):
    """Reverse complement of a nucleotide string (N-safe)."""
    return sequence.translate(COMPLEMENT)[::-1]


# ── Masking out the insertion sequences ─────────────────────────────────────

def mask_intervals(sequence, intervals):
    """Blank out the given 1-based inclusive intervals with N.

    Input:  a contig sequence, plus the intervals to hide - in practice every
            insertion sequence ISEScan called on that contig.
    Output: a new string of the same length with those stretches replaced by N.

    Why this exists: IS elements carry terminal inverted repeats and leave short
    direct repeats of target DNA behind when they transpose, so a region with a
    few IS copies is densely repetitive for reasons that have nothing to do with
    ICE integration. The de novo scan below skips any k-mer containing N, so
    masking removes those decoys before they can outrank the real att site. The
    spec names this as the most likely way to get att search wrong.

    Masking deliberately does NOT shift coordinates - the string keeps its length
    - so every position reported afterwards still refers to the real genome.
    """
    if not intervals:
        return sequence

    characters = list(sequence)
    for start, end in intervals:
        # Clamp into the sequence and convert 1-based inclusive -> 0-based slice.
        first = max(1, min(start, end)) - 1
        last = min(len(sequence), max(start, end))
        for position in range(first, last):
            characters[position] = "N"
    return "".join(characters)


# ── Reading the tRNA annotation ─────────────────────────────────────────────
#
# These four helpers answer one question between them: does a candidate repeat
# copy sit inside an annotated tRNA gene, and if so which one? That is the test
# that separates an integration scar from a pair of ordinary paralogous genes,
# and it is used in two places - search_maximal_repeat, to rank and to reject,
# and count_repeat_copies_outside_trnas, to count credibly.


def trna_covering(copy, trnas):
    """Which annotated tRNA gene, if any, does this att copy sit inside?

    Same test as overlaps_any_trna but returns the FEATURE rather than a boolean,
    because the caller needs to know which amino acid the tRNA carries - see
    same_trna_species below for why that distinction matters.
    """
    start, end = copy[0], copy[1]
    for trna in trnas:
        if start <= trna["end"] and end >= trna["start"]:
            return trna
    return None


def trna_species(trna):
    """The amino acid a tRNA carries: 'tRNA-Phe(gaa)' -> 'phe'.

    Bakta writes tRNA products as tRNA-<AminoAcid>(<anticodon>). Only the amino
    acid is taken, so the several anticodons of one amino acid count as the same
    species - they are the paralogues the guard below is aimed at.
    Returns '' when the product cannot be parsed, which makes the caller fall back
    to the old, stricter behaviour rather than guessing.
    """
    if not trna:
        return ""
    match = re.search(r"tRNA[-_ ]([A-Za-z]{3})", trna.get("name", ""))
    return match.group(1).lower() if match else ""


def same_trna_species(left_trna, right_trna):
    """True when both att copies sit in tRNAs carrying the SAME amino acid.

    THE DISTINCTION THIS DRAWS, and why it was worth adding. The guard in
    search_maximal_repeat rejects a repeat whose two copies both sit inside tRNA
    genes, on the reasoning that those are two paralogous copies of one tRNA
    rather than the scar of an integration event. That reasoning is sound for
    paralogues and WRONG for anything else, because a genome's tRNA genes are not
    all copies of each other.

    Measured on ICEEc2 (GU725392): its real 22 bp att pair sits in tRNA-Phe at
    one end and tRNA-Ser at the other - two DIFFERENT tRNA species, so not
    paralogues at all. The blanket rule threw the correct boundary away, and the
    element was reported 37 kb short of its true extent. Requiring the SAME amino
    acid keeps the guard's power against real paralogues while letting an element
    that integrated between two unlike tRNAs through.

    When either product cannot be parsed the answer is True, which reproduces the
    old behaviour: an unparseable name is not evidence that the two differ.
    """
    left_species = trna_species(left_trna)
    right_species = trna_species(right_trna)
    if not left_species or not right_species:
        return True
    return left_species == right_species


def overlaps_any_trna(copy, trnas):
    """Does this att copy sit inside an annotated tRNA gene?

    Takes in: one (start, end, mismatches) occurrence, and the contig's tRNA
              features as parsed from the Bakta GFF3.
    Returns:  True if the occurrence overlaps any tRNA at all.

    Used by count_repeat_copies_outside_trnas to enforce the one-in-one-out rule:
    a real integration scar has one copy in the reconstituted host tRNA and the
    other in ordinary sequence. Both inside tRNAs (of the same amino acid) means
    we are looking at two paralogous tRNA genes, not an element.
    """
    start, end = copy[0], copy[1]
    for trna in trnas:
        if start <= trna["end"] and end >= trna["start"]:
            return True
    return False


# ── How short a repeat is still worth believing ─────────────────────────────

def minimum_informative_repeat_length(left_flank_bp, right_flank_bp):
    """Shortest repeat length worth believing when comparing two flanks.

    The problem this solves: any two stretches of DNA share short "repeats" by
    pure chance, and the shorter the k-mer the more of them there are. Comparing
    two windows of L1 and L2 bases, the number of chance k-mer matches is roughly

        L1 * L2 / 4**k

    because there are L1 * L2 pairs of positions and each pair agrees over k
    bases with probability 4**-k. For the shipped 50 kb flanks that is ~149
    expected matches at k=12, ~2.3 at k=15, and ~0.04 at k=18.

    At the shipped window this returns 18 bp - and it also returned 18 at the old
    30 kb window, so widening the window did NOT move the floor. Worth knowing
    before anyone tunes DEFAULT_FLANK_WINDOW_BP expecting it to.

    ⚠ THIS IS A FLOOR, NOT A FILTER — READ THIS BEFORE TRUSTING IT.
    The formula above assumes DNA is a random string of four equally likely
    letters. Real chromosomes are nothing of the sort: they carry rRNA operons
    (7 copies in a typical enterobacterium), REP/BIME elements (hundreds of
    copies), prophage remnants and paralogous gene families. Between two large
    windows of a REAL genome an exact 18-25 bp direct repeat is COMMON.

    Measured on the K. pneumoniae positive control chromosome, 300 randomly
    placed non-ICE spans: this length threshold alone let 22% of them return a
    confident boundary, against the ~1.3% the formula predicts — wrong by ~17x.
    That is why the length test is not the only guard. The guard that does the
    work this function was mistakenly assumed to do is the copy count in
    count_repeat_copies_outside_trnas, applied in search_maximal_repeat.
    """
    pairs = max(1, left_flank_bp) * max(1, right_flank_bp)
    # Smallest k with 4**k >= pairs / threshold.
    required_k = math.log(pairs / DENOVO_MAX_EXPECTED_CHANCE_MATCHES, 4)
    return max(DENOVO_ABSOLUTE_MIN_REPEAT_BP, int(math.ceil(required_k)))


# A genuine attL/attR pair is the scar of ONE integration event: the single
# attP/attB crossover leaves attL at one end and attR at the other, so the repeat
# is expected exactly TWICE on the contig and no more.
#
# Anything occurring more often is a repeat FAMILY, and the families are the
# reason this constant exists. A 25 bp stretch of a 16S rRNA gene is also "an
# exact direct repeat shared by the two flanks", and it clears every length
# threshold - but it occurs seven times, because the cell needs seven ribosomal
# RNA operons. REP/BIME elements run to hundreds of copies through
# enterobacterial intergenic space. Counting copies separates a scar from a
# family cleanly, where measuring length cannot: it was the length test alone
# that let 22% of random spans return a confident boundary (see
# minimum_informative_repeat_length).
#
# Read by conjscan_to_ice.py as well, for its assembly-wide repeat-family check.
DENOVO_MAX_CONTIG_COPIES = 2


def count_repeat_copies(sequence, kmer):
    """How many times this exact repeat occurs in one sequence, both strands.

    Takes in: one sequence, and the repeat to count in it. conjscan_to_ice.py is
              the only caller left: it runs this over EVERY contig of the
              assembly and sums the answers, to ask whether a repeat that looks
              unique on its own contig is really a family member somewhere else.
              That matters on a fragmented draft, where the rest of the family
              may simply have landed on other contigs.
    Returns:  the copy count for that one sequence, which the caller totals and
              compares against DENOVO_MAX_CONTIG_COPIES.

    Both strands are counted because a repeat family can sit either way round
    (rRNA operons certainly do), and a family member pointing the other way is
    still evidence that the sequence is repetitive rather than a unique scar.

    str.count is a C-level scan, so this is cheap enough to run per surviving
    candidate. It counts NON-OVERLAPPING occurrences, which is exactly right
    here: an att site that overlaps itself is a tandem repeat, not a scar.
    """
    copies = sequence.count(kmer)
    reverse = reverse_complement(kmer)
    if reverse != kmer:            # a palindrome would otherwise count twice
        copies += sequence.count(reverse)
    return copies


# ── Maximal-repeat search: what the reference tools actually do ─────────────
#
# WHY THIS REPLACED THE FIXED-LENGTH PROBE
#     The earlier design cut a fixed 25 bp probe from a tRNA 3' end and looked
#     for a second copy. It could not work, and the reason is measurable: the
#     ICEKp direct repeat is 17 bp (CCAGTCAGAGGAGCCAA, Lam et al. 2018,
#     https://pmc.ncbi.nlm.nih.gov/articles/PMC6202445/). A 25 bp probe cannot
#     match a 17 bp repeat under any mismatch budget. On two clinical
#     K. pneumoniae isolates it found nothing at 25 bp and pairs appeared at 18.
#
#     Bacterial ICE direct repeats span roughly 10-60 bp, so NO single length is
#     defensible. Every reference tool searches variable-length:
#       ICEfinder2      vmatch -l 15   (script/single.py L286)
#       icefinder-opt   vmatch -l 15   (2025 fork, unchanged)
#       DEPhT           BLASTN of the left flank vs the right flank
#       DBSCAN-SWA      same, 12 bp floor, ranked by bitscore
#
# WHY EXACT MAXIMAL REPEATS RATHER THAN BLAST
#     `vmatch -l N` returns EXACT MAXIMAL repeats with a minimum length - it
#     finds a seed and extends it as far as the sequences agree. That is what is
#     implemented below, in stdlib Python.
#
#     BLAST (the DEPhT route) would additionally tolerate mismatches and gaps,
#     which is strictly more sensitive. It was not chosen because att_search.py
#     is stdlib-only and the rule that calls it has no conda environment, so
#     shelling out to blastn would make the module depend on whatever happens to
#     be on PATH. Exact maximal repeats is not a compromise invented here: it is
#     precisely what ICEfinder2, the reference ICE tool, uses.
#
# Floor of 15 bp, matching ICEfinder2 and icefinder-opt. DBSCAN-SWA uses 12; 15
# is the more conservative of the two published choices.
MIN_ATT_REPEAT_BP = 15

# Cap on seed matches examined per flank pair. A repeat-dense region can produce
# an enormous number of seeds; without a bound the search would crawl there, and
# a region that dense is not one where a clean att pair is going to be found
# anyway.
#
# ⚠ THE CAP IS A TRUNCATION, NOT JUST A BOUND, AND IT IS CURRENTLY SILENT.
# Hitting it breaks out of the whole scan, so seeds after that point are never
# examined and a real repeat sitting among them would be missed. The element then
# reports boundary_method='none', which is indistinguishable from "the flanks
# carry no repeat" - there is no audit row saying the search was cut short.
# find_maximal_repeats does return a `capped` flag for exactly this purpose, but
# search_maximal_repeat currently discards it; surfacing it would mean a new
# audit row in conjscan_to_ice.py, i.e. a change to the output, so it has been
# left for a deliberate change rather than done quietly here.
#
# Why this is a latent gap rather than a live problem: the re-run recorded in
# docs/methods_att_and_small_plasmids.md covered 246 att searches over 52
# benchmark genomes and no window hit either implementation's internal cap.
MAX_SEED_MATCHES = 200_000


def find_maximal_repeats(left_seq, right_seq, min_length=MIN_ATT_REPEAT_BP):
    """Exact maximal repeats shared by two sequences - `vmatch -l` semantics.

    Takes in: the two flank sequences, and the shortest repeat worth returning.
    Does:     indexes every min_length-mer of the left flank, then for each
              matching k-mer in the right flank EXTENDS the match outwards in
              both directions for as long as the two agree. The result is the
              MAXIMAL repeat containing that seed, whatever length that is.
    Returns:  a list of (length, left_start, left_end, right_start, right_end),
              all 0-based inclusive offsets into the two inputs, plus a flag
              saying whether the seed cap was hit.

    Several seeds inside one long repeat all extend to the same maximal repeat,
    so results are de-duplicated on their extended coordinates rather than on
    the seed - otherwise a 60 bp repeat would be reported 46 times.
    """
    hits = []
    if not left_seq or not right_seq or min_length <= 0:
        return hits, False
    if len(left_seq) < min_length or len(right_seq) < min_length:
        return hits, False

    index = {}
    for offset in range(len(left_seq) - min_length + 1):
        seed = left_seq[offset:offset + min_length]
        if "N" in seed:
            continue                     # masked IS sequence, or an assembly gap
        index.setdefault(seed, []).append(offset)

    seen = set()
    examined = 0
    capped = False
    for right_offset in range(len(right_seq) - min_length + 1):
        seed = right_seq[right_offset:right_offset + min_length]
        if "N" in seed:
            continue
        left_offsets = index.get(seed)
        if not left_offsets:
            continue
        for left_offset in left_offsets:
            examined += 1
            if examined > MAX_SEED_MATCHES:
                capped = True
                break

            # The seed matched, so grow it outwards in both directions for as
            # long as the two flanks keep agreeing. That is what makes the
            # reported repeat MAXIMAL rather than just min_length long.

            # Walk right, starting one base past the seed on each side.
            left_after = left_offset + min_length
            right_after = right_offset + min_length
            while (left_after < len(left_seq) and right_after < len(right_seq)
                   and left_seq[left_after] == right_seq[right_after]
                   and left_seq[left_after] != "N"):
                left_after += 1
                right_after += 1

            # Walk left, starting one base before the seed on each side.
            left_before = left_offset - 1
            right_before = right_offset - 1
            while (left_before >= 0 and right_before >= 0
                   and left_seq[left_before] == right_seq[right_before]
                   and left_seq[left_before] != "N"):
                left_before -= 1
                right_before -= 1

            # Both walks stop one base PAST the repeat - that is the base where
            # the sequences stopped agreeing - so step back in to get the repeat
            # itself, 0-based inclusive.
            left_start, left_end = left_before + 1, left_after - 1
            right_start, right_end = right_before + 1, right_after - 1

            key = (left_start, right_start)
            if key in seen:
                continue                 # another seed inside the same repeat
            seen.add(key)
            hits.append((left_end - left_start + 1,
                         left_start, left_end, right_start, right_end))
        if capped:
            break

    # Longest repeat first: field 0 of each tuple is the repeat's length. The
    # caller re-ranks (tRNA anchoring outranks length), but starting from longest
    # keeps the best-evidence candidates at the front.
    hits.sort(key=lambda hit: hit[0], reverse=True)
    return hits, capped



# How many copies of the repeat may sit OUTSIDE a tRNA. Two, not one, and the
# two cases it has to cover are why:
#
#   tRNA-anchored scar : attL is inside the reconstituted host tRNA and is NOT
#                        counted; attR is outside -> 1 copy outside.
#                        Same-species tRNA paralogues are also inside, so a
#                        genome with five tRNA-Asn genes no longer defeats the
#                        test - which is what the plain copy count got wrong.
#   de novo scar       : neither copy is in a tRNA -> 2 copies outside.
#
# Either way a repeat FAMILY is still rejected: rRNA operons and REP/BIME
# elements have many copies outside tRNAs, which is what this guard was written
# to catch and still catches.
MAX_ATT_COPIES_OUTSIDE_TRNA = 2


def count_repeat_copies_outside_trnas(sequence, kmer, trnas):
    """How many copies of this repeat sit OUTSIDE any annotated tRNA.

    See the note at the call site for why this, and not a plain copy count, is
    the right credibility test. Both strands are counted, because a repeat family
    can sit either way round and a family member pointing the other way is still
    evidence of repetitiveness.
    """
    if not kmer:
        return 0

    # Search both strands, except when the repeat is its own reverse complement.
    # A palindrome reads identically either way, so searching for it twice would
    # count every copy of it twice and make a unique scar look like a family.
    reverse = reverse_complement(kmer)
    probes = [kmer] if reverse == kmer else [kmer, reverse]

    positions = []
    for probe in probes:
        start = sequence.find(probe)
        while start >= 0:
            positions.append(start + 1)        # 1-based
            # Step on by one, not by len(probe), so copies that overlap each
            # other are all counted - tandem repeats must not hide from this.
            start = sequence.find(probe, start + 1)

    outside = 0
    for position in positions:
        if not overlaps_any_trna((position, position + len(kmer) - 1, 0), trnas):
            outside += 1
    return outside


def search_maximal_repeat(sequence, element_start, element_end, trnas=(),
                          flank_window_bp=DEFAULT_FLANK_WINDOW_BP,
                          min_element_bp=DEFAULT_MIN_ELEMENT_BP,
                          max_element_bp=DEFAULT_MAX_ELEMENT_BP,
                          min_repeat_bp=MIN_ATT_REPEAT_BP):
    """Find the att pair by variable-length maximal repeat, the way ICEfinder does.

    Takes in: the contig (already masked against IS calls), the machinery span to
              bracket, and the tRNAs on this contig.
    Returns:  a result dict, or None.

    THE RANKING, and why tRNA anchoring is not simply a separate mode.
        The previous design had two modes: a tRNA-anchored search and a blind
        one, with the tRNA mode strictly preferred and run FIRST off its own
        fixed-length probe. That is the wrong shape. ICEfinder divides the labour
        differently and better - the tRNA LOCATES a candidate site, the repeat
        search DELIMITS it - so there is one search here, and the tRNA is a
        property of the candidates it returns rather than a separate way of
        finding them.

        So: find every maximal repeat that brackets the machinery, discard the
        ones that are not credible, then rank what is left by

            (1) is one copy inside an annotated tRNA?   then  (2) how long is it?

        in that strict order. Length says how unlikely the match is by chance;
        sitting in a tRNA says the match is where integration actually happens.
        The second is evidence about the biology and the first is only evidence
        against coincidence, so length is never allowed to outvote the anchor -
        see the note at the ranking itself for the SPI-7 measurement that settled
        this.

        BOTH copies inside tRNAs of the same amino acid means two paralogous tRNA
        genes, not an att pair, and is rejected outright rather than scored down.
    """
    if not sequence:
        return None

    left_region_start = max(1, element_start - flank_window_bp)
    left_region_end = element_start
    right_region_start = element_end
    right_region_end = min(len(sequence), element_end + flank_window_bp)

    left_flank = sequence[left_region_start - 1:left_region_end]
    right_flank = sequence[right_region_start - 1:right_region_end]
    if not left_flank or not right_flank:
        return None

    # The floor is the LARGER of two requirements, and both are needed:
    #   * MIN_ATT_REPEAT_BP (15) - ICEfinder2's published floor, below which a
    #     repeat is too short to be a credible att site whatever the statistics;
    #   * the chance-match floor for THESE window sizes - because two 50 kb
    #     flanks share a 15 bp repeat by luck about 2.3 times over, so a flat 15
    #     would report noise on any large window. Verified: reinstating this is
    #     what stopped random sequence yielding confident boundaries again.
    #
    # At the shipped 50 kb window this resolves to 18 bp - as it also did at the
    # old 30 kb window, so widening the window did not move the floor.
    #
    # ⚠ NOTE THE AWKWARD BIT, because it is easy to misread. The published ICEKp
    # direct repeat is 17 bp, which is BELOW this 18 bp floor - the module's
    # motivating example would not clear its own threshold if the search were
    # looking for that 17 bp core on its own. It is not: the search returns
    # MAXIMAL repeats, so what it actually reports on those genomes is the longer
    # repeat that CONTAINS the published core and extends past it on one or both
    # sides. That is consistent with the measurement recorded above
    # MIN_ATT_REPEAT_BP - nothing was found at 25 bp and pairs appeared at 18.
    effective_min = max(min_repeat_bp,
                        minimum_informative_repeat_length(len(left_flank), len(right_flank)))
    repeats, _capped = find_maximal_repeats(left_flank, right_flank, effective_min)
    if not repeats:
        return None

    best = None
    for length, l_off, l_end_off, r_off, r_end_off in repeats:
        left_start = left_region_start + l_off
        left_end = left_region_start + l_end_off
        right_start = right_region_start + r_off
        right_end = right_region_start + r_end_off

        if right_start <= left_start:
            continue                       # not a bracketing pair
        element_length = right_end - left_start + 1
        if not (min_element_bp <= element_length <= max_element_bp):
            continue

        left_trna = trna_covering((left_start, left_end, 0), trnas)
        right_trna = trna_covering((right_start, right_end, 0), trnas)
        left_in_trna = left_trna is not None
        right_in_trna = right_trna is not None
        # Two copies of the SAME tRNA are paralogues, not an integration scar.
        # Two copies in DIFFERENT tRNA species are a real possibility and used to
        # be discarded here - see same_trna_species for the ICEEc2 measurement.
        if left_in_trna and right_in_trna and same_trna_species(left_trna, right_trna):
            continue

        kmer = sequence[left_start - 1:left_end]
        # CREDIBILITY: count copies OUTSIDE tRNAs, not copies overall.
        #
        # The plain "occurs at most twice" rule is right for a de novo repeat but
        # WRONG for the tRNA-targeting elements that matter most, and it silently
        # rejected the real thing. ICEKp integrates at tRNA-Asn, so its att core
        # is a piece of the tRNA 3' end - and a genome carrying five tRNA-Asn
        # genes therefore contains that sequence five times whether or not any
        # ICE is present. Measured on both clinical isolates: the published ICEKp
        # repeat CCAGTCAGAGGAGCCAA occurs 5x, so the old rule threw it away.
        #
        # What integration actually leaves is asymmetric:
        #   attL  -> inside the reconstituted host tRNA
        #   attR  -> the second copy, out in ordinary sequence
        #   other tRNA paralogues of the same species -> also inside tRNAs
        # so a genuine scar has exactly ONE copy outside any tRNA. Verified on
        # both genomes: 5 copies, 4 inside tRNAs, 1 outside - and that one is
        # attR, giving a 59,517 bp element in each.
        #
        # The guard keeps all its power against the families it was written for:
        # rRNA operons and REP/BIME elements have many copies OUTSIDE tRNAs and
        # are still rejected.
        if count_repeat_copies_outside_trnas(sequence, kmer, trnas) > MAX_ATT_COPIES_OUTSIDE_TRNA:
            continue

        anchored = left_in_trna or right_in_trna

        # RANK BY ANCHORING FIRST, THEN BY LENGTH.
        #
        # This used to rank on a single score, length + a fixed bonus when the
        # repeat was anchored, so that anchoring "broke ties between repeats of
        # comparable length". That is the wrong ordering, because the two
        # properties are not comparable quantities. Repeat length says how
        # unlikely the match is by chance; sitting at a tRNA 3' end says the match
        # is where integration actually happens. The second is evidence about the
        # biology, the first is only evidence against coincidence, so no bonus
        # large enough to be fair in one case is fair in the other.
        #
        # Measured on SPI-7 (AL513382): a 51 bp repeat in ordinary sequence beat
        # the real 24 bp att pair at tRNA-Phe, and the element came out 50 kb
        # short. Under this ordering the tRNA-anchored pair wins and the call
        # lands on the curated interval.
        #
        # Length still decides among anchored candidates, and among unanchored
        # ones - it is only no longer allowed to outvote the anchor itself. The
        # comparison is strictly greater-than, so among equals the first
        # candidate seen wins, and find_maximal_repeats hands them over longest
        # first.
        candidate = (length, anchored, left_start, left_end,
                     right_start, right_end, kmer)
        if best is None or (anchored, length) > (best[1], best[0]):
            best = candidate

    if best is None:
        return None

    # _length is carried only so the ranking above can read it; the reported
    # length is recomputed from the sequence in build_result.
    _length, anchored, left_start, left_end, right_start, right_end, kmer = best

    # Name the tRNA this element integrated into, for the report. Whichever copy
    # is the one sitting in a tRNA, that is the gene to name.
    trna_name = ""
    if anchored:
        for trna in trnas:
            for start, end in ((left_start, left_end), (right_start, right_end)):
                if start <= trna["end"] and end >= trna["start"]:
                    trna_name = trna.get("name", "")
                    break
            if trna_name:
                break

    return build_result(
        # 'tRNA' when one copy sits in a tRNA - the arrangement integration
        # leaves - and 'denovo' otherwise. Downstream only ACTS on tRNA-anchored
        # boundaries, so this label still carries the same weight it always did.
        method="tRNA" if anchored else "denovo",
        left_copy=(left_start, left_end, 0),
        right_copy=(right_start, right_end, 0),
        att_sequence=kmer,
        mismatches=0,
        trna_name=trna_name,
        trna_start=None,
        trna_end=None,
        repeat_orientation="+",
    )


def build_result(method, left_copy, right_copy, att_sequence, mismatches,
                 trna_name, trna_start, trna_end, repeat_orientation):
    """Package one att-site call into the dict conjscan_to_ice.py writes out.

    Keys mirror the spec's Phase 3 output: attL, attR, att_seq, tRNA,
    boundary_method - plus the refined element interval the boundaries imply and
    the evidence (repeat length) behind the call.

    ON `mismatches`, WHICH IS ALWAYS 0. The search is now exact - maximal repeats
    are grown only while the two flanks agree base for base - so every call made
    today passes 0, and the att_mismatches column in the output table is
    constant. The parameter is kept because the column is part of the published
    output schema and because a future mismatch-tolerant search (BLAST-style, see
    the note above MIN_ATT_REPEAT_BP) would need it again. Do not read a 0 in
    that column as "we checked and there were none".
    """
    return {
        "boundary_method": method,                       # tRNA | denovo
        "att_left": "%d..%d" % (left_copy[0], left_copy[1]),
        "att_right": "%d..%d" % (right_copy[0], right_copy[1]),
        "att_sequence": att_sequence,
        "att_length_bp": len(att_sequence),
        "att_mismatches": mismatches,
        "att_orientation": repeat_orientation,
        "trna": trna_name or "NA",
        "trna_start": trna_start,
        "trna_end": trna_end,
        # The element as the boundaries define it: from the START of attL to the
        # END of attR, since both repeats are part of the integrated element.
        "element_start": left_copy[0],
        "element_end": right_copy[1],
        "element_length_bp": right_copy[1] - left_copy[0] + 1,
    }


NO_BOUNDARY_RESULT = {
    "boundary_method": "none",
    "att_left": "NA",
    "att_right": "NA",
    "att_sequence": "NA",
    "att_length_bp": 0,
    "att_mismatches": 0,
    "att_orientation": "NA",
    "trna": "NA",
    "trna_start": None,
    "trna_end": None,
    "element_start": None,
    "element_end": None,
    "element_length_bp": 0,
}


# ── The entry point conjscan_to_ice.py calls ────────────────────────────────

def find_att_sites(sequence, element_start, element_end, trnas=(),
                   mask=(), flank_window_bp=DEFAULT_FLANK_WINDOW_BP,
                   min_element_bp=DEFAULT_MIN_ELEMENT_BP,
                   max_element_bp=DEFAULT_MAX_ELEMENT_BP):
    """Find the element's real ends. This is the function conjscan_to_ice.py calls.

    Input:  the contig sequence; the machinery span (1-based inclusive) to
            bracket; the tRNAs on this contig (from the Bakta GFF3); the
            intervals to mask (the ISEScan calls); and the search/size bounds.
    Output: always a dict of the shape above - NO_BOUNDARY_RESULT when nothing
            was found, never None, so the caller can write the columns
            unconditionally.

    Two steps, and that is the whole function: mask out the insertion sequences,
    then run the one repeat search over what is left. When it finds nothing the
    machinery span stands unchanged and the report says boundary_method='none' -
    the correct answer for an assembly whose contig ends before the element does.
    """
    if not sequence or element_start is None or element_end is None:
        return dict(NO_BOUNDARY_RESULT)
    if element_start > element_end:
        element_start, element_end = element_end, element_start

    # ONE search, not two modes. See search_maximal_repeat for why the previous
    # tRNA-probe-then-fallback shape was wrong: a fixed-length probe cannot find
    # a repeat shorter than itself, and the ICEKp repeat is 17 bp.
    #
    # Masking against the IS calls applies to the WHOLE search now. Under the old
    # design the tRNA-anchored half was deliberately left UNMASKED, because its
    # probe came from an annotated tRNA and was trusted on that basis. There is
    # no privileged probe any more - every candidate is found by the same repeat
    # search - so IS terminal repeats would flood it exactly as they flood a
    # blind scan. The spec calls this the single most likely way to get this
    # wrong.
    masked_sequence = mask_intervals(sequence, mask)
    hit = search_maximal_repeat(
        masked_sequence, element_start, element_end, trnas,
        flank_window_bp, min_element_bp, max_element_bp)
    if hit is not None:
        return hit

    return dict(NO_BOUNDARY_RESULT)


def group_by_contig(records, contig_key="contig"):
    """Bucket a list of dicts by contig, so per-contig lookups are not O(n) each."""
    grouped = {}
    for record in records:
        grouped.setdefault(record[contig_key], []).append(record)
    return grouped


# ── Thin CLI, for checking one genome by hand ───────────────────────────────

def main(argv=None):
    """Run the att search over a TSV of intervals and print what it found.

    Not used by the workflow - conjscan_to_ice.py imports the functions above
    directly - but invaluable for testing a hypothesis against a real genome
    without running Snakemake, and it keeps this module honestly executable.
    """
    parser = argparse.ArgumentParser(
        description="Find attL/attR direct repeats bracketing an element.")
    parser.add_argument("--genome", required=True,
                        help="genome FASTA (the one Bakta annotated)")
    parser.add_argument("--gff3", default="",
                        help="Bakta GFF3, for the tRNA genes; without it no "
                             "call can be labelled tRNA-anchored")
    parser.add_argument("--intervals", required=True,
                        help="TSV with contig/start/end columns - e.g. an ICE "
                             "candidate table")
    parser.add_argument("--is-table", default="",
                        help="IS element table, masked out before the search")
    parser.add_argument("--flank-window", type=int, default=DEFAULT_FLANK_WINDOW_BP)
    parser.add_argument("--min-element", type=int, default=DEFAULT_MIN_ELEMENT_BP)
    parser.add_argument("--max-element", type=int, default=DEFAULT_MAX_ELEMENT_BP)
    args = parser.parse_args(argv)

    sequences = read_fasta(args.genome)
    trnas_by_contig = group_by_contig(parse_trna_features(args.gff3)) if args.gff3 else {}

    mask_by_contig = {}
    if args.is_table:
        with open(args.is_table, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                contig = row.get("contig")
                start = to_int(row.get("start"))
                end = to_int(row.get("end"))
                if contig and start is not None and end is not None:
                    mask_by_contig.setdefault(contig, []).append((start, end))

    with open(args.intervals, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig = row.get("contig")
            start = to_int(row.get("start"))
            end = to_int(row.get("end"))
            if not contig or start is None or end is None:
                continue
            sequence = sequences.get(contig, "")
            result = find_att_sites(
                sequence, start, end,
                trnas=trnas_by_contig.get(contig, []),
                mask=mask_by_contig.get(contig, []),
                flank_window_bp=args.flank_window,
                min_element_bp=args.min_element,
                max_element_bp=args.max_element)
            print("%s\t%s-%s\tmethod=%s\tattL=%s\tattR=%s\ttRNA=%s\tlen=%s"
                  % (contig, start, end, result["boundary_method"],
                     result["att_left"], result["att_right"], result["trna"],
                     result["element_length_bp"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
