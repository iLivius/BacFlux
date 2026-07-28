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
    tRNA on one side and leaves the second copy as the scar on the other. This is
    why Mode A below, which uses a nearby tRNA's 3' end as the probe, is much
    more trustworthy than a blind repeat search: it is looking for the specific
    thing the biology predicts, not for any repeat that happens to be there.

TWO MODES, TRIED IN THAT ORDER
    Mode A - tRNA-anchored (high precision). Take the last ~25 bp of a tRNA near
        the element, then look for a second copy of it on the other side of the
        machinery. A hit is strong evidence: we predicted the sequence in advance
        from the biology and then found it in the right place.
    Mode B - de novo (lower precision, better recall). No usable tRNA, so look for
        ANY direct repeat shared between the left and right flanks, longest
        first. This will occasionally find a repeat that has nothing to do with
        integration, which is why it is only used when Mode A fails and why its
        result is reported as boundary_method='denovo' so a reader can weight it
        accordingly.

    Neither works -> boundary_method='none' and the machinery span is kept, which
    is the honest answer for a fragmented assembly where the flanks simply are
    not present.

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
import sys


# ── Tunables ────────────────────────────────────────────────────────────────

# How far beyond the machinery span to look for the element's real ends. An ICE
# is usually a good deal larger than its tra cluster, but not unboundedly so;
# 30 kb each side comfortably covers the realistic range while keeping the k-mer
# index small enough to be instant.
DEFAULT_FLANK_WINDOW_BP = 30_000

# Length of the probe taken from a tRNA's 3' end in Mode A. att sites are
# typically 15-25 bp, so 25 is the long end: long enough to be specific in a
# megabase genome, short enough to sit inside a real att site.
TRNA_PROBE_BP = 25

# att sites are conserved but not always perfectly - one substitution between
# attL and attR is common, because only one of the two copies is under selection
# to keep the tRNA functional. More than one and we are no longer looking at a
# recombination scar.
MAX_PROBE_MISMATCHES = 1

# Mode B tries progressively shorter repeats, starting here.
DENOVO_MAX_REPEAT_BP = 25

# Absolute floor, whatever the arithmetic below says. Real att sites are ~15-25 bp
# and nothing shorter than this is worth reporting even in a tiny window.
DENOVO_ABSOLUTE_MIN_REPEAT_BP = 12

# How many chance matches we are willing to tolerate when choosing the shortest
# acceptable repeat length. See minimum_informative_repeat_length: comparing two
# 30 kb flanks, a 12 bp "repeat" is expected about FIFTY times purely by chance,
# so a fixed 12 bp floor would report noise as an element boundary on essentially
# every genome. This was caught by the unit tests below, which expected a random
# sequence to yield no boundaries and got a confident de novo call instead.
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


# ── Approximate matching ────────────────────────────────────────────────────

def count_mismatches(left, right, limit):
    """Hamming distance between two equal-length strings, giving up past `limit`.

    Returns the count, or `limit + 1` as soon as it is exceeded - the caller only
    ever asks "is this within tolerance", so there is no point finishing a
    comparison that has already failed. Any N counts as a mismatch: masked or
    ambiguous sequence must never be allowed to satisfy a match.
    """
    mismatches = 0
    for left_base, right_base in zip(left, right):
        if left_base != right_base or left_base == "N":
            mismatches += 1
            if mismatches > limit:
                return mismatches
    return mismatches


def find_probe_occurrences(sequence, probe, region_start, region_end,
                           max_mismatches=MAX_PROBE_MISMATCHES):
    """Every position in a region where `probe` matches within a mismatch budget.

    Input:  the contig sequence, the probe, and the 1-based inclusive region to
            search in.
    Output: list of (start, end, mismatches) 1-based inclusive, in position order.

    A plain sliding-window scan. For the window sizes here (tens of kb against a
    25 bp probe) this is milliseconds, which is the whole reason no specialised
    index or external tool is needed.
    """
    occurrences = []
    probe_length = len(probe)
    if probe_length == 0:
        return occurrences
    if "N" in probe:
        return occurrences  # a probe drawn from masked sequence proves nothing

    first = max(0, region_start - 1)
    last = min(len(sequence), region_end)
    for offset in range(first, last - probe_length + 1):
        window = sequence[offset:offset + probe_length]
        mismatches = count_mismatches(window, probe, max_mismatches)
        if mismatches <= max_mismatches:
            occurrences.append((offset + 1, offset + probe_length, mismatches))
    return occurrences


# ── Mode A: tRNA-anchored ───────────────────────────────────────────────────

def trna_three_prime_probe(sequence, trna, probe_length=TRNA_PROBE_BP):
    """The last `probe_length` bases of a tRNA gene, in the tRNA's own direction.

    Input:  the contig sequence and one tRNA feature.
    Output: the probe string, or "" when it cannot be taken cleanly.

    Strand awareness is the whole point. A tRNA on the + strand is read left to
    right, so its 3' end is the HIGHER coordinate; one on the - strand is read
    right to left, so its 3' end is the LOWER coordinate and the sequence must be
    reverse-complemented to be read in the gene's own direction. Integration
    happens at the 3' end, so taking the wrong end - or the right end in the
    wrong orientation - would search for a sequence that is simply not there.
    """
    start_index = trna["start"] - 1          # 0-based, inclusive
    end_index = trna["end"]                  # 0-based, exclusive
    if start_index < 0 or end_index > len(sequence) or end_index - start_index < probe_length:
        return ""

    gene_sequence = sequence[start_index:end_index]
    if trna["strand"] == "-":
        gene_sequence = reverse_complement(gene_sequence)
    return gene_sequence[-probe_length:]


def search_trna_anchored(sequence, element_start, element_end, trnas,
                         flank_window_bp, min_element_bp, max_element_bp):
    """Mode A - use a nearby tRNA's 3' end as the att probe.

    Input:  contig sequence; the machinery span (1-based inclusive); the tRNAs on
            THIS contig; and the search/size bounds.
    Output: a result dict (see build_result) or None.

    How it works, following the biology directly:
      1. Consider each tRNA lying within the search region - the machinery span
         plus a flank window each side. One of them may be the gene the element
         integrated into.
      2. Take its 3' end as the probe, in the tRNA's own reading direction.
      3. Find every copy of that probe in the search region, in EITHER genomic
         orientation, because a tRNA on the minus strand yields a probe whose
         copies appear reverse-complemented in plus-strand coordinates. Both
         copies must be in the SAME orientation as each other - attL and attR are
         direct repeats.
      4. Keep pairs that bracket the machinery: one copy at or left of the span's
         start, the other at or right of its end. That bracketing is what makes
         the pair a candidate boundary rather than an incidental repeat.
      5. Of those, prefer the tightest element that still satisfies the size
         bounds, and among equals the pair with the fewest mismatches. Tightest
         wins because a longer interval is easy to manufacture by reaching for a
         more distant copy, whereas the true attL/attR are the innermost pair
         that still contain the machinery.
    """
    search_start = max(1, element_start - flank_window_bp)
    search_end = min(len(sequence), element_end + flank_window_bp)

    best = None
    for trna in trnas:
        # Only tRNAs inside the search region can plausibly be the landing site.
        if trna["end"] < search_start or trna["start"] > search_end:
            continue

        probe = trna_three_prime_probe(sequence, trna)
        if not probe:
            continue

        # Look for the probe in both orientations, but pair only like with like.
        for orientation, oriented_probe in (("+", probe), ("-", reverse_complement(probe))):
            occurrences = find_probe_occurrences(
                sequence, oriented_probe, search_start, search_end)
            if len(occurrences) < 2:
                continue

            left_copies = [o for o in occurrences if o[0] <= element_start]
            right_copies = [o for o in occurrences if o[1] >= element_end]
            for left_copy in left_copies:
                for right_copy in right_copies:
                    if right_copy[0] <= left_copy[0]:
                        continue  # not a bracketing pair
                    element_length = right_copy[1] - left_copy[0] + 1
                    if not (min_element_bp <= element_length <= max_element_bp):
                        continue
                    # EXACTLY ONE copy may sit inside an annotated tRNA, and this
                    # is the test that makes Mode A trustworthy rather than merely
                    # tRNA-flavoured.
                    #
                    # What integration actually leaves behind: the element crosses
                    # over into the 3' end of a tRNA, and the recombination
                    # RECONSTITUTES that tRNA at one end of the element while
                    # placing a second copy of the same short sequence at the far
                    # end, out in ordinary DNA. So one copy is in a tRNA and one is
                    # not.
                    #
                    # If BOTH copies are inside annotated tRNAs, nothing was
                    # integrated: these are two paralogous tRNA genes. Isoacceptor
                    # tRNAs share their 3' ends by design - that is what the probe
                    # is made of - and a genome carries dozens of them, so the
                    # pattern is common and looks superb. It produced 9 spurious
                    # "tRNA-anchored" boundaries in 300 random spans of the KPNIH1
                    # chromosome, each labelled with the method this module treats
                    # as its most precise.
                    if overlaps_any_trna(left_copy, trnas) and overlaps_any_trna(right_copy, trnas):
                        continue
                    total_mismatches = left_copy[2] + right_copy[2]
                    candidate = (element_length, total_mismatches, left_copy,
                                 right_copy, oriented_probe, trna, orientation)
                    # Tightest element first, then fewest mismatches.
                    if best is None or (element_length, total_mismatches) < (best[0], best[1]):
                        best = candidate

    if best is None:
        return None

    _length, mismatches, left_copy, right_copy, oriented_probe, trna, orientation = best
    return build_result(
        method="tRNA",
        left_copy=left_copy,
        right_copy=right_copy,
        att_sequence=oriented_probe,
        mismatches=mismatches,
        trna_name=trna["name"],
        trna_start=trna["start"],
        trna_end=trna["end"],
        repeat_orientation=orientation,
    )


def overlaps_any_trna(copy, trnas):
    """Does this att copy sit inside an annotated tRNA gene?

    Takes in: one (start, end, mismatches) occurrence, and the contig's tRNA
              features as parsed from the Bakta GFF3.
    Returns:  True if the occurrence overlaps any tRNA at all.

    Used to enforce the one-in-one-out rule in search_trna_anchored: a real
    integration scar has one copy in the reconstituted host tRNA and the other in
    ordinary sequence. Both inside tRNAs means we are looking at two paralogous
    tRNA genes, not an element.
    """
    start, end = copy[0], copy[1]
    for trna in trnas:
        if start <= trna["end"] and end >= trna["start"]:
            return True
    return False


# ── Mode B: de novo direct repeats ──────────────────────────────────────────

def minimum_informative_repeat_length(left_flank_bp, right_flank_bp):
    """Shortest repeat length worth believing when comparing two flanks.

    The problem this solves: any two stretches of DNA share short "repeats" by
    pure chance, and the shorter the k-mer the more of them there are. Comparing
    two windows of L1 and L2 bases, the number of chance k-mer matches is roughly

        L1 * L2 / 4**k

    because there are L1 * L2 pairs of positions and each pair agrees over k
    bases with probability 4**-k. For two 30 kb flanks that is ~54 expected
    matches at k=12, ~0.2 at k=16, and effectively zero at k=20.

    ⚠ THIS IS A FLOOR, NOT A FILTER — READ THIS BEFORE TRUSTING IT.
    The formula above assumes DNA is a random string of four equally likely
    letters. Real chromosomes are nothing of the sort: they carry rRNA operons
    (7 copies in a typical enterobacterium), REP/BIME elements (hundreds of
    copies), prophage remnants and paralogous gene families. Between two 30 kb
    windows of a REAL genome an exact 18-25 bp direct repeat is COMMON.

    Measured on the KPNIH1 chromosome, 300 randomly placed non-ICE spans:
    this length threshold alone let 22% of them return a confident "denovo"
    boundary, against the ~1.3% the formula predicts — wrong by ~17x. That is
    why the length test is no longer the only guard; see is_credible_att_repeat,
    which does the work this function was mistakenly assumed to do.
    """
    pairs = max(1, left_flank_bp) * max(1, right_flank_bp)
    # Smallest k with 4**k >= pairs / threshold.
    required_k = math.log(pairs / DENOVO_MAX_EXPECTED_CHANCE_MATCHES, 4)
    return max(DENOVO_ABSOLUTE_MIN_REPEAT_BP, int(math.ceil(required_k)))


# A genuine attL/attR pair is the scar of ONE integration event, so the repeat
# occurs exactly TWICE on the contig: once each side of the element. Anything
# occurring more often is a repeat FAMILY (rRNA, REP/BIME, paralogues), which is
# what actually generates the false boundaries measured above.
DENOVO_MAX_CONTIG_COPIES = 2


def count_repeat_copies(sequence, kmer):
    """How many times this exact repeat occurs on the contig, both strands.

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


# How many ranked candidates get the expensive whole-contig copy check before we
# give up on this repeat length. Each check scans the contig, so an unbounded
# loop would make the search crawl on a repeat-rich genome. Candidates are tested
# best-first, so a real att pair is reached within the first few; a span where
# dozens of candidates in a row are all repeat-family members is a span with no
# credible boundary, and stopping early reaches that answer sooner.
DENOVO_MAX_COPY_CHECKS = 25


def is_locally_unique(left_occurrences, right_occurrences):
    """Cheap first half of the credibility test: unique within each flank?

    Free (the k-mer index already holds the counts), so it runs on every shared
    k-mer. Kills tandem repeats and locally repeated sequence before anything
    expensive happens. NOT sufficient on its own — see is_credible_att_repeat.
    """
    return left_occurrences == 1 and right_occurrences == 1


def is_credible_att_repeat(sequence, kmer, left_occurrences, right_occurrences):
    """Is this shared k-mer plausibly an integration scar, or just repetitive DNA?

    This is the guard that does the real work, and the reasoning is biological
    rather than statistical, which is why the chance-match formula above could
    never have replaced it.

    An att site is created ONCE, when the element recombines into the host: the
    single attP/attB crossover leaves attL at one end and attR at the other. So
    the repeat is expected exactly twice on the contig, full stop.

    A 25 bp stretch of a 16S rRNA gene is also "an exact direct repeat shared by
    the two flanks" — and it satisfies every length threshold — but it occurs
    seven times, because the cell needs seven ribosomal RNA operons. The same
    goes for the REP elements scattered in hundreds of copies through
    enterobacterial intergenic space. Counting copies separates the two cases
    cleanly, where measuring length cannot.

    Two tests, cheap one first:
      1. unique within each flank — kills tandem and locally repeated sequence;
      2. no more than DENOVO_MAX_CONTIG_COPIES copies on the whole contig —
         kills the dispersed families. Test 1 alone is NOT enough: two DIFFERENT
         rRNA operons, one in each flank, are each unique in their own window and
         sail through it. That was the exact false positive found on KPNIH1.
    """
    if not is_locally_unique(left_occurrences, right_occurrences):
        return False
    return count_repeat_copies(sequence, kmer) <= DENOVO_MAX_CONTIG_COPIES


def search_denovo(sequence, element_start, element_end, flank_window_bp,
                  min_element_bp, max_element_bp,
                  max_repeat_bp=DENOVO_MAX_REPEAT_BP,
                  min_repeat_bp=None):
    """Mode B - look for any direct repeat shared by the two flanks.

    Input:  the contig sequence, ALREADY MASKED against the IS calls; the
            machinery span; and the search/size bounds.
    Output: a result dict or None.

    Used only when Mode A found nothing. Because it has no prior expectation
    about what the att site should look like, it is more prone to picking up a
    repeat that has nothing to do with integration - hence the guards:

      * only the flanks are searched, never the element interior, so the repeat
        genuinely brackets the machinery;
      * the sequence is masked against IS calls first (the caller does this), so
        transposon repeats cannot dominate;
      * longest repeat wins - a 25 bp exact match between two specific 30 kb
        windows is far less likely by chance than a short one, so k counts down
        from the top and the first length that yields a valid pair is taken;
      * k never counts down past the point where chance matches become likely,
        which minimum_informative_repeat_length works out from the actual window
        sizes rather than guessing a fixed floor. NOTE that this length rule
        guards only against RANDOM background and is nowhere near sufficient on a
        real chromosome;
      * so the repeat must also be CREDIBLE as an integration scar, meaning it
        occurs exactly twice on the contig - see is_credible_att_repeat. This,
        not the length rule, is what keeps rRNA operons and REP elements out;
      * ties are broken by SYMMETRY. A real att pair sits at comparable distances
        either side of the machinery it brackets, because the machinery is
        somewhere in the middle of the element; a lopsided pair is more likely
        coincidence.
    """
    left_region_start = max(1, element_start - flank_window_bp)
    left_region_end = element_start
    right_region_start = element_end
    right_region_end = min(len(sequence), element_end + flank_window_bp)

    left_flank = sequence[left_region_start - 1:left_region_end]
    right_flank = sequence[right_region_start - 1:right_region_end]
    if not left_flank or not right_flank:
        return None

    # How short a repeat may get before it is indistinguishable from background
    # similarity, derived from these flanks' actual sizes.
    if min_repeat_bp is None:
        min_repeat_bp = minimum_informative_repeat_length(len(left_flank), len(right_flank))
    if max_repeat_bp < min_repeat_bp:
        return None  # the window is so large that no repeat we allow is credible

    for repeat_length in range(max_repeat_bp, min_repeat_bp - 1, -1):
        # Index every k-mer of the left flank once, then stream the right flank
        # past it. This is the "plain dict is fast enough" the spec relies on.
        left_index = {}
        for offset in range(len(left_flank) - repeat_length + 1):
            kmer = left_flank[offset:offset + repeat_length]
            if "N" in kmer:
                continue  # masked-out IS sequence, or an assembly gap
            left_index.setdefault(kmer, []).append(offset)

        if not left_index:
            continue

        # Index the right flank the same way, rather than streaming it past the
        # left index. The extra dict is what makes the copy-number test possible:
        # deciding whether a repeat is a scar or a family member needs to know how
        # often it occurs on BOTH sides, which a one-pass stream cannot tell us.
        right_index = {}
        for offset in range(len(right_flank) - repeat_length + 1):
            kmer = right_flank[offset:offset + repeat_length]
            if "N" in kmer:
                continue
            right_index.setdefault(kmer, []).append(offset)

        # Gather every candidate that passes the FREE tests, then rank them, and
        # only then spend a contig scan on the copy-number test. Doing it in that
        # order matters: the copy test is the one that actually separates scars
        # from repeat families, but it is also ~50x more expensive than
        # everything else here, so it must not run on every shared k-mer.
        candidates = []
        for kmer, right_offsets in right_index.items():
            left_offsets = left_index.get(kmer)
            if left_offsets is None:
                continue
            if not is_locally_unique(len(left_offsets), len(right_offsets)):
                continue
            right_start = right_region_start + right_offsets[0]
            right_end = right_start + repeat_length - 1
            left_start = left_region_start + left_offsets[0]
            left_end = left_start + repeat_length - 1
            if right_start <= left_start:
                continue
            element_length = right_end - left_start + 1
            if not (min_element_bp <= element_length <= max_element_bp):
                continue
            # Symmetry: how differently far the two copies sit from the
            # machinery they bracket.
            asymmetry = abs((element_start - left_end) - (right_start - element_end))
            candidates.append((asymmetry, element_length, left_start, left_end,
                               right_start, right_end, kmer))

        # Best-first: most symmetric, then shortest element.
        candidates.sort(key=lambda item: item[:2])

        for candidate in candidates[:DENOVO_MAX_COPY_CHECKS]:
            _asymmetry, _length, left_start, left_end, right_start, right_end, kmer = candidate
            # THE guard against repeat families. Without it, ~22% of arbitrary
            # spans on a real chromosome come back with a confident boundary
            # built out of rRNA or REP sequence.
            if count_repeat_copies(sequence, kmer) > DENOVO_MAX_CONTIG_COPIES:
                continue
            return build_result(
                method="denovo",
                left_copy=(left_start, left_end, 0),
                right_copy=(right_start, right_end, 0),
                att_sequence=kmer,
                mismatches=0,
                trna_name="",
                trna_start=None,
                trna_end=None,
                repeat_orientation="+",
            )

    return None


# ── Result shape ────────────────────────────────────────────────────────────

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
# anyway. Reaching the cap is reported, never silent.
MAX_SEED_MATCHES = 200_000

# Score bonus for a repeat with exactly ONE copy inside an annotated tRNA - the
# arrangement site-specific integration at a tRNA 3' end actually produces.
# Modest on purpose: it breaks ties between comparable repeats rather than
# letting a marginal 15 bp repeat at a tRNA outrank a convincing 40 bp one.
TRNA_ANCHOR_BONUS_BP = 10


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

            # Extend forwards while the two sequences agree.
            a = left_offset + min_length
            b = right_offset + min_length
            while (a < len(left_seq) and b < len(right_seq)
                   and left_seq[a] == right_seq[b] and left_seq[a] != "N"):
                a += 1
                b += 1
            # Extend backwards.
            c = left_offset - 1
            d = right_offset - 1
            while (c >= 0 and d >= 0
                   and left_seq[c] == right_seq[d] and left_seq[c] != "N"):
                c -= 1
                d -= 1

            left_start, left_end = c + 1, a - 1
            right_start, right_end = d + 1, b - 1
            key = (left_start, right_start)
            if key in seen:
                continue                 # another seed inside the same repeat
            seen.add(key)
            hits.append((left_end - left_start + 1,
                         left_start, left_end, right_start, right_end))
        if capped:
            break

    hits.sort(key=lambda h: h[0], reverse=True)
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
    positions = []
    for probe in (kmer, reverse_complement(kmer)):
        if probe != kmer and probe == kmer:
            continue
        start = sequence.find(probe)
        while start >= 0:
            positions.append(start + 1)        # 1-based
            start = sequence.find(probe, start + 1)
        if reverse_complement(kmer) == kmer:
            break                              # palindrome: do not count twice

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

    THE SCORING, and why tRNA proximity is a bonus rather than a separate mode.
        The previous design had two modes: a tRNA-anchored search and a blind one,
        with the tRNA mode strictly preferred. That is the wrong shape. ICEfinder
        divides the labour differently and better - the tRNA LOCATES a candidate
        site, the repeat search DELIMITS it - and DEPhT scores integrase proximity
        as one term among several rather than as a gate.

        So: find every maximal repeat that brackets the machinery, then rank by

            length  +  bonus if exactly one copy sits in a tRNA

        Longer repeats are better evidence (a 40 bp exact repeat between two
        specific windows is far less likely by chance than a 15 bp one), and a
        repeat with one end in a tRNA is what site-specific integration at a tRNA
        3' end actually leaves behind - the recombination reconstitutes the host
        tRNA at one end and puts the second copy in ordinary sequence at the
        other. BOTH copies inside tRNAs means two paralogous tRNA genes, not an
        att pair, and is rejected outright rather than merely scored down.
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
    #   * the chance-match floor for THESE window sizes - because two 30 kb
    #     flanks share a 15 bp repeat by luck roughly six times over, so a flat
    #     15 would report noise on any large window. Verified: reinstating this
    #     is what stopped random sequence yielding confident boundaries again.
    # For 30 kb flanks this resolves to 18 bp; the published ICEKp repeat is 23,
    # so the real signal is comfortably above it.
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

        left_in_trna = overlaps_any_trna((left_start, left_end, 0), trnas)
        right_in_trna = overlaps_any_trna((right_start, right_end, 0), trnas)
        # Two paralogous tRNAs are not an integration scar - see the note above.
        if left_in_trna and right_in_trna:
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
        # The bonus is deliberately modest: it should break ties between repeats
        # of comparable length, not let a marginal 15 bp repeat at a tRNA beat a
        # convincing 40 bp one elsewhere.
        score = length + (TRNA_ANCHOR_BONUS_BP if anchored else 0)

        candidate = (score, length, anchored, left_start, left_end,
                     right_start, right_end, kmer)
        if best is None or candidate[:2] > best[:2]:
            best = candidate

    if best is None:
        return None

    _score, length, anchored, left_start, left_end, right_start, right_end, kmer = best
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
    the evidence (repeat length, mismatches) behind the call.
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
    """Find the element's real ends, trying the precise method before the loose one.

    Input:  the contig sequence; the machinery span (1-based inclusive) to
            bracket; the tRNAs on this contig; the intervals to mask (the IS
            calls); and the search/size bounds.
    Output: always a dict of the shape above - NO_BOUNDARY_RESULT when nothing
            was found, never None, so the caller can write the columns
            unconditionally.

    Order matters. Mode A is tried first and its answer preferred whenever it
    exists, because a tRNA-anchored hit tested a specific prediction made in
    advance from the biology, while a de novo hit merely found the best repeat
    available. When neither succeeds the machinery span stands unchanged and the
    report says boundary_method='none' - which is the correct answer for a
    short-read assembly whose contig simply ends before the element does.

    Masking is applied only to the DE NOVO search. Mode A's probe comes from an
    annotated tRNA and is verified against the biology, so it should still be
    found even if it happens to lie near an IS; blanking that region would throw
    away the most reliable evidence available.
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
    # design Mode A was deliberately unmasked, because its probe came from an
    # annotated tRNA and was trusted on that basis. There is no privileged probe
    # any more - every candidate is found by the same repeat search - so IS
    # terminal repeats would flood it exactly as they flood a blind scan. The
    # spec calls this the single most likely way to get this wrong.
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
                        help="Bakta GFF3, for tRNA-anchored (Mode A) search")
    parser.add_argument("--intervals", required=True,
                        help="TSV with contig/start/end columns - e.g. an ICE "
                             "candidate table")
    parser.add_argument("--is-table", default="",
                        help="IS element table to mask before the de novo scan")
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
