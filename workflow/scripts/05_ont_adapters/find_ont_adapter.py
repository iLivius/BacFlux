#!/usr/bin/env python3
r"""Find the Oxford Nanopore adapter, and the barcode after it, at the start of ONT reads; check trimmed reads.

Why this is needed
------------------
An ONT read starts with the sequencing adapter that brought the DNA molecule to the pore. In a
barcoded library the barcode and its short flanking sequences follow the adapter. Dorado removes
them when it is run by hand with the kit named, but reads basecalled live, or by older
basecallers, often keep them, and nothing in the FASTQ file says which. Left in the reads, they
end up in the assembly wherever many reads end at the same place: at the ends of linear contigs.
In the BacFlux hybrid assembly of Klebsiella pneumoniae TUM24772 (ONT reads SRR30870266, Rapid
Barcoding Kit SQK-RBK004, R9.4.1), the linear contig_2 ended in about 110 bp of adapter, rapid
barcode RB11 and its flanks, carried unchanged through Flye, dnaapler, Medaka and Polypolish.

fastplong removes such a sequence well when it is told what the sequence is. Left to guess, it
looks only at the first reads of the file and takes the most common sequence at the read ends,
which can be genome (on fungal reads it took the telomere repeat). This script therefore works
the sequence out itself, from all the reads, and gives fastplong the result.

How the sequence is found (subcommand scan)
-------------------------------------------
1. The anchor. The first 150 bp of every read are searched for the 20 bp core of the ONT adapter
   top strand, TTCGTTCAGTTACGTATTGC. The whole adapter is CCTGTACTTCGTTCAGTTACGTATTGCT; the
   ligation, native and rapid adapters of the kit 14 chemistry share it (ONT Chemistry Technical
   Document, section 15). A read holding the core exactly is "anchored". Reads whose adapter
   carries a basecalling error in the core are missed, so the share of anchored reads is a lower
   bound. Older kits use other adapters: the RBK004 rapid adapter ends in TTCGTTCAGTTACG and
   then the rapid flank, without the core's last six bases, so its reads are never anchored and
   step 5 finds them.
2. What follows the adapter. In a barcoded library every anchored read continues with a barcode
   block (native barcoding: AAGGTTAA, the 24 bp barcode, CAGCACCT: 40 bp). In a library without
   barcodes it continues straight into the genome, which differs from read to read. For block
   lengths of 8, 12, ... 120 bp, the script finds the most common sequence that follows the
   adapter, and keeps the longest length whose most common sequence is shared by at least 60 % of
   the anchored reads. On the public Botrytis cinerea B05.10 reads (SRR40741277) the 40 bp native
   block was shared by 87 % of the anchored reads, and the most common 44 bp sequence by 1 %.
3. ONT's published barcodes. The sequence after the adapter is also compared with the native
   barcodes NB01-96 and the rapid barcodes RB01-96 that ONT publishes in its Chemistry Technical
   Document (section 14; the rule passes BacFlux's copy, workflow/data/adapters/ont_chtd.fasta).
   A read counts for a barcode when the 24 bases where the barcode belongs differ from it at 3
   positions or fewer. This names the barcode of the sample and finds every barcode of a file that
   holds several. The block trimmed for each barcode is the one found in its own reads (step 2 on
   those reads alone), unless that block is shorter than the published one: step 2 keeps only a
   stretch that 60 % of the reads share base for base, and with about 1 % basecalling errors per
   base fewer than 60 % of the reads carry a whole 90 bp rapid block without one, so the stretch
   ends inside the Mu end (or inside the 3' flank of a native block). The published block is then
   trimmed instead; fastplong matches it despite the errors. A kit or barcode that is not in the
   list is trimmed all the same, as found in the reads.
4. The decision. fastplong gets the most common sequence with -s, on one line of options, and
   every sequence to trim, in both orientations, in a FASTA file it reads with -a:
   * fewer than 1 % of reads anchored, and fewer than 1 % starting with a barcode block either:
     the reads are already trimmed. fastplong runs with -A, which copies the reads unchanged.
   * at least 1 % anchored, with published barcodes in the reads (each in at least 1 % of the
     anchored reads, and at least 10 reads): the adapter followed by each barcode's block.
     fastplong also searches the read END for the reverse complement of every sequence, which is
     how the barcode and adapter sit at the 3' end of a native-barcoded read.
   * at least 1 % anchored, no published barcode: the adapter followed by the block of step 2,
     when there is one, or the adapter alone (a library without barcodes).
   * fewer than 1 % anchored, but at least 1 % of reads start with a barcode block (step 5): each
     barcode found in at least 1 % of all reads and at least 10 reads.
   * a sequence written in the configuration (parameters.<mode>.ont_adapter): that sequence.
   Two or more barcodes that each hold at least 10 % of the barcoded reads print a NOTE, because
   such a file may hold several samples. A barcode below 10 % is usually a few reads of a
   neighbouring sample on the same flow cell assigned to this one ("bleed-through"); it is
   trimmed without a note, and the second most common barcode is written down with its share
   either way, so the report shows it.
   fastplong also splits a read where it finds a sequence in the middle (two molecules joined
   into one read). It allows the same share of mismatches there as at the ends, so a short
   sequence matches genome by chance: with the 28 bp adapter alone it cut 1,041 B05.10 reads,
   1,001 of them in error, and with the 68 bp adapter + barcode it cut 73, all of them real
   joins. Splitting is therefore switched off (--disable_chimera_splitting) when any sequence
   to trim is shorter than 50 bp.
5. Barcode blocks without the adapter core. A read that is not anchored is searched for a whole
   barcode block in its first 150 bp, in two ways:
   * the native block (AAGGTTAA, 24 bp, CAGCACCT) or the rapid block (GCTTGGGTGTTTAACC, 24 bp,
     and the first 16 bp of the 50 bp sequence after the rapid barcode), with exact flanks;
   * for rapid barcodes, the first 16 bp of that 50 bp sequence alone, GTTTTCGCATTTATCG, with
     the 24 bp in front of it compared with ONT's rapid barcodes (up to 3 differences, as in
     step 3). The 50 bp sequence is the right end of phage Mu, identical to the last 50 bp of
     the Mu genome (NC_000929.1), the end that a Mu transposase joins to DNA. This second search
     exists because the first misses whole runs: on R9.4.1 reads the 16 bases right after the
     adapter are often called wrong. In SRR30870266 the RB11 block sat at the start of 18,284 of
     the 19,091 reads (BLAST), but its 16 bp flank was exact in 37 of them, while RB11 followed
     by the start of the Mu end was exact in 3,180.
   Each barcode found this way in at least 1 % of all reads, and at least 10 reads, is trimmed:
   for a published barcode its published block (fastplong tolerates the basecalling errors and
   removes the unrecognised adapter in front of the block with it), for a block that ONT does not
   publish the block as found. Genome reads that hold a Mu prophage end are far too few to reach
   1 %, and their 24 bp in front of the Mu end are not a barcode.
6. Safety checks. A block that is a tandem repeat of a short unit (10 bp or less, such as a
   homopolymer or a microsatellite) is never added. When the anchored reads continue with two or
   more different sequences that are NOT published barcodes, only the part they share is
   trimmed, with a NOTE: those barcodes then stay in the reads, and the contig adapter screen
   reports them. A block with the flanks of a native or rapid barcode around a barcode that is
   not in the list is trimmed as found, with a NOTE.

The check (subcommand check)
----------------------------
The reads Filtlong kept from the trimmed reads, the reads the assembler gets, are searched again:
for the adapter core in their first 150 bp and, when a block was trimmed, for the last 20 bp of
the most common trimmed sequence. More than 1 % of reads still starting with either prints a
NOTE; it does not stop the run. The second count matters for older kits, whose adapter lacks
the core: there, it is the only sign that trimming missed something.

Usage (rules find_ont_adapter and check_ont_trimming, workflow/rules/shared/05_ont_adapters.smk):
    python find_ont_adapter.py scan --reads FASTQ --setting detect|SEQUENCE --sample NAME \
        [--barcode-list FASTA] --out-scan TSV --out-choice TSV --out-options TXT \
        --out-sequence FASTA
    python find_ont_adapter.py check --reads TRIMMED_FASTQ --choice TSV --sample NAME --out TSV

  --reads          the sample's raw ONT reads (scan), or the reads Filtlong kept after fastplong
                   (check); FASTQ, plain or gzip-compressed (.gz)
  --setting        detect, or the sequence written in parameters.<mode>.ont_adapter
  --barcode-list   ONT's published barcodes (workflow/data/adapters/ont_chtd.fasta); without it,
                   no barcode is named, reads with the adapter get only the block they share (so
                   a file with several barcodes is not trimmed of each), and the Mu-end search of
                   step 5 finds nothing; barcode blocks with exact flanks are still trimmed
  --out-scan       one row per 10,000 reads: how many start with the adapter, how many with a
                   barcode block and no adapter core, and how many with the Mu end and no
                   adapter core; a batch of trimmed reads in an otherwise untrimmed file shows up
                   here
  --out-choice     the decision and the numbers behind it, one "name<TAB>value" line each
  --out-options    the fastplong options, one line: "-s <most common sequence>", with
                   "--disable_chimera_splitting" added when a sequence is shorter than 50 bp, or
                   "-A" when nothing is trimmed
  --out-sequence   every sequence to trim, each in both orientations (fastplong reads it with -a,
                   and -a does not add the reverse complement itself), or an empty file; the
                   contig adapter screen adds it to its list of adapters
  --choice         (check) the --out-choice file of the scan
  --out            (check) "name<TAB>value" lines: reads, and how many still carry the adapter or
                   the end of the trimmed sequence

Lines starting with [BacFlux] are printed for the console (the rules copy them from the log).
Exit status 0; 1 with a message on stderr when a file cannot be read.

00_common.smk imports this file too, for is_tandem_repeat(), clean_sequence(), ADAPTER_CORE and
MIN_SPLIT_LENGTH, to check a sequence given in the configuration while Snakemake reads the workflow.
The method was first written for FunFlux, BacFlux's sister workflow for fungal genomes; the Mu-end
search of step 5, the bleed-through record and the second NOTE of the check were added here. Uses
only the Python standard library. Tests: test_find_ont_adapter.py (same folder).
"""

import argparse
import collections
import gzip
import random
import re
import sys

# The ONT adapter top strand and its 20 bp core (ONT Chemistry Technical Document, section 15).
ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
ADAPTER_CORE = "TTCGTTCAGTTACGTATTGC"
# Where the adapter ends, counted from the start of the core: the core is followed by one more
# adapter base (T) before the block or the genome begins.
ADAPTER_TAIL_AFTER_CORE = len(ADAPTER) - (ADAPTER.index(ADAPTER_CORE) + len(ADAPTER_CORE))

# The flanks around the barcodes (ONT Chemistry Technical Document, section 14). Native:
# AAGGTTAA, 24 bp barcode, CAGCACCT (40 bp). Rapid: GCTTGGGTGTTTAACC, 24 bp barcode, then a 50 bp
# constant sequence (90 bp). The 50 bp sequence (RAPID_TAIL) is the right end of phage Mu: it is
# identical to the last 50 bp of the Mu genome, NC_000929.1. The two patterns find a barcode block
# with no adapter in front; for the rapid block only the first 16 bp of the 50 bp sequence are
# required, and the whole 50 bp are added to the sequence to trim.
NATIVE_BLOCK = re.compile(r"AAGGTTAA[ACGT]{24}CAGCACCT")
RAPID_TAIL = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
RAPID_BLOCK = re.compile(r"GCTTGGGTGTTTAACC[ACGT]{24}" + RAPID_TAIL[:16])
# The Mu-end anchor of step 5: the first 16 bp of RAPID_TAIL. The 24 bp in front of it are where
# a rapid barcode sits.
MU_END_ANCHOR = RAPID_TAIL[:16]
# Where the 24 bp barcode starts inside each kind of block: after the 8 bp native flank, after
# the 16 bp rapid flank.
BARCODE_LENGTH = 24
BARCODE_START = {"native": 8, "rapid": 16}

START_WINDOW = 150            # bp at the read start searched for the adapter or a block
BLOCK_LENGTHS = range(8, 121, 4)
MIN_READ_SHARE = 0.01         # below this share of reads, the adapter (or block) counts as absent
MIN_BLOCK_SHARE = 0.60        # share of anchored reads one block must reach to be trimmed
SECOND_BLOCK_SHARE = 0.10     # a second barcode this common means several barcodes in one file
MIN_BARCODE_SHARE = 0.01      # share of anchored reads a published barcode needs to be trimmed
MIN_BARCODE_READS = 10        # ... and at least this many reads
MAX_BARCODE_MISMATCHES = 3    # positions in which a read's barcode may differ from the published one
MIN_SPLIT_LENGTH = 50         # shorter sequences are trimmed at the ends only
MAX_TANDEM_UNIT = 10          # a block repeating a unit this short or shorter is refused
STORED_FOLLOWERS = 100_000    # how many "what follows the adapter" sequences are kept for steps 2-3
READS_PER_ROW = 10_000        # reads per row of the --out-scan table
CHECK_TAIL_LENGTH = 20        # bp at the end of the trimmed sequence searched for by the check
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def clean_sequence(text):
    """Upper-case a sequence and remove spaces; return None if anything but A, C, G, T is left.

    Takes in:  text, as written in the configuration.
    Used by:   scan() for --setting, and 00_common.smk (parameters.<mode>.ont_adapter).
    """
    sequence = re.sub(r"\s+", "", str(text)).upper()
    if not sequence or not re.fullmatch(r"[ACGT]+", sequence):
        return None
    return sequence


def is_tandem_repeat(sequence, max_unit=MAX_TANDEM_UNIT):
    """True if the sequence is a short unit (max_unit bp or less) repeated end to end.

    TTAGGGTTAGGGTTAGGG is TTAGGG three times and returns True; an adapter or a barcode does not.
    A partial copy at the end counts too (TTAGGGTTAGGGTT), but the unit must be there at least
    twice in full: otherwise any sequence of up to max_unit bp would be "its own unit once".
    Sequences of 8 bp or less are never called a repeat, so that a short block is judged by the
    other rules.
    Used by:   choose_block() and 00_common.smk.
    """
    if len(sequence) <= 8:
        return False
    for unit_length in range(1, max_unit + 1):
        if len(sequence) < 2 * unit_length:
            break
        unit = sequence[:unit_length]
        repeated = unit * (len(sequence) // unit_length + 1)
        if repeated[:len(sequence)] == sequence:
            return True
    return False


def read_sequences(path):
    """Yield the sequence line of every FASTQ record (the second of its four lines).

    Takes in:  path, a FASTQ file; gzip-compressed when the name ends in .gz.
    Does:      reads the file line by line, so a file of any size fits in memory.
    Used by:   scan() and check().
    """
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt")
    else:
        handle = open(path)
    with handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 1:
                yield line.strip().upper()


# ======================================================================
# ONT'S PUBLISHED BARCODES (step 3)
# ======================================================================

def read_barcode_list(path):
    """Read the native and rapid barcode blocks of ONT's published list.

    Takes in:  path, workflow/data/adapters/ont_chtd.fasta (or "" for no list). Its records are
               named ONT_native_barcode_NB01 ... NB96 (AAGGTTAA, barcode, CAGCACCT) and
               ONT_rapid_barcode_RB01 ... RB96 (GCTTGGGTGTTTAACC, barcode, 50 bp sequence); the
               adapter records in the same file are skipped.
    Produces:  a list of dicts, one per barcode: name (NB24), kind (native or rapid), block (the
               barcode with its flanks, the sequence that follows the adapter in a read) and
               barcode (the 24 bp alone).
    Used by:   scan().
    """
    barcodes = []
    if not path:
        return barcodes
    name = None
    parts = []
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))

    for record_name, block in records:
        match = re.fullmatch(r"ONT_(native|rapid)_barcode_([A-Z]+\d+)", record_name)
        if not match:
            continue
        kind, barcode_name = match.group(1), match.group(2)
        start = BARCODE_START[kind]
        barcodes.append({
            "name": barcode_name,
            "kind": kind,
            "block": block,
            "barcode": block[start:start + BARCODE_LENGTH],
        })
    return barcodes


def mismatches(first, second):
    """The number of positions at which two sequences of the same length differ."""
    return sum(1 for a, b in zip(first, second) if a != b)


def build_piece_index(barcodes):
    """Index the published barcodes by four 6 bp pieces, for closest_barcode().

    Takes in:  the list of read_barcode_list().
    Does:      cuts every 24 bp barcode into four pieces of 6 bp and files the barcode under
               (kind, offset, piece). A read barcode with up to 3 differences keeps at least one
               of its four pieces unchanged, so looking up its four pieces finds every published
               barcode within the tolerance, without comparing it with all 192.
    Produces:  a dict: (kind, offset, piece) -> list of barcodes.
    Used by:   count_barcodes() and scan() (the Mu-end search of step 5).
    """
    piece_length = BARCODE_LENGTH // 4
    by_piece = collections.defaultdict(list)
    for barcode in barcodes:
        sequence = barcode["barcode"]
        for offset in range(0, BARCODE_LENGTH, piece_length):
            by_piece[(barcode["kind"], offset, sequence[offset:offset + piece_length])].append(barcode)
    return by_piece


def closest_barcode(piece, kind, by_piece, best=None, best_mismatches=MAX_BARCODE_MISMATCHES + 1):
    """The published barcode of one kind closest to a 24 bp read piece, within the tolerance.

    Takes in:  piece, the 24 bases where a barcode of this kind would sit in a read; kind (native
               or rapid); the index of build_piece_index(); and, optionally, the best barcode
               found so far with its number of differences, so that a caller can compare both
               kinds and keep the overall closest.
    Produces:  (barcode, differences): the closest barcode with at most MAX_BARCODE_MISMATCHES
               differences and fewer than best_mismatches, or the best passed in.
    Used by:   count_barcodes() and scan().
    """
    piece_length = BARCODE_LENGTH // 4
    if len(piece) < BARCODE_LENGTH:
        return best, best_mismatches
    for offset in range(0, BARCODE_LENGTH, piece_length):
        for barcode in by_piece.get((kind, offset, piece[offset:offset + piece_length]), []):
            differences = mismatches(piece, barcode["barcode"])
            if differences < best_mismatches:
                best = barcode
                best_mismatches = differences
    return best, best_mismatches


def count_barcodes(followers, barcodes):
    """Sort the followers by the published barcode their block holds.

    Takes in:  followers, the sequences that follow the adapter in a sample of the anchored reads;
               barcodes, the list of read_barcode_list().
    Does:      looks at the 24 bases where a native barcode would sit (after 8 bases of flank) and
               where a rapid barcode would sit (after 16), and finds the published barcode that
               differs from them at MAX_BARCODE_MISMATCHES positions or fewer (the closest one;
               on a tie the native barcode, which is looked at first, is kept).
    Produces:  a dict: barcode name -> the followers that carry it (scan() counts them and finds
               each barcode's block in them). A follower that matches no barcode is left out.
    Used by:   scan().
    """
    by_piece = build_piece_index(barcodes)
    followers_of = collections.defaultdict(list)
    for follower in followers:
        best = None
        best_mismatches = MAX_BARCODE_MISMATCHES + 1
        for kind, start in BARCODE_START.items():
            piece = follower[start:start + BARCODE_LENGTH]
            best, best_mismatches = closest_barcode(piece, kind, by_piece, best, best_mismatches)
        if best is not None:
            followers_of[best["name"]].append(follower)
    return followers_of


def name_block(block, barcodes):
    """Name the published barcode a block holds, for example "native barcode NB24", or return "".

    Takes in:  a barcode block as found in the reads (native: AAGGTTAA, barcode, CAGCACCT; rapid:
               GCTTGGGTGTTTAACC, barcode, 50 bp sequence) and the list of read_barcode_list().
    Does:      tells the kind of block from its flanks first, and then compares the 24 bases where
               the barcode of that kind sits with the barcodes of that kind only (up to
               MAX_BARCODE_MISMATCHES differences). The kind must come first: ONT uses the same 24 bp
               sequences for many native and rapid barcodes (RB13-RB96 equal NB13-NB96), and only
               the flanks tell them apart.
    Used by:   scan(), for the barcode blocks found with no adapter in front.
    """
    if NATIVE_BLOCK.match(block):
        kind = "native"
    elif RAPID_BLOCK.match(block):
        kind = "rapid"
    else:
        return ""
    start = BARCODE_START[kind]
    piece = block[start:start + BARCODE_LENGTH]
    for barcode in barcodes:
        if barcode["kind"] == kind and mismatches(piece, barcode["barcode"]) <= MAX_BARCODE_MISMATCHES:
            return f"{kind} barcode {barcode['name']}"
    return ""


def looks_like_barcode_block(block):
    """True if the block carries the flanks of a native or rapid barcode around 24 bases."""
    return bool(NATIVE_BLOCK.search(block) or RAPID_BLOCK.search(block))


# ======================================================================
# THE BLOCK THAT FOLLOWS THE ADAPTER (step 2)
# ======================================================================

def choose_block(followers):
    """Find the block that follows the adapter in most anchored reads (step 2 of the method).

    Takes in:  followers, the sequences that follow the adapter (up to 120 bp each) in a sample of
               the anchored reads.
    Does:      for each length in BLOCK_LENGTHS, counts the most common sequence of that length
               and its share of all followers (a follower shorter than the length cannot match).
               The share can only fall as the length grows, because a longer sequence is an
               extension of a shorter one, so the search stops at the first length below
               MIN_BLOCK_SHARE. That next length also tells whether the reads continue with two
               or more different sequences of at least SECOND_BLOCK_SHARE each (several barcodes).
    Produces:  (block, share, several_blocks): the block ("" when none reaches the threshold),
               the share of followers that carry it, and True when several different blocks
               continue after it. A block that is a tandem repeat is returned as "".
    Used by:   scan().
    """
    block = ""
    block_share = 0.0
    several_blocks = False
    if not followers:
        return block, block_share, several_blocks
    for length in BLOCK_LENGTHS:
        counts = collections.Counter(f[:length] for f in followers if len(f) >= length)
        if not counts:
            break
        ranked = counts.most_common(2)
        top_share = ranked[0][1] / len(followers)
        if top_share >= MIN_BLOCK_SHARE:
            block = ranked[0][0]
            block_share = top_share
            continue
        # The first length below the threshold: do two different sequences continue here?
        if len(ranked) == 2 and ranked[1][1] / len(followers) >= SECOND_BLOCK_SHARE:
            several_blocks = True
        break
    # The lengths go in steps of 4 bp, so the block can end up to 3 bp short of its true end (a
    # 90 bp rapid block would stop at 88). Try the next three lengths one base at a time.
    if block:
        for length in range(len(block) + 1, len(block) + 4):
            counts = collections.Counter(f[:length] for f in followers if len(f) >= length)
            if not counts:
                break
            longer, count = counts.most_common(1)[0]
            if count / len(followers) < MIN_BLOCK_SHARE:
                break
            block, block_share = longer, count / len(followers)
    if block and is_tandem_repeat(block):
        return "", 0.0, several_blocks
    return block, block_share, several_blocks


# ======================================================================
# FILES
# ======================================================================

def write_key_values(path, pairs):
    """Write "name<TAB>value" lines, in the given order."""
    with open(path, "w") as out:
        for name, value in pairs:
            out.write(f"{name}\t{value}\n")


def read_key_values(path):
    """Read the "name<TAB>value" lines that write_key_values() wrote, as a dict."""
    values = {}
    with open(path) as handle:
        for line in handle:
            name, _, value = line.rstrip("\n").partition("\t")
            values[name] = value
    return values


def write_sequences(path, sample, sequences):
    """Write the sequences to trim, each in both orientations, as FASTA (an empty file for none).

    Takes in:  sequences, a list of (sequence, label) pairs, the most common first.
    Does:      writes <sample>_trim_<n> and <sample>_trim_<n>_reverse_complement for each one; the
               label (for example "adapter + native barcode NB24") follows the name in the header.
    Used by:   scan().
    """
    with open(path, "w") as out:
        for number, (sequence, label) in enumerate(sequences, start=1):
            out.write(f">{sample}_trim_{number} {label}\n{sequence}\n")
            out.write(f">{sample}_trim_{number}_reverse_complement {label}\n{reverse_complement(sequence)}\n")


def share_text(name, share, suffix=""):
    """One barcode and its share as the report writes it, for example "NB24 92.9 %"."""
    return f"{name} {100 * share:.1f} %{suffix}"


# ======================================================================
# THE TWO SUBCOMMANDS
# ======================================================================

def read_the_reads(path):
    """Read every read once and collect what the decision needs (steps 1, 2 and 5).

    Takes in:  path, the raw ONT reads.
    Does:      for every read, looks for the adapter core in the first START_WINDOW bp (an
               anchored read) and, if it is not there, for a native or rapid barcode block with
               exact flanks and, failing that, for the start of the Mu end (step 5). The sequence
               after the adapter of an anchored read goes into a sample of at most
               STORED_FOLLOWERS sequences drawn evenly from the whole file (reservoir sampling
               with a fixed seed, so that every run gives the same result); the barcode blocks
               without adapter are counted per exact sequence, and the 24 bp in front of the Mu
               end per exact 24 bp sequence.
    Produces:  a dict: reads, anchored, block_only, mu_end_only (counts), followers (the sample),
               block_counts (Counter of exact barcode blocks without adapter), mu_pieces (Counter
               of the 24 bp in front of the Mu end) and rows (the --out-scan rows).
    Used by:   scan().
    """
    sampler = random.Random(20261109)
    found = {"reads": 0, "anchored": 0, "block_only": 0, "mu_end_only": 0, "followers": [],
             "block_counts": collections.Counter(), "mu_pieces": collections.Counter(), "rows": []}
    row_reads = row_anchored = row_block_only = row_mu_end_only = 0

    for sequence in read_sequences(path):
        found["reads"] += 1
        row_reads += 1
        start = sequence[:START_WINDOW]
        core_position = start.find(ADAPTER_CORE)
        if core_position != -1:
            found["anchored"] += 1
            row_anchored += 1
            adapter_end = core_position + len(ADAPTER_CORE) + ADAPTER_TAIL_AFTER_CORE
            follower = sequence[adapter_end:adapter_end + max(BLOCK_LENGTHS)]
            # Reservoir sampling: keep the first STORED_FOLLOWERS, then replace a kept one with
            # falling probability, so that every anchored read had the same chance to be kept.
            if len(found["followers"]) < STORED_FOLLOWERS:
                found["followers"].append(follower)
            else:
                slot = sampler.randrange(found["anchored"])
                if slot < STORED_FOLLOWERS:
                    found["followers"][slot] = follower
        else:
            native = NATIVE_BLOCK.search(start)
            rapid = RAPID_BLOCK.search(start)
            if native:
                found["block_only"] += 1
                row_block_only += 1
                found["block_counts"][native.group(0)] += 1
            elif rapid:
                found["block_only"] += 1
                row_block_only += 1
                found["block_counts"][rapid.group(0)[:-16] + RAPID_TAIL] += 1
            else:
                # Step 5, second search: the start of the Mu end, with room for a barcode in
                # front of it. Whether those 24 bp are a published rapid barcode is decided in
                # scan(), which has the list.
                anchor = start.find(MU_END_ANCHOR)
                if anchor >= BARCODE_LENGTH:
                    found["mu_end_only"] += 1
                    row_mu_end_only += 1
                    found["mu_pieces"][start[anchor - BARCODE_LENGTH:anchor]] += 1
        if row_reads == READS_PER_ROW:
            found["rows"].append((found["reads"] - row_reads + 1, found["reads"], row_reads,
                                  row_anchored, row_block_only, row_mu_end_only))
            row_reads = row_anchored = row_block_only = row_mu_end_only = 0
    if row_reads:
        found["rows"].append((found["reads"] - row_reads + 1, found["reads"], row_reads,
                              row_anchored, row_block_only, row_mu_end_only))
    return found


def barcodes_without_adapter(found, barcodes):
    """Step 5: the barcodes of the reads without the adapter core, counted per barcode.

    Takes in:  found, from read_the_reads(); barcodes, the list of read_barcode_list().
    Does:      names every exact barcode block with name_block(), and every 24 bp piece in front
               of the Mu end with closest_barcode() among the rapid barcodes. Reads of the same
               published barcode are added up over both searches. An exact block that ONT does
               not publish keeps its own entry, as found; a piece in front of the Mu end that is
               no published rapid barcode is left out, because nothing says what it is.
    Produces:  a dict: label -> {"reads": n, "kind": native/rapid/"", "name": NB24 or "",
               "block": the sequence to trim}. For a published barcode the block is its
               published block; for an unpublished one, the block as found in the reads.
    Used by:   scan().
    """
    by_name = {barcode["name"]: barcode for barcode in barcodes}
    entries = {}
    for found_block, count in found["block_counts"].items():
        label = name_block(found_block, barcodes)
        if label:
            name = label.split()[-1]
            entry = entries.setdefault(label, {"reads": 0, "kind": by_name[name]["kind"], "name": name,
                                               "block": by_name[name]["block"]})
        else:
            label = "a barcode not in ONT's published list: " + found_block
            entry = entries.setdefault(label, {"reads": 0, "kind": "", "name": "", "block": found_block})
        entry["reads"] += count
    by_piece = build_piece_index(barcodes)
    for piece, count in found["mu_pieces"].items():
        barcode, _ = closest_barcode(piece, "rapid", by_piece)
        if barcode is None:
            continue
        label = f"rapid barcode {barcode['name']}"
        entry = entries.setdefault(label, {"reads": 0, "kind": "rapid", "name": barcode["name"],
                                           "block": barcode["block"]})
        entry["reads"] += count
    return entries


def scan(args):
    """Subcommand scan: read every read, decide what fastplong removes, write the four files.

    Takes in:  args.reads (the raw ONT reads), args.setting (detect, or a given sequence),
               args.barcode_list (ONT's published barcodes, or ""), args.sample, and the four
               output paths.
    Does:      read_the_reads(), then the decision of step 4 of the method.
    Produces:  the scan table, the choice file, the options line and the sequence FASTA; one
               [BacFlux] summary line and, if needed, [BacFlux] NOTE lines on standard output.
    Used by:   main().
    """
    given = None
    if args.setting.strip().lower() != "detect":
        given = clean_sequence(args.setting)
        if given is None:
            sys.exit(f"[BacFlux] find_ont_adapter {args.sample}: the given sequence is not A, C, G, T only: {args.setting!r}")

    barcodes = read_barcode_list(args.barcode_list)
    found = read_the_reads(args.reads)
    reads = found["reads"]
    followers = found["followers"]
    anchored_share = found["anchored"] / reads if reads else 0.0
    block, block_share, several_blocks = choose_block(followers)

    # Step 3: the published barcodes among the anchored reads, most common first, and those
    # common enough to be trimmed.
    barcode_by_name = {barcode["name"]: barcode for barcode in barcodes}
    followers_of = count_barcodes(followers, barcodes)
    ranked_anchored = sorted(followers_of, key=lambda name: len(followers_of[name]), reverse=True)
    present = []
    for name in ranked_anchored:
        count = len(followers_of[name])
        share = count / len(followers)
        if count >= MIN_BARCODE_READS and share >= MIN_BARCODE_SHARE:
            present.append((name, share))
    barcodes_text = ", ".join(share_text(name, share) for name, share in present) or "none found"

    # Step 5: the barcodes of the reads without the adapter core, most common first.
    without_adapter = barcodes_without_adapter(found, barcodes)
    ranked_without = sorted(without_adapter, key=lambda label: without_adapter[label]["reads"], reverse=True)
    barcoded_reads = sum(entry["reads"] for entry in without_adapter.values())

    # Step 4: the sequences to trim, as (sequence, label) pairs, the most common first, the NOTEs
    # for the console, and the second most common barcode of the path taken (bleed-through).
    to_trim = []
    notes = []
    second_barcode = "none"
    if given:
        to_trim = [(given, "the sequence given in the configuration")]
    elif anchored_share >= MIN_READ_SHARE and present:
        # Published barcodes in the reads: the adapter followed by each barcode's block, found in
        # the reads of that barcode alone (step 2 on them). When those reads share no block that
        # holds the barcode, the published block is used instead.
        for name, share in present:
            barcode = barcode_by_name[name]
            own_block = choose_block(followers_of[name])[0]
            # Keep the block of the reads when it holds this barcode where the barcode belongs,
            # allowing the same few differences as when the reads were counted, and when it is
            # at least as long as the published block. A shorter one stopped where basecalling
            # errors made fewer than 60 % of the reads agree base for base (step 3 of the
            # method): trimming it would leave the rest of the block on the reads.
            start = BARCODE_START[barcode["kind"]]
            own_barcode = own_block[start:start + BARCODE_LENGTH]
            if len(own_barcode) < BARCODE_LENGTH or mismatches(own_barcode, barcode["barcode"]) > MAX_BARCODE_MISMATCHES:
                own_block = barcode["block"]
            elif len(own_block) < len(barcode["block"]):
                own_block = barcode["block"]
            to_trim.append((ADAPTER + own_block, f"adapter + {barcode['kind']} barcode {name}"))
        common = [name for name, share in present if share >= SECOND_BLOCK_SHARE]
        if len(common) >= 2:
            notes.append(
                f"the reads carry {len(common)} barcodes in at least {SECOND_BLOCK_SHARE:.0%} of the reads "
                f"with the adapter ({barcodes_text}); all of them are trimmed. If they are different samples, "
                "demultiplex the reads and run each sample on its own."
            )
        if len(ranked_anchored) >= 2:
            name = ranked_anchored[1]
            second_barcode = share_text(name, len(followers_of[name]) / len(followers), " of the reads with the adapter")
    elif anchored_share >= MIN_READ_SHARE:
        # No published barcode: the adapter followed by the common block, if there is one.
        if block:
            to_trim = [(ADAPTER + block, f"adapter + the {len(block)} bp block after it")]
            # Barcode flanks around a barcode that the list does not hold (only said when the
            # list was given: without it, every barcode would look unknown).
            if barcodes and looks_like_barcode_block(block):
                notes.append(
                    "the barcode after the adapter is not in ONT's published list; it is trimmed as "
                    "found in the reads."
                )
        else:
            to_trim = [(ADAPTER, "adapter (no common block after it)")]
        if several_blocks:
            notes.append(
                "the reads continue with several different sequences after the adapter that are not "
                "published ONT barcodes; only the part they share is trimmed, so those sequences stay in "
                "the reads and the contig adapter screen will report them. Write the whole sequence into "
                "parameters.<mode>.ont_adapter to trim it."
            )
    elif reads and barcoded_reads / reads >= MIN_READ_SHARE:
        # Step 5: barcode blocks without the adapter core in front, each barcode found in at least
        # 1 % of all reads (and 10 reads). The barcodes reported are these, as shares of all reads.
        named = []
        for label in ranked_without:
            entry = without_adapter[label]
            count = entry["reads"]
            if count < MIN_BARCODE_READS or count / reads < MIN_READ_SHARE:
                continue
            if entry["name"]:
                named.append(share_text(entry["name"], count / reads, " of all reads"))
                text = f"{label} without the adapter core in front"
            else:
                text = "a barcode not in ONT's published list, without the adapter core in front"
            to_trim.append((entry["block"], text))
        if named:
            barcodes_text = ", ".join(named)
        common = [label for label in ranked_without
                  if without_adapter[label]["reads"] / barcoded_reads >= SECOND_BLOCK_SHARE
                  and without_adapter[label]["reads"] >= MIN_BARCODE_READS]
        if len(common) >= 2:
            listed = ", ".join(share_text(without_adapter[label]["name"] or "unlisted block",
                                          without_adapter[label]["reads"] / barcoded_reads)
                               for label in common)
            notes.append(
                f"the reads carry {len(common)} barcodes in at least {SECOND_BLOCK_SHARE:.0%} of the reads "
                f"that start with a barcode ({listed}); all of them that reach 1 % of the reads are trimmed. "
                "If they are different samples, demultiplex the reads and run each sample on its own."
            )
        # The second most common PUBLISHED barcode: an unlisted block is usually a read or two
        # with basecalling errors in a flank, not a neighbouring sample.
        published_ranked = [label for label in ranked_without if without_adapter[label]["name"]]
        if len(published_ranked) >= 2:
            entry = without_adapter[published_ranked[1]]
            second_barcode = share_text(entry["name"], entry["reads"] / barcoded_reads,
                                        " of the reads that start with a barcode")

    # The fastplong options: -s with the most common sequence (fastplong searches the read end
    # for its reverse complement), and --disable_chimera_splitting when any sequence to trim is
    # short; -A when there is nothing to trim. The rule always adds -a with the sequence FASTA.
    if to_trim:
        main_sequence, main_label = to_trim[0]
        shortest = min(len(sequence) for sequence, label in to_trim)
        if shortest < MIN_SPLIT_LENGTH:
            options = f"-s {main_sequence} --disable_chimera_splitting"
            splitting = "off"
        else:
            options = f"-s {main_sequence}"
            splitting = "on"
        if len(to_trim) == 1:
            summary = f"{main_label} ({len(main_sequence)} bp)"
        else:
            summary = f"{len(to_trim)} sequences: " + "; ".join(label for sequence, label in to_trim)
    else:
        main_sequence = ""
        options = "-A"
        splitting = "off"
        if anchored_share < MIN_READ_SHARE and (not reads or barcoded_reads / reads < MIN_READ_SHARE):
            summary = "nothing (fewer than 1 % of the reads start with the adapter or a barcode)"
        else:
            summary = "nothing (no single barcode reaches 1 % of the reads and 10 reads)"

    with open(args.out_scan, "w") as out:
        out.write("first_read\tlast_read\treads\treads_with_adapter\treads_with_barcode_block_only"
                  "\treads_with_mu_end_only\n")
        for row in found["rows"]:
            out.write("\t".join(str(value) for value in row) + "\n")
    write_key_values(args.out_choice, [
        ("sample", args.sample),
        ("setting", "given sequence" if given else "detect"),
        ("reads", reads),
        ("reads_with_adapter", found["anchored"]),
        ("percent_with_adapter", f"{100 * anchored_share:.2f}"),
        ("reads_with_barcode_block_only", found["block_only"]),
        ("reads_with_mu_end_only", found["mu_end_only"]),
        ("reads_with_barcode_without_adapter", barcoded_reads),
        ("percent_with_barcode_without_adapter", f"{100 * barcoded_reads / reads:.2f}" if reads else "0.00"),
        ("block", block or "none"),
        ("published_barcodes", barcodes_text),
        ("second_barcode", second_barcode),
        ("trimmed_sequence", main_sequence or "none"),
        ("trimmed_length", len(main_sequence)),
        ("trimmed_sequences", len(to_trim)),
        ("chimera_splitting", splitting),
        ("summary", summary),
        ("note", " ".join(notes) or "none"),
    ])
    with open(args.out_options, "w") as out:
        out.write(options + "\n")
    write_sequences(args.out_sequence, args.sample, to_trim)

    print(
        f"[BacFlux] find_ont_adapter {args.sample}: {found['anchored']:,} of {reads:,} reads "
        f"({100 * anchored_share:.1f} %) start with the ONT adapter, {barcoded_reads:,} with a barcode and "
        f"no adapter core; fastplong removes: {summary}."
    )
    for note in notes:
        print(f"[BacFlux] NOTE find_ont_adapter {args.sample}: {note}")


def check(args):
    """Subcommand check: how many trimmed reads still start with the adapter (or the block)?

    Takes in:  args.reads (the reads Filtlong kept after fastplong), args.choice (the scan's
               choice file), args.sample, args.out.
    Does:      counts the reads with the adapter core in their first START_WINDOW bp and, when
               the most common trimmed sequence was longer than the adapter, the reads with its
               last CHECK_TAIL_LENGTH bp there (barcode left behind).
    Produces:  args.out, "name<TAB>value" lines; a [BacFlux] NOTE when more than MIN_READ_SHARE
               of the reads still start with the adapter, or with the end of the trimmed sequence.
    Used by:   main().
    """
    choice = read_key_values(args.choice)
    construct = choice.get("trimmed_sequence", "none")
    tail = ""
    if construct not in ("none", "") and len(construct) > len(ADAPTER):
        tail = construct[-CHECK_TAIL_LENGTH:]

    reads = 0
    with_adapter = 0
    with_tail = 0
    for sequence in read_sequences(args.reads):
        reads += 1
        start = sequence[:START_WINDOW]
        if ADAPTER_CORE in start:
            with_adapter += 1
        if tail and tail in start:
            with_tail += 1

    adapter_share = with_adapter / reads if reads else 0.0
    tail_share = with_tail / reads if reads else 0.0
    write_key_values(args.out, [
        ("sample", args.sample),
        ("reads", reads),
        ("reads_with_adapter", with_adapter),
        ("percent_with_adapter", f"{100 * adapter_share:.2f}"),
        ("reads_with_end_of_trimmed_sequence", with_tail if tail else "NA"),
        ("percent_with_end_of_trimmed_sequence", f"{100 * tail_share:.2f}" if tail else "NA"),
    ])
    tail_text = f", {with_tail:,} ({100 * tail_share:.2f} %) with the end of the trimmed sequence" if tail else ""
    print(
        f"[BacFlux] check_ont_trimming {args.sample}: after fastplong and Filtlong, {with_adapter:,} of "
        f"{reads:,} reads ({100 * adapter_share:.2f} %) still start with the ONT adapter{tail_text}."
    )
    if adapter_share > MIN_READ_SHARE:
        print(
            f"[BacFlux] NOTE check_ont_trimming {args.sample}: more than {MIN_READ_SHARE:.0%} of the reads "
            "kept after trimming still start with the ONT adapter. The adapter may differ from the one BacFlux searched "
            "for; the contig adapter screen shows whether it reached the assembly."
        )
    if tail_share > MIN_READ_SHARE:
        print(
            f"[BacFlux] NOTE check_ont_trimming {args.sample}: more than {MIN_READ_SHARE:.0%} of the reads "
            "kept after trimming still start with the end of the trimmed sequence, so part of the barcode block "
            "stayed in the reads. The contig adapter screen shows whether it reached the assembly."
        )


def build_parser():
    parser = argparse.ArgumentParser(description="Find the ONT adapter in reads, or check trimmed reads.")
    subcommands = parser.add_subparsers(dest="command", required=True)

    scan_parser = subcommands.add_parser("scan", help="decide what fastplong removes")
    scan_parser.add_argument("--reads", required=True)
    scan_parser.add_argument("--setting", required=True, help="detect, or the sequence to trim")
    scan_parser.add_argument("--sample", required=True)
    scan_parser.add_argument("--barcode-list", default="", help="ONT's published barcodes (FASTA)")
    scan_parser.add_argument("--out-scan", required=True)
    scan_parser.add_argument("--out-choice", required=True)
    scan_parser.add_argument("--out-options", required=True)
    scan_parser.add_argument("--out-sequence", required=True)

    check_parser = subcommands.add_parser("check", help="count the adapters left after fastplong")
    check_parser.add_argument("--reads", required=True)
    check_parser.add_argument("--choice", required=True)
    check_parser.add_argument("--sample", required=True)
    check_parser.add_argument("--out", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "scan":
            scan(args)
        else:
            check(args)
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        sys.exit(f"[BacFlux] {args.command} {args.sample}: cannot read an input file: {error}")


if __name__ == "__main__":
    main()
