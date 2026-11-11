# ONT adapter trimming and the adapter screen

An Oxford Nanopore (ONT) read starts with the sequencing adapter that brought the DNA molecule
to the pore, and in a barcoded library the barcode and its flanks follow it. Left in the reads,
both can end up in the assembly. BacFlux therefore removes them from the ONT reads before
assembly, in `nanopore` and `hybrid` mode, and searches every finished genome for adapter
sequence afterwards, in all four modes.

!!! tip "Practical guidance"

    - The default, `ont_adapter: detect`, works out from the reads whether anything is left
      to trim. Reads that are already trimmed pass through unchanged, so leaving the default
      on costs a few minutes and nothing else.
    - Set `ont_adapter: off` to keep an assembly made by an earlier BacFlux unchanged: with
      trimming on, a long-read output folder made before v2.2.0 assembles again (see
      [Upgrading an output folder](#upgrading-an-output-folder)).
    - Read the **Adapters** table of the MultiQC report: it says what was trimmed from the
      reads, which barcodes they carried, and whether the finished genome still holds adapter.
    - The screen covers the common adapters only. Before a GenBank submission, also run
      NCBI's own screen, FCS-adaptor, on the genome.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **adapter** | the short synthetic DNA that brings a molecule to the nanopore; the kit 14 ligation, native and rapid adapters share the top strand `CCTGTACTTCGTTCAGTTACGTATTGCT` |
    | **barcode block** | in a barcoded library, the barcode with its two flanks: for native barcoding `AAGGTTAA`, the 24 bp barcode and `CAGCACCT` (40 bp); for rapid barcoding `GCTTGGGTGTTTAACC`, the 24 bp barcode and a 50 bp constant sequence (90 bp) |
    | **Mu end** | that 50 bp constant sequence: the right end of phage Mu, identical to the last 50 bp of the Mu genome (NC_000929.1) |
    | **trimmed sequence** | what fastplong removes from the reads of one sample |
    | **chimeric read** | one read that holds two molecules joined end to end, with an adapter between them |
    | **end hit** | adapter sequence within 100 bp of the end of a linear contig |
    | **bleed-through** | reads of one sample on a multiplexed flow cell that are assigned to another sample's barcode |

## Why

Dorado can remove adapters and barcodes when it basecalls, but reads made by older basecallers,
or with other settings, can keep them, and a FASTQ file does not say which. In the reads, the adapter and barcode sit at the read ends. In the assembly,
they end up wherever many reads end at the same place, which is at the ends of linear contigs.
A circular contig has no ends, so a complete bacterial genome is rarely affected; an
incomplete one, or a linear replicon, is.

In the BacFlux hybrid assembly of *Klebsiella pneumoniae* TUM24772 made without trimming, the
linear contig `contig_2` (16,144 bp) ended in 112 bp of library construct: part of the rapid
adapter, the rapid flank, rapid barcode RB11 and the Mu end. The ONT reads (SRR30870266) were
made with the Rapid Barcoding Kit SQK-RBK004, and the construct sat at the start of most of
them. It was already at the end of Flye's contig, and it was still there after dnaapler
(which leaves a linear contig unrotated), Medaka and Polypolish. NCBI screens genome
submissions for contamination, and FCS-adaptor, the adaptor screen NCBI provides for genomes
before submission to GenBank, trims adaptor found within 100 bp of a sequence end
([FCS-adaptor output][fcs-output]): it would have flagged that contig end.

Trimming the reads with the exact sequence also lets fastplong split chimeric reads at the
adapter between their two molecules.

## Settings

The key goes in the block of the mode, `parameters.nanopore` or `parameters.hybrid`:

| `ont_adapter` | Effect |
|---|---|
| `detect` (the default; also empty or missing) | BacFlux works the sequence out from the reads, or copies reads that are already trimmed unchanged |
| a sequence of 20 to 200 bases of A, C, G and T | that sequence is trimmed, as written |
| `off` (also an unquoted `off`, `no` or `false`) | no trimming: Filtlong reads the raw reads |

An unquoted `on`, `yes` or `true` names none of these and stops the run. A sequence that is a
short unit repeated end to end, such as a homopolymer or a microsatellite, is refused when the
workflow starts, because trimming it would cut genome from the read ends. A sequence that lacks
the adapter core, or is shorter than 50 bp, gets a NOTE. The adapter screen has no settings
apart from where UniVec is downloaded from (`links.univec_link` and `links.univec_sha256`, see
[Configuration](reference/configuration.md)).

## How the sequence is found

Rule `find_ont_adapter` reads every read of the sample once, in seconds to a minute for an
isolate.

1. **The adapter.** The first 150 bp of every read are searched for the 20 bp core of the
   adapter top strand, `TTCGTTCAGTTACGTATTGC`. The core must match exactly, so reads whose
   adapter carries a basecalling error there are not counted, and the share of reads with the
   adapter is a lower bound. Older kits use other adapters: the RBK004 rapid adapter ends in
   `TTCGTTCAGTTACG` and then the rapid flank, without the core's last six bases, so its reads
   are found by step 4 instead.
2. **The block after the adapter.** In a barcoded library most reads continue after the
   adapter with the same barcode block; in a library without barcodes they continue into the
   genome, which differs from read to read. BacFlux takes the longest block, up to 120 bp, that
   at least 60 % of the reads with the adapter share.
3. **ONT's published barcodes.** The sequence after the adapter is compared with the 96 native
   and 96 rapid barcodes ONT publishes in its Chemistry Technical Document; a read counts for a
   barcode when at most 3 of the 24 barcode bases differ. The list names the barcode and sorts
   the reads of a file that holds several barcodes. The sequence trimmed for each barcode is the
   block found in its own reads, unless that block is shorter than the published one: with about
   1 % basecalling errors per base, fewer than 60 % of the reads carry a whole 90 bp rapid block
   without one, so the shared stretch ends inside the Mu end, and the published block is trimmed
   instead. A kit or barcode the list lacks is trimmed as found.
4. **Barcode blocks without the adapter core.** A read without the core is searched for a whole
   barcode block in its first 150 bp, in two ways: a native or rapid block with exact flanks,
   and, for rapid barcodes, the start of the Mu end (`GTTTTCGCATTTATCG`) with a published rapid
   barcode in the 24 bp in front of it. The second search exists because the first misses whole
   runs: on R9.4.1 reads the 16 bases right after the adapter are often called wrong. In the
   TUM24772 reads, BLAST found the RB11 block at the start of 18,284 of the 19,091 reads, but
   its 16 bp flank was exact in 37 of them; RB11 followed by the start of the Mu end was exact
   in 3,180.
5. **The decision.**
    - When fewer than 1 % of the reads start with the adapter, and fewer than 1 % with a
      barcode block, the reads count as already trimmed and pass through unchanged.
    - When published barcodes follow the adapter, each in at least 1 % of the reads with the
      adapter (and at least 10 reads), the adapter followed by each barcode's block is trimmed.
    - When the adapter is followed by no published barcode, the adapter and the common block
      after it are trimmed, or the adapter alone in a library without barcodes.
    - When the adapter core is missing but at least 1 % of all reads start with a barcode block
      (step 4), each such barcode is trimmed: for a published barcode its published block,
      which fastplong matches despite basecalling errors and removes together with whatever
      unrecognised adapter sits in front of it.

    Two or more barcodes that each hold at least 10 % of the barcoded reads get a NOTE,
    because such a file may hold several samples. A barcode below 10 % is usually bleed-through
    from a neighbouring sample on the same flow cell: it is trimmed when it reaches 1 % of the
    reads, without a note, and the second most common barcode is reported with its share either
    way. A block that is a short tandem repeat is never added. When the adapter is followed by
    several different blocks that are not published barcodes, only the part they share is
    trimmed, with a NOTE.

The file `01.reads/{sample}/ont/adapters/{sample}_adapter_scan.tsv` counts the reads with the
adapter, and with a barcode block and no adapter core, in every block of 10,000 reads. A batch
of already-trimmed reads in an otherwise untrimmed file shows up there.

## What fastplong does with it

Rule `trim_ont_adapters` runs fastplong 0.7.1 with the most common sequence (`-s`) and a FASTA
file of every sequence to trim, in both orientations (`-a`). fastplong searches the start of
every read for the `-s` sequence and the end for its reverse complement; the `-a` sequences are
searched at both ends as written, which is how every barcode of a file holding several is
trimmed, and how a rapid block copied to the read end in its forward orientation is found.
fastplong allows a quarter of the bases to differ and removes the sequence together with
everything on the outer side of it. BacFlux adds:

| Option | Why |
|---|---|
| `--trimming_extension 0` | fastplong's default cuts 10 more bases past the sequence, which after a whole adapter and barcode are genome |
| `-Q -L` | fastplong's quality and length filters stay off: Filtlong, the next step, is the read filter of the workflow |
| `--disable_chimera_splitting`, when a sequence to trim is shorter than 50 bp | fastplong also cuts a read where it finds the sequence in the middle (a chimeric read), with the same tolerance as at the ends, so a short sequence matches genome by chance. On public reads of the fungus *Botrytis cinerea* B05.10 (SRR40741277), the 28 bp adapter alone split 1,041 reads, 1,001 of them in error, while the 68 bp adapter and barcode split 73, all genuine joins |

The two pieces of a split read are named `chimera-segment-1-<read>` and
`chimera-segment-2-<read>`. Read headers are otherwise kept as they are, including the
basecaller tags the Medaka model check reads. When there is nothing to trim, fastplong runs with
`-A` and copies the reads byte for byte. fastplong's own guess of the adapter is never used: it
looks only at the first reads of a file, which can be a batch of already-trimmed reads, and takes
the most common sequence at the read ends, which can be genome.

Rule `check_ont_trimming` then counts, among the reads Filtlong kept (the reads Flye
assembles), those that still start with the adapter core and those that still start with the
end of the trimmed sequence. More than 1 % for either prints a NOTE; the run goes on. The
second count matters for older kits, whose adapter lacks the core.

## The adapter screen

Rule `blast_adapters` searches the finished genome, `contigs_final.fasta`, of every sample in
every mode, with BLAST (`blastn -task blastn-short`) for:

| Sequences | Source |
|---|---|
| ONT adapters, the 96 native barcodes and the 96 rapid barcodes with their flanks | ONT Chemistry Technical Document (shipped with BacFlux) |
| PacBio and MGI adapters | NCBI UniVec (shipped with BacFlux) |
| Illumina adapters and primers | NCBI UniVec build 10.1, downloaded into the run, because Illumina's sequences are copyrighted |
| the sequence trimmed from the sample's reads | rule `find_ont_adapter`, when BacFlux trimmed the reads |

The genome as delivered is searched, not Flye's draft: it is the genome every later step
annotates and the one a user submits, it exists in all four modes, and dnaapler rotates only
circular contigs. A linear contig keeps Flye's ends through dnaapler and polishing.

Rule `screen_adapters` keeps a hit when at least 90 % of its aligned bases are identical and it
is at least 20 bp long within 100 bp of the end of a linear contig, or at least 28 bp long, a
whole adapter, inside a contig. The end rule does not ask how much of an adapter record is
covered, because adapter is often cut short at the very end of a contig. Overlapping hits form
one locus. Two cases need care in bacteria:

- **Circular contigs have no ends.** Their two ends are an arbitrary cut, so every hit on a
  circular contig is judged by the internal rule. To see an adapter lying across that cut, the
  search runs on a copy of each circular contig extended by its own first 200 bp. Topology comes
  from Flye, through the Bakta replicon table of the long-read modes; in `illumina` and
  `contigs` mode it is not known, and every contig is treated as linear.
- **The Mu end is genome in some bacteria.** A genome that carries a Mu prophage holds the same
  50 bp as the end of every rapid barcode. A locus made only of the Mu end is therefore listed
  but not counted when it lies inside a contig, where a rapid-kit remnant would bring its
  barcode along; at the end of a linear contig it counts.

Each locus gets a `counted` (yes/no) and a `reason` column: `adapter_at_contig_end`,
`mu_end_at_contig_end`, `adapter_inside_contig` or `mu_end_inside_contig`. The screen changes no
contig.

!!! note "Before a GenBank submission"

    The lists cover the common adapters only. NCBI provides its own adaptor screen,
    FCS-adaptor ([FCS][fcs]), distributed as a container, to clean genomes before submission
    to GenBank; run it on the genome before you submit.

## Output

| File | Content |
|---|---|
| `01.reads/{sample}/ont/adapters/{sample}_adapter_scan.tsv` | reads with the adapter, with a barcode block and no adapter core, and with the Mu end and no adapter core, per 10,000 reads |
| `01.reads/{sample}/ont/adapters/{sample}_adapter_choice.tsv` | the decision and the numbers behind it: reads with the adapter, the barcodes found, the second barcode, what was trimmed, any NOTE |
| `01.reads/{sample}/ont/adapters/{sample}_fastplong_options.txt`, `{sample}_adapter_sequence.fasta` | what fastplong was given |
| `01.reads/{sample}/ont/adapters/{sample}_fastplong.json`, `.html` | fastplong's own reports: reads and bases in and out, reads trimmed, reads split |
| `01.reads/{sample}/ont/adapters/{sample}_adapter_check.tsv` | reads still starting with the adapter, or with the end of the trimmed sequence, after Filtlong |
| `02.assembly/adapter_lists/univec_illumina.fasta`, `PROVENANCE.txt` | UniVec's Illumina records, and where and when they were downloaded |
| `02.assembly/{sample}/eval/adapter_screen/{sample}_adapter_query.fasta` | every adapter the genome was searched for |
| `02.assembly/{sample}/eval/adapter_screen/{sample}_adapter_blast.tsv` | every BLAST hit, including those that do not pass the rules above |
| `02.assembly/{sample}/eval/adapter_screen/{sample}_adapter_loci.tsv` | one row per locus, with `counted` and `reason` |
| **Adapters** table of the MultiQC report | one row per sample: the read numbers before and after trimming, the barcodes, what was trimmed, the contig ends with adapter, the loci inside contigs and a verdict |

The files in `01.reads/{sample}/ont/adapters/` exist only for ONT samples that BacFlux trimmed.

## Public runs scanned

Before trimming was added, the scan of step 1 to 5 was run on eight public bacterial ONT
runs:

| Run | Organism | Kit | Reads | What the scan found | Decision |
|---|---|---|---|---|---|
| SRR30870266 | *K. pneumoniae* TUM24772 | SQK-RBK004, R9.4.1 | 19,091 | RB11 in front of the Mu end in 19.9 % of the reads (step 4); no adapter core | trim the 90 bp RB11 block |
| ERR16717241 | *Bacillus cereus* | SQK-RBK114 | 70,147 | RB10 without adapter core in 31.8 % of the reads, nearly all of them shorter than 1 kb; RB04 in 2.4 % of the barcoded reads | trim the 90 bp RB10 block |
| SRR41339580 | *K. pneumoniae* | SQK-NBD114-96, basecalled live in MinKNOW 25.09 | 56,128 | adapter in 0.02 %, its barcode NB23 in 0.2 % | nothing to trim |
| SRR41339581 | *Escherichia coli* | as above, same flow cell | 88,966 | adapter in 0.05 %, NB22 in 0.35 % | nothing to trim |
| SRR41339582 | *Citrobacter freundii* | as above, same flow cell | 54,032 | adapter in 0.04 %, NB21 in 0.3 % | nothing to trim |
| SRR41610995 | *Pectobacterium punjabense* | SQK-NBD114 | 81,463 | adapter in 0.87 % | nothing to trim |
| SRR41617180 | *Methylobacterium ajmalii* | SQK-LSK114 with EXP-NBD114, trimmed by its submitters | 95,398 | barcode block in 0.26 % | nothing to trim |
| SRR36766647 | *K. pneumoniae* | not stated | 103,071 | adapter in 0.31 % | nothing to trim |

Six of the eight runs were already trimmed, the three runs basecalled live in MinKNOW among
them. In those three runs of one flow cell, the few reads still starting with a barcode block
carried mostly the sample's own barcode (NB23, NB22 and NB21: 120, 311 and 161 reads); no other
barcode appeared in more than 5 reads. In the *B. cereus* run the rapid construct stayed mainly
on short fragments: 62.7 % of the reads shorter than 1 kb start with it, against 0.1 to 0.2 %
of the longer ones, and Filtlong drops reads shorter than 1 kb. The TUM24772 run carried it at
the start of nearly every read, long and short alike (18,284 of 19,091 by BLAST). When the reads of TUM24772 and of
*B. cereus* were put into one file, as a mislabelled file would be, the scan trimmed both
barcodes and printed the several-barcodes NOTE.

## What it did on *K. pneumoniae* TUM24772

The TUM24772 reads (SRR30870266) went through `nanopore` mode twice, with the same settings,
once with `ont_adapter: detect` and once with `off`:

| | `off` | `detect` |
|---|---|---|
| reads fastplong trimmed | — | 18,636 of 19,091, about 113 bp each: the RB11 block and the adapter in front of it |
| chimeric reads split | — | 13, into 25 pieces |
| reads Flye assembled that still start with the adapter core, or with the end of the block | — | 0, and 30 of 12,116 (0.25 %) |
| Flye's contigs | 6: the same five circles, and a linear contig of 15,344 bp | 5, all circular: the chromosome and four plasmids |
| contig ends with adapter | 1: the end of the linear contig, 88 bp of the RB11 block 24 bp from the end | none |

The linear contig of the untrimmed reads is a near-copy of part of the 23.8 kb plasmid (99.4 %
identity over 11.6 kb), and it ends in the construct. With the reads trimmed, Flye assembles no
such contig; its sequence lies inside the circular plasmid, whose coverage rises from 87× to
106×, about what the linear contig held (17×). The hybrid assembly of the same untrimmed reads
had the same kind of contig, the 16,144 bp `contig_2` described under [Why](#why). The `off`
arm was run as far as Flye, and its contigs were searched with the same screen.

## Upgrading an output folder

With trimming on, Filtlong reads a new file, the trimmed reads, so an output folder made before
v2.2.0 runs Filtlong, Flye and every later step again for every long-read sample, even when the
scan finds nothing to trim: the copy fastplong then makes is still a new file, and Snakemake
compares file times, not contents. In `hybrid` mode the Illumina steps run again too, because
their trimmed reads are temporary files that the later steps need. `ont_adapter: off` keeps
Filtlong's earlier input, and with it the earlier assembly; the adapter screen then runs on it
on its own, as it does in `illumina` and `contigs` mode, where nothing else changes.

## Limits

- The barcode check sees only what is still in the reads. Reads demultiplexed and trimmed at
  basecalling carry nothing to count, so the check then says nothing about bleed-through or a
  mislabelled file.
- Below 1 % of the reads, the adapter or a barcode counts as absent, and nothing is trimmed. A
  few untrimmed reads ending at a linear contig end could still leave adapter there; the screen
  reports it.
- The screen searches for the adapters on its lists. An adapter of another platform or kit is
  not found unless it is on them, or unless BacFlux trimmed it from the reads.

## References

- Oxford Nanopore Technologies. Chemistry Technical Document, version CHTD_500_v1_revAV_24July2026,
  sections 14 and 15 (adapter and barcode sequences).
  <https://nanoporetech.com/document/chemistry-technical-document>
- Chen S (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and
  deduplication using fastp. *iMeta* 2:e107. <https://doi.org/10.1002/imt2.107>
- Chen S, Zhou Y, Chen Y, Gu J (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor.
  *Bioinformatics* 34:i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>
  (fastplong, <https://github.com/OpenGene/fastplong>, asks to be cited through these two
  papers.)
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009). BLAST+:
  architecture and applications. *BMC Bioinformatics* 10:421.
  <https://doi.org/10.1186/1471-2105-10-421>
- National Center for Biotechnology Information. The UniVec Database.
  <https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/> (files at
  <https://ftp.ncbi.nlm.nih.gov/pub/UniVec/>)
- NCBI Foreign Contamination Screen (FCS), FCS-adaptor. [FCS][fcs]; the trimming rules in
  [FCS-adaptor output][fcs-output].
- Enterobacteria phage Mu, complete genome. NCBI RefSeq NC_000929.1, derived from GenBank
  AF083977.1, the sequence of Morgan GJ, Hatfull GF, Casjens S, Hendrix RW (2002).
  Bacteriophage Mu genome sequence: analysis and comparison with Mu-like prophages in
  Haemophilus, Neisseria and Deinococcus. *Journal of Molecular Biology* 317:337–359.
  <https://doi.org/10.1006/jmbi.2002.5437>

[fcs]: https://github.com/ncbi/fcs
[fcs-output]: https://github.com/ncbi/fcs/wiki/FCS-adaptor-output
