# Method reference: hybrid mode on five closed genomes

Five bacterial isolates with published, high-accuracy reference genomes, run end to end
through `hybrid` mode with every optional layer enabled. The genomes come from the
Autocycler benchmark, where each was sequenced deeply on ONT and Illumina and assembled
to a curated reference — which makes them the rare case where BacFlux's output can be
checked against a known answer rather than against another prediction.

Two things came out of it: a measurement of when short-read correction of long reads
helps and when it does harm, and an end-to-end picture of what the workflow delivers on
genomes whose true content is known.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies, the long-read platform |
    | **IS** | insertion sequence — the smallest mobile element, and the main cause of contig breaks |
    | **AMR** | antimicrobial resistance |
    | **ICE** | integrative and conjugative element — carries its own conjugation machinery |
    | **replicon** | one complete DNA molecule: a chromosome or a plasmid |
    | **concatemer** | a circular sequence assembled as two or more tandem copies of itself |

## The isolates

| species | reference | replicons | IS per Mb | Illumina depth |
|---|--:|--:|--:|--:|
| *Providencia rettgeri* | 4,465,806 bp | 2 | 1.8 | 102× |
| *Listeria innocua* | 2,972,545 bp | 2 | 3.7 | 223× |
| *Enterobacter hormaechei* | 5,384,747 bp | 4 | 7.4 | 92× |
| *Klebsiella pneumoniae* | 5,990,196 bp | 5 | 15.7 | 82× |
| *Shigella flexneri* | 4,828,487 bp | 8 | 81.2 | 67× |

IS density is counted by the mobilome module on each finished assembly, so it is a
number the workflow produces rather than one taken from the literature.

!!! warning "The published depth is not a realistic depth"

    The ONT sets published with these isolates run to **306×–1,070×**, because they were
    sequenced for an assembler benchmark. A sequencing facility multiplexing bacterial
    isolates on one flow cell delivers a small fraction of that, and nobody outside a
    methods study will see these depths.

    That difference is not cosmetic. `hybrid` mode applies **no coverage cap** to its ONT
    reads — the short reads are expected to guide the assembly instead — and at these
    depths Flye assembles nothing at all, reporting `Assembled 0 disjointigs` on all five
    genomes. The runs described here therefore use the **50× subsets** published
    alongside the full sets, which is also the depth the original benchmark assembled.

    If you are running `hybrid` on unusually deep ONT data, subsample it first.

## Short-read correction of long reads: it can help, and it can lose a replicon

In `hybrid` mode Filtlong is given the Illumina reads (`-1`/`-2`) and scores each long
read by how well it agrees with them, trimming (`--trim`) and splitting (`--split`)
where agreement fails. `nanopore` mode has no short reads to offer, so it filters on
length and quality alone and caps coverage with `--target_bases`.

Three arms were assembled with Flye from the same 50× reads, differing only in the
filtering before them:

| arm | filtering |
|---|---|
| **guided** | Filtlong with `-1`/`-2`, `--trim`, `--split` — what `hybrid` ships |
| **unfiltered** | none |
| **unguided** | Filtlong with `--min_length`, `--keep_percent`, `--target_bases` — the `nanopore` settings |

| genome | reference | guided | unfiltered | unguided |
|---|--:|--:|--:|--:|
| *Providencia rettgeri* | 2 | 2c, +25 | 2c, +0 | 2c, +0 |
| *Listeria innocua* | 2 | **2c, −1,134** | 3c, +6,020 | 3c, +6,019 |
| *Enterobacter hormaechei* | 4 | **3c, −91,386** | 4c, −1,134 | 4c, +1,418 |
| *Klebsiella pneumoniae* | 5 | 4c, +2,274 | 4c, +2,281 | 4c, −1,234 |
| *Shigella flexneri* | 8 | **40c, −35,210** | 6c, +22,974 | 6c, +21,340 |

The unguided arm tracks the unfiltered arm on every genome. Filtlong's own length and
quality filtering is therefore harmless; **the entire effect comes from the short-read
guidance**, and the mechanism is visible in the flags. `--split` severs a read wherever
1,000 consecutive bases lack short-read support, and the regions where support thins are
repeats and IS elements — the same regions whose spanning reads hold an assembly
together. The severed fragments still pass `--min_length`, so mean read length and N50
both *rise* while contiguity collapses.

**On *Shigella*, with 81 IS per Mb and the shallowest short-read coverage of the set,
guidance turned a 6-contig assembly into 40.** On *Enterobacter* it cost a whole
replicon and 91 kb. On *Listeria* — the least repetitive genome, with by far the deepest
short reads — it was the only arm to recover both replicons.

!!! note "What can be said in advance, and what cannot"

    Five genomes cannot support a threshold. What they do show is a direction: guidance
    was harmful where IS density was high or short-read coverage was thin, and helpful on
    the one genome that was neither. Both quantities are measurable before assembly for
    the short reads, and only after assembly for IS density.

    The practical reading: **treat short-read guidance as unsafe on repeat-rich genomes
    and on shallow short-read data.** A fragmented ONT assembly next to a much less
    fragmented short-read assembly of the same isolate is the signal that it has bitten.

## What the workflow delivered

Assemblies, against the published references, using the shipped configuration:

| genome | contigs | assembled | difference |
|---|--:|--:|--:|
| *Providencia rettgeri* | 2 | 4,465,802 | **−4 bp** |
| *Listeria innocua* | 2 | 2,971,407 | −1,138 bp |
| *Klebsiella pneumoniae* | 4 | 5,992,470 | +2,274 bp |
| *Enterobacter hormaechei* | 3 | 5,293,335 | −91,412 bp |
| *Shigella flexneri* | 40 | 4,793,237 | −35,250 bp |

*Providencia* is four bases from a reference built with Trycycler, Medaka, Polypolish and
Pypolca. Two failure modes account for the rest, and both are worth recognising in your
own data.

**Small plasmids go missing.** *Klebsiella*'s reference carries a 1,240 bp plasmid and
*Enterobacter*'s an 80,745 bp one; neither survived. This is the documented small-plasmid
loss of long-read assembly, described in
[att sites and small plasmids](methods_att_and_small_plasmids.md).

**Small plasmids also arrive doubled.** On *Shigella*, three of the six small plasmids
assembled as exact tandem concatemers — 2 × 6,790, 3 × 3,181 and 4 × 3,835 bp, matching
whole multiples to within 9 bp. A read longer than a small circle wraps it more than
once, and the assembler represents that as a tandem repeat. The consequence is not
cosmetic: the replicon looks several times its true size and its genes appear more than
once, so gene content and copy number read from it are wrong. Four of the eight replicons
were recovered exactly, three as concatemers, two lost.

### Mobilome

Stage `08.mobilome` with all four optional layers enabled:

| genome | AMR genes | IS elements | ICE/IME candidates | named elements | mobility tiers |
|---|--:|--:|--:|--:|---|
| *Listeria innocua* | 2 | 11 | 1 | 0 | all tier 1 |
| *Providencia rettgeri* | 5 | 8 | 1 | 0 | 3 × tier 1, 2 × tier 6 |
| *Shigella flexneri* | 39 | 392 | 0 | 0 | tiers 1, 2, 3, 5 |
| *Enterobacter hormaechei* | 44 | 40 | 4 | 2 | 10 × tier 1, **15 × tier 4**, 19 × tier 6 |
| *Klebsiella pneumoniae* | 55 | 94 | 3 | 2 | 11 × tier 1, 44 × tier 6 |

The distribution matches what is known about these organisms. *Listeria innocua* is a
non-pathogen and returns two AMR genes, both intrinsic candidates. *Klebsiella* returns
44 genes as predicted self-transmissible. *Shigella*'s 392 IS elements are 4× the next
genome's and are the direct cause of its fragmented assembly.

*Enterobacter* produced **15 tier 4 rows**, all fifteen genes inside a single curated
transposon, Tn*SMR478*, on the chromosome — a stress- and metal-resistance cluster. Tier
4 requires the opt-in TnCentral layer, a curated element covering the gene, and no
plasmid or ICE claiming it first ([the mobility ladder](mobilome/mobility-ladder.md)),
which is why it is rare; a single chromosomal transposon carrying fifteen genes satisfies
all three at once.

Read that 15 as one element rather than fifteen findings. The other genome where tier 4
has been measured, ATCC BAA-2146, produced exactly one row from seven named elements,
because six of the seven sat on plasmids and were claimed by tier 5 or 6 first
([Worked example](mobilome/worked-example.md)). The two results bracket the behaviour:
the row count tracks genes, and how many rows an element yields is a property of the
element.

## Acknowledgement

The isolates, the sequencing, the reference genomes and the read subsets are the work of
Ryan Wick, Benjamin Howden and Timothy Stinear, published with the Autocycler paper and
released on Figshare. BacFlux assembles with Flye alone and is not a competitor to a
consensus assembler; where an assembly matters more than the analysis that follows it,
Autocycler is the better tool, and its output can be fed to BacFlux through `contigs`
mode.

> Wick RR, Howden BP, Stinear TP (2025). *Autocycler: long-read consensus assembly for
> bacterial genomes.* **Bioinformatics** 41(9), btaf474.
> <https://doi.org/10.1093/bioinformatics/btaf474>

Data: <https://figshare.unimelb.edu.au/projects/Autocycler/247142>, released **CC BY-NC
4.0**. BacFlux ships none of it; the runs described here were made from a local download.
