# Method reference: hybrid mode on five closed genomes

Five bacterial isolates with published, high-accuracy reference genomes, run end to end
through `hybrid` mode with every optional layer enabled. The genomes come from the
Autocycler benchmark, where each was sequenced deeply on ONT and Illumina and assembled
to a curated reference — providing an external reference for assessing assembly content
and sequence.

Two things came out of it: a measurement of when short-read correction of long reads
helps and when it does harm, and an end-to-end picture of what the workflow delivers on
genomes with curated reference assemblies.

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
number the workflow produces rather than one taken from the literature. It counts every
ISEScan call, including the one or two look-alike calls described under Mobilome below.

!!! note "The experiment used the published 50× ONT subsets"

    The full ONT sets published with these isolates reach **306×–1,070×**, far deeper
    than routine isolate sequencing. The full-depth runs had short-read guidance on, and
    with guidance on `hybrid` passes Filtlong no `--target_bases` cap: Flye reported
    `Assembled 0 disjointigs` on all five genomes. The runs described here therefore use
    the **50× subsets** published alongside the full sets, which is also the depth the
    original benchmark assembled.

    With guidance off, the default, Filtlong caps the ONT reads at 500 Mb with
    `--target_bases`. With guidance on, subsample unusually deep ONT data before the run.

## Short-read correction of long reads: it can help, and it can lose a replicon

`hybrid` mode can give Filtlong the Illumina reads (`-1`/`-2`), so that it scores each
long read by how well it agrees with them and trims (`--trim`) or splits (`--split`) it
where agreement fails. That is what `parameters.hybrid.short_read_guidance` controls.
Without it, Filtlong sees only the long reads: it filters on length and quality and caps
coverage with `--target_bases`, which is also what `nanopore` mode does.

Three arms were assembled with Flye from the same 50× reads, differing only in the
filtering before them:

| arm | filtering |
|---|---|
| **guided** | Filtlong with `-1`/`-2`, `--trim`, `--split` — `short_read_guidance: true` |
| **unfiltered** | none |
| **unguided** | Filtlong with `--min_length`, `--keep_percent`, `--target_bases` — `short_read_guidance: false` |

| genome | reference | guided | unfiltered | unguided |
|---|--:|--:|--:|--:|
| *Providencia rettgeri* | 2 | 2c, +25 | 2c, +0 | 2c, +0 |
| *Listeria innocua* | 2 | **2c, −1,134** | 3c, +6,020 | 3c, +6,019 |
| *Enterobacter hormaechei* | 4 | **3c, −91,386** | 4c, −1,134 | 4c, +1,418 |
| *Klebsiella pneumoniae* | 5 | 4c, +2,274 | 4c, +2,281 | 4c, −1,234 |
| *Shigella flexneri* | 8 | **40c, −35,210** | 6c, +22,974 | 6c, +21,340 |

The unguided and unfiltered arms produced the same contig counts on all five
genomes. This comparison points to the guidance options as the main source of the
observed differences; it does not establish that length and quality filtering are
harmless on other datasets. `--split` severs a read wherever
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

    The practical reading: short-read guidance warrants particular scrutiny on
    repeat-rich genomes or shallow Illumina data. It is disabled by default. A
    fragmented ONT assembly should prompt comparison of the filtering arms and read
    support; fragmentation alone does not identify its cause.

## What the workflow delivered

Assemblies, against the published references. These runs had short-read guidance
**on**. The table measures the delivered genomes, after polishing, whereas the arm
comparison above measured the Flye assemblies, so total lengths differ from the
`guided` column by up to 40 bp. The mobilome results that follow were produced from
these delivered genomes.

| genome | contigs | assembled | difference |
|---|--:|--:|--:|
| *Providencia rettgeri* | 2 | 4,465,802 | **−4 bp** |
| *Listeria innocua* | 2 | 2,971,407 | −1,138 bp |
| *Klebsiella pneumoniae* | 4 | 5,992,470 | +2,274 bp |
| *Enterobacter hormaechei* | 3 | 5,293,335 | −91,412 bp |
| *Shigella flexneri* | 40 | 4,793,237 | −35,250 bp |

*Providencia* differs in total assembly length by four bases from a reference built
with Trycycler, Medaka, Polypolish and Pypolca. Equal or similar total length does
not establish sequence accuracy. Two failure modes account for the rest, and both are
worth recognising in your own data.

**Small plasmids go missing.** *Klebsiella*'s reference carries a 1,240 bp plasmid and
*Enterobacter*'s an 80,745 bp one; neither survived. This is the documented small-plasmid
loss of long-read assembly, described in
[att sites and small plasmids](methods_att_and_small_plasmids.md).

**Small plasmids also arrive doubled.** In the unguided arm, three of *Shigella*'s six
small plasmids assembled as tandem concatemers — 2 × 6,790, 3 × 3,181 and 4 × 3,835 bp,
matching whole multiples to within 9 bp. A read longer than a small circle wraps it more
than once, and the assembler represents that as a tandem repeat. The replicon then
appears several times its true size and its genes more than once, so gene content and
copy number read from it are wrong. Aligned against the reference, that assembly
recovered three of the eight replicons as single contigs (the chromosome, the 222 kb
plasmid and the 5,057 bp plasmid), three as concatemers, and lost the two smallest
plasmids (2,088 and 1,240 bp). Total assembly length cannot separate these outcomes;
only a per-replicon alignment can.

### Mobilome

Stage `08.mobilome` with all four optional layers enabled:

| genome | AMR genes | IS elements | ICE/IME candidates | named elements | mobility tiers |
|---|--:|--:|--:|--:|---|
| *Listeria innocua* | 2 | 11 | 1 | 0 | all tier 1 |
| *Providencia rettgeri* | 5 | 8 | 1 | 0 | 3 × tier 1, 2 × tier 6 |
| *Shigella flexneri* | 39 | 392 | 0 | 0 | tiers 1, 2, 3, 5 |
| *Enterobacter hormaechei* | 44 | 40 | 4 | 2 | 10 × tier 1, **15 × tier 4**, 19 × tier 6 |
| *Klebsiella pneumoniae* | 55 | 94 | 3 | 2 | 11 × tier 1, 44 × tier 6 |

The IS counts include ISEScan's look-alike calls over *dnaA* and over the chromosome's
helicase loader, which the module removes before anything reads the IS table (see
[False IS21 calls over replication genes](mobilome/index.md#false-is21-calls-over-replication-genes)).
On these genomes that is one call fewer for *Providencia* and two fewer for each of the
other four. No AMR gene's tier and no ICE/IME candidate depends on those calls.

In these runs, *L. innocua* returned two tier-1 intrinsic candidates; this label
does not establish species-wide intrinsic resistance. *Klebsiella* returned 44 genes
with predicted self-transmissible context. *Shigella* had 392 IS calls, roughly four
times the next-highest count. Repeat density is consistent with its assembly difficulty,
but the count alone does not establish the cause of each contig break.

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
a consensus assembly approach such as Autocycler may be appropriate. Its output
can be analysed by BacFlux in `contigs` mode.

> Wick RR, Howden BP, Stinear TP (2025). *Autocycler: long-read consensus assembly for
> bacterial genomes.* **Bioinformatics** 41(9), btaf474.
> <https://doi.org/10.1093/bioinformatics/btaf474>

Data: <https://figshare.unimelb.edu.au/projects/Autocycler/247142>, released **CC BY-NC
4.0**. BacFlux ships none of it; the runs described here were made from a local download.
