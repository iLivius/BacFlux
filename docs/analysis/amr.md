# Antimicrobial resistance

BacFlux looks for resistance genes three times, with three methods that fail in
different ways. Each one sees things the others cannot.

| Leg | Input | Database | Present when | Modes |
|---|---|---|---|---|
| **ABRicate** | the finished contigs | eight curated sets, one run each | ≥80% identity over ≥70% of the reference gene | all four |
| **BBMap → CARD** | the trimmed read pairs | CARD protein homolog model | ≥70% of the reference gene's **length** covered by reads | `illumina`, `hybrid` |
| **AMRFinderPlus** | Bakta's proteins + the contigs | NCBI's AMRFinderPlus database | AMRFinderPlus's own per-gene curated cutoffs | all four, *optional* |

The first two are stage `05.amr` and always run. The third runs inside the
mobilome module, so it appears only when `mobilome.run: true`, and its output
lands in `08.mobilome` rather than `05.amr` — see
[the mobilome module](../mobilome/index.md).

**Why three.** Assembly collapses repeats, so a gene present in several copies —
or one sitting on a repeat-rich mobile element that broke the assembly — can be
under-represented or missing in the contigs while sitting plainly in the reads.
Mapping side-steps the assembly entirely. ABRicate gives breadth across eight
independently curated gene sets and is the leg carrying the EFSA reporting
thresholds. AMRFinderPlus gives each call a coordinate, a method label and a
drug class, which is what the mobility analysis needs to intersect AMR genes with
mobile elements, and it is the only leg that can report resistance-conferring
**point mutations**.

None of the three replaces another, and their counts are not meant to agree.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **CARD** | Comprehensive Antibiotic Resistance Database |
    | **EFSA** | European Food Safety Authority |
    | **GTDB** | Genome Taxonomy Database, the reference taxonomy GTDB-Tk classifies against |

## ABRicate on the contigs

*Rules `amr_contigs`, `AMR_summary`. ABRicate v1.2.0.*

ABRicate BLASTs the delivered genome (`contigs_final.fasta`) against one curated
resistance or virulence database and reports every hit clearing the cutoffs. It
consumes nothing but the finished assembly, so it behaves identically in all four
modes, including [`contigs`](../modes/contigs.md).

The same contigs are screened against each database separately and the eight
tables are kept side by side rather than merged:

```
argannot  card  ecoh  ecoli_vf  megares  ncbi  resfinder  vfdb
```

Different databases cover different gene sets, and merging them would hide which
resource made a call. `abricate --summary` then collapses one sample's eight
tables into a presence/absence matrix, gene by database — the readable overview of
the whole screen for that isolate. Its column order follows the database order
above.

```bash
abricate --db card contigs_final.fasta --minid 80 --mincov 70 --nopath --quiet
```

`--minid 80` and `--mincov 70` are the EFSA reporting thresholds: a hit must show
at least 80% identity over at least 70% of the reference gene's length. They
belong to **this leg only**. The CARD read leg below gates on reference-length
coverage instead, and AMRFinderPlus applies its own curated per-gene cutoffs; a
blanket 80/70 across all three would fight that curation.

!!! warning "A hit just below 80% is not an absent gene"

    The 80% cutoff is brittle in exactly the place it is most often quoted. EFSA's
    own *Bacillus* catalogue (EFSA 2024b) measured *rphB* in *B. subtilis* at
    **85.3%** of genomes when queried as nucleotide and **1.3%** when queried as
    protein — because most genomes sit at 79–79.9% identity, just under the line.
    *rphC* and *satA* behave the same way. The same sequence reads as
    "near-universal" or "essentially absent" depending on which side of the cutoff
    it lands.

    Read a borderline call from the identity value in the table, not from
    presence or absence. This is also why AMRFinderPlus's per-gene cutoffs are
    left alone rather than overridden with 80/70.

## CARD read mapping

*Rules `download_amr_db` (or `download_amr_db_local`), `map_amr_db`,
`card_mapping_report`. BBMap v39.33, tested against CARD v4.0.1.*

Short-read modes only — `illumina` and `hybrid`. The gate is a capability flag,
"does this mode produce short reads", not a mode name, so `nanopore` and
`contigs` get no CARD leg and never build the BBMap environment.

The fastp-trimmed pairs are mapped straight onto CARD's
`nucleotide_fasta_protein_homolog_model.fasta` — the CARD model that holds
acquired resistance genes, as opposed to its mutation-based models. BBMap's
`covstats` output then says, for every reference sequence, how much of it the
reads actually covered. A gene is called present when **at least 70% of its
length** is covered.

Length coverage, not read count: a short conserved domain can attract many
reads without the gene being present.

Set `links.card_link` to CARD's `broadstreet-*.tar.bz2` archive, or point
`directories.card_db` at an extracted copy you already hold to skip the download
— see [reference databases](../getting-started/databases.md). BBMap is capped at
24 threads and asks the JVM for `min(RAM, 32)` GB.

### What the identity filter actually does

!!! warning "The 0.99 identity claim was wrong. Do not quote it."

    Until v2.0.0 this rule passed BBMap `idfilter=0.99`, and the code and the
    README both described the leg as screening reads at 99% identity. **It never
    did.** `idfilter` does not filter the primary alignment of a properly-paired
    read: in BBMap 39.33's own source,
    `align2/AbstractMapThread.java::processIDFilter` clears the mapping only
    `if(!r.paired() && identity < IDFILTER)`, and the loop that filters the
    remaining sites runs `for(int i=sites.size()-1; i>0; i--)` — it stops before
    index 0, so the top site is never reached.

    Measured on real BacFlux reads (400k pairs) under the old setting: **1,240
    alignments retained, 1,238 of them below 99% identity, the lowest at 48.59%**.
    The same reads run as single-end kept 1. So every CARD result this workflow
    has ever produced was screened at BBMap's default of **0.76**.

    The rule now passes `minid=0.76`, which does apply. That is not a tightening
    and it changes no output — it states what has always been running, using the
    flag that enforces it.

The identity floor was deliberately **not** raised to 0.99. At a real 0.99 the
same reads keep 2 alignments instead of 1,240, which would leave almost nothing
in the one AMR leg that is immune to assembly collapse.

The leg is a sensitive screen whose specificity comes from the ≥70%
reference-length gate, not from per-read identity. A methods section can say
"reads were mapped with BBMap at its default minimum identity (0.76) and a CARD
sequence was reported when reads covered ≥70% of its length". It cannot claim a
99% identity screen.

Two further BBMap settings shape the numbers. `ambiguous=best` gives a read that
maps equally well to several CARD entries one best site, and `secondary=f`
suppresses secondary alignments — so near-identical alleles of the same gene
family share the reads between them rather than each appearing fully covered.

### The report and its five categories

A raw CARD hit list overstates a resistome, and on an environmental
Gram-negative it overstates it badly. CARD's protein homolog model is not a list
of acquired resistance genes: it also holds the subunits of multi-component
efflux pumps, the transcriptional regulators of those pumps, and entries where
resistance comes from the gene being **absent**. A single RND efflux system
contributes an inner-membrane transporter, a periplasmic adaptor and an
outer-membrane channel — three separate CARD entries, all core chromosomal
machinery, all counted as "AMR genes found" by a flat total.

CARD already classifies every entry in `aro_index.tsv`. `card_mapping_report.py`
joins that classification onto the hits, so the inflation is visible instead of
silent. The category is a **label, never a filter**: coverage decides what reaches
the table, and nothing is then dropped for what kind of entry it turned out to be.

| `category` | What the entry is | Examples |
|---|---|---|
| `resistance_determinant` | Inactivating enzymes, target protection, target alteration, target replacement — what a reader normally means by "an AMR gene" | CTX-M-15, OXA-58 |
| `efflux_other` | Efflux, but not one of the multi-subunit families below: a single-protein pump, or a family the list does not name | — |
| `efflux_component` | One subunit of a multi-subunit efflux system (RND, ABC, MFS, SMR, MATE, outer-membrane porin). One system contributes several rows | MexF, TriC |
| `regulator` | Regulates a resistance system rather than conferring resistance itself | CpxR, MexR |
| `presence_indicates_susceptibility` | CARD's mechanism for the entry is "resistance by absence". Detecting it **by presence** indicates the opposite of resistance | mgrB, OmpK35, OmpK36, OprD, LamB, carO |

The last category is a sign error, not a borderline call, which is why it gets
its own label and a `note` column spelling out what the presence means.

The `regulator` label uses the long-standing convention that a name ending in R
is the regulator of the operon named without it (MexR regulates *mexAB*). The
convention is not a guarantee, so it is applied **only** where CARD's own
mechanism for the entry is efflux — CARD's regulator entries are almost all
efflux-pump regulators. Without that guard the suffix rule misfires on genuine
determinants that happen to end in R: measured on a *Bacillus* isolate from this
project's own test set, *vmlR* — an ABC-F ribosomal protection protein, CARD
mechanism "antibiotic target protection", a real resistance gene — was demoted to
`regulator`.

Three rows from a real *Pseudomonas* isolate, in the order the report puts them:

| `aro_accession` | `aro_name` | `covered_percent` | `category` | `resistance_mechanism` | `reference_organism` |
|---|---|--:|---|---|---|
| ARO:3000804 | MexF | 100.00 | `efflux_component` | antibiotic efflux | Pseudomonas aeruginosa PAO1 |
| ARO:3003681 | TriC | 95.77 | `efflux_component` | antibiotic efflux | Pseudomonas aeruginosa PAO1 |
| ARO:3000829 | CpxR | 99.00 | `regulator` | antibiotic efflux | Escherichia coli |

Rows are sorted so the ones to look at first come first — determinants above
efflux machinery — then by coverage, which is why CpxR at 99% sits below TriC at
95.77%. The rule's log prints the per-category counts, so the shape of the table
is visible without opening it:

```text
Sample <name>: <total> CARD sequences at >= 70% covered length.
  resistance_determinant: <n>
  efflux_other: <n>
  efflux_component: <n>
  regulator: <n>
  presence_indicates_susceptibility: <n>
```

Only the categories actually present are printed, and a
`presence_indicates_susceptibility` count is followed by a line saying in words
that those rows are not evidence of resistance. A bare total is the number this
report exists to stop people quoting.

`reference_organism` is the organism CARD's **reference sequence** came from,
never a statement about your sample. In the rows above, CARD attributes its CpxR
sequence to *Escherichia coli* in a report on a *Pseudomonas* isolate: that is
CARD's provenance, not a second organism in the sample. The full column list —
including `amr_gene_family` and `drug_class` — is in
[output files](../reference/output.md).

### Three files, one to read

| File | What it is |
|---|---|
| `{sample}_CARD_report.tsv` | **The file to open.** Coverage joined to CARD's classification, categorised and sorted. |
| `{sample}_covstats.tsv` | BBMap's raw per-reference coverage table, re-sorted by descending covered percent. |
| `{sample}_AMR_legend.tsv` | The `aro_index.tsv` rows for every feature covered ≥70%, kept as raw evidence. |

An isolate carrying nothing above 70% is a normal result, and
`map_amr_db` runs without `set -e` on purpose so that it does not become a
pipeline failure. The report script draws the opposite line: a `covstats` file
with no `#ID` header means BBMap did not finish writing it, and the script stops
rather than turning an empty table into a clean report saying the genome has no
AMR genes.

## AMRFinderPlus, inside the mobilome module

*Rules `amrfinder_organism`, `amrfinderplus` — **optional**, gated on
`mobilome.run`. AMRFinderPlus 4.2.7, shipped inside the Bakta environment.*

Bakta already runs AMRFinderPlus internally but surfaces only the gene name and
product. Running it directly costs no new environment and no new database — the
tool is in the Bakta conda environment and its database is inside the Bakta
database — and the full report adds five things the mobility analysis needs:

1. the `Method` column (`EXACTX`, `BLASTX`, `PARTIALX`, `HMM`, `POINTX`…), which
   is a confidence tier straight from the tool;
2. the `PARTIAL*` methods, which flag a hit truncated at a contig end — the
   honest signal for a fragmented assembly;
3. element type and subtype (AMR / STRESS / VIRULENCE);
4. drug class and subclass;
5. resistance-conferring **point mutations**, when the organism is known.

Point mutations are the one thing neither of the other legs can reach. ABRicate
asks "is this gene present?" and answers by BLAST — exactly right for an acquired
gene, and structurally blind to a mutated one. Query a susceptible genome for
*gyrA* and you get a hit at ~99.9% identity; query a resistant one and you get a
hit at ~99.9% identity. The substitution that decides the phenotype is averaged
away inside the percentage, and no threshold separates them, because the
difference is one residue at one coordinate. AMRFinderPlus with `--organism`
loads a curated list of known resistance-conferring substitutions for that
species and inspects the coordinate itself.

### `--organism`, and why it is usually absent

`--organism` only engages for the organisms AMRFinderPlus curates mutations for —
**31** in version 4.2.7, overwhelmingly clinical. `amrfinder_organism` translates
this sample's GTDB-Tk species call into one of those names, or into nothing, and
writes an audit file recording the decision and its reason either way.

The translation is not a string match. GTDB splits genera and marks the split-off
lineages with a suffix (`Pseudomonas_E`, `Klebsiella_A`), and gives unnamed
genomes a placeholder epithet (`sp024807945`). A naive match would turn
`s__Pseudomonas_E sp010095445` into `--organism Pseudomonas_aeruginosa` and score
that genome against *P. aeruginosa*'s curated mutation list — worse than making no
call at all.

No match is the normal case for environmental isolates. AMRFinderPlus curates
around thirty organisms, nearly all of them clinical, so a soil or plant-associated
*Bacillus*, *Paenibacillus* or *Arthrobacter* matches none of them and the audit file
records the same reason each time: *no curated organism for this taxon*.

### GTDB to AMRFinderPlus organism

`--organism` only works for the names AMRFinderPlus curates, and those are **NCBI**
names. BacFlux calls taxonomy with GTDB-Tk, which frequently disagrees: GTDB splits
genera (`Campylobacter_D coli`, `Enterococcus_B faecium`) and suffixes epithets
(`Haemophilus influenzae_E`). Handed one of those, AMRFinderPlus matches nothing.

Two committed lookup tables bridge the gap, both generated per GTDB release by
`generate_gtdb_organism_table.py` from GTDB's full metadata — never hand-edited:

- **`gtdb_organism_genus_rules.tsv`** — may this organism be matched from the genus
  alone? Measured, not assumed: *Salmonella* scores 100%, *Escherichia* 98%, while
  *Campylobacter* scores **0**, because GTDB places them all in `Campylobacter_D`.
- **`gtdb_organism_equivalences.tsv`** — the per-species overrides, each row carrying
  the percentage of genomes on which GTDB and NCBI agree and how many genomes that
  rests on.

!!! note "Rebuilding the tables is a manual step"

    `generate_gtdb_organism_table.py` is not invoked by any rule. It is run by hand,
    and both tables are committed to the repository; at run time
    `gtdb_amrfinder_organism.py` only reads them from disk, and needs no network
    access. They need rebuilding when GTDB or AMRFinderPlus changes, not per sample,
    and nothing in the workflow warns that they have gone stale. A stale table maps
    fewer species; an unmapped species silently loses its `--organism` flag, which
    means no point-mutation screening for that isolate. `check_gtdb_organism_table.py`
    validates a rebuild.

At run time `amrfinder_organism` looks this sample's classification up in those tables
and writes either an organism name or nothing. Nothing is the safe default and the
usual outcome: a wrong `--organism` would produce confidently wrong point-mutation
calls, so the rule refuses to guess — including when a hybrid sample's two assemblies
disagree. Every decision, and its reason, lands in
`{sample}_amrfinder_organism_audit.tsv`.

!!! warning "An empty mutations file means 'not assessed', never 'no mutations found'"

    When no organism matches, `{sample}_amrfinderplus_mutations.tsv` is written
    with a single line saying point-mutation screening was not available. Writing
    it either way keeps the output set the same whether or not an organism
    matched, and keeps the refusal and its reason in the output instead of leaving
    a gap that reads like a negative result.

## What BacFlux does not decide: intrinsic or acquired

None of the three legs classifies a gene as intrinsic or acquired, and the
mobilome module does not either. EFSA's definitions are population-level:

> **Intrinsic AMR gene** — "Gene inherent to strains of a bacterial species … An
> AMR gene is considered 'intrinsic' when it is shared by the vast majority of
> wild type strains of the same species (or subspecies) and is restricted to those
> located on the chromosome."
>
> **Acquired AMR gene** — "A resistance gene novel for the strain under
> assessment, acquired through horizontal transfer … Acquired AMR genes could be
> integrated in the bacterial chromosome **or** harboured on a separate genetic
> element."
>
> — EFSA Scientific Committee (2025), Glossary

**You cannot establish "intrinsic" from one genome** — the claim is about the
species, and settling it means comparing gene presence across many genomes of
that taxon. BacFlux performs no population comparison of any kind; every
analysis here is per-genome. And **chromosomal does not mean
intrinsic**: an acquired gene may well sit on the chromosome. That is why the
mobility ladder's tier 1 is labelled `intrinsic_candidate` and not `intrinsic`,
and why the CARD report's mechanism category is a hint rather than a
determination — efflux is enriched for chromosomal core systems but includes
plenty of mobile ones, and "antibiotic inactivation" spans both the intrinsic
chromosomal AmpC and the acquired CTX-M.

**BacFlux produces supporting evidence for the intrinsic/acquired judgement. It
does not produce the judgement.** The confirmatory work is a species-wide
distribution analysis plus phenotypic testing.

Where the three legs do line up with the guidance: EFSA (2024a) requires the AMR
search to run against at least two maintained, curated databases at ≥80% identity
and ≥70% length coverage. BacFlux runs ten database queries across the three legs
— eight through ABRicate, plus CARD by read mapping and AMRFinderPlus — with those
thresholds applied on the ABRicate leg.

For what happens to an AMR call once the mobilome module has it, see
[the mobility ladder](../mobilome/mobility-ladder.md).

## Where the files land

| Path under `output_dir` | Modes |
|---|---|
| `05.amr/abricate/{sample}/{db}.tsv` (×8) | all four |
| `05.amr/abricate/{sample}/AMR_summary.txt` | all four |
| `05.amr/mapping/{sample}/{sample}_CARD_report.tsv` | `illumina`, `hybrid` |
| `05.amr/mapping/{sample}/{sample}_covstats.tsv` | `illumina`, `hybrid` |
| `05.amr/mapping/{sample}/{sample}_AMR_legend.tsv` | `illumina`, `hybrid` |
| `08.mobilome/{sample}/{sample}_amrfinderplus.tsv` | all four, `mobilome.run: true` |
| `08.mobilome/{sample}/{sample}_amrfinderplus_mutations.tsv` | all four, `mobilome.run: true` |
| `08.mobilome/{sample}/{sample}_amrfinder_organism_audit.tsv` | all four, `mobilome.run: true` |

## References

- EFSA Scientific Committee (2025). *Guidance on the characterisation of
  microorganisms in support of the risk assessment of products used in the food
  chain.* EFSA Journal 23(11):e9705. <https://doi.org/10.2903/j.efsa.2025.9705>
  — source of the definitions quoted above (Glossary) and the AMR decision tree
  (§3.2.1).
- EFSA BIOHAZ Panel (2023). *Statement on how to interpret the QPS qualification
  on 'acquired antimicrobial resistance genes'.* EFSA Journal 21(10):8323.
  <https://doi.org/10.2903/j.efsa.2023.8323>
- EFSA (2024a). *Statement on the requirements for whole genome sequence analysis
  of microorganisms intentionally used in the food chain.* EFSA Journal
  22(8):e8912. <https://doi.org/10.2903/j.efsa.2024.8912> — source of the
  ≥2-databases requirement and the 80%/70% thresholds.
- EFSA (2024b). *Catalogue of antimicrobial resistance genes in species of
  Bacillus used to produce food enzymes and feed additives.* EFSA Supporting
  Publication 2024:EN-8931. <https://doi.org/10.2903/sp.efsa.2024.EN-8931> — a
  technical report, not guidance; source of the *rphB* / *rphC* / *satA* figures
  above.
- Alcock, B. P. *et al.* (2023). CARD 2023: expanded curation, support for machine
  learning, and resistome prediction at the Comprehensive Antibiotic Resistance
  Database. *Nucleic Acids Research* 51(D1):D690–D699.
  <https://doi.org/10.1093/nar/gkac920>
- Feldgarden, M. *et al.* (2021). AMRFinderPlus and the Reference Gene Catalog.
  *Scientific Reports* 11:12728. <https://doi.org/10.1038/s41598-021-91456-0>
- Seemann, T. ABRicate. <https://github.com/tseemann/abricate>
- Bushnell, B. BBMap. <https://sourceforge.net/projects/bbmap/>

The claim that a BLAST-style presence search is *structurally* unable to resolve
point mutations is inference from what the tools compute, not a statement any of
the sources above makes. The `idfilter` behaviour was read from BBMap's own source
and confirmed by mapping a test isolate's reads with the filter set both ways.
