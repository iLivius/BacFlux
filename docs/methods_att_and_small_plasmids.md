# Method reference: att-site detection and small-plasmid recovery

*Destined for the MkDocs site, and intended to be quotable in a methods section.
Every claim here carries its source. Two of them are NEGATIVE results — places
where the literature has no answer — and those are flagged explicitly, because
they are the ones worth being careful about in a paper.*

Established 2026-07-28 by reading the primary sources and the tools' own source
code, prompted by two concrete failures on clinical *Klebsiella pneumoniae*
hybrid data (isolates TUM24772 / PRJNA1168299 and K3 / PRJNA1291976).

---

## Part 1 — How att sites are actually detected

### The failure that prompted this

BacFlux's att search took a **fixed 25 bp probe** from the 3′ end of a tRNA and
looked for a second copy bracketing the candidate element. On two clinical
ICE*Kp* elements it found nothing at 25 bp, but candidate pairs appeared at
18 bp.

The reason is measured, not inferred: **the ICE*Kp* direct repeat is 17 bp**
(`CCAGTCAGAGGAGCCAA`), reported by Lam *et al.* 2018.
A 25 bp probe cannot match a 17 bp repeat under any mismatch budget. This was a
structural defect, not a tuning problem.

> Lam MMC *et al.* (2018) *Genetic diversity, mobilisation and spread of the
> yersiniabactin-encoding mobile element ICEKp in Klebsiella pneumoniae
> populations.* Microbial Genomics.
> <https://pmc.ncbi.nlm.nih.gov/articles/PMC6202445/>

### What the reference tools do — none uses a fixed length

| tool | att method | length rule | source |
|---|---|---|---|
| **ICEfinder** (Ou lab, the reference ICE tool) | ARAGORN locates tRNA/tmRNA 3′ termini; **Vmatch** finds the direct repeats marking the tRNA-distal boundary | maximal exact repeats — variable by construction | [ICEberg 2.0 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6323972/) |
| **ICEfinder2** | same, in code | **`vmatch -l 15`** — floor 15 bp, length read from the output | [`script/single.py` L286](https://github.com/EBI-Metagenomics/icefinder2/blob/main/script/single.py) |
| **icefinder-opt** (maintained 2025 fork) | unchanged | still `vmatch -l 15` | [repo](https://github.com/guogenglin/icefinder-opt) |
| **DEPhT** | **BLASTN of the left flank against the right flank**, scoring *all* resulting pairs | whatever the local alignment returns | [DEPhT paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9303363/) |
| **DBSCAN-SWA** | same flank-vs-flank BLAST pattern | explicit **12 bp** minimum, ranked by bitscore | [`bin/dbscan-swa.py`](https://github.com/HIT-ImmunologyLab/DBSCAN-SWA/blob/master/bin/dbscan-swa.py) |
| **Islander / TIGER** | sidesteps repeats entirely — integration splits a tDNA and the island restores it, so it hunts the **displaced tDNA fragment** by sensitive BLASTN | whatever the alignment returns | [Islander paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4383910/) |
| **IntegronFinder** | covariance model (Infernal) for *attC* | bounded variable length, defaults 40–200 bp | [NAR 2016](https://academic.oup.com/nar/article/44/10/4539/2516972) |
| IslandViewer components (IslandPath-DIMOB, SIGI-HMM), alien_hunter, PAIDB | **do not detect att sites at all** — compositional / mobility-gene methods | n/a | [IslandViewer 4](https://academic.oup.com/nar/article/45/W1/W30/3787837) |

IslandViewer integrates Islander *specifically because* its own compositional
predictors have poor boundary accuracy — so they offer no precedent for
boundary calling.

### Why a fixed length cannot work

- Bacterial ICE direct repeats span roughly **10–60 bp**
  ([FEMS Microbiol Rev](https://academic.oup.com/femsre/article/41/4/512/3089980)),
  so any single value sits arbitrarily inside the real distribution.
- Integrase family predicts the range: **tyrosine**-integrase elements at tRNA
  sites (which ICE*Kp* is) sit at **20–45 bp** cores, **serine**-integrase
  elements go as low as **3–12 bp** ([DEPhT](https://pmc.ncbi.nlm.nih.gov/articles/PMC9303363/)).

### ⚠ NEGATIVE RESULT: no precedent for "step the probe down"

The obvious repair — try 25, then 24, 23 … and take the longest that works — has
**no published precedent** in att-site or genomic-island detection. Nobody needs
it: a maximal-repeat search (Vmatch) or a local alignment (BLAST) returns the
longest repeat in a single pass. Exact k-mer stepping is also strictly *weaker*,
because it cannot absorb the indels and mismatches an aligner handles for free,
and it yields no bitscore to weight confidence with.

### ⚠ NEGATIVE RESULT: no published false-positive rate for boundary calling

There is **no published false-positive rate or statistical treatment specific to
att-site detection**. Existing benchmarks measure genomic-island *calling*, not
*boundary* calling. Islander's paper admits only "the few false positives", with
no number.

This matters for writing up: our own measured **22% spurious de-novo boundary
rate** (300 randomly placed non-ICE spans on the KPNIH1 chromosome, with the
real IS mask applied) is better characterised than anything in the published
literature, and can be stated as such.

### What BacFlux does as a result

Flank-vs-flank **BLASTN**, following DEPhT and DBSCAN-SWA:

```
blastn -query <left flank> -subject <right flank> -task blastn-short -dust no
```

- **`-task blastn-short`** — tuned for query lengths under ~50 bp
  ([BLAST+ manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/))
- **`-dust no`** is load-bearing: DUST masking hides low-complexity att cores,
  which are common
- **minimum 15 bp**, matching ICEfinder2 and icefinder-opt; DBSCAN-SWA uses 12
- ranked by **bitscore**, as DEPhT does
- **tRNA proximity is a scoring bonus, not a separate mode** — the same shape as
  DEPhT's integrase-proximity bonus, and as ICEfinder's division of labour
  (tRNA locates, Vmatch delimits)

No new dependency: NCBI BLAST+ is already a BacFlux dependency, is **public
domain**, and carries no commercial restriction.

### Why not Vmatch

Spec §3.3 previously rejected Vmatch as "not on bioconda and licence-restricted".
**The bioconda half is wrong** — `bioconda/vmatch 2.3.1` exists for linux-64 and
osx-64. The rejection stands for the other reason: the recipe declares
`license: Unknown / OTHER` and vmatch.de was unreachable, so its terms cannot be
verified — which under the project's §11 rule is still a blocker for an
MIT-licensed workflow. The spec has been corrected to give the right reason.

### Licensing note

ICEfinder2 is **CC BY-NC-SA 4.0**. Its source was read only to establish the
algorithm and its parameters, which spec §11 explicitly permits. **No code was
copied**, and none of the above requires it — the flank-vs-flank BLAST recipe
comes from DEPhT and DBSCAN-SWA.

---

## Part 2 — Why small plasmids go missing, and what to do

### It is not the assembler

Flye is one of the **better** long-read assemblers for small plasmids, not an
offender. Recovery of plasmids **<10 kb**
([Microbial Genomics 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10272865/)):

| assembler | <10 kb recovery |
|---|---|
| Unicycler (hybrid) | 100% |
| Canu | 100% |
| **Flye `--nano-raw`** | **79%** |
| Flye `--meta` | 73% |
| **Flye `--nano-hq`** | **67%** |
| Miniasm/Minipolish | 64% |
| Raven | 39% |

Wick & Holt independently report NECAT "failed to assemble many plasmids",
NextDenovo "performed poorly on plasmid assembly", and Raven "wasn't good with
small plasmids"
([[Genome Biol]](https://pmc.ncbi.nlm.nih.gov/articles/PMC6966772/)).
**Switching assemblers would mostly make this worse.**

> **Actionable side-finding:** `--nano-hq` is the *worst* of Flye's three modes
> for small plasmids (67% vs 79% for `--nano-raw`), and BacFlux auto-selects the
> mode. Worth revisiting if small replicons matter more than consensus accuracy.
> Note also that Flye's `--plasmids` flag **no longer exists** — added in 2.4,
> removed in 2.9.

### The two causes that do explain it

**1. ONT ligation library prep.** Plasmids <20 kb are under-represented in
ligation read sets by a mean factor of **~4**, and up to **>100-fold** for the
smallest (one 2.4 kb plasmid). Unfragmented circular plasmids never receive a
blunt-end adapter. **Rapid (transposase) prep shows no such bias.**
This is unfixable *in silico*.
([Microbial Genomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC8549360/))

**2. Length-weighted read QC — ours.** Wick's Feb 2026 read-QC benchmark (5
genomes, 11 assemblers) found length-weighted filtlong **wiped out sub-10 kb
plasmids in 3 of the 5 genomes that had them**, and consequently produced *more*
structural errors. <https://rrwick.github.io/2026/02/05/read_qc_testing.html>

BacFlux had `--length_weight 10`, ten times filtlong's default of 1. Filtlong
scores reads as

```
(Length^lw × MeanQ^mqw)^(1/(lw+mqw)) × WindowQ
```

so at `lw=10` the ranking is dominated by length, and `--keep_percent` then
deletes the bottom of it. **A plasmid cannot produce reads longer than itself**,
so its reads sit at the bottom by construction.

### Measured on our own data

TUM24772, whose 5,596 bp Col plasmid was missing from the assembly while sitting
complete in the Illumina data:

| ONT read set | reads mapping to the plasmid |
|---|---|
| raw | **614** |
| `length_weight 10` (old default) | **93** — 85% destroyed |
| `length_weight 1`, no `keep_percent` | **602** — 98% recovered |

Read-length distribution, same sample: 47.0% of raw reads were ≤6 kb; after
filtering, **3.1%**, with the **1–3 kb bin emptied entirely** despite every read
in it being above `min_length`. Total bases barely moved (208 → 227 Mb), so this
was never about depth — only about *which* reads survived.

**Fix applied:** `parameters.long_read_qc.length_weight` now defaults to **1**
and `keep_percent` to **95**, both configurable.

### The doubled-circle artifact

Our 23,928 bp contig is the 11,970 bp Col plasmid **duplicated within one
contig** (2 × 11,970 = 23,940). This is documented: start/end overlap in small
plasmids produces ~200% contiguity, and the Plassembler paper names Flye
specifically as having "multiplicated many small plasmids" on real data.

**dnaapler does not fix this.** dnaapler only *rotates* contigs, which cures a
few duplicated bases at a circular junction, not a whole-plasmid concatemer.
`autocycler trim` does handle it, but cannot distinguish artifact from genuine
duplication and expects to run inside the Autocycler pipeline.

### Options, ranked by effort-to-benefit

| option | effort | verdict |
|---|---|---|
| **`length_weight` 10 → 1** | one line | **done.** Recovers 98% of the plasmid's reads |
| **SPAdes rescue** | ~½ day, no new deps | Published recommendation, not a hack: small plasmids "usually appear as circular contigs" in a short-read graph ([Wick/Judd/Holt](https://pmc.ncbi.nlm.nih.gov/articles/PMC9980784/)). We already write `assembly_graph_with_scaffolds.gfa` |
| **Plassembler** | ~½ day + a database | **MIT, bioconda 1.8.3**, purpose-built. Pools reads that do *not* map to the Flye contigs and hybrid-assembles them with Unicycler — exactly where a Flye-absent plasmid lands. Can **reuse our existing Flye assembly** (`--flye_directory`), and writes empty outputs when it finds nothing, so the DAG never breaks. Cost: PLSDB is mandatory (`-d`, no skip flag) and states **no licence at all** — handle exactly as TnCentral/ICEberg per spec §5.5 |
| **Hybracter** | large | MIT, bioconda, wraps Plassembler — but replaces the *entire* assembly stage and nests Snakemake inside Snakemake, with two schedulers competing for cores |
| **Autocycler** | large | **Would probably not have helped.** Aimed at chromosome consensus accuracy; its own paper reports the smallest plasmid tested (2.5 kb) was occasionally missed and needed manual curation. Its clustering discards any cluster seen in too few input assemblies — precisely the small-plasmid case. Its own docs say: *if small plasmids matter, add Plassembler* |
| **Trycycler** | n/a | Unusable unattended; superseded by Autocycler |
| **Union of multiple long-read assemblers** | medium | Would probably not have helped — the failure is shared across assemblers and driven by the read set |

### ⚠ NEGATIVE RESULT: newest assemblers unbenchmarked here

No 2025/2026 head-to-head benchmark scores small-plasmid recovery numerically
for Myloasm, metaMDBG, LJA or hifiasm. Recent accuracy benchmarks exist; none
report per-assembler small-plasmid recovery rates.

---

## Part 3 — The compounding trap in our own pipeline

Independent of everything above, BacFlux has a second route to losing a plasmid,
**upstream of the assembler**:

```
decontamination selects contigs by genus
        ↓
only reads mapping to SELECTED contigs become SEL_R1/SEL_R2
        ↓
filtlong uses SEL_R1/SEL_R2 as its short-read reference
        ↓
ONT reads for unselected sequence score as low quality and are discarded
        ↓
Flye never sees the plasmid
```

This bites when a plasmid's best BLAST hit is a **different genus from the
host** — which is not exotic, since plasmids cross genus boundaries constantly.
On KPNIH1, `auto` mode dropped a genuine plasmid because *E. coli* database
entries outnumbered *Klebsiella* ones 58,709 to 16,253: BLAST bestsum follows
database composition, not biology.

**It did not cause the TUM24772 loss** — that plasmid survived decontamination
(verified: present in `contigs_filt.fasta`) and was lost purely to
`length_weight`. Both mechanisms must be ruled out separately.

**How to check, then fix:** see the annotated `decontamination:` block in
`config/config_v2.yaml`. In short: read
`contig_taxonomy_decisions.tsv` for a discarded plasmid-sized contig whose genus
differs from the sample's; then re-run with `mode: include` naming both genera,
or `mode: off`, or `discard_no_hit: false`. Snakemake redoes only what changed.
