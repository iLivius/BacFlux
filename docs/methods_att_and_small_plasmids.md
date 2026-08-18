# Method reference: att-site detection and small-plasmid recovery

Two problems that look unrelated and are not. **Where does a mobile element stop?**
and **why did a small plasmid vanish from the assembly?** Both come down to short
repeated sequences, and both were found the same way — by two concrete failures on
real data.

Every claim carries its source. Two are NEGATIVE results — questions the published
literature simply does not answer — and they are flagged as such, because those are
the ones to be careful with in a paper.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | ***att* site** | the short repeated sequence left at each end of an element when it integrates into a genome. Finding the pair is how an element's true edges are established |
    | **attL / attR** | the left and right copies of that repeat |
    | **ICE** | integrative and conjugative element — a mobile element that sits in the chromosome and can move itself to another cell |
    | **direct repeat** | two copies of a sequence in the same orientation |
    | **HMM** | hidden Markov model, a statistical profile used to recognise a protein family |

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

This matters for writing up: our own measured spurious de-novo boundary rate is
better characterised than anything in the published literature, and can be stated
as such. Measured on 300 randomly placed 15 kb non-ICE spans of a clinical
*K. pneumoniae* chromosome, with the real IS mask applied: **22%** of them
returned a confident de novo "boundary", falling to **16%** once the
repeat-family guard was added, and to **~1%** reaching the AMR table once de novo
repeats stopped being applied at all. Quote **16%** as the current de novo
false-boundary rate — 22% is the *before* figure, and the two are easy to confuse
because both appear in the code comments.

*(The source comment in `att_search.py` names this chromosome "KPNIH1" and
`mobilome_worked_example.md` attributes it to ATCC BAA-2146. Those are two
different genomes — see the note in the worked example — and which one carried
this particular measurement has not been re-established, so no strain is claimed
here. The rate itself is unaffected either way.)*

This is also why de novo repeats are **reported but not applied**: at a
one-in-six error rate a de novo pair is a reasonable lead for a human to follow
and an unacceptable basis for silently redefining an element, since widening an
interval turns every gene inside it into predicted cargo. tRNA-anchored pairs
*are* applied, because they start from a position integrases are known to target.

### What BacFlux does as a result

**Exact maximal repeats between the two flanks, computed in standard-library
Python** — `vmatch -l N` semantics, which is what ICEfinder and ICEfinder2 use.
The search is in `workflow/scripts/80_mobilome/att_search.py`
(`find_maximal_repeats` / `search_maximal_repeat`).

> **Correction, 2026-07-31.** This section previously described a flank-vs-flank
> `blastn -task blastn-short -dust no` call, following DEPhT and DBSCAN-SWA. That
> was the design at the time of writing; it is **not** what shipped. Since commit
> `9c05c3e` the search calls no external program at all. The reason is stated in
> the code: `att_search.py` is stdlib-only and the Snakemake rule that calls it
> has no conda environment of its own, so shelling out to `blastn` would make the
> module depend on whatever happened to be on `PATH`. The BLAST recipe is
> genuinely more sensitive — it absorbs mismatches and gaps, which an exact
> repeat search cannot — so this is a deliberate trade of sensitivity for a
> dependency-free rule, not an equivalence.

> **Correction, 2026-07-31 (second pass).** Three further things in this section
> had gone stale against commit `4a93d89`, and are fixed below rather than
> quietly rewritten:
> 1. the flank window is **50 kb**, not 30 kb — the bullets used 30 kb as their
>    worked example throughout;
> 2. the ranking test was described as "**exactly** one copy in a tRNA". The code
>    tests *at least* one, which is not a detail: the paralogue guard that
>    removes the both-in-a-tRNA case now only fires when the two tRNAs carry the
>    **same amino acid**, so a genuine element sitting between two unlike tRNAs
>    reaches the ranking with both copies anchored;
> 3. "two 30 kb flanks share a 15 bp repeat by luck roughly six times over" was
>    simply wrong arithmetic — the formula the code uses gives ~0.8 at 30 kb and
>    ~2.3 at 50 kb. The conclusion it supported (a flat 15 bp floor is not safe
>    on a large window) survives; the number did not.
>
> The licensing note at the end of this part needed the same treatment, for a
> reason that matters more — see the correction there.

How it works, and where each choice comes from:

- **Exact maximal repeats**: index every *k*-mer of the left flank, find matches
  in the right flank, then extend each match outwards for as long as the two
  sequences agree. One pass returns the longest repeat, so nothing has to "step
  the probe down".
- **A 50 kb flank window each side** (`DEFAULT_FLANK_WINDOW_BP`), raised from
  30 kb. Measured on SPI-7 (`AL513382`): the correct attR sits in a tRNA-Phe
  5,800 bp *outside* a 30 kb window, so the element was reported 50 kb short for
  want of anywhere to look. Widening to 50 kb recovers it, and 80 kb, 120 kb and
  200 kb change no element's answer on the benchmark — so this is where the curve
  flattens, not simply a bigger round number. The cost is more candidate repeats
  to rank, which is why the ranking rule below matters more at 50 kb than at 30.
- **Floor of 15 bp**, matching ICEfinder2 and icefinder-opt; DBSCAN-SWA uses 12,
  so 15 is the more conservative of the two published choices. The floor actually
  applied is the **larger** of 15 and a chance-match threshold computed from the
  two window sizes: at 15 bp, two 50 kb flanks are expected to share **~2.3**
  repeats by luck alone, so a flat 15 would report noise on any large window. At
  the shipped 50 kb window this resolves to **18 bp** — as it also did at 30 kb,
  so widening the window did not move the floor.
  Note this threshold assumes DNA is a random string of four equally likely
  letters, which real chromosomes are not; it is a floor against random
  background, not a filter against repetitive sequence. The 16% figure above is
  what happens when it is asked to do more than that.
- **One search, not two modes.** tRNA proximity is a **ranking term inside the
  single search**, not a separate search: every maximal repeat that brackets the
  machinery is found first, and only then are they ranked. This follows
  ICEfinder's division of labour — the tRNA *locates* a candidate site, the
  repeat search *delimits* it — and DEPhT's treatment of integrase proximity as
  one term among several rather than as a gate.
- **Anchoring outranks length.** Candidates are ordered by whether *at least one*
  copy sits inside an annotated tRNA **first**, and by repeat length only within
  each group. The two are not comparable quantities: length says how unlikely the
  match is by chance, a tRNA 3′ end says the match is where integration actually
  happens. Measured on SPI-7 (`AL513382`), a 51 bp repeat in ordinary sequence
  beat the real 24 bp *att* pair at tRNA-Phe and the element came out 50 kb short;
  under this ordering the tRNA-anchored pair wins and the call lands on the
  curated interval. (A small `TRNA_ANCHOR_BONUS_BP = 10` also exists in the score,
  from the earlier tie-break design.)
- **Both copies inside tRNAs of the *same amino acid* is rejected outright** —
  that is two paralogous tRNA genes, not an integration scar. Two copies in
  tRNAs of *different* amino acids is allowed, and the distinction is not
  academic: ICE*Ec2* (`GU725392`) has its real 22 bp *att* pair in tRNA-Phe at
  one end and tRNA-Ser at the other, and the earlier blanket rule threw that
  boundary away, reporting the element 37 kb short. This is why the ordering
  above tests "at least one copy in a tRNA" rather than "exactly one": by the
  time ranking happens, the paralogue case has already been removed.
- The label written to `boundary_method` is therefore an **outcome**, not a mode:
  `tRNA` when the winning repeat has a copy in a tRNA, `denovo` when it does not,
  `none` when nothing survived.

**No new dependency, and no external program**: the search is pure Python
standard library.

### Why not Vmatch

Spec §3.3 previously rejected Vmatch as "not on bioconda and licence-restricted".
**The bioconda half is wrong** — `bioconda/vmatch 2.3.1` exists for linux-64 and
osx-64. The rejection stands for the other reason: the recipe declares
`license: Unknown / OTHER` and vmatch.de was unreachable, so its terms cannot be
verified — which under the project's §11 rule is still a blocker for an
MIT-licensed workflow. The spec has been corrected to give the right reason.

There is now also a **measured** reason, which matters more than the licence one.
Because the search above already computes vmatch's own semantics, a second
implementation was built on vmatch's actual data structure — prefix-doubling
suffix array, Kasai LCP array, cross-flank MEM enumeration — and run against ours
on every real flank window in the benchmark.

The script is
`BacFlux_v2_validation/phase7_benchmark/verify_att_equivalence.py`
and it re-runs in about ten minutes. It does not generate windows of its own: it
replays the actual caller (`conjscan_to_ice.py`) over the benchmark genomes with
`find_maximal_repeats` wrapped in a recorder, so every window tested is one the
pipeline genuinely searched — after IS masking, at the real element coordinates,
at the real computed floor.

**Result, re-run 2026-07-31 at commit `4a93d89`:** 52 benchmark genomes, 246 att
searches, **180 distinct flank windows, 1,033 repeats, exact set equality on
every window**. No window hit either implementation's internal cap. So the
conclusion holds: **a vmatch-based search would find nothing that is not already
found**, and the licence blocker costs no sensitivity.

> **Correction, 2026-07-31.** The figure previously quoted here — "**35 of 35
> windows, 122 repeats**" — **does not reproduce**, and it should not be cited.
> The equality result reproduces and is now measured over roughly five times as
> many windows, but the counts do not match under any scoping that could be
> reconstructed: not all searches (246), not distinct windows (180), not the
> 18-genome ICE pilot alone (137 searches / 98 windows), and not the subset of
> windows that returned at least one repeat (142 and 62 respectively). Since the
> original script was never kept, what it was run over cannot now be recovered —
> which is exactly the failure mode that prompted retaining this one. Quote the
> numbers above, which are reproducible by running the script.

One tempting follow-up idea is disposed of by the same measurement: shipping both
searches and treating their agreement as a confidence signal would be worthless.
Two implementations of identical semantics are *expected* to agree exactly, so
their agreement carries no information about whether a boundary is real. It tests
the code, once, and that is all it is for.

### Licensing note

ICEfinder2 is **CC BY-NC-SA 4.0**. Its source was read only to establish the
algorithm and its parameters, which spec §11 explicitly permits. **No code was
copied.**

> **Correction, 2026-07-31.** This paragraph used to end "*none of the above
> requires it — the flank-vs-flank BLAST recipe comes from DEPhT and
> DBSCAN-SWA*". That sentence was the licence argument, and it stopped being
> true when the implementation changed: BacFlux no longer uses the DEPhT recipe,
> it computes the same maximal exact repeats ICEfinder2 gets from Vmatch. So the
> separation has to be stated properly rather than by pointing at a different
> tool.

What was taken from ICEfinder2 is **which algorithm to use and with what
parameters** — maximal exact repeats between the flanks, floor 15 bp — and spec
§11 is explicit that design decisions and thresholds are facts, not expression.
The algorithm itself is not ICEfinder2's to license: maximal exact repeat search
is standard published string processing (Vmatch, and REPuter before it), and
`find_maximal_repeats` is an independent stdlib Python implementation written
from that definition. Vmatch's own terms are a separate question and are the
reason it is not installed — see above.

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

**2. Our own read QC.** BacFlux hard-coded `--keep_percent 90 --length_weight 10`;
filtlong's own default for the second is 1. Filtlong scores reads as

```
(Length^lw × MeanQ^mqw)^(1/(lw+mqw)) × WindowQ
```

and `--keep_percent` deletes the bottom of that ranking. **A plasmid cannot
produce reads longer than itself**, so its reads sit at the bottom by
construction, and the more length dominates the score the further down they sit.

Wick's Feb 2026 read-QC benchmark (5 genomes, 11 assemblers) found the same
failure and is worth reading carefully, because it does **not** say what we first
took it to say. His `Filtlong-defaults` arm ran filtlong at its **own defaults**
(`--target_bases` = 100 × genome size, so `length_weight` 1 — the value we now
ship), and that arm *"removed all reads below 10 kbp, essentially erasing these
plasmids"* in the 3 of 5 genomes that had sub-10 kb plasmids, with more
structural errors as a result. So it is not independent confirmation of our
`length_weight 10` finding: it is a warning that filtlong's default length
weighting is already enough to erase small plasmids once you cull hard to a
target. His recommendation is to go **below** the default —
`--length_weight 0 --window_q_weight 0`. Note this bears directly on `nanopore`
mode, which culls with `--target_bases` exactly as his arm did.
<https://rrwick.github.io/2026/02/05/read_qc_testing.html>

### Measured on our own data

*Klebsiella pneumoniae* **TUM24772** — BioProject PRJNA1168299, closed genome
GCA_043950115.1, chromosome CP171785.1. A clinical test isolate, and **not** one
of the BAA-2146 / KPNIH1 mobilome positive controls; this project has conflated
those before, so the distinction is worth stating. Its 5,596 bp Col2 plasmid
(pMTY24772_Col2, CP171790.1) was missing from our assembly while sitting complete
in the Illumina data.

Reads mapping to that plasmid, for all four combinations of the two keys. Counts
are `minimap2 -x map-ont` against `CP171790.1`, keeping reads with an alignment
block of at least 1 kb; on that measure the **raw** read set holds **603**:

| | `length_weight 10` | `length_weight 1` |
|---|---|---|
| **`keep_percent 90`** — the old value | **93** | 413 |
| **`keep_percent 95`** — what ships | 500 | **502** |

**Both keys matter, and `keep_percent` is the one to reach for.** At the 95 that
now ships, changing `length_weight` is worth two reads; at the old 90 it is worth
320. The plasmid needed 90 and 10 acting together to disappear, and either change
on its own brings most of it back.

> **Correction, 2026-08-16.** This section previously showed two arms only and
> concluded the plasmid was lost to `length_weight`. Filling in the missing arms
> of the 2×2 — they were run for this correction; the original experiment had no
> `keep_percent 95` arm at all — reverses that conclusion, and takes two figures
> with it:
>
> - The **602** quoted here and in three other places was measured with **no
>   `--keep_percent` at all** (arm `C_no_keeppct`), which is not a configuration
>   BacFlux can produce — both filtlong rules always pass the flag. The shipped
>   defaults give **502**. The **614** raw figure quoted alongside it does not
>   reproduce either; on the measure above the raw set holds 603.
> - "93 — too few for Flye to assemble it" implied that some other setting did
>   assemble it. **None did.** Re-checked against the six Flye runs still on
>   disk: no arm, raw unfiltered reads included, produced a 5,596 bp contig, and
>   every arm misses that same replicon. (The arms are not otherwise identical —
>   `D_permissive` also fails to close the 87.9 kb replicon — but no arm trades
>   one of those for the Col2 plasmid.)

### So why is `length_weight` still 1?

Because the two keys protect different size classes, which only shows in the
read-length distribution. Reads per band, same sample — the raw ONT set holds
2,954 reads of 1–3 kb and 3,418 of 3–6 kb:

| | 1–3 kb | 3–6 kb |
|---|---|---|
| `keep_percent 90`, `lw 10` | **0** | 318 |
| `keep_percent 95`, `lw 10` | **0** | 2,666 |
| `keep_percent 90`, `lw 1` | 123 | 2,331 |
| `keep_percent 95`, `lw 1` | 605 | 2,935 |

`length_weight 10` empties the 1–3 kb band at **both** `keep_percent` values.
Raising `keep_percent` rescues the 3–6 kb band, which is why it rescued this
5.6 kb plasmid, but only lowering `length_weight` rescues 1–3 kb. Put plainly:
**`keep_percent` protects plasmids of a few kb; `length_weight` protects plasmids
under ~3 kb**, and this one happened to sit in `keep_percent`'s range.

Total bases barely moved across the arms (208 → 227 Mb), so none of this was ever
about depth — only about *which* reads survived. At the old 90/10, 3.1% of kept
reads were ≤6 kb, against 47.0% of the raw set.

**Fix applied:** `parameters.long_read_qc.keep_percent` now defaults to **95**
(was 90) and `length_weight` to **1** (was 10), both configurable. Note that
`length_weight` reaches filtlong in **hybrid mode only** — the `nanopore` rule
never passes `--length_weight`, so in two of the four modes the key does nothing.

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
| **`keep_percent` 90 → 95** (and `length_weight` 10 → 1) | two lines | **done.** Recovers 502 of the plasmid's 603 reads — but *not* the assembled plasmid |
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
On *K. pneumoniae* **ATCC BAA-2146**, `auto` mode dropped a genuine plasmid —
pMYS, `NZ_CP006660.1`, 2,014 bp — because *E. coli* database entries outnumbered
*Klebsiella* ones 58,709 to 16,253: BLAST bestsum follows database composition,
not biology. (This page previously attributed that run to KPNIH1; `NZ_CP006660.1`
is a BAA-2146 replicon, and KPNIH1 is `CP008827.1`.)

**It did not cause the TUM24772 loss** — that plasmid survived decontamination
(verified: present in `contigs_filt.fasta`), and was lost at the read-filtering
step instead, to `keep_percent 90` and `length_weight 10` acting together. Both
mechanisms must be ruled out separately.

**How to check, then fix:** see the annotated `decontamination:` block in
`config/config.yaml`. In short: read
`contig_taxonomy_decisions.tsv` for a discarded plasmid-sized contig whose genus
differs from the sample's; then re-run with `mode: include` naming both genera,
or `mode: off`, or `discard_no_hit: false`. Snakemake redoes only what changed.
