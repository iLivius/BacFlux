# Draft assemblies

Most BacFlux runs are short-read, and every validation this module originally had was
done on closed genomes. This page is what happens in between.

**The one-sentence version:** on a fragmented assembly the module keeps telling you *what
kind of element* it found and demotes honestly when the evidence breaks, but it stops
telling you *how big* the element is. Read `mge_class` and `confidence` alongside
`spans_contigs` and `at_contig_boundary`, read `boundary_method` before any coordinate,
and expect more low-confidence calls than a closed genome would give.

## How this was measured

40 benchmark genomes — 30 carrying curated ICE or IME elements, 12 negative controls —
were cut into contigs at three contiguities, **re-annotated with Bakta** and re-run
through the whole chain with the same commands and parameters as the closed runs. 120
assemblies, 191 draft calls.

| Arm | Target N50 | Achieved (median) | Contigs (median) | Curated elements still on one contig |
|---|---|---|---|---|
| good | 150 kb | 135 kb | 32 | 21/30 |
| typical | 50 kb | 48 kb | 82 | 14/30 |
| poor | 20 kb | 20 kb | 180 | 8/30 |

Breaks were placed at IS elements first, then rRNA operons, then at random — because
that is what breaks a real assembly — and sequence was deleted at each break. Contigs
were length-sorted and about half reverse-complemented, so cross-contig clustering could
not be flattered by contigs happening to sit in genome order.

**For calibration:** a survey of 28 real BacFlux Illumina assemblies gives a median N50
of **307 kb**, with 26 of the 28 at or above 150 kb. The sweep is deliberately
pessimistic — a failure at "good" is a serious result, a failure at "poor" may simply be
the assembly.

## What degrades, and by how much

| | closed | good (150 kb) | typical (50 kb) | poor (20 kb) |
|---|--:|--:|--:|--:|
| ICE detected (of 18 curated) | 15 | **15** | **15** | 12 |
| IME detected (of 12 curated) | 5 | 4 | 3 | 3 |
| median called ÷ **true** length | 0.55 | 0.45 | 0.34 | **0.22** |
| median called ÷ **largest surviving piece** | 0.55 | 0.57 | 0.77 | 0.72 |
| median call precision (overlap ÷ call length) | 1.0 | 1.0 | 1.0 | 1.0 |
| `spans_contigs` = TRUE | 0 | 13 | 27 | 32 |
| `at_contig_boundary` = TRUE | 3 | 7 | 21 | 23 |
| confidence high | 24% | 19% | 10% | 3% |
| confidence low | 10% | 34% | 63% | 65% |
| `boundary_method` = `none` | 44% | 77% | 82% | 90% |
| calls on the 12 negative controls | 5 | 4 | 3 | 2 |

Three things there matter more than the rest.

**Detection is flat down to 50 kb N50.** ICE recall is *identical* to the closed genomes
at both 150 kb and 50 kb, and only breaks at 20 kb. The IME arm is weak everywhere,
closed genomes included: most IME misses are elements below the module's size floor
(522 bp, 1.2 kb, 5.1 kb) that were already missed before any cutting, so fragmentation
costs 2 IMEs, not 9.

**The length collapses, but the caller is not what lost it.** Against the *true* element
the median recovered length falls from 0.55 to 0.22. Against the **largest surviving
contig piece** — the ceiling any per-contig caller can reach — it is 0.55 → 0.57 → 0.77 →
0.72, flat or better. The caller keeps recovering the same share of what is still
visible. The missing bases are the assembly's.

**Nothing is invented.** Across all 191 draft calls: no fabricated *att* sequence (every
pair reported on a draft occurs at least twice in the corresponding closed genome), no
class promoted at high confidence, no new false positive on a negative control, and not
one violation of the caller's own confidence rules in 163 checks.

## Where the line is

Not between "good" and "poor" — between **classification** and **extent**.

| N50 | What to do with the calls |
|---|---|
| **≥ 150 kb** | usable as-is. ICE detection identical to closed genomes; every high-confidence call was corroborated by the closed run |
| **≈ 50 kb** | the **class** is still trustworthy (ICE 15/18, and 91% of calls land on some curated entry for these accessions). The **extent** is not: median 0.34× the true length, 82% of calls with `boundary_method=none`. *Read the class, ignore the length* |
| **≈ 20 kb** | detection itself starts to go (12/18), 65% of calls are low confidence, and 59% of MacSyFinder's systems join genes from different contigs. Not *wrong*, but it has stopped saying much |

Each sample's own contig count and N50 are written as the **first row** of its
`{sample}_ice_discarded.tsv` (`action = assembly_qc`), with one sentence saying which
band it falls in. Read that before the calls.

## Which columns to read

These are the columns of `{sample}_ice_candidates.tsv`, the element table. The AMR
table beside it repeats the two contig flags under slightly different names —
`is_at_contig_boundary` for `at_contig_boundary`, and `spans_contigs` unchanged — so
check the header of whichever file you have open
(`head -1 file | tr '\t' '\n' | nl`).

**These survive fragmentation and never over-claim:**

| Column | Why |
|---|---|
| `mge_class` | ICE / IME / AICE / passive island / conjugative region. Demotes honestly when evidence is lost — never promotes |
| `mobility` | the sentence that goes with the class. Always "**predicted** self-transmissible" |
| `confidence` | every high call in the benchmark satisfied the documented rule and was corroborated by the closed genome |
| `spans_contigs` | TRUE = the evidence is on more than one contig. **Absolute cap at low** — 72 TRUE calls in the sweep, 72 at low, no exceptions |
| `at_contig_boundary` | TRUE = the element runs into a contig end. Also an absolute cap: 54 TRUE, 54 at low |
| `boundary_method` | **read before any coordinate.** `tRNA` = real boundaries; `denovo` = a lead, never acted on; `none` = the interval is the machinery span only |
| `machinery_intact`, `n_anchor_classes` | what the class actually rests on |
| `mpf_from_other_contig` | TRUE = this is an ICE only because the mating-pair genes are on a *different* contig |

**Treat as a floor, never a measurement, whenever `boundary_method = none`:** `start`,
`end`, `length_bp`. At the typical arm, calls on elements fragmentation had genuinely
broken reported a **median 19% of the element**.

**Do not do on a draft:** compare element lengths between samples, or do arithmetic on
`length_bp` — total mobilome burden, percent of genome mobile, and so on. Those numbers
are governed by your assembly, not by the biology.

## What a low-confidence call means

**That the assembly could not support a stronger claim — not that the call is wrong.** Low
confidence is set by exactly the guards above, and on a fragmented assembly those
conditions are the normal state of affairs, which is why low-confidence calls go from 10%
on closed genomes to 63% at 50 kb N50.

Evidence that the low tier is not simply noise: on the 12 negative controls the number of
calls went **down** under fragmentation (5 → 4 → 3 → 2), none at high confidence, and not
one at a locus the closed genome had not already called. Contig ends did not manufacture
false positives; they removed evidence.

So a low-confidence call on a draft is a lead worth following — look at `spans_contigs`,
the missing components and the audit file to see *what* was lost — not a result to report,
and not one to discard silently.

## The direction that stays safe, and the one exemption

Per call, back-translated onto the closed genome's coordinates:

| | good | typical | poor |
|---|--:|--:|--:|
| **demoted** vs the closed call | 8 | 26 | 34 |
| … ICE → IME (*tra* operon lost) | 1 | 5 | 5 |
| … ICE → conjugative region (integrase lost) | 7 | 16 | 20 |
| **promoted** vs the closed call | 0 | 1 | 0 |
| … at high confidence | 0 | **0** | 0 |
| calls created where the closed genome had nothing | 1 | 1 | 1 |
| … at high confidence | 0 | **0** | 0 |

The dangerous direction is essentially empty, and the single promotion carries
`spans_contigs=TRUE` and low confidence — the guard caught it. Total ICE-class calls fall
37 → 32 → 24 → 16: fragmentation makes the module *more* conservative about tier 6, not
less.

!!! warning "`mpf_from_other_contig = TRUE` — the row that reads as self-contradictory"

    On a draft the *tra* operon routinely lands on a different contig from the relaxase,
    so a cluster with no mating-pair gene of its own is still called an ICE when the typed
    system's other hits are elsewhere. Without that exemption real ICEs get demoted the
    moment an assembly breaks. The cost is a row saying `mobility = predicted
    self-transmissible` next to `missing_components = mating-pair apparatus`.

    The flag fired on 9 of 191 draft calls and on **0** of 68 closed-genome calls. All 9
    also carry `spans_contigs = TRUE`, so all 9 are already at low confidence. Read such a
    row as *conjugation machinery is present in this genome and may belong to this
    element* — the evidence for the mating bridge is not on the contig you are looking at,
    and cannot be checked there.

## One worked call, closed and broken

ICE*Msp*.M1D, a 199,376 bp chromosomal ICE in a *Mesorhizobium* sp. (`CP034444.1`).
Values taken verbatim from `{sample}_ice_candidates.tsv` in each run.

| | closed | good | typical | poor |
|---|---|---|---|---|
| element on how many contigs | 1 | 3 | 9 | 13 |
| largest surviving piece | 100% | 79% | 60% | 14% |
| `mge_class` | **ice** | conjugative_region | conjugative_region | conjugative_region |
| `length_bp` | 39,437 | 16,002 | 16,002 | 13,132 |
| `has_integrase` | TRUE | **FALSE** | **FALSE** | FALSE |
| `boundary_method` | denovo | **none** | **none** | **none** |
| `spans_contigs` | FALSE | FALSE | FALSE | **TRUE** |
| `confidence` | **high** | **high** | **high** | **low** |

The class demotes for a stated reason: the integrase is lost to a contig break, and
without one the caller refuses to say "ICE". It does not guess.

**The uncomfortable part** is the *good* column. That call is high confidence with both
warning flags FALSE, and both flags are literally true — the machinery cluster really is
entirely on one contig, well away from its ends — while the element continues on contigs
the caller never looked at. `confidence` answers *is this conjugation machinery, and of
what class*. It does **not** answer *does this interval delimit the element*.
`boundary_method` is the only column that answers the second question, which is why it
has to be read first. This is not a one-off: at *good*, 5 of the 7 detected-but-broken
elements carry no warning flag at all while reporting a median 19% of the element.

**What to conclude from that row:** this genome carries conjugation machinery of MPF type
T on one contig, with a relaxase and an intact apparatus, and the region is a real
finding worth following up. **What not to conclude:** that the element is 16 kb; that it
is not an ICE (the integrase is elsewhere in the assembly, not absent from the genome);
or that the 16 kb interval defines its cargo.

## What this measurement could not establish

These are simulated breaks, not real short-read assemblies.

- **The cuts are cleaner than an assembler's.** Repeat-driven breakpoints are 68% of all
  breaks at the *good* arm but only 33% at *typical* and 14% at *poor* — the genomes
  simply do not contain enough IS and rRNA to reach 20 kb N50 — so the two harsher arms
  cut closer to uniform than reality, which under-breaks elements and makes those arms
  **optimistic**.
- **No misassemblies.** A real draft contains contigs that chimerically join
  non-contiguous sequence *within* one contig. That is exactly what per-contig clustering
  cannot see, and nothing here tests it.
- **No read-level effects** — coverage variation, collapsed repeats and the IS
  copy-number collapse are not represented at all. The optional
  [ISOSDB layer](optional-layers.md) is what puts a number on the last of those.

The IME arm is underpowered (3–5 of 12) and the negative-control arm is small (2–5 calls),
so those rows carry wide uncertainty. And several of ICEberg's own element boundaries are
themselves predictions, so "agrees with the curated length" means agreement with a
curation, not with truth.

Full working, including the two guards added after this measurement:
[`mobilome_draft_assemblies.md`](https://github.com/iLivius/BacFlux/blob/main/docs/mobilome_draft_assemblies.md).
