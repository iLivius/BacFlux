# Reading mobilome ICE/IME calls on a draft assembly

**Who this is for:** anyone running the mobilome module on short-read (Illumina)
assemblies rather than closed genomes — which is most BacFlux runs.

**The one-sentence version:** on a fragmented assembly the module keeps telling you
*what kind of element* it found and demotes honestly when the evidence breaks, but
it stops telling you *how big the element is*; read `mge_class` + `confidence`
alongside `spans_contigs` and `at_contig_boundary`, read `boundary_method` before
you read any coordinate, and expect more low-confidence calls than you would get
off a closed genome.

---

## 1. Why this document exists

Every validation this module ever had was run on **closed, single-contig genomes**.
On the whole 63-call closed benchmark, `spans_contigs` was `TRUE` on **zero** calls —
so the guards that exist specifically for fragmentation had never actually fired in
a measurement. The module's own spec (§2.3) calls short-read fragmentation "the
central constraint", but nothing had ever measured what it costs.

That has now been measured. 40 benchmark genomes (30 carrying curated ICE/IME
elements, 12 negative controls) were cut into contigs at three contiguities,
**re-annotated with Bakta** and re-run through the whole chain — ISEScan, CONJscan,
ICEscan, `conjscan_to_ice.py` — with the same commands and parameters as the closed
runs. 120 assemblies, 191 draft calls.

| arm | target N50 | achieved N50 (median) | contigs (median) | curated elements still on one contig |
|---|---|---|---|---|
| good | 150 kb | 135 kb | 32 | 21/30 |
| typical | 50 kb | 48 kb | 82 | 14/30 |
| poor | 20 kb | 20 kb | 180 | 8/30 |

Breaks were placed at IS elements first, then rRNA operons, then at random —
because those are what actually break a real assembly — and sequence was deleted at
each break. Contigs were length-sorted and about half reverse-complemented, so that
MacSyFinder's cross-contig clustering could not be flattered by contigs happening to
be in genome order.

**For calibration:** a survey of 28 real BacFlux Illumina assemblies on this host
gives a median N50 of **307 kb**, with 26 of the 28 at or above 150 kb. The sweep is
therefore deliberately pessimistic: a failure at "good" is a serious result, a
failure at "poor" may simply be the assembly.

---

## 2. What degrades, and by how much

| metric | closed | good (150 kb) | typical (50 kb) | poor (20 kb) |
|---|---|---|---|---|
| **ICE detected** (of 18 curated) | 15 | **15** | **15** | 12 |
| **IME detected** (of 12 curated) | 5 | 4 | 3 | 3 |
| recall, elements ≥50 kb (of 14) | 13 | 12 | 12 | 10 |
| median called / **true** length | 0.553 | 0.446 | 0.338 | **0.219** |
| median called / **largest surviving piece** | 0.553 | 0.574 | 0.771 | 0.721 |
| median call precision (overlap ÷ call length) | 1.0 | 1.0 | 1.0 | 1.0 |
| calls | 68 | 64 | 67 | 60 |
| `spans_contigs` = TRUE | 0 | 13 | 27 | 32 |
| `at_contig_boundary` = TRUE | 3 | 7 | 21 | 23 |
| confidence **high** | 16 (24%) | 12 (19%) | 7 (10%) | 2 (3%) |
| confidence **low** | 7 (10%) | 22 (34%) | 42 (63%) | 39 (65%) |
| `boundary_method` = `none` | 30 (44%) | 49 (77%) | 55 (82%) | 54 (90%) |
| calls on the 12 negative controls | 5 | 4 | 3 | 2 |
| MacSyFinder systems spanning >1 contig | 0 | 39/140 (28%) | 69/139 (50%) | 87/148 (59%) |

Three things in that table matter more than the rest.

**Detection is flat down to 50 kb N50.** ICE recall is *identical* to the closed
genomes at both 150 kb and 50 kb. It only breaks at 20 kb. The IME arm is weak
everywhere including on closed genomes — most IME misses are elements below the
module's 2 kb size floor (522 bp, 1.2 kb, 5.1 kb) that were already missed before any
cutting, so fragmentation costs 2 IMEs, not 9.

**The length collapses, but the caller is not to blame for it.** Against the *true*
element the median recovered length falls from 0.55 to 0.22. Against the **largest
surviving contig piece** — the ceiling any per-contig caller can possibly reach — it
is 0.55 → 0.57 → 0.77 → 0.72, i.e. flat or better. The caller keeps recovering the
same share of what is still visible; the missing bases are the assembly's.

**Nothing is invented.** Across all 191 draft calls: no fabricated att sequence
(every att pair reported on a draft occurs at least twice in the corresponding closed
genome), no class promoted at high confidence, no new false positive on a negative
control, and not one violation of the caller's own confidence rules in 163 checks.
Every high-confidence draft call sits where the closed genome also called something,
and 20 of 21 were high confidence there too.

---

## 3. Where the line actually is

**Not between "good" and "poor" — between *classification* and *extent*.**

- **N50 ≥ 150 kb** — usable as-is. ICE detection identical to closed genomes, every
  high-confidence call corroborated by the closed run, guards firing correctly.
- **N50 ≈ 50 kb** — the **class** is still trustworthy (ICE 15/18, and 91% of calls
  land on some entry in ICEberg's curation of these accessions, up from 83% on the
  closed genomes). The **extent** is not: median 0.34× the true length, and
  82% of calls have `boundary_method=none`. *Read the class, ignore the length.*
- **N50 ≈ 20 kb** — detection itself starts to go (ICE 12/18), 65% of calls are low
  confidence, and 59% of MacSyFinder's "systems" join genes from different contigs.
  Still not *wrong*, but it has stopped saying much: 2 high-confidence calls in 60.

Each sample's own contig count and N50 are written as the first row of its
`{sample}_ice_discarded.tsv` audit file (`action = assembly_qc`), together with one
sentence saying which of these three bands it falls in. Read that before the calls.

---

## 4. Which columns to read, and which to ignore

**Read these. They survive fragmentation and they never over-claim.**

| column | why |
|---|---|
| `mge_class` | ICE / IME / AICE / passive island / conjugative region. Demotes honestly when evidence is lost — never promotes. |
| `mobility` | the sentence that goes with the class. Always "**predicted** self-transmissible". |
| `confidence` | high / medium / low. Every high call in the benchmark satisfied the documented rule and was corroborated by the closed genome. |
| `spans_contigs` | TRUE = this system's hits are on more than one contig. **Absolute cap at low confidence** — 72 TRUE calls in the sweep, 72 at low, no exceptions. |
| `at_contig_boundary` | TRUE = the element runs into a contig end. Also an absolute cap: 54 TRUE, 54 at low. |
| `boundary_method` | **read this before any coordinate.** `tRNA` = real boundaries. `denovo` = a repeat reported as a lead, never used. `none` = the interval is the machinery span only. |
| `machinery_intact`, `n_anchor_classes` | what the class actually rests on. |
| `mpf_from_other_contig` | TRUE = this is an ICE only because the mating-pair genes are on a *different* contig (see §6). |

**Treat these as a floor, never a measurement, whenever `boundary_method = none`:**
`start`, `end`, `length_bp`. At the typical arm, calls on elements that fragmentation
had genuinely broken reported a **median 19% of the element**.

**Ignore on a draft:** any attempt to compare element lengths between samples, and
any arithmetic on `length_bp` (total mobilome burden, % of genome mobile, and so on).
Those numbers are governed by your assembly, not by the biology.

---

## 5. What a low-confidence call means on a draft

**It means the assembly could not support a stronger claim — not that the call is
wrong.** Low confidence is set by exactly the guards above: hits on more than one
contig, or an element running into a contig end, or incomplete machinery, or fewer
than three anchor classes. On a fragmented assembly those conditions are the *normal*
state of affairs, which is why low-confidence calls go from 10% of calls on closed
genomes to 63% at 50 kb N50.

Evidence that the low tier is not just noise: on the 12 negative-control genomes the
number of calls went **down** under fragmentation (5 → 4 → 3 → 2), none at high
confidence, and **not one** was at a locus the closed genome had not already called.
Contig ends did not manufacture false positives; they removed evidence.

So: a low-confidence ICE call on a draft is a lead worth following — look at
`missing_components`, `spans_contigs` and the audit file to see *what* was missing —
not a result to report and not a result to discard silently.

---

## 6. What fragmentation does to the class, and the one place it can still overclaim

Per call, back-translated onto the closed genome's coordinates:

| direction | good | typical | poor |
|---|---|---|---|
| **demoted** vs the closed call | 8 | 26 | 34 |
| … ICE → IME (tra operon lost) | 1 | 5 | 5 |
| … ICE → conjugative_region (integrase lost) | 7 | 16 | 20 |
| **promoted** vs the closed call | 0 | 1 | 0 |
| … at high confidence | 0 | **0** | 0 |
| calls **created** where the closed genome had nothing | 1 | 1 | 1 |
| … at high confidence | 0 | **0** | 0 |

The dangerous direction is essentially empty. The single promotion carries
`spans_contigs=TRUE` and low confidence — the guard caught it. Total ICE-class calls
fall 37 → 32 → 24 → 16: fragmentation makes the module *more* conservative about
tier 6, not less.

**The one exemption that can still buy a tier-6 label.** On a draft, the *tra* operon
routinely lands on a different contig from the relaxase, so a cluster with no
mating-pair gene of its own is still called an ICE when the typed system's other hits
are elsewhere. Without that exemption real ICEs get demoted the moment an assembly
breaks (a named positive control did exactly that in an earlier sweep). The cost is a
row that reads as self-contradictory: `mobility = predicted self-transmissible` next
to `missing_components = mating-pair apparatus`.

**`mpf_from_other_contig = TRUE` is the column that explains it.** It fired on 9 of
191 draft calls (2 good / 1 typical / 6 poor) and on **0** of 68 closed-genome calls.
All 9 also carry `spans_contigs = TRUE`, so all 9 are already at low confidence. Two
of them are on genuinely curated elements. Treat such a row as "conjugation machinery
is present in this genome and *may* belong to this element" — the evidence for the
mating bridge is not on the contig you are looking at, and cannot be checked there.

**One thing this cannot see.** MacSyFinder runs with `--db-type ordered_replicon`, so
on a draft it treats the whole assembly as one replicon and can build "systems" from
genes at opposite ends of the chromosome: 28% / 50% / 59% of its systems span more
than one contig, and in 87–93% of those the joined sequence was not contiguous in the
real genome (median gap bridged: 1.36 Mb at the good arm). `conjscan_to_ice.py`
clusters per contig, so **no reported interval ever joins non-contiguous sequence** —
the damage is confined to the inherited system label, and `spans_contigs` plus
`mpf_from_other_contig` are what make it visible.

---

## 7. Worked example: the same element, closed and fragmented

**ICE*Msp*.M1D**, a 199,376 bp chromosomal ICE in a *Mesorhizobium* sp.
(`CP034444.1`), one of the Phase 7 curated pilot elements. Values are taken verbatim
from `{sample}_ice_candidates.tsv` in each run.

| column | **closed** | **good** (N50 148 kb) | **typical** (N50 49 kb) | **poor** (N50 20 kb) |
|---|---|---|---|---|
| element on how many contigs | 1 | 3 | 9 | 13 |
| largest surviving piece | 100% | 79% | 60% | 14% |
| `mge_class` | **ice** | conjugative_region | conjugative_region | conjugative_region |
| `mobility` | predicted self-transmissible | conjugative region, boundaries not established | conjugative region, boundaries not established | conjugation machinery without a relaxase … |
| `start`–`end` | 1,854,246–1,893,682 | 16,921–32,922 | 2,328–18,329 | 85–13,216 |
| `length_bp` | 39,437 | 16,002 | 16,002 | 13,132 |
| `n_anchor_classes` | 4 | 3 | 3 | 2 |
| `has_integrase` | TRUE | **FALSE** | **FALSE** | FALSE |
| `has_relaxase` | TRUE | TRUE | TRUE | **FALSE** |
| `boundary_method` | denovo | **none** | **none** | **none** |
| `spans_contigs` | FALSE | FALSE | FALSE | **TRUE** |
| `at_contig_boundary` | FALSE | FALSE | FALSE | **TRUE** |
| `confidence` | **high** | **high** | **high** | **low** |

**How to read that, line by line.**

- The class demotes for a stated reason: the **integrase is lost to a contig break**,
  and without an integrase the caller refuses to say "ICE" (spec §8 Phase 4: relaxase
  + T4SS and no integrase is "report, don't call ICE"). It does not guess. That is
  the demotion working, not a bug.
- The reported length is 16 kb at *good* and *typical*. The element is **199 kb**.
  The call is the *machinery span* and nothing more — and `boundary_method = none` is
  the column that says so. `length_bp` is a floor.
- **The uncomfortable part:** at *good* and *typical* this call is `high` confidence
  with `spans_contigs = FALSE` and `at_contig_boundary = FALSE`. Both flags are
  literally true — the machinery cluster the caller saw really is entirely on one
  contig, well away from its ends — but the element continues on contigs it never
  looked at. Confidence answers "is this conjugation machinery, and of what class",
  **not** "does this interval delimit the element". `boundary_method` is the only
  column that answers the second question, which is why it must be read first.
  This is not a one-off: at *good*, 5 of the 7 detected-but-broken elements carry no
  warning flag at all while reporting a median 19% of the element, 3 of them at high
  confidence.
- At *poor* the guards finally fire — `spans_contigs`, `at_contig_boundary`, low
  confidence — and the row correctly reads "I am not sure about this."

**What you should conclude from the *good* row:** this genome carries conjugation
machinery of MPF type T on contig NODE_14, with a relaxase and an intact apparatus,
and the *Mesorhizobium* symbiosis island region is a real finding worth following up.

**What you must not conclude:** that the element is 16 kb; that it is not an ICE
(the integrase is elsewhere in the assembly, not absent from the genome); or that the
16 kb interval defines its cargo. On a closed genome all three of those questions
would have answers. On this draft they do not.

---

## 8. The two guards added after this measurement

Both were added because the measurement showed a real leak. Neither changes any
result on a closed genome — re-running the caller with the guards in place over the
closed benchmark reproduced its **68 calls with zero differences in any column**, and
re-running it over all 120 fragmented assemblies changed nothing except the att
columns and the new flag: same classes, same confidences, same coordinates.

**1. De novo att repeats are now counted across the whole assembly, not just the
contig.** A real integration scar exists in exactly two copies (attL and attR), so
`att_search` already refuses a repeat with three or more. But it counted copies on the
contig it was given — and on a closed genome the contig *is* the assembly, so nobody
noticed the two questions were different. On a draft they are not: **all 35** de novo
repeats reported across the three arms had exactly 2 copies on their own contig, but
**10 of them had 3–30 copies across the assembly.** The worst was an 89 bp repeat with
30 genome-wide copies attached to a *high-confidence* call. Those 10 boundary claims
are now withdrawn (`boundary_method` returns to `none`, att columns cleared, audit
reason `denovo_att_is_a_repeat_family`). No coordinate moves — a de novo repeat never
widened an element in the first place.

**2. `mpf_from_other_contig`** — a new output column, described in §6. Informative
only: it changes no class and no confidence, because the class is what stops real
ICEs being demoted by a contig break and the confidence is already capped at low by
`spans_contigs`.

**A guard that was considered and rejected on the evidence.** "Refuse to widen an
element onto an att pair when either copy sits within *N* bp of a contig end" sounds
prudent, and it is wrong here. Four tRNA-anchored att pairs on drafts sit within 2 kb
of a contig end; **three of them reproduce the closed genome's own answer exactly**
(identical element lengths of 20,571 bp and 25,404 bp), and the fourth is ICE*Ec2* on
`GU725392`, a named Phase 7 positive control recovered at 99.9% precision — whose att
pair sits **26 bp** from the contig end on the closed genome too, because the element
*is* the deposited record. The guard would have cost a positive control and caught
nothing: of the 43 att pairs reported across the three draft arms, **zero** are
fabricated — every one is sequence that occurs at least twice in the corresponding
closed genome. It was not added.

---

## 9. What this measurement could not establish

These are simulated breaks, not real short-read assemblies. Three consequences that
remain untested:

- **The cuts are cleaner than a real assembler's.** Repeat-driven breakpoints are 68%
  of all breaks at the *good* arm but only 33% at *typical* and 14% at *poor* — the
  genomes simply do not contain enough IS and rRNA to reach 20 kb N50 — so the two
  harsher arms cut closer to uniform than reality, which under-breaks elements and
  makes those arms **optimistic**.
- **No misassemblies.** A real draft contains contigs that chimerically join
  non-contiguous sequence *within* one contig. That is exactly the failure the
  per-contig clustering guard cannot see, and nothing here tests it.
- **No read-level effects** — coverage variation, collapsed repeats, and the IS
  copy-number collapse the spec warns about are not represented at all.

Also: the IME arm is underpowered (3–5 detected of 12) and the negative-control arm is
small (2–5 calls), so those two rows carry wide uncertainty; and several of ICEberg's
own element boundaries are themselves predictions, so "agrees with the curated length"
means agreement with a curation, not with truth.

---

## 10. Where the numbers come from

- Fragmented assemblies, coordinate translation and per-call scoring:
  `<validation-root>/draft_validation/`
  (`assembly_stats.tsv`, `element_coordinates.tsv`, `element_summary.tsv`,
  `calls_all.tsv`, `scoring/*.tsv`).
- Closed-genome baseline: `<validation-root>/phase7_benchmark/`
  (read-only; not modified by this work).
- Re-run of the caller with the guards in place, used for §8:
  `draft_validation/reverify/`.
- The closed-genome benchmark itself is written up in `docs/mobilome_tuning_guide.md`
  and `docs/methods_att_and_small_plasmids.md`; the design reasoning behind the
  fragmentation constraint is `docs/mobilome_module_SPEC.md` §2.3.
