# Method reference: intrinsic versus acquired AMR, and what BacFlux can actually say

A regulator asking about a resistance gene wants to know one thing: did this strain
always have it, or did it acquire it? This page explains why that question is harder
than it sounds, what BacFlux can answer, and what it deliberately refuses to answer.

Every claim carries its source, and the final section separates peer-reviewed
guidance from technical reports and from inference, so the page can be quoted in a
methods section. For the short version, see
[The mobility ladder](mobilome/mobility-ladder.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **intrinsic** | resistance shared by essentially all strains of a species — a property of the species, not of this isolate |
    | **acquired** | resistance gained horizontally, from a plasmid or mobile element |
    | **EFSA** | European Food Safety Authority, whose guidance drives the intrinsic/acquired question |
    | **FEEDAP** | the EFSA panel that issues the guidance for feed additives |
    | **QPS** | Qualified Presumption of Safety, EFSA's list of species accepted as safe |
    | **MIC** | minimum inhibitory concentration — the lab measurement of resistance |
    | **ANI** | average nucleotide identity, a genome-to-genome similarity measure |
    | **IS** | insertion sequence, the smallest kind of mobile element |

Three questions this page answers, all of which arise when an `08.mobilome` table
is read next to a regulatory checklist:

1. AMRFinderPlus reports point mutations. How is that different from a BLAST
   search like ABRicate's, and why does it need `--organism`?
2. EFSA defines an "intrinsic" AMR gene as one shared by the vast majority of
   strains of a species. Does that mean BacFlux would have to compare against
   every genome of that taxon in NCBI?
3. EFSA publishes a catalogue of AMR genes in *Bacillus*. Is that a new AMR
   database? Should BacFlux use it?

Short answers: it is a fundamentally different question, not a stricter search;
yes, and that is why BacFlux does not attempt it; and no, it is not a database
and BacFlux does not use it.

---

## 1. Two kinds of resistance, and why one is invisible to BLAST

A bacterium becomes resistant in one of two ways.

**It acquires a gene it did not have.** A *bla*CTX-M or a *tet*(M) arrives by
horizontal transfer. Before, the gene was absent; after, it is present.

**It mutates a gene it already had.** A single amino-acid substitution in *gyrA*
(fluoroquinolones), *rpoB* (rifampicin) or *rpsL* (streptomycin). These are
essential genes, present in every strain of the species whether resistant or
susceptible. One residue differs.

ABRicate answers *"is this gene present?"* by BLAST, and reports identity and
coverage. That question is exactly right for the first mechanism and structurally
blind to the second: query a susceptible genome for *gyrA* and you get a hit at
~99.9% identity; query a resistant one and you get a hit at ~99.9% identity. The
substitution that decides the phenotype is averaged away inside the percentage.
No threshold separates them, because the difference is not one of degree — it is
one residue at one coordinate. Lowering the cutoff only adds noise.

AMRFinderPlus with `--organism` asks a different question. It loads a curated set
of *known resistance-conferring substitutions* for that species, aligns the query
to the reference allele, and inspects the specific coordinate: not "is *gyrA*
here" but "is residue 83 of GyrA a serine or a leucine". Hits are reported with
`Method: POINTX` / `POINTP`, naming the substitution.

Two further differences matter even where no mutation is involved:

- **Per-gene curated cutoffs** rather than one blanket threshold. Gene families
  differ in natural diversity, so a single cutoff is either too permissive for
  conserved families or too strict for variable ones. This is why BacFlux applies
  the EFSA 80%/70% thresholds to the **ABRicate leg only** — imposing them on
  AMRFinderPlus would replace curation with a blunter rule.
- **Structured output**: element type (AMR / STRESS / VIRULENCE), class and
  subclass, and a partial-at-contig-end flag. The mobilome module consumes these.

### Why this matters for the intrinsic/acquired question

Point mutations *are* the chromosomal, non-transferable mechanism. So this is the
one place a single genome can speak to the intrinsic side. Everything else the
module does — replicon, mobile-element context, mobility tier — addresses the
acquired and mobile side.

### The caveat, measured on real data

`--organism` only engages for organisms AMRFinderPlus curates mutations for:
roughly three dozen, overwhelmingly clinical — *Escherichia*, *Salmonella*,
*Klebsiella*, *Staphylococcus aureus*, *Campylobacter*, *Pseudomonas aeruginosa*,
*Enterococcus*, *Vibrio*, *Neisseria* and similar.

On a batch of plant- and environment-associated strains — the Bacillaceae and
soil Pseudomonadota that make up most environmental collections — `--organism`
was applied to **none of them**. Every
`{sample}_amrfinder_organism_audit.tsv` recorded the same reason: *no curated
organism for this taxon*. Every `{sample}_amrfinderplus_mutations.tsv` was
therefore header-only.

**Read those empty files as "not assessed", never as "no mutations found."**
BacFlux writes the refusal and its reason to the audit file so the distinction
survives.

---

## 2. Why BacFlux does not classify genes as intrinsic or acquired

EFSA's definitions are population-level, and the wording is explicit:

> **Intrinsic AMR gene** — "Gene inherent to strains of a bacterial species …
> An AMR gene is considered 'intrinsic' when it is shared by the vast majority of
> wild type strains of the same species (or subspecies) and is restricted to those
> located on the chromosome."
>
> **Acquired AMR gene** — "A resistance gene novel for the strain under
> assessment, acquired through horizontal transfer … Acquired AMR genes could be
> integrated in the bacterial chromosome **or** harboured on a separate genetic
> element."
>
> — EFSA Scientific Committee (2025), Glossary

Two consequences follow, and both constrain what this workflow may claim.

**You cannot establish "intrinsic" from one genome.** The claim is about the
species, not the isolate. Determining it means comparing gene presence across many
genomes of that taxon — which is why EFSA publishes a separate method
(EFSA BIOHAZ Panel, 2023a) and a separate tool for it, described as a *"Pipeline
for the automated analysis of gene distribution in microbial species"*.
BacFlux implements no population comparison of any kind: every analysis is
per-genome.

**Chromosomal does not mean intrinsic.** An acquired gene may sit on the
chromosome. This is why mobility tier 1 is labelled **`intrinsic_candidate`** and
not `intrinsic`. The qualifier is doing real work and should not be dropped: the
module observed that a gene is chromosomal with no mobile-element context, which
is *consistent with* intrinsic and does not establish it.

So the honest framing, which the documentation uses throughout: **BacFlux produces supporting
evidence for the intrinsic/acquired judgement. It does not produce the
judgement.** The confirmatory work is a species-wide distribution analysis, plus
phenotypic testing (MIC against the relevant antimicrobial), per the decision tree
in EFSA (2025) §3.2.1.

Where BacFlux does align well: EFSA (2024a) requires the AMR search to run against
**at least two** maintained/curated databases, at ≥80% identity and ≥70% length
coverage. BacFlux runs ten database queries across three independent legs —
ABRicate over eight databases (`argannot`, `card`, `ecoh`, `ecoli_vf`, `megares`,
`ncbi`, `resfinder`, `vfdb`), AMRFinderPlus, and a BBMap→CARD read-mapping leg
that is immune to assembly collapse — with those thresholds applied on the
ABRicate leg. That alignment was not designed against the citation; it predates it.

---

## 3. The EFSA *Bacillus* catalogue — what it is, and why we do not use it

**It is not an AMR gene database.** It is a **prevalence catalogue**: how often
each of a set of AMR gene sequences occurs among published genomes of five
*Bacillus* species. It answers the population question BacFlux structurally
cannot.

### How it was built

*(Read from the report itself, EFSA 2024b, §3.)*

**Genomes.** NCBI **RefSeq**, restricted to assembly level *complete genome*,
snapshot **30 November 2023**.

**The query set is the unusual part.** The AMR sequences searched for were **not**
taken from CARD, ResFinder or NCBI as a database. They were **extracted from
dossiers submitted to EFSA by applicants** up to December 2023 — the AMR hits
those applicants had themselves reported in their production or active-agent
strains. Identifiers arrived in four different currencies (NCBI protein ID, NCBI
nucleotide accession, **CARD ARO ID**, **UniProt** accession) and were all
normalised back to NCBI protein IDs. Hits that resolved to no database record were
dropped.

So CARD and UniProt appear only as *lookup services* for resolving identifiers,
not as the source of the gene set. The set is bounded by what applicants happened
to report, not by AMR biology.

**Species selection.** Nineteen *Bacillus* and related species from the QPS list
were screened; those with **≥30 complete genomes** were analysed. Five qualified:

| Species | Screened | Confirmed (ANI ≥95%) |
|---|--:|--:|
| *B. velezensis* | 329 | 328 (99.7%) |
| *B. subtilis* | 299 | 299 (100%) |
| *B. amyloliquefaciens* | 76 | 57 (75%) |
| *B. licheniformis* | 40 | 40 (100%) |
| *B. paralicheniformis* | 26 | 26 (100%) |

*B. paralicheniformis* was included below the threshold because many dossiers
concerned it. *Priestia megaterium* had 40 genomes and was excluded anyway, because
no dossier had ever been submitted for it. The catalogue is shaped by the
regulatory caseload, not by the genus.

**Pipeline.** fastANI v1.32 (one-to-many) against each species' reference/type
strain; genomes at **ANI ≥95%** retained; ANIclustermap for trees. A local BLAST
database built from the retained genomes, with **plasmid sequences separated from
chromosomes** by a custom script so plasmid hits are reported separately. Searches
with **tblastn** (protein queries) and **blastn** (nucleotide queries), E-value
0.05. Then: deduplicate within a genome keeping the highest identity, apply
**≥80% identity / ≥70% coverage**, compute the frequency of strains matching each
query, and the median identity and coverage. Python 3 with R 4.3.3 for plots.

**Controls, which are worth copying.** Every table carries internal positive
controls — *rpoB*, *gyrA*, *gyrB* and 16S rRNA from the species' own reference
genome, which should be near-universal — and a negative control, *fimA* from
*Salmonella*, which should be absent. In the published tables the housekeeping
genes return ~100% of genomes and *fimA* returns 0% every time. The pipeline
demonstrates its own validity in each result table.

**Scale of the output** (Annex A): 358 raw matches across the five species,
192 valid and non-duplicated, and **95 hits above threshold** actually analysed.

### The finding most relevant to BacFlux

The catalogue is a clear demonstration of **threshold brittleness at 80%
identity**:

- *rphB* in *B. subtilis*: found above threshold in **85.3%** of genomes when
  queried as nucleotide, but **1.3%** when queried as protein — because most
  genomes sit at **79–79.9%** identity, just under the cutoff.
- *rphC* behaves the same way in *B. amyloliquefaciens* (79.8% median) and
  *B. velezensis* (79.8%): absent by the rule, near-universal in reality.
- *satA* in *B. amyloliquefaciens*: ~91% of genomes by nucleotide, ~12–21% by
  protein, with a median identity of 79.2%.

The same sequence is "present in nearly all strains" or "essentially absent"
depending on whether it is queried as DNA or as protein. Since BacFlux applies the
same 80%/70% EFSA thresholds on its ABRicate leg, the lesson transfers directly:
**a gene just under the cutoff is not absent, and a borderline call should be read
from the actual identity value, not from presence/absence.** This is a further
argument for keeping AMRFinderPlus's per-gene curated cutoffs unmodified.

Also instructive: *tet*L sits in ~63% of *B. subtilis* chromosomes — neither
"vast majority" nor rare, i.e. genuinely undecidable by the intrinsic rule. Genes
originating in other species (ANT(4′)-Ib from an *S. aureus* plasmid, *erm*B from
*E. faecium*, *cat* from a *Streptococcus* phage) appear in 1–4% of genomes, which
is what a real acquired gene looks like.

### Why it is not built into BacFlux

**Coverage.** It spans five species. On an environmental batch of the kind this
document was written from, fewer than a fifth of the isolates fell inside those
five; roughly another third were *Bacillus* or *Priestia* species outside them,
and the rest were other genera entirely. For most of a typical environmental
collection it would return nothing.

**It is a snapshot.** RefSeq at 30 November 2023, dossiers to December 2023. It
carries no update mechanism.

**It answers a question BacFlux does not ask.** Wiring it in as another ABRicate
database would be a category error — it is not a detection resource. The only
sensible use is the reverse: after AMRFinderPlus reports a gene in one of the five
species, look up whether EFSA already regards that gene as widespread there. That
is a manual cross-check, and it is how we treat it.

**Recommended use:** consult it by hand for *B. subtilis*, *B. velezensis*,
*B. amyloliquefaciens*, *B. licheniformis* and *B. paralicheniformis* isolates
when a regulatory question is actually in play. Note that the general-purpose
pipeline behind it is published separately and can be run for *any* species with
enough genomes — that, not the catalogue, is the reusable artefact.

---

## Sourcing: standing, licence, and stated limitations

Not all of these carry the same weight, and the difference matters if any of this
is quoted in a dossier.

**Peer-reviewed EFSA guidance (highest standing)**

- **EFSA Scientific Committee (2025).** *Guidance on the characterisation of
  microorganisms in support of the risk assessment of products used in the food
  chain.* EFSA Journal 23(11):e9705. doi:10.2903/j.efsa.2025.9705. Adopted
  24 September 2025. Open access (CC BY-ND). **This is the current guidance**, and
  its Appendix A lists the parts of the 2018 FEEDAP guidance it supersedes.
  Definitions quoted above are from its Glossary; the AMR decision tree is §3.2.1.
- **EFSA BIOHAZ Panel (2023a).** *Statement on how to interpret the QPS
  qualification on 'acquired antimicrobial resistance genes'.* EFSA Journal
  21(10):8323. doi:10.2903/j.efsa.2023.8323. The method for discriminating
  intrinsic from acquired.
- **EFSA (2024a).** *EFSA statement on the requirements for whole genome sequence
  analysis of microorganisms intentionally used in the food chain.* EFSA Journal
  22(8):e8912. doi:10.2903/j.efsa.2024.8912. Source of the ≥2-databases
  requirement and the 80%/70% thresholds.

**Technical report (lower standing — not a peer-reviewed panel opinion)**

- **EFSA (2024b).** *Catalogue of antimicrobial resistance genes in species of
  Bacillus used to produce food enzymes and feed additives.* EFSA Supporting
  Publication 2024:EN-8931, 33 pp. doi:10.2903/sp.efsa.2024.EN-8931. Approved
  3 July 2024. ISSN 2397-8325. *"Reproduction is authorised provided the source is
  acknowledged."*

  It is an EFSA **Technical Report**, not an EFSA Journal guidance or opinion —
  a commissioned technical output. Its standing comes from being cited by the
  2025 Guidance (as EFSA, 2024b) as an available implementation. Cite it as a
  supporting technical report, not as guidance.

  **Limitations the report states about itself:** 95% ANI cannot separate
  *B. velezensis* from *B. amyloliquefaciens*, so their discrimination "only
  relies on the annotation by RefSeq, which constitutes a limitation in the
  analysis"; several NCBI accessions used were already superseded at publication.

  **Limitation we add:** the query set is dossier-derived. It reflects what
  applicants reported, so absence from the catalogue is not evidence that a gene
  is absent from the species — only that no dossier raised it.

**Non-EFSA**

- **Siguier P, Gourbeyre E, Varani A, Ton-Hoang B, Chandler M (2015).**
  *Everyman's guide to bacterial insertion sequences.* Microbiology Spectrum
  3(2):MDNA3-0030-2014. doi:10.1128/microbiolspec.mdna3-0030-2014. Cited for IS
  size range (0.7–2.5 kb) and abundance.

**Verified from software, not literature**

- The AMRFinderPlus organism list, the `--organism` refusals and the empty
  mutation tables were read from this repository's own run outputs and audit
  files, not from documentation.

**Inference, not quoted from any source**

- That ABRicate-style BLAST is *structurally* unable to resolve point mutations.
  This follows from what the tools compute and is not a claim any of the cited
  documents makes.
- That the threshold-brittleness seen in the catalogue transfers to BacFlux's own
  ABRicate leg. The mechanism is identical (same thresholds, same kind of search),
  but no source states it about BacFlux.
