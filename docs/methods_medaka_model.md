# The Medaka polishing model

This page explains how BacFlux selects a Medaka model, why the bacterial model can
be appropriate for canonical basecalls from native bacterial DNA, and what the cited
benchmarks establish.

!!! tip "Practical guidance"

    - `medaka_model: auto`, the default, asks Medaka for the bacterial consensus
      model compatible with the basecaller recorded in the reads. If none resolves,
      the run stops before assembly.
    - The cited ten-genome comparison reported 26 residual errors with the bacterial
      model and 100 with the standard model. The benefit depends on the dataset.
    - Check assembly structure before polishing: missing replicons and unresolved
      circular overlaps can cause misleading read alignments.
    - For amplified DNA, set `medaka_model` to the standard model matching the
      basecaller (see [Overriding the model](#overriding-the-model)). This rests on
      the bacterial model's stated scope; its effect on amplified libraries has not
      been measured.

The sections below describe the implementation and evidence separately.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies, the long-read sequencing platform |

---

## What BacFlux does

In long-read modes (`nanopore`, `hybrid`) with `medaka_model: auto`, BacFlux asks
Medaka to resolve the model from the reads themselves, requesting the *bacterial*
variant — `workflow/scripts/12_medaka_check/medaka_model_check.py`:

```
medaka tools resolve_model --auto_model consensus_bacteria <reads>
```

Medaka reads the `basecall_model_version_id` tag out of the FASTQ headers, maps it
to the matching consensus model, and then substitutes the bacterial model — but
only if the detected basecaller is on its compatibility list. On a typical R10.4.1
`sup` run this resolves to `r1041_e82_400bps_bacterial_methylation`.

That command is character-for-character what `medaka_consensus` itself runs for its
`--bacteria` flag, so `auto` really does get the bacterial model — this is verified
against the installed Medaka source, not inferred. If resolution fails (no
basecaller tag, or a tag Medaka's lookup does not know) the check exits 1 with a
suggestion table; it never falls back to a non-bacterial model. Model choice is
worth that strictness: across ten genomes it was the difference between 26 and 100
residual errors ([Wick 2025](https://rrwick.github.io/2025/02/07/dorado-polish.html)).

This happens in `check_medaka_model` (`workflow/rules/shared/12_medaka_check.smk`),
which runs *before* assembly so a bad or incompatible model fails in seconds rather
than after the assembler has spent an hour.

The same check puts the model's weights on disk. Medaka ships the weights of only a
few models and downloads the rest from GitHub the first time they are needed:
`resolve_model` above does that download as part of resolving the model, and an
explicit name is followed by `medaka tools download_models --models <name>`. The
check is a local rule, so on a cluster it runs on the computer that runs Snakemake,
and the polishing job itself needs no internet.

!!! note "Naming drift: `dorado polish --bacteria`"

    Wick's 2026 posts recommend `dorado polish --bacteria` in place of Medaka. As of
    his 2025-02-07 measurement the two shipped identical model weights, though he
    hedged that with "At the time of writing" and nobody has re-checked since. Noted
    only so the tool names in those posts do not read as a contradiction of this
    page; BacFlux stays on Medaka.

## Basecalling and polishing models { #why-this-is-not-a-mismatch }

The objection assumes "methylation model" means the same thing at both stages. It
does not — the two are different tools doing different jobs:

| | Dorado modified-base model | Medaka bacterial model |
|---|---|---|
| Stage | basecalling | consensus polishing |
| Job | **reports** methylation (`MM`/`ML` tags, per-base 5mC/6mA) | **corrects errors caused by** methylation |
| Effect on the A/C/G/T sequence | none | changes it |
| Output | an epigenetic annotation | a more accurate consensus |

The Medaka bacterial model is not the polishing counterpart of a Dorado
modified-base model. It is a consensus model **for ordinary canonical basecalls**,
specialised for the case where the template DNA is methylated.

**The model-matching rule is still satisfied.** Medaka models must match the error
profile of the basecaller that produced the reads — ONT states this directly:

> "it is important to specify the correct inference model, according to the
> basecaller used."
> — Medaka README (nanoporetech/medaka, v2.2.1)

The bacterial model declares compatibility with nine **canonical** basecaller
models (hac and sup, v4.2.0 → v6.0.0; `medaka/options.py`,
`bact_methyl_compatible_models`). There are no modified-base models on that list.
If it required methylation-aware basecalling, it would name those instead. Medaka
enforces this: an incompatible basecaller raises `RuntimeError` rather than
silently polishing with the wrong model (`medaka/models.py`).

## Methylation-related consensus errors { #the-point-that-resolves-the-paradox }

The methylation errors exist because the basecaller was *not*
methylation-aware. Bacteria methylate their own genomes (6mA, 5mC, 4mC, largely via
restriction–modification systems). A modified base shifts the raw nanopore current
away from the canonical signal, so a canonical basecaller — which can only emit
A/C/G/T — misreads those positions systematically, at the same sequence motifs
every time. Because the errors are systematic rather than random, **more coverage
does not remove them.**

Those errors are already in the basecalls. The bacterial polishing model is
trained to correct exactly that residue, which is why ONT scopes it to native
material:

> "For native data with bacterial modifications, such as bacterial isolates,
> metagenomic samples, or plasmids expressed in bacteria, there is a research model
> that shows improved consensus accuracy."
> — Medaka README (nanoporetech/medaka, v2.2.1)

ONT's own words: a **research model**, whose compatibility is deliberately broad —
*"compatible with several basecaller versions for the R10 chemistries"*.

In Dorado (checked against version 2.1.1 on 2026-08-21), modified-base calling writes
modification probabilities beside the canonical calls and does not change which base
is called. Basecalling with a modified-base model therefore leaves the same
methylation-induced errors, and the same Medaka bacterial model still applies. The
evidence and its limits are listed at the end of this page.

## Does it actually help?

The peer-reviewed measurement is Nagy *et al.* / NEKSUS Consortium (2026):
92 clinical Enterobacterales, R10.4.1 chemistry, Dorado `sup` v5.0.0 basecalls,
assembled with Flye, Hybracter and Autocycler and polished with Medaka v2.

> "Medaka polishing with un-subsampled long-reads resulted in small improvements in
> indels, but not SNVs for both Flye and Autocycler assemblies."
>
> "Seven-locus multi-locus sequence type, antimicrobial resistance, virulence and
> stress gene annotation was equivalent across assembler/polisher combinations."
>
> — Nagy *et al.* / NEKSUS Consortium (2026), *Microbial Genomics* 12(2):001631.
> https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001631

Of the benchmarks cited here, it is the only one run on Flye assemblies, which is what
BacFlux polishes: Wick's polishing comparisons are on Autocycler assemblies or
deliberately degraded Raven drafts. In that dataset, polishing gave small improvements
in indels and none in SNVs, and the assessed annotations were equivalent. These results
do not guarantee unchanged typing, AMR or virulence calls for another dataset.

A smaller and earlier observation, on two bacterial genomes, compared Medaka v2 with
the bacterial model against the previous release:

| Genome | Errors before | Errors after |
|---|---|---|
| *Campylobacter lari* | 18 | 2 |
| *Enterobacter cloacae* | 11 | 10 |

> Wick RR (2024) *Medaka v2: progress and potential pitfalls.*
> https://rrwick.github.io/2024/10/17/medaka-v2.html

Two genomes is a small sample, and the gain is variable — large on one, marginal on
the other. The result does not establish a fixed benefit across organisms or preparations.

The largest tally is 147 assemblies, where the bacterial model improved 129 and made
6 worse ([Wick 2026](https://rrwick.github.io/2026/06/19/dorado-v2-polishing.html)).
The reported benefit was not universal. Assembly structure is one possible source
of polishing errors, as described below; these counts alone do not establish the
cause of every worsened assembly.

!!! warning "An older recommendation said the opposite, and is still findable"

    In October 2023, Wick tested Medaka on `sup` basecalls and concluded: *"it
    performed erratically with the sup/res assemblies, making things worse more often
    than better"*, recommending *"Use Medaka if you are assembling hac reads but skip
    it for sup/res reads"*
    ([post](https://rrwick.github.io/2023/10/24/ont-only-accuracy-update.html)).

    That test used `r1041_e82_400bps_sup_v4.2.0` — the **standard** model. The bacterial
    model did not exist yet. Measured against each other on the same ten genomes, the
    standard model left 100 residual errors and the bacterial model 26, which is the gap
    that reverses the advice. The comparison supports the bacterial model for the native
    bacterial data tested. It does not establish a rule for every chemistry, basecaller
    or DNA preparation.

## Assembly structure and polishing { #a-separate-pitfall-worth-knowing }

Unrelated to model choice: **polish only structurally sound assemblies.** In the 2024
post above, where small plasmids were absent from the assembly, their reads
misaligned elsewhere and Medaka introduced over 100 erroneous changes. A missing
replicon can therefore affect the sequence of other assembled replicons during polishing.

Wick broadened the warning in a 2025-02-25 revision of his assembly guide:

> "Structural errors in your assembly (such as missing plasmids or start-end overlap
> of circular contigs) can sometimes lead to Medaka increa[s]ing the number of
> errors."
> — https://rrwick.github.io/2020/10/30/guide-to-bacterial-genome-assembly.html

So a missing plasmid is one case of a wider problem: an unresolved start–end overlap on
a circular contig can also cause polishing errors. BacFlux's reorientation,
plasmid-recovery and replicon-audit steps run before polishing partly for this reason.

## Overriding the model

The bacterial model assumes **native, unamplified DNA**. Amplification — ONT's PCR
barcoding kit (SQK-RPB114.24) or a whole-genome amplification protocol — copies the
template with unmodified bases, so the methylation the model corrects for is no longer
there. If a library was amplified, pin the matched standard model instead of `auto`:

```yaml
parameters:
  nanopore:            # or: hybrid
    medaka_model: r1041_e82_400bps_sup_v4.2.0
```

This is an example of a standard model; select one compatible with the actual
basecaller. The 100-versus-26 comparison above was not an amplified-DNA benchmark
and cannot quantify the cost or benefit of this choice on amplified libraries.

`check_medaka_model` validates any explicit name against the installed Medaka's model
list, fetches its weights, and fails early with a suggestion table if the name is
wrong. Amplification status cannot
be read from the reads; it has to come from the lab record.

## Sourcing: what is documented, what is inference

Stated so the claims above can be defended or challenged individually.

**Verified from primary sources**

- The `consensus_bacteria` call, and the early-validation rule — read from this
  repository (`medaka_model_check.py`, `12_medaka_check.smk`).
- The nine compatible canonical basecaller models, the substitution logic, and the
  `RuntimeError` on incompatibility — read from the installed Medaka 2.2.2 source
  (`options.py`, `models.py`).
- "native data with bacterial modifications… improved consensus accuracy",
  "research model", "compatible with several basecaller versions", and the
  model-matching requirement — quoted from ONT's Medaka README (v2.2.1).
- That `--auto_model consensus_bacteria` is the same call `medaka_consensus` makes
  for its `--bacteria` flag, and the exit-1-on-failure behaviour — read from the
  installed Medaka source and from `medaka_model_check.py`.
- The 92-isolate Flye/Hybracter/Autocycler result and the annotation-equivalence
  statement — quoted verbatim from Nagy *et al.* (2026), cited above.
- The error counts for *C. lari* and *E. cloacae*, and the small-plasmid pitfall —
  Wick (2024), cited above.
- The broadened structural-error warning, and the 26-versus-100 residual-error
  figure for model choice — Wick (2025), cited above.

**Inference, not quoted from documentation**

- The biological mechanism (6mA/5mC/4mC, restriction–modification systems, signal
  deviation at modified bases, motif-systematic errors that coverage cannot
  average out). Standard domain knowledge, but not taken from the sources above.
- **That amplified DNA is a contraindication.** ONT scopes the model to "native
  data" but does not explicitly warn against amplified input. The caution above is
  read off that scoping, not quoted from a warning.
- That methylation-induced errors are predominantly indels, and that frameshifted
  or pseudogene counts in the annotation are therefore a good readout for testing
  the model's benefit.

**Not re-checked:** whether `dorado polish --bacteria` still ships the same model
weights as Medaka's bacterial model. Wick reported that in February 2025 and hedged
it with "At the time of writing"; nothing here verifies it since.

**Verified against Dorado 2.1.1, 2026-08-21.** Modified-base calling does not feed back
into which base is called, so it cannot repair a methylation-induced miscall and does
not remove the need for this consensus model. Two independent lines establish it, and
the SAM specification is consistent with it:

- Dorado's basecall decoder emits only four symbols — `alphabet = {'A','C','G','T'}`
  in `basecall/decode/beam_search.cpp`. There is no symbol for a modified base. The
  modbase node reads the called sequence and writes only modification probabilities
  beside it; it never assigns to the sequence.
- ONT's maintainers state the invariant outright, and treat it as a design constraint:
  *"modification calling has no effect on the canonical basecall — you should get
  exactly the same sequence whether you do modbase or not. This is intentional"*
  ([dorado#1494](https://github.com/nanoporetech/dorado/issues/1494#issuecomment-3302392321),
  given as the reason for refusing a feature that would have let a modification
  probability alter the emitted base).
- The SAM specification describes how modification annotations relate to the primary sequence:
  base modifications are *"a series of edits from the primary unmodified sequence as
  originally reported by the sequencing instrument"*
  ([SAMtags](https://raw.githubusercontent.com/samtools/hts-specs/master/SAMtags.tex)),
  counted against the primary `SEQ` field. This format describes the output
  representation; it does not by itself prove how a basecaller uses modification
  information internally.

ONT's own documentation does not say this in words — it describes modified basecalling
only as "an extension to the normal `simplex` and `duplex` basecalling subcommands"
(<https://software-docs.nanoporetech.com/dorado/latest/basecaller/mods/>).

Two limits on the claim. It says no information flows back into base determination, not
that output is byte-identical: enabling modification models changes GPU memory use and
can shift the auto-selected batch size, and Dorado is not bit-reproducible across batch
sizes — a difference that is non-determinism rather than feedback, and cannot
systematically correct motif-specific errors. And it is specific to Dorado's split
architecture; it does not hold for older basecallers where basecalling and modification
calling happened inside one model.
