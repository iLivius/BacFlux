# The Medaka polishing model

Why BacFlux polishes an ONT assembly with a bacterial methylation model when Dorado
basecalled it with a canonical one. Written to be quotable in a methods section: every
claim carries its source, and the last section states which come from primary
documentation and which are inference.

Established 2026-08-04 by reading Medaka's own source and ONT's documentation,
prompted by an objection: **if Dorado basecalled without a methylation
model, why does BacFlux polish with a methylation-aware one?**

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

This happens in `check_medaka_model` (`workflow/rules/shared/12_medaka_check.smk`),
which runs *before* assembly so a bad or incompatible model fails in seconds rather
than after the assembler has spent an hour.

## Why this is not a mismatch

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

## The point that resolves the paradox

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

**Corollary.** Had the run been basecalled *with* `--modified-bases`, the A/C/G/T
sequence would be unchanged: modified-base calling adds an annotation layer on top
of the canonical calls, it does not re-derive them. You would gain methylation
calls and keep exactly the same methylation-induced errors — so the same Medaka
model would still be the right choice. (See the sourcing note below: this
corollary is inference, not quoted documentation.)

## Does it actually help?

Independently measured on two bacterial genomes, comparing Medaka v2 with the
bacterial model against the previous release:

| Genome | Errors before | Errors after |
|---|---|---|
| *Campylobacter lari* | 18 | 2 |
| *Enterobacter cloacae* | 11 | 10 |

> Wick RR (2024) *Medaka v2: progress and potential pitfalls.*
> https://rrwick.github.io/2024/10/17/medaka-v2.html

The gain is real but variable — large on one genome, marginal on the other. Expect
it to scale with how heavily methylated the organism is, not to be a fixed benefit.

## A separate pitfall worth knowing

From the same source, and unrelated to model choice: **polish only structurally
sound assemblies.** Where small plasmids were absent from the assembly, their reads
misaligned elsewhere and Medaka introduced over 100 erroneous changes. A missing
replicon is therefore not merely an omission: it corrupts the sequence that is
present. BacFlux's plasmid-recovery and replicon-audit steps run before polishing
partly for this reason.

## When this is the wrong choice, and how to override

The model assumes **native, unamplified DNA**. PCR or whole-genome amplification
erases methylation, so on amplified libraries the model is asked to correct signal
distortions that are not there.

If any library was amplified, pin the matched standard model explicitly instead of
using `auto` — for R10.4.1 `sup` v4.2.0 basecalls that is:

```yaml
parameters:
  nanopore:            # or: hybrid
    medaka_model: r1041_e82_400bps_sup_v4.2.0
```

`check_medaka_model` validates any explicit name against the installed Medaka's
model list and fails early with a suggestion table if it is wrong.

Library chemistry cannot be read reliably from FASTQ headers — the `sample_id`
field is free text typed by the operator, not authoritative protocol metadata. If
amplification status matters for a dataset, it has to come from the lab record.

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
- The error counts for *C. lari* and *E. cloacae*, and the small-plasmid pitfall —
  Wick (2024), cited above.

**Inference, not quoted from documentation**

- The biological mechanism (6mA/5mC/4mC, restriction–modification systems, signal
  deviation at modified bases, motif-systematic errors that coverage cannot
  average out). Standard domain knowledge, but not taken from the sources above.
- That Dorado modified-base calling leaves the canonical A/C/G/T sequence
  unchanged, and therefore that modbase basecalling would not remove the need for
  this model. Consistent with how Dorado is architected, but not verified here
  against Dorado's own documentation.
- **That amplified DNA is a contraindication.** ONT scopes the model to "native
  data" but does not explicitly warn against amplified input. The caution above is
  read off that scoping, not quoted from a warning.
- That methylation-induced errors are predominantly indels, and that frameshifted
  or pseudogene counts in the annotation are therefore a good readout for testing
  the model's benefit.

**Not consulted:** Dorado's own documentation; any primary paper on modified-base
effects on nanopore basecalling accuracy. A search also surfaced a BMC Genomics
study reporting that the bacterial methylation-aware model performed best among
those tested, but the full text was not read and it is not relied on above.
