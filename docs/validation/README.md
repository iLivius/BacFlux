# Example mobilome output — read this before using these files

These are **real output files from one real run**, kept so you can see the shape
of the mobilome module's tables without running it yourself. They are examples of
*format*, not a current statement of what the module calls.

## What genome this is

***Klebsiella pneumoniae* ATCC BAA-2146** (`GCF_000364385.3`) — the "NDM-1
superbug": one chromosome (`CP006659.2`) plus four plasmids, pan-resistant, and
described in the published literature, which is what makes it usable as a
positive control. The walkthrough is in
[`../mobilome_worked_example.md`](../mobilome_worked_example.md).

> ⚠ **The `sample` column inside these files says `KPNIH1`. That is wrong, and it
> has been left as it is on purpose.**
>
> `CP006659.2` is ATCC BAA-2146. **KPNIH1 is a different genome, `CP008827.1`** —
> another carbapenem-resistant ST258 clinical isolate, with different resistance
> content and its own curated ICEberg element (`ICEKpnKPNIH1-1`). The run that
> produced these files was simply configured with the wrong sample name.
>
> The file *names* have been corrected to `BAA-2146_*`. The file *contents* have
> not been edited, because they are the literal output of that run and rewriting
> them would hide the mistake rather than record it. Treat the `sample` column
> here as a label typed by a human, which is all it ever is.

## When these were produced, and what has changed since

Produced **2026-07-27**, several commits before `4a93d89`. Two kinds of drift
have accumulated since, and both matter if you compare these files against your
own output.

### 1. The schema has grown

| File | Columns here | Columns now | Added since |
|---|---:|---:|---|
| `BAA-2146_ice_candidates.tsv` | 44 | **51** | `mobility_tier`, `mobility_tier_reason`, `evidence_sources`, `icescan_model_class`, `evidence_level`, `missing_components`, `machinery_gap_bp` |
| `BAA-2146_amr_mobility.tsv` | 44 | **46** | `named_element`, `named_element_type` |

Nothing was removed or renamed, so every column here still exists — the files are
short, not wrong, in this respect.

### 2. The calls themselves have changed, and one changed class

This is the important one. On the chromosome, `BAA-2146_ice_candidates.tsv`
reports:

```
NZ_CP006659.2|conjugative_region-4610740:4644558   33,819 bp   boundary_method=none
```

The module **no longer calls that locus a `conjugative_region`**. Re-run at
`4a93d89`, the same locus comes out as:

```
CP006659.2|ice-4604073:4644558   4,603,807–4,658,749   54,943 bp
  boundary_method = tRNA (43 bp repeat at tRNA-Phe(gaa))
  anchor_classes  = integrase,relaxase,t4cp,t4ss      confidence = high
```

Two separate fixes are responsible. An integrase Bakta annotates as
`DNA integration/recombination/inversion protein` was being missed by the product
regex, which is what demoted the element from `ice` to `conjugative_region`; and
the reworked *att* search then found the tRNA-anchored repeat that fixes its ends.
Against ICEberg's curated `ICEKpnATCCBAA-2146-1` (58,048 bp) the current call
recovers **0.946** of the element, versus roughly 0.70 for the interval shown in
this file.

The current run also finds a **third** chromosomal element these files do not
contain, `CP006659.2|ime-589955:596310`.

**So: do not read `conjugative_region` here as the module's verdict on this
locus.** Shipping an example where the module declines to call an ICE, when it
now calls one confidently, would mislead — which is the whole reason this note
exists.

### What has *not* changed

`BAA-2146_replicon_calls.tsv`, `BAA-2146_plasmid_concordance.tsv` and
`BAA-2146_contig_taxonomy_decisions.tsv` are replicon- and contig-level outputs
that the ICE and *att* work did not touch. They remain representative.

## A caveat on the AMR table

`BAA-2146_amr_mobility.tsv` has 66 rows reaching tiers 1, 2, 3, 5 and 6 — but
**every `mge_name` is `NA`**, because this run did not have the optional TnCentral
naming layer switched on. Tier 4 ("inside a *named* transposon or integron") is
therefore not exercised anywhere in this example. See the README's *What the
benchmark does not show* for why that gap is worth knowing about.

## Files

| File | What it is |
|---|---|
| `BAA-2146_amr_mobility.tsv` | The deliverable: one row per AMR/stress gene, with its mobile-element context, mobility tier and confidence |
| `BAA-2146_ice_candidates.tsv` | ICE / IME / conjugative-region calls — **see the correction above** |
| `BAA-2146_replicon_calls.tsv` | Chromosome or plasmid per contig, and plasmid mobility class |
| `BAA-2146_plasmid_concordance.tsv` | Platon versus geNomad, where both were run |
| `BAA-2146_contig_taxonomy_decisions.tsv` | The decontamination audit: every contig kept or dropped, with its reason |
| `baa2146_positive_control_config.yaml` | The config this run used |

*Licensing note:* this run had `mobilome.run: true` and `phage.caller: genomad`,
both of which pull in non-commercially-licensed components (CONJscan models,
CC BY-NC-SA 4.0; geNomad, academic use only). See
[`../mobilome_worked_example.md`](../mobilome_worked_example.md) and
[`../../CITATIONS.md`](../../CITATIONS.md).
