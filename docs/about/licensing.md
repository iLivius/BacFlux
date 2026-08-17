# Licensing

BacFlux is MIT licensed. Some of the tools, model sets and databases it can
invoke are not, and this page is the one place their terms are written down.

It records what each licensor publishes, and draws no conclusion about whether a
particular use is permitted — that depends on who you are and what you are
doing, and it is yours to establish with the licensor.

## BacFlux's own code

[MIT](https://github.com/iLivius/BacFlux/blob/main/LICENSE), © 2024 Livio
Antonielli. Nothing you switch on below changes that.

The repository contains no third-party database, model set or sequence file. It
ships **code and URLs**. Every database a run needs is either one you already
hold and point the config at, or one the workflow fetches for you, on your
machine, from the publisher — the arrangement `bakta_db`, `gtdbtk_db`,
`platon_db` and the CARD link have always had, and the one every optional
mobilome layer uses too. See
[Databases](../getting-started/databases.md) for what that means in practice.

## Data and code are two different problems

| | A database or model set | Source code |
|---|---|---|
| **What BacFlux distributes** | a URL in `config/config.yaml` | the file itself, in the repository |
| **Who obtains it** | you, under your own agreement with the licensor | anyone who clones the repository |
| **Does a config switch settle it?** | yes — nothing is shipped in the first place | **no** |

!!! note "Why default-off settles one and not the other"

    A config switch controls what *runs*. It does not control what is
    *distributed*. Share-alike terms — the "if you pass it on, pass it on under
    the same terms" clause in CC BY-NC-SA — attach to source that ships in the
    repository, so a CC BY-NC-SA file would carry its terms into BacFlux whether
    or not any switch ever ran it. The only two remedies are to not copy it, or
    to relicense the whole workflow. This project took the first.

    Running such a tool is a different act from redistributing its source, and
    does not affect BacFlux's own licence.

## The default path

VirSorter2 and not geNomad is the default virus caller, and the reason is the
licence rather than the science
(`workflow/rules/shared/00_common.smk:104-111`). The licence types a default run
depends on:

| Component | Licence type |
|---|---|
| VirSorter2 — default virus caller | GPLv2 |
| Platon — plasmid caller, always on | GPLv3 |
| CheckV — grades whatever the virus caller found | LBNL BSD |

Those three are the components this repository states a licence type for. Every
other tool BacFlux invokes — on the default path or inside the mobilome module —
is distributed under its own licence by its own authors;
[CITATIONS.md](https://github.com/iLivius/BacFlux/blob/main/CITATIONS.md) lists
them all, and each tool's repository carries the licence itself. The databases
every run needs — Bakta, GTDB, CARD and the rest — are fetched the same way,
under their publishers' terms;
[Databases](../getting-started/databases.md) says what to download and from where.

## Components you switch on yourself

Each of these is off in the shipped configuration and stays off until you set
the key in the middle column.

| Component | Turned on by | Licence as published | Read from |
|---|---|---|---|
| **geNomad** — alternative virus caller, and the second opinion in the plasmid stage | `phage.caller: genomad` | Berkeley Lab licence, **academic / non-commercial use only** — *"User must be an accredited academic institution"*; commercial use requires a separate LBNL licence | the raw `LICENSE` file, checked 2026-07-22 |
| **CONJScan models** — the HMM profiles and system definitions that detect conjugation machinery | `mobilome.run: true` (fetched whenever the module runs) | **CC BY-NC-SA 4.0**, Institut Pasteur / CNRS | the package's own `metadata.yml` |
| **MacSyFinder** — the engine those models run in | same | GPLv3 | — |
| **ICEscan models** — the optional second model set, a fork of CONJScan distributed inside ICEfinder2's bundle | `mobilome.icescan.run: true` | **CC BY-NC-SA 4.0**, Institut Pasteur / CNRS — the same terms as the models it forks | the package's own `LICENSE` file |
| **TnCentral** — names transposons and integrons | `mobilome.tncentral.url` or `.dir` | an **"All Rights Reserved"** notice, and **no terms page at all** | the publisher's site, checked 2026-07-28 |
| **ICEberg 3.0** — names the ICE and IME candidates | `mobilome.iceberg.urls` or `.dir` | **no licence, terms or reuse statement anywhere**; every page carries only *"Copyright © 2023 All Rights Reserved by Microbial Bioinformatics Group in MML, SJTU."* | the publisher's site, checked 2026-07-28 |
| **ISOSDB** — counts the IS copies the assembly collapsed | `mobilome.isosdb.fasta_url` + `family_map_url`, or `.dir` | **MIT**, in the pseudoR repository | the repository |

Two things here are easy to get backwards. **MacSyFinder is not
CONJScan**: the engine is GPLv3, the models it loads are CC BY-NC-SA 4.0, and
it is the models that BacFlux fetches. And **geNomad's bioconda recipe is tagged
`BSD-4-Clause`, which does not match the raw `LICENSE` file** — the file is what
was read for the row above.

!!! warning "No published terms is not the same as permissive"

    TnCentral and ICEberg publish no terms of use. That is an absence of a
    statement, not a grant. Both layers are opt-in, neither is fetched unless
    you configure a source, and BacFlux redistributes neither.

## ISfinder is not wired in

ISfinder is the reference catalogue of bacterial insertion sequences. Its terms are
stricter than any other resource BacFlux touches:

> *"It is not permitted to download the ISfinder database without written
> authorization. Moreover, it is also not permitted at any time to distribute the
> database to third parties either individually or as part of web-based software."*
>
> — ISfinder terms of use, read 2026-08-17

**So BacFlux never offers you a way to fetch it.** Every other optional database in the
mobilome module works the same way: the workflow ships a URL, you download under the
publisher's terms. For ISfinder there is no URL anywhere in the configuration, because
supplying one would route around an authorisation only you can obtain.

**This does not cost you anything.** ISEScan, which is how BacFlux actually finds
insertion sequences, uses its own bundled profile HMMs and needs no ISfinder data.

**Where it could still reach you.** TnCentral serves several endpoints, and two of them
bundle ISfinder content. Only `nc/tn`, the transposon set on its own, appears in the
configuration; `tn_in_is` and `prot/tn_is` are named in the config comments **solely to
record that they are deliberately left out** (`config/config.yaml`, the comment above `mobilome.tncentral`). Point
`mobilome.tncentral.url` at one of those two and ISfinder's terms become yours to
satisfy — which is a decision to make deliberately, not one to make by copying a URL out
of a comment.

## What a run tells you

The workflow says this out loud rather than leaving it in a document:

- Selecting geNomad prints a notice at parse time, before any job starts
  (`workflow/rules/shared/00_common.smk:123-127`).
- `conjscan_models` and `icescan_models` write the CC BY-NC-SA notice to the top
  of their own logs
  (`workflow/rules/shared/80_mobilome.smk:596`, `:730`).
- Every opt-in download comes from an unversioned address, so each writes a
  `PROVENANCE.txt` beside the data. `tncentral_db`, `iceberg_db` and `isosdb_db`
  record the source, the fetch date, the checksum of the FASTA that actually
  arrived, the sequence count and the licence statement for that source;
  `icescan_models` records the package's own version line and profile count in
  place of a checksum. That file is the only record of which release you
  analysed — see [Optional layers](../mobilome/optional-layers.md).

## What this repository deliberately does not contain

**No CC BY-NC-SA code.** EBI's `mobilome-annotation-pipeline` is Apache-2.0, but
three of its scripts — `bin/ice_boundary_refinement.py`,
`bin/map_tools/icefinder_process.py` and `bin/prescan_to_fasta.py` — derive from
ICEfinder2 and carry CC BY-NC-SA headers. None of it is copied here. What was
taken from that pipeline is **decisions**: that ICEscan exists and is worth
running, the Sequence Ontology terms for the mobilome GFF, the
`contig_id|mge_type-start:end` element ID format, the discard-with-reason file
pattern, and the rule that a gene counts as inside an element only when at least
0.9 of it lies within the boundaries. Those are facts and design choices, not
expression. Reading source to understand an algorithm and then writing your own
is not copying, and that is the line kept here.

Consuming that pipeline's `mobilome.gff.gz` **output** carries no licence
consequence, which is why it is offered as an escape hatch for deeper mobilome
work.

**No vmatch.** The *att*-site search computes maximal exact repeats in plain
Python. bioconda's vmatch package declares `license: Unknown / OTHER` and
vmatch.de is unreachable, so its terms cannot be established — which is weaker
ground than a known-restrictive licence. That the substitute returns the same
repeats was measured rather than assumed; see
[Validation](../mobilome/validation.md).

## Citing

Licensing and citation are separate obligations: a permissive licence still
expects the citation. [Citation](citation.md) covers how to cite BacFlux and the
tools that do the work.

---

!!! note "Standing caveat"

    This page records what each licensor published, as read on the dates given.
    Licences change and pages move. Before relying on any of it, read the current
    terms at the source. Nothing here is legal advice.
