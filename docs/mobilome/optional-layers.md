# Optional layers

`mobilome.run: true` is the whole switch for the module itself
([Turning it on](enabling.md)). It is **not** the switch for the four layers
described here. Each of them stays off until you give it a source of its own, and
turning the module on changes none of them.

## The four layers

| Layer | What it adds | Switched on by | Modes | Rules |
|---|---|---|---|:--:|
| **ICEscan models** | a second MacSyFinder model set run beside CONJscan: the AICE class, roughly twice as many IMEs, and integrase profiles CONJscan does not carry | `mobilome.icescan.run: true` | all | 2 |
| **TnCentral** | curated transposon and integron **names** — the only thing that makes ladder tier 4 reachable | `mobilome.tncentral.url` or `.dir` | all | 3 |
| **ICEberg 3.0** | curated **names** for the elements the ICE/IME caller already found. Changes no gene's tier | `mobilome.iceberg.urls` or `.dir` | all | 3 |
| **ISOSDB copy number** | how many IS copies the assembler collapsed, from read depth. Changes no gene's tier | `mobilome.isosdb.fasta_url` + `family_map_url`, or `.dir` | `illumina`, `hybrid` | 4 |

Twelve rules in total, against the ten the module runs on its own. Every layer is
independent: any combination is valid, and a layer you leave alone is skipped in
full — its rules are never defined, nothing is downloaded, and none of its outputs
appear.

!!! note "Where the data comes from"

    All four fetch third-party data at run time. BacFlux ships none of it: the
    workflow distributes a URL and you download under your own agreement with the
    licensor, the same arrangement `bakta_db`, `gtdbtk_db` and the CARD link have
    always had. Licence types are listed once, on
    [Licensing](../about/licensing.md).

## Half-configured is an error, not a silent skip

The gating is deliberately blunt: a layer is on when it has a source and off when
it does not. Two configurations are ambiguous rather than off, and both stop the
run at parse time — before any job starts — rather than three hours in:

```text
[BacFlux] mobilome.icescan.run is true but neither mobilome.icescan.url nor
mobilome.icescan.dir is set, so there is no model package to use.
```

```text
[BacFlux] mobilome.isosdb.fasta_url is set but mobilome.isosdb.family_map_url is
empty, so the IS family map cannot be downloaded and every element would be
reported as 'unassigned'.
```

Each message goes on to name the ways out: set the missing key, point `dir:` at a
copy you already hold, or switch the layer off. The reasoning is the same in both
cases — asking for a layer and getting today's results back, with nothing to say
the layer never existed, is worse than a stop.

## ICEscan models

A second set of machinery models, run over the same Bakta proteins as CONJscan.
The two hit tables are merged before the ICE/IME caller sees them.

It is not a second opinion. ICEscan is a **fork of CONJScan** by the same
Institut Pasteur authors — its own `metadata.yml` still calls itself CONJScan —
and it is one minor version behind the release BacFlux installs (2.0.1 against
2.1.0). It is therefore added *alongside* CONJscan and never instead of it: the
fork drops the MOB relaxase models, the decayed-machinery models and the whole
plasmid set, so a swap would lose elements BacFlux currently detects.

**What it adds.** Two system models CONJscan 2.1.0 does not carry — `IME` and
`AICE` — and 21 profile HMMs: eight integrase families, the Gram-positive and
IME relaxases, and the AICE machinery. The integrases are the substantive part,
because CONJscan carries none at all: conjugation is its subject and integration
is not. Four of the eight are used as element anchors; the other four are
integron or chromosomal housekeeping recombinases, deliberately ignored because
they would anchor an element onto the wrong gene.

The measured effect is on IMEs: over the 28-genome benchmark, 10 `ime` rows
without the layer and 21 with it. Turning it off roughly halves the IME rows; it
does not remove them, because BacFlux's own rules call an IME whenever the anchors
are an integrase plus a relaxase with no mating apparatus, and CONJscan's MOB
models supply plenty of those. `aice` really is ICEscan-only. The rest of the
numbers, including what the layer does *not* improve, are in
[Validation](validation.md).

**What it does not add: boundaries.** MacSyFinder reports gene ordinals — "the
41st protein in the file" — never base pairs, so every coordinate BacFlux prints
still comes from its own Bakta GFF3 join, its own clustering and its own *att*
search. ICEscan's own spans are never used as element edges.

```yaml
mobilome:
  icescan:
    run: true        # the only key you must change — a URL ships filled in
    sha256: ""       # optional, recommended: pins which release you analysed with
    dir: ""          # optional: a MacSyFinder models directory holding ICEscan/
```

This is the one layer whose `url` is already populated in the shipped config (the
ICEfinder2 database bundle, ~61 MB, of which only `macsydata/ICEscan` is kept and
the rest discarded). So `run: true` is genuinely all it takes.

!!! note "Both machinery searches run at the same stringency, and must"

    CONJscan and ICEscan are both given `mobilome.coverage_profile` (default
    `0.5`, MacSyFinder's own default). Their hit tables are merged, so different
    stringencies would make an element's class depend on which model set happened
    to be more permissive — untanglable downstream. What lowering it buys and
    costs, measured, is on [Tuning](tuning.md).

**Outputs.** `08.mobilome/icescan_models/` once per run, and
`08.mobilome/{sample}/icescan/` per sample. No new column appears in the mobility
table: the layer works by changing what the ICE/IME caller can see.

## TnCentral

Curated transposons and integrons, matched by `blastn` against the whole
assembly, then merged into element copies and named.

Tier 3 is an **inference**: two IS copies of one family, the right distance
apart, a gene between them, so we call it a composite transposon. The pattern
only stays right because of hand-written exceptions to it — IS*26* flanks its
cargo in *direct* orientation, so it is exempted from the same-orientation test
every other family has to pass, and without that exemption the most important
architecture in clinical AMR would be the one BacFlux missed. A TnCentral hit
needs no such exception. It matches an element somebody characterised, named and
deposited, whose architecture is already known. Where both fire on one gene the
curated hit wins, and the composite call it displaced is written to the audit.
This is the only route to [tier 4](mobility-ladder.md).

```yaml
mobilome:
  tncentral:
    url: "https://tncentral.ncc.unesp.br/api/download_blast/nc/tn"
    sha256: ""                    # optional, recommended — the endpoint is unversioned
    dir: ""                       # optional: a directory already holding tncentral.fa
    min_identity: 90.0            # percent identity before a hit may confer a name
    min_reference_coverage: 0.8   # fraction of the CURATED element that must be present
```

`min_reference_coverage` is measured against the reference, not against your
contig, and that is the point: a 7 kb transposon inside a 300 kb contig covers 2%
of the contig and 100% of itself. The question is whether the whole of a known
element is here — a fragment of a transposon is not that transposon. On a draft
assembly a low value usually means the element is split across contigs. Both
thresholds are naming conventions rather than biological boundaries, and every
hit refused a name reaches the discard audit with the measured identity and
coverage that refused it.

The `nc/tn` endpoint above is the transposon set on its own. The variants that
also serve ISfinder content bring ISfinder's terms with them — written
authorisation to download, no redistribution — so they are deliberately not wired
in.

!!! warning "Tier 4 is a working code path, not a measured one"

    Reaching it needs a TnCentral source **and** at least `min_reference_coverage`
    of a curated element present. On the clinical runs that had the layer on,
    every Tn*Ecp1.1* candidate was refused at 12% and 49% coverage of a 3,417 bp
    reference. **No run kept on disk contains a tier-4 row.** Treat the layer as
    what it is — the naming cascade works, and this particular outcome has not
    been observed. See [Validation](validation.md).

**Outputs.** `08.mobilome/tncentral_db/` once per run;
`{sample}_tncentral_blast.tsv`, `{sample}_named_elements.tsv` and
`{sample}_named_elements_discarded.tsv` per sample. The named elements are handed
to the co-localisation step as a third element table, alongside the IS calls and
the ICE candidates.

## ICEberg 3.0

Curated names for the candidates the ICE/IME caller already found — ICEs, IMEs,
AICEs and genomic islands. A conjugative region never takes a name, because it
has no integrase and calling it an ICE is exactly what the classifier refused to
do.

This layer **labels**; it never decides. The caller has already said what is an
ICE and what class it is, and turning this on cannot change any gene's tier. What
it changes is "predicted self-transmissible element" into a name like
`ICEKpnATCCBAA-2146-1`, which is what lets you go and read about the thing.

```yaml
mobilome:
  iceberg:
    urls:
      - https://tool2-mml.sjtu.edu.cn/ICEberg3/data/download/ICE_seq_all.fas
      - https://tool2-mml.sjtu.edu.cn/ICEberg3/data/download/IME_seq_all.fas
    dir: ""                     # optional: a directory of .fas files you already hold
    min_identity: 80.0          # identity floor for a name
    min_overlap_fraction: 0.5   # how much of OUR candidate the curated element must cover
```

Note the host: the older `bioinfo-mml.sjtu.edu.cn` path 404s for ICEberg 3.0,
because it still serves ICEberg 2.0 — which is exactly why a stale URL kept
looking plausible. `ICE_seq_all.fas` is around 100 MB and the server is slow, so
the download rule is allowed to resume a part-finished transfer.

`min_overlap_fraction` is lenient next to the transposon cascade on purpose. ICEs
are mosaic and their cargo varies between strains, so demanding near-completeness
would refuse to name exactly the divergent elements a name helps with most. A
separate rule handles the other direction: when less than **80%** of the
*reference* is present, the name is suffixed `-like`, because we have part of a
known element rather than the whole of it. That 80% is not a config key — it is
`EXACT_NAME_REFERENCE_COVERAGE` in
`workflow/scripts/80_mobilome/name_ice_elements.py`, named here so there is
something to grep for.

!!! note "The real reason to switch this layer on"

    It measures how far the module's boundaries fall short. On the
    *K. pneumoniae* ATCC BAA-2146 positive control the ICE call spans 54,943 bp
    against ICEberg's 58,048 bp for the same element: 0.946 of it recovered,
    stopping 3,138 bp inside its far end, and the naming audit says so in as many
    words. Any AMR gene in that last 3 kb is scored as though it were outside the
    ICE. That is the good case — a closed genome where the *att* search found a
    tRNA-anchored repeat. With `boundary_method = none` the interval is only the
    machinery span and the shortfall is far larger
    ([Draft assemblies](draft-assemblies.md)).

    A second caution that comes with any name from this layer: ICEs of one
    species are near-identical across strains, so one real element matches dozens
    of curated entries. The name is a group label, not a unique identification.

**Outputs.** `08.mobilome/iceberg_db/` once per run;
`{sample}_iceberg_blast.tsv`, `{sample}_ice_candidates_named.tsv` and
`{sample}_ice_naming.tsv` per sample. Whenever this layer is on, the
co-localisation step reads the named candidate table instead of the unnamed one:
the same rows, `mge_name` filled in, and six columns added carrying the evidence
behind each name — identity, how much of our interval the curated element covers,
how much of the curated element is present, and how many other entries fit about
as well.

## ISOSDB copy number

The module says repeatedly that on a fragmented assembly the located IS count is
a **floor, not a count**. True, and unquantified: you cannot tell whether the
floor is 2 short or 40 short. This layer puts a number on it.

Reads are immune to assembly collapse — every copy of an IS contributes its own
reads whether or not the assembler kept them apart — so an IS present in five
copies attracts about five times the read depth of the single-copy chromosome:

```text
copy number  ≈  depth over the IS  /  depth over the genome
```

Both depths are measured with BBMap at identical settings, one against ISOSDB and
one against the sample's own assembly, so their ratio means something.

```yaml
mobilome:
  isosdb:
    fasta_url: "https://raw.githubusercontent.com/joshuakirsch/pseudoR/main/ISOSDB.V3.fna.zip"
    family_map_url: "https://raw.githubusercontent.com/joshuakirsch/pseudoR/main/IS_fam_annot.txt"
    dir: ""                    # optional: a directory already holding both files
    min_covered_percent: 90.0  # a database entry must be covered end to end to be believed
    min_copy_number: 0.5       # below this multiple of the genome baseline, treated as absent
```

Set both URLs or neither. A partly covered database entry is usually a conserved
domain shared with another family, and averaging it in would inflate every
estimate, which is what `min_covered_percent` refuses. `min_copy_number` sits
below 1.0 on purpose: a real single-copy IS lands near 1×, and sampling noise
plus mapping loss routinely drag it to 0.6–0.8×.

`illumina` and `hybrid` only. In `nanopore` and `contigs` mode there are no short
reads to map and the layer is simply absent, whatever the config says.

**Outputs.** `08.mobilome/isosdb_db/` once per run;
`{sample}_isosdb_covstats.tsv`, `{sample}_assembly_covstats.tsv`,
`{sample}_is_copy_number.tsv` and `{sample}_is_copy_number_audit.tsv` per sample.
The deliverable is one row per IS family:

| Column | What it is |
|---|---|
| `located_copies` | copies ISEScan found on the contigs |
| `estimated_copies` | copies the read depth implies |
| `collapse_delta` | `estimated − located` — the assembler's collapse, measured |
| `db_informative` | `TRUE` only when ISOSDB detected the family *and* the estimate reached the located count |
| `n_db_entries_detected`, `max_entry_copy_number`, `genome_baseline_depth` | the working behind the estimate |

Reported per **family** rather than per entry because ISOSDB is dereplicated at
95% identity and families remain similar: the family sum is robust to which
near-identical entry an ambiguous read happened to land on, and the per-entry
split is not.

!!! warning "This table changes no AMR gene's tier, and must not"

    It says how many copies exist, never *where* they are, so it cannot place a
    gene inside anything. Nothing downstream reads it — it is a quality metric on
    the IS inventory, read by a person, and the honest companion to the "the
    count is a floor" warning.

## `PROVENANCE.txt`

Three of these four addresses carry **no version at all**, and the fourth is a
GitHub `main` branch that can move under you. Fetch `nc/tn` today and again in six
months and you may get different data from the same URL, with nothing in the file
to say which is which — and no way, later, to say which release a result came
from.

So each download rule writes a `PROVENANCE.txt` beside the data it fetched,
recording what actually arrived:

```text
source:      https://tncentral.ncc.unesp.br/api/download_blast/nc/tn (sha256 8f3c…)
fetched:     2026-08-16T09:41:07Z
fasta_sha256: 5b90…
sequences:   533
```

| File | Written by | Records |
|---|---|---|
| `08.mobilome/icescan_models/PROVENANCE.txt` | `icescan_models` | source, fetch date, the package's own version line, number of profiles |
| `08.mobilome/tncentral_db/PROVENANCE.txt` | `tncentral_db` | source, fetch date, checksum of the repaired FASTA, sequence count |
| `08.mobilome/iceberg_db/PROVENANCE.txt` | `iceberg_db` | source, fetch date, checksum of the concatenated FASTA, sequence count |
| `08.mobilome/isosdb_db/PROVENANCE.txt` | `isosdb_db` | source, fetch date, checksum of `ISOSDB.V3.fna`, sequence count |

Under those fields each file carries a few lines on why it exists and on the
terms the data came under.

Two practical uses. **Keep the file with the results** — it is the line a methods
section needs, and it is the only record that exists. And **check it when two runs
disagree**: a different checksum on the same URL is the first thing to rule out.

`conjscan_models`, the module's always-on model fetch, writes none. It installs a
*versioned* package through `msf_data install`, so the version is already
recoverable.

!!! note "Only two of the four can be pinned in advance"

    `mobilome.icescan.sha256` and `mobilome.tncentral.sha256` are checked against
    the download and fail the rule on a mismatch, with a message saying the
    upstream file changed. Set both once you have decided which release you are
    working with.

    ICEberg and ISOSDB have **no `sha256` key**. ICEberg serves two files rather
    than one archive and the large one has to be resumable over a slow link, so a
    single pinned digest is the wrong instrument; ISOSDB comes from GitHub raw
    paths, where the source repository's own history is the version record. For
    both, the observed checksum is recorded after the fact — you can say later
    what you used, but you cannot demand it up front.

## What the download rules do besides downloading

- **They refuse an error page.** A bot block or a redirect arrives with a 200 and
  would otherwise be indexed as an empty database that silently names nothing —
  a wrong answer that looks exactly like a real one. Each rule tests the archive
  or the FASTA and fails loudly instead, printing the first bytes it received.
  TnCentral's server rejects curl's default user-agent, so the rule identifies as
  a browser; if that stops working the rule fails rather than producing an empty
  database.
- **They repair the FASTA before indexing it.** Both curated nucleotide sets ship
  with deflines glued onto the *end* of a sequence line instead of starting their
  own — 21 of 533 records in the TnCentral release checked, 1 of 1,774 in
  ICEberg's. The damage runs both ways: those records are invisible to
  `makeblastdb` (for TnCentral that included Tn*7* itself and eleven integrons,
  the very class tier 4 exists to name), and the records they were glued to
  become chimeric. `In_Tn6162` measured 41,492 bp instead of 8,911 — 4.7× its
  real length — and since the naming cascade tests coverage as
  `alignment / reference length`, an inflated reference makes an element almost
  impossible to name, silently. After repair, each rule asserts that every `>`
  begins a line and stops if it does not.
- **They rebuild the BLAST index as v5.** TnCentral ships a v4 index built from
  the *unrepaired* FASTA, so it is missing the same 21 records. It is deleted
  rather than reused.

## Practical notes

**The databases live in the output directory**, under `08.mobilome/`, and are
fetched once and shared by every sample in that run. A second project therefore
downloads them again. That is what the `dir:` key on every layer is for: point it
at a copy you already hold and nothing is fetched. `dir:` always takes precedence
over the URLs, exactly as `directories:` beats `links:` elsewhere in the config.

**Network and environments.** `icescan_models` uses `wget` and `tar` from the
environment you launched Snakemake in, with no conda environment of its own — the
same arrangement as the CARD download ([Installation](../getting-started/installation.md)).
The other three download rules, and the BLAST searches that follow them, share one
small conda environment carrying BLAST, `curl` and `unzip`, built on first use.
The ISOSDB read mapping uses the same BBMap environment the CARD leg already
builds.

**Do not mix runs.** Switching a layer on or off mid-project changes what the
mobility table can report, and nothing in the table says which layers were on when
it was written — only the `PROVENANCE.txt` of each layer that ran. Re-run
everything rather than comparing a half-and-half result set; a mixed one is very
hard to unpick later.

## Related pages

| Page | Read it for |
|---|---|
| [Turning it on](enabling.md) | the master switch, and what the module does on the default configuration |
| [Reading the output](output.md) | what these layers change in `{sample}_amr_mobility.tsv`, column by column |
| [Tuning](tuning.md) | the naming and coverage thresholds, with the measurements behind them |
| [Validation](validation.md) | what each layer was measured to add, and what it did not |
| [Licensing](../about/licensing.md) | the terms attached to each of the four data sources |
| [Configuration](../reference/configuration.md) | every key in the `mobilome` block, in the file's own order |
