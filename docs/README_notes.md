# README notes — pending clarifications

Working list of things discussed and agreed that still need to be written into the real
`README.md` (or the consolidated v2.0.0 README, per `unification_migration_plan.md`), but
haven't been yet. This file is a to-do list, not documentation — nothing here should be
treated as the final wording. Delete each entry once it's actually folded into the README.

---

## 1. Snakemake launch flags: `--cores` / `--resources cpus=N` / `--jobs` — urgent

> **RESOLVED 2026-07-22 — the underlying cause is FIXED. Read this before writing the README.**
> Every cpu-bound rule now declares Snakemake's built-in `threads:` instead of the custom
> `resources: cpus`. Thread NUMBERS are unchanged (`capped_cpus(N)` returns the same value);
> only the mechanism changed, and the built-in one is enforced by `--cores` automatically.
>
> **So the README must document the SIMPLE form:**
> ```bash
> snakemake --sdm conda --cores 24 --configfile config/config_v2.yaml
> ```
> `--cores N` alone. No `--resources cpus=N`. Still never `--jobs`/`-j`.
>
> Verified after the change: with `--cores 24` the per-rule caps are preserved
> (`virsorter2_db`=4, ONT QC/Medaka=8, `trim_adapters`=16, the rest 24); with `--cores 8`
> Snakemake caps every request at 8 on its own, except `virsorter2_db` whose lower cap of 4
> correctly wins. That is the behaviour the old `resources: cpus` form could never give.
>
> Everything below this line is the HISTORICAL diagnosis, kept because the `--jobs` half of it
> is still true and still needs saying. **Do not copy the `--resources cpus=N` workaround into
> the README** — it is obsolete.
>
> One upgrade note worth a line in the release notes: because Snakemake hashes rule definitions,
> this change makes it want to re-run outputs produced by a PRE-conversion v2. Irrelevant for
> released v1.3.1 users; it only affects anyone who ran a v2 development build.

**Status: the CURRENT README's documented launch commands are actively wrong, not just
under-explained.**

### The problem

`README.md`, "Running BacFlux" section, currently tells every user to launch with:

```bash
snakemake --sdm conda --jobs 4 --cores 12
```

and FastaFlux with:

```bash
snakemake --sdm conda --snakefile workflow/FastaFlux --jobs 2 --cores 12
```

Both combine `--jobs`/`-j` with `--cores`. We verified empirically (not just from the docs)
that for local execution this combination can silently allow real CPU oversubscription: two
rules each declaring `threads: 8`, launched with `--jobs 2 --cores 8`, both received the
full declared 8 threads and ran fully concurrently — 16 real threads against a declared
8-core budget. The identical rules launched with `--cores 8` alone (no `--jobs`) correctly
ran one at a time. Snakemake's own docs hint at this — `--jobs` is described as *"an alias
for `--cores`"* for local execution, recommending `--cores` instead — but don't spell out
that combining the two, rather than being merely redundant, can actively break the
accounting.

**Separately, and just as important:** every rule in BacFlux (checked — zero exceptions)
declares its CPU request via a custom `resources: cpus = ...`, not Snakemake's built-in
`threads:` directive. Snakemake does not automatically enforce a custom-named resource
against `--cores`. This is intentional, documented behavior — *"If no limits are given, the
resources are ignored in local execution"* — because `resources:` is a deliberately
open-ended mechanism for arbitrary user-defined constraints (many cluster-oriented
workflows use a `cpus`-like resource purely to forward a value into a SLURM/PBS submission
template, where the *cluster's* scheduler does the real enforcement — not a mistake there,
just a different execution mode than a local run). The practical upshot for a local BacFlux
run: `--cores N` alone gives **no real protection** against oversubscription, because
nothing in BacFlux speaks `--cores`'s native currency (`threads:`).

### What the corrected README needs to say

1. **Fix the launch commands.** Replace `--jobs 4 --cores 12` with:
   ```bash
   snakemake --sdm conda --cores 12 --resources cpus=12
   ```
   using whatever number `resources.threads` is set to in `config.yaml` (the shipped
   example is `24`, so a real run would read `--cores 24 --resources cpus=24`). Same fix
   for the FastaFlux command. **Never pass `--jobs`/`-j` alongside `--cores`.**

2. **Add a short, plain-language explanation**, near the `resources:` config subsection
   and/or "Running BacFlux":
   - `--cores N` is the one real total CPU budget for a local run.
   - `--resources cpus=N` (matching `resources.threads` in `config.yaml`) is what actually
     makes BacFlux respect that budget, because BacFlux's rules use a custom `cpus`
     resource name rather than Snakemake's built-in `threads:` mechanism.
   - Don't add `--jobs`/`-j` for local runs at all — it is not an independent "N parallel
     jobs, M cores each" setting (a natural but incorrect assumption); it's a synonym for
     `--cores`, and combining the two can silently reintroduce the oversubscription the
     ceiling exists to prevent.
   - Different rules already request different amounts of CPU automatically (e.g. BLAST
     vs. a lightweight step) — nothing needs configuring per rule for this; each rule caps
     itself via `min(CPUS, N)` against the one shared config value, which is the design
     already agreed for v2 (config sets overall machine scale; each rule encodes its own
     tool-appropriate cap).

3. ~~**Longer-term, not urgent**: switch to `threads:`.~~ **DONE 2026-07-22** — see the
   RESOLVED box at the top. `--cores` alone is now sufficient and safe. Item 1 of this
   file is therefore reduced to a single instruction: document `--cores N`, and warn
   against `--jobs`.

### Where this lands in the current README
- "Running BacFlux" (~line 334–352): fix both example launch commands.
- "Configuration" → `resources` subsection (~line 298–303): expand with the explanation
  above.

---

## 2. Sample-name character policy — v1→v2 behavior change, needs stating plainly

**Status: informational, not urgent — but v1's README documented a stricter rule than v2
actually enforces, so leaving it unstated invites confusion.**

The current README (`input_dir` bullet, "Configuration") says sample names *"cannot contain
underscores"* — that was true for three of the four v1 workflows, but v2 standardizes on
**allowing** underscore everywhere (see `00_common.smk`'s `BAD_CHARS`), because it's common
in real sample names and already appears inside the `{sample}_R1`/`{sample}_ont` discovery
pattern itself.

The exact v2 rule, plain and complete:
- **Rejected characters** (anywhere in a sample name): `* # @ % ^ / ! (space) ? & : ; | < >`.
  A name containing any of these is refused at parse time, before anything runs.
- **Everything else is accepted** — including underscore (`_`) and dot (`.`). A dotted name
  like `my.genome.fasta` in contigs mode is read as sample `my.genome` (the *last* dot
  splits off the file extension; earlier dots stay part of the sample name), and processes
  normally — filenames, locus tags, and report labels all handle it safely.

### What the corrected README needs to say
Replace the "cannot contain underscores" line with the actual v2 rule above — ideally as an
explicit list of disallowed characters (so a user can check their own sample names against
it directly) rather than calling out only one or two examples. Note it applies identically
to all four modes (illumina/nanopore/hybrid/contigs), since the check lives once in
`00_common.smk`, not per front end.

### Where this lands in the current README
- "Configuration" → `directories` → `input_dir` bullet (~line 275): replace the "Strain
  names cannot contain underscores" sentence.

---

## 3. Download rules assume `wget` (and tar/sha256sum) in the launch environment — deferred decision

**Status: informational + one deferred design decision. Not urgent.**

Several rules that fetch reference databases run WITHOUT their own conda env, i.e. in
whatever environment Snakemake was launched from, and call `wget` / `tar` / `sha256sum` /
`awk` directly. Confirmed in Stage 2a for `cazyme_db_download` (dbCAN); the same pattern
applies to the front-end download rules still to be ported (`download_phix`,
`download_amr_db`/CARD, `viral_db`). This is faithful to v1, but `wget` in particular is not
part of a default miniconda/micromamba base env and is absent on some minimal/HPC login
shells — so a DB download can fail at run time with a bare `wget: command not found`, after
the DAG has already started.

**Decision to make (either is defensible — this is the tool author's call):**
- (a) Keep them env-less (v1-faithful) and DOCUMENT in the README that the launch/base
  environment must provide `wget`, `tar`, `sha256sum`, `awk`; or
- (b) Give the download rules a tiny pinned bioconda env (e.g. `workflow/envs/download.yaml`
  with `wget`) so their tools are pinned like every other rule, or switch `wget` → `curl -fSL`
  (more commonly present, `-f` makes HTTP errors fail the rule).

Until decided, Stage 2a keeps `cazyme_db_download` env-less per v1.

### Where this lands in the current README
- "Running BacFlux" or "Installation": if option (a), state the launch-env tool requirement.

---

## 6. The mobilome spans stages 06 + 07 + 08 — explain this, don't rename anything

**Status: DECIDED 2026-07-22 (option 3 — document, keep the layout as-is).**

Plasmids and prophages ARE mobile genetic elements, so the output layout can read as if
`08.mobilome/` were "the mobilome" and `06.plasmids/`/`07.phages/` were something else. They
are not. The mobilome is spread across three stages, split by *when it runs and what question
it answers*, not by whether the element is mobile:

| Stage | Question | Runs |
|---|---|---|
| `06.plasmids/`, `07.phages/` | *What replicons and prophages are in this genome?* (detection/inventory) | always — standard WGS characterisation, as in v1 |
| `08.mobilome/` | *Is each AMR gene embedded in a mobile element, and how transferable?* (integration/interpretation) + IS/transposon/integron/ICE detection | opt-in (`mobilome.run: false`), heavier, needs extra and partly licence-encumbered DBs |

The mobilome module **consumes** 06/07 rather than re-detecting them — see
`mobilome_module_SPEC.md` §7 (its co-localisation inputs include the "Platon/MOB-suite replicon
call"), §9 (its output table carries `replicon(chromosome|plasmid_id)` and lists `plasmid` as
one of the `mge_context` values), and §2.5 (mobility ladder tiers 5-6 are plasmid-based).

**README wording to add** (Output section): a short paragraph making the above explicit — that
06/07 are the always-on MGE *detection* stages, 08 is the opt-in *interpretation* layer that
builds on them, and together they are the mobilome. Renaming `08.mobilome` was considered and
deliberately rejected in favour of documenting it.

---

## 5. Phage/plasmid: VirSorter2 is default; geNomad is opt-in and NON-COMMERCIAL

**Status: DECIDED 2026-07-22. README must state the geNomad licence clearly.**

BacFlux v2's default virus/prophage caller is **VirSorter2** (GPLv2) and the default plasmid
deliverable is **Platon** (GPLv3) — both allow commercial use, as does **CheckV** (LBNL
permissive BSD, runs in every phage path). The default pipeline is therefore fully
commercial-usable.

**geNomad is opt-in only** (`config.phage.caller: genomad`), because it is licensed
**ACADEMIC / NON-COMMERCIAL-USE-ONLY** (Berkeley Lab — "User must be an accredited academic
institution"; commercial use reserved, separate LBNL licence required). Selecting it turns on
BOTH geNomad as the virus caller AND the Platon+geNomad plasmid concordance (D9). On the
default, geNomad never runs.

**README wording to add** (Configuration → `phage`, and/or a licensing note):
- State that geNomad is optional and academic/non-commercial-licensed; commercial users must
  either stay on the default (virsorter2) or obtain a commercial geNomad licence from LBNL.
- Note that enabling geNomad also enables the plasmid concordance; the default plasmid output
  is Platon's `verified_plasmids.txt`.
- Correct any earlier claim that geNomad is BSD-licensed (the bioconda tag is wrong).

---

## 4. dbCAN is deliberately pinned to v5.1.2 (+ the Zenodo-hosted DB) — reframe, don't apologize

**Status: DECIDED 2026-07-21 (evidence-based). Reframe the README wording; no code change.**

BacFlux pins `dbcan=5.1.2` and fetches its database from a version-pinned copy on the
author's Zenodo (`zenodo.org/records/18622157`). This is a **deliberate reproducibility
pin**, not a workaround — and it was validated by actually testing the current release
(5.2.9) against the 5.1.2 baseline on a real Bakta protein set:

- **The official run_dbcan AWS S3 mirror only hosts the CURRENT database** (`db_v5-2_...`,
  for 5.2.x). There is no official copy of the v5.1.2 DB — so moving off the personal
  Zenodo is inseparable from bumping the tool to 5.2.x.
- **Bumping to 5.2.9 is NOT drop-in:** the `substrate_prediction` subcommand changed its
  required arguments (now needs `--input_raw_data` and `--mode`); output schema shifted
  (extra `Substrate` column, renamed files); and `easy_substrate` writes to a different
  path layout. An unresolved empty-`CGC.faa`/PUL-blast degradation appeared in testing.
- **It changes results:** 5.2.9 called 205 CAZyme genes vs 315 in 5.1.2 on the identical
  input. The ~110 dropped genes were almost all weak single-tool (`dbCAN_sub`-only,
  non-recommended) hits — 5.2.9 correctly filters them via the e-value-threshold fix
  (5.2.7) that 5.1.2 was ignoring. So 5.2.9 is arguably *more* correct, but the numbers
  differ materially, and high-confidence (2–3 tool) calls are preserved either way.

**README wording to add** (Configuration / Installation, near the dbCAN link):
- State that dbCAN is pinned to 5.1.2 for reproducibility, and that the DB is a
  version-pinned Zenodo copy (a persistent CERN-backed DOI), because run_dbcan's official
  DB server has been unreliable and its commands/output change between releases in ways
  that alter results.
- Note for power users: 5.1.2 over-reports weak `dbCAN_sub`-only hits (the e-value bug);
  the high-confidence calls are unaffected.
- Migration to 5.2.x is a deliberate future task (resolve the CGC issue, switch to
  `easy_substrate`, re-baseline counts, update parsing), not a config-URL swap.

---

## 7. `medaka_model` and `flye_input_mode` are COUPLED in auto mode — document this

**Status: real trap, confirmed against the BacFluxL baseline (2026-07-22).**

When `flye_input_mode: auto` (the default), the Flye read mode is chosen from the Medaka
model name: an explicit model whose name contains **"fast"** switches Flye to `--nano-raw`;
anything else gives `--nano-hq`. The logic is in `00_common.smk` (~line 949) and is inherited
unchanged from v1, so this is documentation-only, not a code change.

The consequence users will not expect: **setting `medaka_model` changes the ASSEMBLER, not
just the polisher.** The published BacFluxL baseline is itself an example — its
`config_custom.yaml` has

```yaml
flye_input_mode: auto
medaka_model: r941_min_fast_g507     # <- contains "fast"
```

and its `flye.log` confirms it therefore assembled with `--nano-raw`, not `--nano-hq`.

### What the README needs to say
- State the coupling explicitly under `parameters.nanopore` / `parameters.hybrid`.
- To pin the Medaka model WITHOUT touching the assembler, set `flye_input_mode` explicitly
  (`nano-hq` or `nano-raw`) instead of leaving it on `auto`.
- Note that the published long-read baselines were produced with `--nano-raw` via this route.

---

## 8. `medaka_model: auto` only works if the FASTQ headers carry a basecaller tag

**Status: confirmed empirically (2026-07-22).**

`auto` runs `medaka tools resolve_model --auto_model consensus_bacteria <reads>`, which reads
the basecaller model out of the FASTQ headers. If the headers do not contain exactly one
model reference, it fails with:

```
ValueError: Input file did not contain precisely 1 basecaller model reference.
```

The CDRTa11 ONT test data has **no basecaller tag at all** — verified on both the raw FASTQ
and the filtlong output, so this is a property of the data, not of which file BacFlux feeds
in. `auto` simply cannot work for such a dataset, in v1 or v2. That is exactly why the
BacFluxL baseline pins the model explicitly.

### What the README needs to say
State that `auto` requires ONT reads basecalled by a version that stamps the model into the
headers (Guppy/Dorado do; older or re-headered/public data often does not), and that the
remedy is to name the model explicitly — while remembering item 7 above, since an explicit
"fast" model also flips the assembler.

---

## 9. Comparing against a published baseline: check WHICH config it used

**Status: process lesson, learned the hard way (2026-07-22).**

The BacFluxL / BacFluxL+ baselines in `*_test/output_dir` were produced from
`config/config_custom.yaml`, **not** `config/config.yaml`. The two differ in ways that change
results:

| | `config.yaml` | `config_custom.yaml` (what the baselines used) |
|---|---|---|
| `threads` | 24 | **56** |
| `medaka_model` | `auto` | **`r941_min_fast_g507`** (and so, via item 7, `--nano-raw`) |

A v2 validation set up from `config.yaml` therefore ran `--nano-hq` at 24 threads against a
baseline built with `--nano-raw` at 56 threads, and the resulting 1 bp assembly difference
looked like tool non-determinism when it was simply a different command. Always read the
baseline's own `flye.log` / `params.json` for the command that actually ran, rather than
inferring it from the repo default config.
