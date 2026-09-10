# Running BacFlux

One command runs all four modes. Which pipeline you get is decided by `mode:`
inside the config file, not on the command line, so the launch line looks the
same for an Illumina batch and for a folder of finished assemblies.

```bash
conda activate snakemake
cd /path/to/BacFlux

# check what would run, without running it
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml -n

# launch
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

`config/config_custom.yaml` is your own copy of the shipped example
(`cp config/config.yaml config/config_custom.yaml`), git-ignored so your edits
survive a `git pull`. See [Quick start](../getting-started/quick-start.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **DAG** | directed acyclic graph — the job graph Snakemake builds from the rules, which decides what runs and in what order |
    | **RAM** | the machine's memory, as opposed to disk |
    | **CPU** | a processor core, the unit `--cores` counts |
    | **OOM** | out of memory — a job killed by the system for using more than the machine had |
    | **QC** | quality control — the read and assembly checks |
    | **YAML** | the config file's plain-text format; indentation is significant and tabs are refused |
    | **GTDB** | Genome Taxonomy Database, the reference GTDB-Tk places each genome in |

## Launch from the repository root

Three things are anchored to the directory you launch from, so make it the
repository root every time:

- **`workflow/Snakefile`** is found automatically from there, which is why it
  never has to be named on the command line.
- **`.snakemake/`** is created in the launch directory. It holds the per-rule
  conda environments, the working-directory lock, and the metadata that decides
  what a re-run rebuilds. Launch from somewhere else and Snakemake builds a
  second copy of all 31 environments from scratch.
- **Relative paths in the config** resolve against the launch directory. BacFlux
  deliberately does not use Snakemake's `workdir:` directive, so the process
  never changes directory; `directories.output_dir` is resolved to an absolute
  path once at parse time and every output path is built off it
  (`workflow/rules/shared/00_common.smk:160`). Input directories are used as
  given — absolute paths are safest.

!!! warning "The config must be named on the command line"

    `workflow/Snakefile` declares no default `configfile:`. That is deliberate:
    Snakemake **merges** a hard-coded config with the one passed on the command
    line rather than replacing it, so a default would quietly fill in any key
    your copy omits — with the shipped example's empty placeholders
    (`bakta_db: ""` and friends) — and the run would carry on against a path
    that does not exist. Without `--configfile` the run stops at once:

    ```text
    [BacFlux] No configuration supplied. Pass one explicitly, e.g.:
      snakemake --sdm conda --cores N --configfile config/config.yaml
    ```

## Always dry-run first

`-n` (or `--dry-run`) builds the job list and prints it without running anything
or writing a file. It is worth the few seconds every time, because BacFlux does
nearly all of its validation while the workflow is being *read* — before the
first job starts. A dry run therefore catches, among others:

| Checked at parse time | Where the check lives |
|---|---|
| `mode` missing, or not one of the four | `00_common.smk:83` |
| `phage.caller` that is neither `genomad` nor `virsorter2` | `00_common.smk:121` |
| an optional `directories.*` database path that is not a directory, or holds something other than that database | `00_common.smk:242`, `244` |
| the input directory this mode needs, unset or not a directory | `00_common.smk:1665`, `1667` |
| a sample name containing one of <code>* # @ % ^ / ! ? &amp; : ; &#124; &lt; &gt;</code> or a space | `00_common.smk:1676` |
| an Illumina R1 with no R2 beside it | `00_common.smk:1690` |
| in hybrid mode, a sample that has reads of only one kind | `00_common.smk:1734` |
| `links.card_link` or `links.phix_link` unset in a mode that needs them | `00_common.smk:1528`, `1542` |
| `parameters.eggnog.dbmem: true` with `resources.ram_gb` below 42 | `00_common.smk:1287` |

The same parse step prints a header naming every choice the config made. Reading
it costs nothing and catches a config you did not mean to run:

```text
Mode: contigs — Downstream analysis of bacterial whole-genome assemblies.
Phage caller: virsorter2 (default). Plasmid stage: Platon-only (geNomad off).
Decontamination mode: auto (discard_no_hit=true).
Using CheckV database from link: 'https://zenodo.org/records/21510554/files/checkv-db-v1.5.tar.gz' (db_id='checkv-db-v1.5').
Sample BAA-2146 will be processed.
```

Optional parts announce themselves too: the mobilome module adds a line, and
eggNOG's `--dbmem` prints the exact flag to add (see [Memory](#memory)).

## `--sdm conda`

`--sdm` is short for `--software-deployment-method`, and replaces the deprecated
`--use-conda`. It tells Snakemake to build each rule's own environment from the
31 YAML files in `workflow/envs/` and run that rule inside it. Only Snakemake
itself is installed by hand; SPAdes, Flye, Bakta, GTDB-Tk and the rest are
installed by the workflow on first use.

The environments are built once and reused. By default they land in
`.snakemake/conda/` under the launch directory. `--conda-prefix` moves them:

```bash
snakemake --sdm conda --conda-prefix ~/.bacflux-envs --cores 24 \
          --configfile config/config_custom.yaml
```

That is worth doing when several checkouts or several projects should share one
set rather than each building its own — together they are tens of gigabytes, and
building them takes a while.

!!! tip "Build the environments ahead of time"

    `--conda-create-envs-only` builds every environment and stops. Useful before
    a long run, and on machines where only the login node reaches the network.

## Pin files: the exact environments

Beside every environment file in `workflow/envs/` sits a `<name>.linux-64.pin.txt`.
It lists every package in that environment, down to the exact build, as it was
when the workflow was validated. On a linux-64 machine Snakemake finds it
automatically and installs from it instead of solving the yaml, so two installs
made months apart on different machines get identical environments. Nothing
needs configuring; the run log shows `Using pinnings from ...` for each one.

On any other platform the pin file is ignored and the yaml is solved fresh. The
Python version and the tool version are pinned there too, so the environment is
still close to the validated one, but the other packages may move.

!!! note "One environment has no pin file, on purpose"

    `dbcan` is installed through pip, and a pin file is a conda package list: it
    cannot record pip packages, and an install made from a pin file skips the
    yaml's pip section entirely. That environment therefore always builds from
    its yaml. The comment at the top of `dbcan.yaml` says the same, so nobody
    regenerates a pin file for it.

To regenerate a pin file after deliberately changing an environment, validate
the rebuilt environment first, then:

```bash
conda list --explicit -p <path to the built environment> > workflow/envs/<name>.linux-64.pin.txt
```

## `--cores` is the CPU budget

`--cores N` is the one CPU knob, and it is enough on its own. Every rule that
wants more than one core declares Snakemake's built-in `threads:`, which is the
directive `--cores` enforces automatically; the two ABRicate rules declare none
and count as one core each. Rules whose tool stops scaling ask for less, through
`capped_cpus(n)`, which is `min(resources.threads, n)` against the config value
(`00_common.smk:421-429`):

| Thread request | Rules |
|---|---|
| `CPUS` (the full `resources.threads`) | `map_phix`, `illumina_assembly` (SPAdes), `ont_assembly` (Flye), `map_contigs`, `map_sel_contigs` |
| `capped_cpus(24)` | BLAST screening, QUAST, CheckM, Qualimap, GTDB-Tk, Bakta, eggNOG, antiSMASH, dbCAN, BBMap→CARD, Platon, VirSorter2, geNomad, CheckV, Polypolish, Snippy, `fix_start`, Medaka **in nanopore mode** |
| `capped_cpus(16)` | `trim_adapters` (fastp); the mobilome depth and ISOSDB mapping legs |
| `capped_cpus(8)` | NanoPlot raw and filtered QC; Medaka **in hybrid mode** (the two long-read modes carry different values, kept as they were until timings say otherwise); the mobilome tool rules (AMRFinderPlus, ISEScan, CONJScan, ICEscan, the two BLAST naming legs); the two CheckV database rules |
| `capped_cpus(4)` | `virsorter2_db` |

So `--cores 24` gives each rule its declared cap, and `--cores 8` caps every
request at 8 by itself — except `virsorter2_db`, whose own lower cap of 4
correctly wins. Lowering `--cores` is how you run fewer heavy jobs side by side.

!!! warning "Never add `--jobs` alongside `--cores`"

    For a local run `--jobs`/`-j` is an **alias** for `--cores`, not a second,
    independent "N jobs of M cores each" setting — a natural but wrong
    assumption. Passing both switches the core budget off rather than adding a
    cap on top of it. Measured here: two rules each declaring `threads: 8`,
    launched with `--jobs 2 --cores 8`, both received their full 8 threads and
    ran together — 16 real threads against a declared 8-core budget. The same
    rules with `--cores 8` alone ran one at a time. `--jobs` belongs on cluster
    runs, where it caps how many jobs are queued with the scheduler at once.

!!! note "No priority scheme — the DAG orders the rules"

    A rule runs when its inputs exist and never before, which is the only
    ordering BacFlux relies on. It used to carry a `priority:` on 94 rules; all
    94 were removed in v2.0.0. Snakemake's default scheduler maximises the *sum*
    of the priorities of the jobs it starts, unweighted by how many cores each
    takes, so several small per-sample jobs together outscored the one big job the
    numbers were meant to start first — the scheme did the opposite of its intent
    whenever more than a couple of jobs were ready at once.

## Memory

`resources.threads` and `resources.ram_gb` in the config describe your machine.
Snakemake enforces only the CPU side. `ram_gb` is handed to the tools that accept
a memory ceiling, so two large jobs running at once can each take it:

| Rule | What `ram_gb` does |
|---|---|
| `illumina_assembly` | SPAdes `-m`, a hard ceiling it aborts rather than exceed (`illumina/20_assembly.smk:70`, `hybrid/20_assembly.smk:64`) |
| `map_evaluation` | Qualimap JVM heap, `min(ram_gb, 64)` GB (`shared/20_qc.smk:281`) |
| `map_amr_db`, and the mobilome depth and ISOSDB legs | BBMap `-Xmx`, `min(ram_gb, 32)` GB (`shared/50_amr.smk:345`; `shared/80_mobilome.smk:1602`, `1643`) |

Three rules also declare a memory figure as a named Snakemake `resource` — `ram`
on SPAdes, `java_mem` on Qualimap, `mem_gb` on eggNOG — and Snakemake schedules
against a named resource **only when the launch line passes it**, as in
`--resources mem_gb=64`. All three are gigabyte figures, not core counts.

### eggNOG `--dbmem`

eggNOG-mapper's annotation phase does random-access lookups into the 39 GB
`eggnog.db` once per seed ortholog, which is why `functional_annotation` is
routinely the last rule still running. Setting `parameters.eggnog.dbmem: true`
loads that database wholly into RAM instead, at roughly **42 GB per concurrent
eggNOG job** (`00_common.smk:1281`). It is off by default.

Turning it on prints the flag to add, with your own `ram_gb` filled in:

```text
eggNOG --dbmem is ON: each functional_annotation job loads the 39 GB eggnog.db into RAM
(~42 GB/job; 1 fit in ram_gb=64). To have Snakemake cap concurrent eggNOG jobs to that
many, add '--resources mem_gb=64' to your launch command.
```

```bash
snakemake --sdm conda --cores 24 --resources mem_gb=64 \
          --configfile config/config_custom.yaml
```

Without the flag the run still works; eggNOG concurrency then falls back to the
`--cores` bound instead of the RAM one, which is how a machine with 64 GB ends up
running two 42 GB jobs at once. Setting `dbmem: true` with `ram_gb` below 42 is
refused at parse time rather than allowed to fail on the first job.

## Resuming an interrupted run

Snakemake rebuilds only what is missing or out of date, so an interrupted run is
relaunched with the same command. Everything that finished cleanly is kept.

An interruption leaves two kinds of debris.

- **Part-written outputs.** A job killed mid-write leaves a file Snakemake marks
  incomplete. `--rerun-incomplete` (`--ri`) rebuilds those instead of trusting
  them.
- **A locked working directory.** A run killed outright — `Ctrl-C` twice, a node
  failure, an OOM kill — can leave `.snakemake/` locked. Snakemake says so on the
  next launch; `snakemake --unlock` clears it.

Changing the config can be enough to mark a step out of date on its own.
Snakemake compares more than file timestamps: a changed rule parameter, script or
conda environment counts too (`--rerun-triggers` defaults to
`code input mtime params software-env`). A dry run answers "what would re-run?"
directly, and faster than reasoning about it.

!!! warning "`--ignore-incomplete` and `--keep-incomplete` are not resume flags"

    They are the opposite of `--rerun-incomplete`: `--keep-incomplete` leaves a
    failed job's part-written output on disk, and `--ignore-incomplete` stops
    Snakemake checking for it, so that file can be accepted as a finished result.
    They are useful when you want to *inspect* what a failing rule produced.
    After an interruption, use `--rerun-incomplete`.

    `--keep-going` (`-k`) is different and safe: it lets independent jobs carry on
    when one fails, so a whole batch is not lost to one bad sample. It applies
    only to runtime failures, never to a parse or DAG error.

## When a QC tool fails

Four tools produce reports that nothing else in the workflow reads: QUAST
(assembly metrics), Qualimap (read-mapping quality), NanoPlot (long-read QC) and
CheckV (prophage grading). MultiQC gathers what they wrote; no analysis step
depends on them. So if one of them crashes on a sample, the rule finishes
anyway: it creates its output directory first, lets the tool fail, and writes a
note beginning `NOTE:` into that sample's log under `logs/`. The run carries on,
and that sample simply has no section for that tool in the MultiQC report.

The reverse is deliberate too. CheckM, the genome-staging step and the
assembly-depth step look like QC but feed later rules — GTDB-Tk, QUAST itself,
the IS copy-number layer — so a failure there stops the run rather than letting
a bad input flow downstream.

For everything else, `--keep-going` (`-k`) is the flag to reach for: when one
sample's job fails, the other samples keep running to completion instead of the
whole batch stopping. Fix the one sample, then resume with `--rerun-incomplete`.

## What survives between runs

| Kept where | Contents | Cost of deleting it |
|---|---|---|
| `.snakemake/conda/` in the launch directory | the 31 per-rule environments | rebuilt from scratch, needs network |
| `.snakemake/` metadata | what is up to date, and the lock | the next run recomputes more than it needs to |
| under `output_dir` | the databases the workflow fetches and keeps — antiSMASH, dbCAN, CheckV, the phage caller's own database, and the optional mobilome databases | re-downloaded |
| `output_dir/logs/` | one log per rule and sample (`00_common.smk:539`) | nothing recomputes; you lose the record |

PhiX and CARD are the two that do not persist: both are declared `temp()`
(`illumina/10_reads.smk:67`, `shared/50_amr.smk:201`), so Snakemake deletes them
as soon as the last rule reading them finishes, and every run fetches them again.

Reusing the same launch directory and the same `output_dir` is what keeps the
rest from being rebuilt. Databases you already hold elsewhere can be pointed at
instead of downloaded — see
[Reference databases](../getting-started/databases.md).

## Cleaning up a finished project

`workflow/scripts/clean_workdir.sh` strips a finished output directory down to
results, reports and audit trails, deleting bulky intermediates (trimmed and
filtered reads, SPAdes and Flye working directories, Medaka scratch, CheckM
scratch) and every database the workflow keeps in the output tree. Raw input
reads are never touched — they are read in place, never copied into the output
tree.

```bash
# preview — this is what --run would delete
workflow/scripts/clean_workdir.sh /path/to/output_dir

# actually delete
workflow/scripts/clean_workdir.sh --run /path/to/output_dir
```

It is **dry-run by default**; only `--run` deletes, and the listing printed is
identical in both modes, so the preview is exactly what you get. It recognises
the current layout and the four older ones it replaced, so output directories
made before v2.0.0 can still be cleaned — see
[Coming from v1](../getting-started/from-v1.md).

!!! warning "This is an archive step, not a between-runs tidy"

    It is deliberately aggressive. Re-running the workflow on a cleaned directory
    recomputes a great deal, assembly and polishing included. Run it once the
    results are known good and you are ready to archive or hand over.

`--include-snakemake` additionally deletes a `.snakemake/` directory found
*inside* the output directory — which only exists if you launched from in there
rather than from the repository root. It is opt-in because it is categorically
heavier than the rest: that directory holds the conda environments, routinely the
largest single item in a finished project. When one is present but not opted
into, the script prints its size rather than quietly leaving it.

!!! note "`miscellaneous/clean_workdir.sh` still works"

    It is a wrapper that forwards every argument to the real script, kept because
    older notes and habits point at that path.

## For contributors

Every script a rule runs lives in a directory named after the rule file that runs
it — `workflow/scripts/50_amr/` is run by `workflow/rules/shared/50_amr.smk`, and
so on — with each test beside its script. A script sitting loose in
`workflow/scripts/` is **not** part of a run; it is a maintenance tool you run by
hand, as `clean_workdir.sh` is. See `workflow/scripts/README.md`.

```bash
# the test suite
pytest workflow/scripts -q

# fail on any rule name, path, config key or constant that a comment
# mentions but that no longer exists
python miscellaneous/check_comment_references.py
```

### Before tagging a release

Run the whole workflow, not a dry run. A dry run resolves the graph and checks
inputs; it never builds an environment and never executes a tool. Building the
environments without running them is not enough either: a package can install
cleanly and still break at runtime.

The check that catches these: a fresh clone, a fresh `--conda-prefix`, one small
public genome, every mode. It takes an hour or two, and it is the only test that
exercises what a new user actually gets.

The version number is declared in five places, and they have to move together:

| File | What to change |
|---|---|
| `workflow/Snakefile` | the header comment and the banner it prints |
| `README.md` | the label under the ASCII banner |
| `docs/index.md` | "This site documents version ..." |
| `.zenodo.json` | `version` |
| `CITATION.cff` | `version`, and `date-released` for the tag day |

Leave every other mention of an old version alone: most of them are statements
about what that release changed, and they stay true.

After Zenodo archives the release, add the version DOI it mints to the
`identifiers` block in `CITATION.cff`, replacing the previous version entry.

## See also

- [Configuration](configuration.md) — every key named here, in the file's own order.
- [Output files](output.md) — what a run writes, and where.
- [Troubleshooting](../troubleshooting.md) — what to do when a run stops.
