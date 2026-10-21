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
  second copy of every environment from scratch.
- **Relative paths in the config** resolve against the launch directory. BacFlux
  deliberately does not use Snakemake's `workdir:` directive, so the process
  never changes directory; `directories.output_dir` is resolved to an absolute
  path once at parse time and every output path is built off it (`OUT` in
  `workflow/rules/shared/00_common.smk`). Input directories are used as
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

| Checked at parse time | Where the check lives in `00_common.smk` |
|---|---|
| `mode` missing, or not one of the four | where `MODE` is set |
| `phage.caller` that is neither `genomad` nor `virsorter2` | where `PHAGE_CALLER` is set |
| an optional `directories.*` database path that is not a directory, or holds something other than that database | `_resolve_optional_db_dir()` |
| `resources.threads`, `ram_gb` or `runtime_hours` that is not a whole number of 1 or more | `_whole_number()` |
| `parameters.eggnog.dbmem: true` with `resources.ram_gb` below 42 | `EGGNOG_DBMEM_GB` |
| `parameters.decontamination.mode` or `auto_vote` that is not one of its options | `_decontam_settings()` |
| `links.card_link` or `links.phix_link` unset in a mode that needs them | where `CARD_LINK` and `PHIX_LINK` are set |
| the input directory this mode needs, unset or not a directory | `_require_input_dir()` |
| a sample name containing one of <code>* # @ % ^ / ! ? &amp; : ; &#124; &lt; &gt;</code> or a space | `_check_sample_names()` |
| an Illumina R1 with no R2 beside it | `_require_r2_mates()` |
| in hybrid mode, a sample that has reads of only one kind | the hybrid sample discovery (`_missing_long`, `_missing_short`) |

Each of these stops the run before any job starts, with exit status 1, so a
wrapper script, a CI job or a cluster submission sees the failure as a failure.
The message begins with `[BacFlux]`, except for a missing `mode`, which Python
reports as a `KeyError`.

The same parse step prints a header naming every choice the config made. Reading
it costs nothing and catches a config you did not mean to run:

```text
Mode: contigs — Downstream analysis of bacterial whole-genome assemblies.
Phage caller: virsorter2 (default). Plasmid stage: Platon-only (geNomad off).
Decontamination mode: auto (auto_vote=bp, discard_no_hit=true).
Using CheckV database from link: 'https://zenodo.org/records/21510554/files/checkv-db-v1.5.tar.gz' (db_id='checkv-db-v1.5').
Sample BAA-2146 will be processed.
```

Optional parts announce themselves too: the mobilome module adds a line, and
eggNOG's `--dbmem` prints the exact flag to add (see [Memory](#memory)).

## `--sdm conda`

`--sdm` is short for `--software-deployment-method`, and replaces the deprecated
`--use-conda`. It tells Snakemake to build each rule's own environment from the
YAML files in `workflow/envs/` and run that rule inside it. Only Snakemake
itself is installed by hand; SPAdes, Flye, Bakta, GTDB-Tk and the rest are
installed by the workflow on first use. The small steps get pinned environments
too: BacFlux's own Python scripts run in `python.yaml`, and the download steps
take `wget`, `curl`, `tar` and the unpackers from `download.yaml` or from the
environment of the tool they serve, so the computer's own Python and download
tools never matter.

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
    a long run, and on machines where only the login node reaches the network
    ([On a cluster](#on-a-cluster)).

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
(`capped_cpus()` in `00_common.smk`):

| Thread request | Rules |
|---|---|
| `CPUS` (the full `resources.threads`) | `map_phix`, `illumina_assembly` (SPAdes), `ont_assembly` (Flye), `map_contigs`, `map_sel_contigs`, `map_final_contigs` |
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
    cap on top of it. Measured: two rules each declaring `threads: 8`,
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
On a single computer Snakemake enforces only the CPU side. `ram_gb` is handed to
the tools that accept a memory ceiling, so two large jobs running at once can each
take it:

| Rule | What `ram_gb` does |
|---|---|
| `illumina_assembly` | SPAdes `-m`, a hard ceiling it aborts rather than exceed (the rule's `ram` resource, in `illumina/20_assembly.smk` and `hybrid/20_assembly.smk`) |
| `map_evaluation` | Qualimap's Java heap, the most memory Java may use: `ram_gb` minus 4 GB, at most 64 GB and at least 1 GB (`QUALIMAP_HEAP_GB` in `00_common.smk`). The 4 GB are left for what Java needs beyond its heap |
| `map_amr_db`, and the mobilome depth and ISOSDB legs | BBMap `-Xmx`, `min(ram_gb, 32)` GB (the `max_ram` value in each rule's `params:`) |

Three rules also declare a memory figure as a named Snakemake `resource` — `ram`
on SPAdes, `java_mem` on Qualimap, `mem_gb` on eggNOG — and Snakemake schedules
against a named resource **only when the launch line passes it**, as in
`--resources mem_gb=64`. All three are gigabyte figures, not core counts.
`--resources java_mem=N` also lowers each Qualimap heap to N GB, which is the
way to give Qualimap less memory without lowering `ram_gb` for everything else.

The heavy rules also declare `mem_mb` and `runtime`, the memory and time a job
asks a cluster scheduler for; on a cluster they are how `ram_gb` and
`resources.runtime_hours` reach the scheduler — see [On a cluster](#on-a-cluster).
A single computer ignores `runtime`, and uses `mem_mb` only when the launch line
passes `--resources mem_mb=N`: each heavy step then counts as `ram_gb` × 1024 MB
against N, so with N at or below that figure only one heavy step runs at a time.

### eggNOG `--dbmem`

eggNOG-mapper's annotation phase does random-access lookups into the 39 GB
`eggnog.db` once per seed ortholog, which is why `functional_annotation` is
routinely the last rule still running. Setting `parameters.eggnog.dbmem: true`
loads that database wholly into RAM instead, at roughly **42 GB per concurrent
eggNOG job** (`EGGNOG_DBMEM_GB` in `00_common.smk`). It is off by default.

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

The same comparison applies after updating BacFlux. A new release changes rule
code and environments, so the first run in an existing output directory plans
to redo those steps and everything that depends on them, which can reach back to
the assembly. Run `-n` first to see the plan. To keep the results you have and
rebuild only what is missing or older than its inputs, add
`--rerun-triggers mtime` to every run in that directory, or finish the project
with the release it was started with.

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
| `.snakemake/conda/` in the launch directory, or the `--conda-prefix` folder | the per-rule environments, and the Medaka model weights stored inside the medaka one | rebuilt from scratch, needs network |
| `.snakemake/` metadata | what is up to date, and the lock | the next run recomputes more than it needs to |
| under `output_dir` | the databases the workflow fetches and keeps — the PhiX genome, antiSMASH, dbCAN, CheckV, the phage caller's own database, and the optional mobilome databases | re-downloaded |
| `output_dir/logs/` | one log per rule and sample (`LOGS` in `00_common.smk`) | nothing recomputes; you lose the record |

CARD is the one download that does not persist: rule `download_amr_db` declares it
`temp()`, so Snakemake deletes it as soon as the last rule reading it finishes,
and the next run that needs it fetches it again.

Reusing the same launch directory and the same `output_dir` is what keeps the
rest from being rebuilt. Databases you already hold elsewhere can be pointed at
instead of downloaded — see
[Reference databases](../getting-started/databases.md).

## On a cluster

BacFlux runs on a Slurm cluster in one of two ways:

| | A. Snakemake sends each step to Slurm | B. The whole run in one batch job |
|---|---|---|
| Where Snakemake runs | on the login node, or another node with internet, for the whole run | inside one batch job, on one compute node |
| Where the steps run | each step is a Slurm job of its own, on any compute node; the download steps and the Medaka model check run where Snakemake runs | all on that node |
| Command | `--executor slurm --jobs N` ([Run A](#cluster-run-a)) | `--cores N`, as on a single computer ([Run B](#cluster-run-b)) |
| Memory and time | each heavy step asks Slurm for `resources.ram_gb` and `resources.runtime_hours` | the batch job asks once, for the whole run |
| Internet on the compute nodes | not needed | not needed once the [preparation](#cluster-preparation) has run; the long-read modes need [one more step](#cluster-medaka) for the Medaka model weights |
| Needs | the Slurm executor plugin; a login node on which Snakemake may run for many hours | one node with the cores, memory and time of the whole run |
| Checked | the local rules by the test suite, and every job's requests by dry runs (`-n` prints them); not yet run on a Slurm cluster | as on a single computer |

In short: A spreads the samples and steps over many nodes and keeps every
download on the login node; B is a single-computer run inside one job. Both start
with the same preparation, and both need the same [shared storage](#cluster-storage).

### Preparation, on a computer with internet { #cluster-preparation }

Run these commands from the BacFlux folder on a computer that has internet and
sees the same storage as the compute nodes: the login node, or a data-transfer
node.

1. Copy the template and fill it in. Put every path on
   [shared storage](#cluster-storage), and set `resources.threads`, `ram_gb` and
   `runtime_hours` from the limits of the cluster ([table below](#cluster-resources)).

    ```bash
    conda activate snakemake
    cd /path/to/BacFlux
    cp config/config.yaml config/config_custom.yaml
    ```

2. Build the conda environments of every step of the run. Some of them fetch data
   while they are built — CheckM's reference data, for one — so this step needs
   internet too. conda first downloads every package into its package cache, by
   default inside the conda installation; if that sits in a home folder with a
   small quota, move the cache to project space first, as in the first line:

    ```bash
    export CONDA_PKGS_DIRS=/path/to/project_space/conda_pkgs
    snakemake --sdm conda --conda-prefix /path/to/project_space/bacflux_envs \
        --conda-create-envs-only --configfile config/config_custom.yaml
    ```

3. Fetch everything BacFlux downloads itself, and nothing else. Add `-n` first to
   see the plan.

    ```bash
    snakemake download_all --sdm conda --conda-prefix /path/to/project_space/bacflux_envs \
        --cores 4 --configfile config/config_custom.yaml
    ```

4. Check. The dry run must print
   `Nothing to be done (all requested files are present and up to date).`

    ```bash
    snakemake download_all -n --sdm conda --conda-prefix /path/to/project_space/bacflux_envs \
        --configfile config/config_custom.yaml
    ```

| Rule | Why |
|---|---|
| Write `download_all` right after `snakemake`. | `--configfile` and `--sdm` each take one or more values, so at the end of the line `download_all` is read as one of them and Snakemake stops with an error. |
| Run steps 2 to 4 and the real run from the same folder, with the same configuration and the same `--conda-prefix`. | Relative paths in the configuration are read from the launch folder, and Snakemake keeps its records in `.snakemake/` there. Without `--conda-prefix` the environments are in that `.snakemake/conda/` too. Snakemake also records where each step's environment was, so with another prefix the real run builds the environments again and runs every download step again (reason `Software environment definition has changed since last execution`). |
| Run steps 2 to 4 again after any change to the configuration. | A short-read mode adds PhiX and CARD, the other phage caller adds its own database, and each mobilome layer switched on adds its own. |

Behind a proxy that re-signs HTTPS traffic, point both `SSL_CERT_FILE` and
`REQUESTS_CA_BUNDLE` at your site's certificate file before starting Snakemake. The
download tools check servers against the certificate list in their conda
environment, not the computer's: `wget` and `curl` read `SSL_CERT_FILE`, and the
Python downloaders, among them Medaka's download of its model weights, read
`REQUESTS_CA_BUNDLE`.

What `download_all` fetches, according to the configuration:

| Download | Rule | When | Where it lands, under `output_dir` |
|---|---|---|---|
| PhiX genome | `download_phix` | `illumina`, `hybrid` | `01.reads/phix/` |
| CARD | `download_amr_db` | `illumina`, `hybrid`, unless `directories.card_db` names your copy | `05.amr/card_db/`; temporary, so the run deletes it after its last reader |
| dbCAN database | `cazyme_db_download` | unless `directories.dbcan_db` names your copy | `04.annotation/dbcan/<release>/` |
| antiSMASH databases | `secondary_metabolites_db` | unless `directories.antismash_db` names your copy | `04.annotation/antismash/databases/` |
| VirSorter2 database | `virsorter2_db` | `phage.caller: virsorter2`, unless `directories.vs2_db` names your copy | `07.phages/vs2_db/` |
| geNomad database | `genomad_db` | `phage.caller: genomad`, unless `directories.genomad_db` names your copy | `07.phages/genomad_db/` |
| CheckV database | `checkv_db` | unless `directories.checkv_db` names your copy | `07.phages/checkv_db/` |
| CONJScan models | `conjscan_models` | `mobilome.run` | `08.mobilome/conjscan_models/` |
| ICEscan models | `icescan_models` | `mobilome.icescan.run` | `08.mobilome/icescan_models/` |
| TnCentral | `tncentral_db` | the `tncentral` layer is on | `08.mobilome/tncentral_db/` |
| ICEberg | `iceberg_db` | the `iceberg` layer is on | `08.mobilome/iceberg_db/` |
| ISOSDB | `isosdb_db` | the `isosdb` layer is on, in `illumina` and `hybrid` | `08.mobilome/isosdb_db/` |

`download_all` leaves out two things. The five databases you provide yourself
(`bakta_db`, `blast_db`, `eggnog_db`, `gtdbtk_db`, `platon_db`) are named by path
and never downloaded by the workflow
([Reference databases](../getting-started/databases.md)). And the Medaka model
weights are fetched during the run, by the model check — see
[below](#cluster-medaka). The list itself is `download_targets()` in
`00_common.smk`.

### Run A: Snakemake sends each step to Slurm { #cluster-run-a }

Install the Slurm executor plugin once into the launcher environment:

```bash
conda activate snakemake
conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

Then, on the login node, start Snakemake inside `screen` or `tmux`, so that it
keeps running when you log out:

```bash
snakemake --executor slurm --jobs 20 --latency-wait 60 --sdm conda \
    --conda-prefix /path/to/project_space/bacflux_envs \
    --default-resources slurm_account=<account> slurm_partition=<partition> mem_mb=8000 runtime=240 \
    --configfile config/config_custom.yaml
```

| Option | What it does |
|---|---|
| `--executor slurm` | submits every step to Slurm with `sbatch`, except the local rules: the download steps and `check_medaka_model` run on the login node |
| `--jobs 20` | the most jobs Snakemake keeps queued or running at once; it refuses to start without `--jobs` |
| `--latency-wait 60` | after a job ends, Snakemake waits up to 60 seconds for its files to appear on the login node. Shared storage can show a file written on a compute node only some seconds later, and the default wait, 5 seconds, is often too short |
| `--default-resources` | the values for every job that does not set its own: the account and partition for every job (leave them out if your cluster needs none), memory (MB) and time (minutes) for the light steps. A rule's own value wins over it |

`--cores` is not needed: `resources.threads` already caps the CPUs of each job.
Given together with `--executor`, it does not cap the jobs' total; it lowers each
job's threads to at most that number, and the job then asks Slurm for that many.

The local rules run on the login node with their own thread counts, and Snakemake
9.14.6 has no setting that limits them there. The heaviest is the CheckV download,
which builds a DIAMOND index with up to 8 threads. Running the
[preparation](#cluster-preparation) first means every download is done before the
run starts, so during the run the only local step is the Medaka model check.

Each Slurm job writes a log of its own under `.snakemake/slurm_logs/` in the
launch folder, and when a job fails the `Error in rule` block names that file.
The messages of the tool itself are in `logs/` of the output folder, as on a
single computer.

#### How the settings reach Slurm { #cluster-resources }

| Setting | Each job asks Slurm for | Must not exceed | Shown by |
|---|---|---|---|
| `resources.threads` | `--cpus-per-task`: the threads of its step, at most this number (the thread table under `--cores` above) | the CPUs of one compute node | `sinfo -o "%n %c"` |
| `resources.ram_gb` | `--mem`: `ram_gb` × 1024 MB for every heavy step (table below) | the memory Slurm offers one compute node | `sinfo -o "%n %m"`, in MB: divide by 1024 and round down |
| `resources.runtime_hours` | `-t`: `runtime_hours` × 60 minutes for every heavy and every long step | the time limit of the partition | `sinfo -o "%P %l"` |
| `--default-resources mem_mb=... runtime=...` | `--mem` and `-t` for the light steps | the same limits | the same commands |

| Steps | Memory asked for | Time asked for |
|---|---|---|
| the heavy steps: SPAdes (`illumina_assembly`), Flye (`ont_assembly`), Medaka (`long_read_consensus`), Polypolish (`short_read_correction`), the read mappings of decontamination (`map_contigs` in the read-based modes, `map_final_contigs`), the BLAST screens (`blast_contigs`, `blast_final_contigs`), CheckM (`completeness_and_contamination`), GTDB-Tk (`taxonomic_assignment`), Bakta (`annotation`), eggNOG-mapper (`functional_annotation`), antiSMASH (`secondary_metabolites_analysis`), BBMap against CARD (`map_amr_db`), geNomad (`genomad_end_to_end`), VirSorter2 (`viral_identification_virsorter2`), and the two read mappings of the mobilome module (`assembly_depth`, `isosdb_map`) | `ram_gb` | `runtime_hours` |
| Qualimap (`map_evaluation`) | its Java heap plus 4 GB, never more than `ram_gb` ([Memory](#memory)) | `runtime_hours` |
| the long steps that hold little in memory: the PhiX filter (`map_phix`) and the read selection of `hybrid` (`map_sel_contigs`) | nothing | `runtime_hours` |
| every other step | nothing | nothing |

The figures come from `HEAVY_MEM_MB`, `QUALIMAP_JOB_MB` and `LONG_RUNTIME_MINUTES`
in `00_common.smk`, and a dry run prints each job's requests.

In short: set `threads`, `ram_gb` and `runtime_hours` so that one job fits one
node and one partition. GTDB-Tk alone needs about 55 GB, so give `ram_gb` at least
about 60. A `ram_gb` well below the node's memory lets several heavy jobs share a
node.

- A step that asks for nothing gets Slurm's defaults, often a GB or two per CPU
  and the partition's default time. If those are smaller than a few GB or an
  hour, give `--default-resources mem_mb=8000 runtime=240` as above, and raise it
  if a light step is stopped for memory or time.
- Changing `threads` or `runtime_hours` re-runs nothing. Changing `ram_gb`
  re-runs the BBMap steps when `min(ram_gb, 32)` changes with it, because that
  number is one of their parameters. In `illumina` and `hybrid` that reaches much
  further: the trimmed reads those steps read are temporary, so Snakemake trims
  the reads again, and everything that reads them runs again too, the assembly
  included. Settle `ram_gb` before the first run.

### Run B: the whole run in one batch job { #cluster-run-b }

| Step | What to do |
|---|---|
| Request | one node with N cores and enough memory: at least `resources.ram_gb`, which SPAdes takes as its ceiling and Qualimap sizes its heap from, and at least about 60 GB, because GTDB-Tk needs 55 GB |
| Configure | `resources.threads: N`, the number of cores you request |
| Prepare | when the node has no internet: the [preparation](#cluster-preparation) above, on the login node, with the folder and `--conda-prefix` of the batch script; in the long-read modes also the [Medaka step](#cluster-medaka) |
| Launch | inside the job, from the same folder, the usual command with `--cores N` ([script below](#cluster-batch-script)) |
| Time limit | the whole run, plus the building of the environments and the downloads when they were not prepared. After a job is stopped at its time limit, submit the same script again ([Resuming an interrupted run](#resuming-an-interrupted-run)); if Snakemake asks for `--unlock`, run the command once with it first |

In short: treat the node like a single computer. On one node Snakemake counts
only threads: `resources.ram_gb` is then only the memory ceiling of SPAdes, the
size of the Qualimap heap and the cap of BBMap, and `runtime_hours` is ignored.

#### A batch script { #cluster-batch-script }

For Slurm; other schedulers take the same commands with their own resource
lines. Save it as `bacflux.sbatch` and submit it with `sbatch bacflux.sbatch`.

```bash
#!/bin/bash
#SBATCH --job-name=bacflux
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24      # N: the same number as resources.threads
#SBATCH --mem=128G              # at least resources.ram_gb (Request row above)
#SBATCH --time=24:00:00         # Time limit row above

# Make "conda activate" work inside a batch job, then activate Snakemake.
eval "$(conda shell.bash hook)"
conda activate snakemake

cd /path/to/BacFlux
# --rerun-incomplete: rebuild what a time limit cut short, when the job is resubmitted.
snakemake --sdm conda --conda-prefix /path/to/project_space/bacflux_envs \
    --cores "$SLURM_CPUS_PER_TASK" --rerun-incomplete \
    --configfile config/config_custom.yaml
```

`$SLURM_CPUS_PER_TASK` holds the number of cores Slurm gave the job, so `--cores`
always matches the request.

### Medaka model weights without internet { #cluster-medaka }

In `nanopore` and `hybrid` mode the Medaka model check, rule
`check_medaka_model`, finds the model and fetches its weights during the run. It
is a local rule: in Run A it runs on the login node, where there is internet, and
nothing more is needed. `download_all` cannot fetch the weights ahead of the run,
because the check reads the filtered long reads to identify the model, and those
exist only once the run has filtered them.

So in Run B, on a node without internet, run the workflow once up to the check,
on the login node, from the same folder and with the same `--conda-prefix`:

```bash
snakemake --until check_medaka_model --sdm conda \
    --conda-prefix /path/to/project_space/bacflux_envs \
    --cores 4 --configfile config/config_custom.yaml
```

That filters the long reads there (Filtlong, one core per sample; in `hybrid`
with `short_read_guidance: true` the Illumina steps before it run too), checks the
model and stores its weights in the medaka environment, where the batch job finds
them. The other way is to switch polishing off, with `medaka_model: false`.

### Shared storage { #cluster-storage }

| Must be on storage every compute node sees, under the same path | Why |
|---|---|
| the BacFlux folder and the launch folder, with its `.snakemake/` | every job reads the workflow and Snakemake's records |
| the conda installation, with the launcher environment in it | in A, each job starts Snakemake again with the same Python, and Snakemake switches to each step's environment with that installation's `activate` script; in B, the batch script activates it |
| the `--conda-prefix` folder | every step runs in its environment from there |
| the output folder | every step reads and writes it |
| every path in the configuration: the input folders, the five databases you provide, any `directories.*` database copy, the decontamination files | the steps read them; in A every job also runs the start-up checks again, which look at these paths |
| your home folder, when Medaka could not write into its environment | Medaka then keeps the weights in `~/.medaka/data` |

In short: on most clusters this means project space plus the home folder.

- Put the environments on project space with one `--conda-prefix`, not on a home
  folder with a small quota: together they are tens of gigabytes.
  `export CONDA_PKGS_DIRS=/path/to/project_space/conda_pkgs` moves conda's
  package cache there too.
- Launch every run of a project from the same folder: the environments, and the
  records of what is up to date, are tied to the launch folder and the prefix.

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
# the test suite: the scripts' unit tests, and the guard tests on the rule files
pytest workflow/scripts miscellaneous -q

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

The oldest Snakemake BacFlux accepts is a separate number: `min_version()` near
the top of `workflow/Snakefile` stops an older Snakemake at once. When it moves,
move the README badge and the installation page with it.

After Zenodo archives the release, add the version DOI it mints to the
`identifiers` block in `CITATION.cff`, replacing the previous version entry.

## See also

- [Configuration](configuration.md) — every key named here, in the file's own order.
- [Output files](output.md) — what a run writes, and where.
- [Troubleshooting](../troubleshooting.md) — what to do when a run stops.
