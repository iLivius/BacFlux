# Installation

To run BacFlux, obtain the repository, create a Snakemake environment and prepare
the required reference databases. The repository contains the workflow rules,
helper scripts and configuration template. Snakemake creates separate conda
environments for the analysis tools used by the selected mode.

The instructions below cover the repository and software setup. Prepare the
[Reference databases](databases.md) before starting an analysis.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **BAM** | the binary alignment file a read mapper writes |

## 1. Get the repository

```bash
git clone https://github.com/iLivius/BacFlux.git
cd BacFlux
```

A release archive from the repository's Releases page unpacks to the same tree.
The clone is easier to keep current.

Four places in it matter:

| Path | What is in it |
|------|---------------|
| `config/config.yaml` | The configuration template: every setting of a run, commented in place. Copy it and edit the copy — never this file. |
| `workflow/Snakefile` | The entry point. It reads `mode` from your config and includes that mode's rule files plus the shared ones. |
| `workflow/rules/` | The rules themselves, in `shared/` and one directory per mode (`illumina/`, `nanopore/`, `hybrid/`, `contigs/`). |
| `workflow/envs/` | 33 conda environment specifications, one per tool stack. |

A fifth directory, `.snakemake/`, appears in whatever directory you launch from:
it holds the conda environments Snakemake builds and its own bookkeeping. It is
git-ignored, so it never ends up in a commit.

## 2. Create the launcher environment

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=9.14.6
conda activate snakemake
```

Snakemake **9.14.6 or newer**. That is the version BacFlux is tested with, and
the Snakefile declares it as the minimum: an older Snakemake stops at once, before
reading a single rule, with a message naming the version it needs. The pin in the
command above gives you exactly the tested version.

This environment only runs Snakemake, which reads the configuration, plans the run
and starts jobs; it holds no analysis software. Snakemake builds each tool's
environment from `workflow/envs/` itself, and a rule that has its own environment
does not see tools added to this one.

## 3. What `--sdm conda` does

Every BacFlux run is started the same way, from the repository root:

```bash
# config_custom.yaml is your own copy of the shipped config/config.yaml
snakemake --sdm conda --configfile config/config_custom.yaml --cores 24
```

`--sdm` is short for `--software-deployment-method`, the modern spelling of
`--use-conda`, which has been deprecated since Snakemake 8.4.7.

Each rule names an environment file in `workflow/envs/` — the Bakta rules point
at `bakta.yaml`, the Flye rule at `flye.yaml`, and so on. Snakemake works out
which of them this run needs and builds them all before the first job starts,
then reuses them on every later run. They are created under `.snakemake/conda`
in the launch directory.

Two practical consequences:

- **The first run starts slowly.** Solving and downloading a dozen or more conda
  environments takes minutes. After that they are on disk and the workflow starts
  straight away.
- **Tools never have to coexist.** GTDB-Tk, antiSMASH, VirSorter2 and MacSyFinder
  each get their own environment, which reduces dependency conflicts between tools. The
  environments are pinned to different degrees: some to an exact version because a
  parser downstream depends on the output columns (ISEScan 1.7.3, MacSyFinder
  2.1.6), others only loosely.

!!! note "Only the environments your run actually needs are built"

    The Snakefile globs `workflow/rules/<mode>/*.smk` for the mode in your config,
    so an `illumina` run never solves the Flye or Medaka environment and a
    `nanopore` run never solves SPAdes. Several shared rules are gated the same
    way, on whether the run has short reads, long reads or reads at all.

| Run | Environments built |
|-----|--------------------|
| Any mode | blast, blobtools, quast, checkm, gtdbtk, bakta, antismash, dbcan, eggnog-mapper, abricate, platon, checkv, multiqc, python, download, and either virsorter or genomad |
| `illumina` adds | fastp, bowtie, spades, bbmap, qualimap |
| `nanopore` adds | nanoplot, filtlong, flye, dnaapler, medaka, minimap, qualimap |
| `hybrid` adds | everything in the two rows above, plus polypolish and snippy. Bowtie2 builds the first screen's coverage track; minimap maps the long reads for the second screen, over the delivered genome |
| `contigs` adds | minimap — its own front-end rule is plain `awk`, but BlobTools still wants a BAM, so the contigs are mapped against themselves |
| `mobilome.run: true` adds | isescan, macsyfinder, and tncentral for the optional naming and copy-number layers |

## 4. The tools BacFlux expects to find, and does not install

Each analysis tool runs in a conda environment that Snakemake builds (section 3).
Three things come from the computer instead, and BacFlux never installs them:

- **conda and the launcher environment** with Snakemake (section 2);
- **bash**, which Snakemake uses to run the commands of every rule;
- **the standard Unix command-line utilities.** Rules call them for the small jobs
  around the analysis tools — creating folders, copying and linking files, checking
  a download — and, unless their own environment happens to include one, take them
  from the system `PATH`. The complete list is `awk`, `basename`, `cat`, `chmod`,
  `cp`, `cut`, `date`, `dirname`, `find`, `grep`, `head`, `ln`, `ls`, `md5sum`,
  `mkdir`, `mv`, `rm`, `sed`, `sha256sum`, `sort`, `tail`, `tee`, `tr` and `wc`.
  They belong to the coreutils, findutils, grep, sed and awk packages of every Linux
  installation; only a stripped-down system, such as a minimal container image, may
  lack some of them.

**Download tools are the exception.** A minimal Linux system can lack them — a
minimal Debian or Ubuntu system has no `wget`, and a minimal Red Hat-family system
no `tar` — so BacFlux provides them:

- `wget`, `tar`, `gzip` and `bzip2` come from `workflow/envs/download.yaml` (PhiX,
  CARD, dbCAN, the ICEscan models) or from the tool's own environment (geNomad,
  VirSorter2, CheckV);
- `curl` and `unzip` come from `workflow/envs/tncentral.yaml`, for the optional
  mobilome layers (TnCentral, ICEberg, ISOSDB);
- the CheckM package downloads its reference data (about 275 MB) while conda builds
  its environment, using the `wget`, `tar` and `gzip` in `workflow/envs/checkm.yaml`;
  building that environment therefore needs internet.

The download recipes on [Reference databases](databases.md), which you run yourself,
use the `wget` and `tar` of your own system.

**Rules without an environment.** Eight rules declare no conda environment, and two
more, `all` and `download_all`, are targets that run no command of their own. The
eight run in the environment you launched Snakemake from and use only utilities from
the list above:

| Group | Rules | Uses |
|-------|-------|------|
| Filter, copy or stage a FASTA | `filter_contigs`, `finalize_contigs`, `stage_qc_genomes` | `awk`, `grep`, `head`, `cat`, `cp`, `mkdir`, `ls`, `dirname`, `tee` |
| Build a symlink view of a database you already hold | `download_amr_db_local`, `cazyme_db_local`, `secondary_metabolites_db_local`, `genomad_db_local`, `virsorter2_db_local` | `mkdir`, `ln`, `basename` |

**BacFlux's own helper scripts.** Twelve rules — `select_contigs`,
`screen_final_contigs`, `build_replicons` and nine rules of the mobilome module — run
their script in `workflow/envs/python.yaml`, a pinned Python with nothing else in it.
The scripts import only Python's standard library, so there is no `pip install` step.

The helper scripts' unit tests need nothing but Python and `pytest`, so they run
from the launcher environment. `pytest` is not in the environment created above, so
add it first:

```bash
conda install -n snakemake -c conda-forge pytest
pytest workflow/scripts miscellaneous -q    # three tests skip unless BACFLUX_TEST_* variables name real GTDB files
```

`miscellaneous/` holds the workflow's own guard tests. They read the rule files and
check, for example, that every helper rule names its Python environment and that
every download rule is a local rule.

## 5. Reference databases

Five databases must be downloaded by hand before a run: Bakta, NCBI `core_nt`,
eggNOG, GTDB and Platon. A further six may be pointed at a copy you already hold,
and the rest BacFlux fetches itself on the first run that needs them, or all at
once with the `download_all` target. Sizes, download recipes and the config keys
that hold the paths are all on [Reference databases](databases.md).

Required database keys must exist and point to prepared databases. A missing
key stops configuration parsing, whereas an empty or invalid path may fail later,
when the database is used. For the six optional local-database paths, BacFlux
also checks for an existing directory and database-specific files during planning.
A successful dry run does not establish that every database is complete.

## Next

- [Quick start](quick-start.md) — clone to first result, in one page.
- [Choosing a mode](choosing-a-mode.md) — which of the four front ends your data belongs in.
- [Reference databases](databases.md) — what to download, and where to put it.
- [Running BacFlux](../reference/running.md) — the command line, CPU budget and restarting.
