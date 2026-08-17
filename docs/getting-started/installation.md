# Installation

BacFlux is a Snakemake workflow, not a program you install. The clone holds the
rule files, the helper scripts and the configuration template, and nothing is
written to a system directory. Only Snakemake has to be present before the first
run: every analysis tool the workflow calls — SPAdes, Flye, Bakta, GTDB-Tk,
Platon, VirSorter2 and the rest — is installed by conda while the run is already
under way.

Three things to put in place:

1. the clone,
2. a conda environment holding Snakemake, which is what you launch from,
3. the reference databases BacFlux does not fetch for you — those have their own
   page, [Reference databases](databases.md).

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
| `workflow/envs/` | 31 conda environment specifications, one per tool stack. |

A fifth directory, `.snakemake/`, appears in whatever directory you launch from:
it holds the conda environments Snakemake builds and its own bookkeeping. It is
git-ignored, so it never ends up in a commit.

## 2. Create the launcher environment

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

Snakemake **9.10.1 or newer**, the minimum BacFlux declares. This documentation
was written against 9.14.6.

The environment is a launcher and nothing more. It reads the config, plans the
run and starts jobs; it contains no analysis software. Adding Bakta or SPAdes to
it does not help — a rule runs inside its own conda environment and will not see
anything you put here.

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

- **The first run spends a while before it does any science.** Solving and
  downloading a dozen or more conda environments takes minutes. After that they
  are on disk and the workflow starts straight away.
- **Tools never have to coexist.** GTDB-Tk, antiSMASH, VirSorter2 and MacSyFinder
  each get their own environment, so an upgrade to one cannot break another. The
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
| Any mode | blast, blobtools, quast, checkm, gtdbtk, bakta, antismash, dbcan, eggnog-mapper, abricate, platon, checkv, multiqc, and either virsorter or genomad |
| `illumina` adds | fastp, bowtie, spades, bbmap, qualimap |
| `nanopore` adds | nanoplot, filtlong, flye, dnaapler, medaka, minimap, qualimap |
| `hybrid` adds | everything in the two rows above except minimap (a hybrid run has short reads, so it maps with Bowtie2), plus polypolish and snippy |
| `contigs` adds | minimap — its own front-end rule is plain `awk`, but BlobTools still wants a BAM, so the contigs are mapped against themselves |
| `mobilome.run: true` adds | isescan, macsyfinder, and tncentral for the optional naming and copy-number layers |

## 4. The tools BacFlux expects to find, and does not install

Twenty-three of the workflow's rules declare no conda environment. They run in
the environment you launched Snakemake from, and use whatever is on its `PATH`.
**BacFlux does not install these tools and does not check for them in advance** —
if one is missing, the rule fails with a bare `command not found`, and it fails
after the run has already started.

They fall into four groups.

| Group | Rules | Needs |
|-------|-------|-------|
| Fetch and unpack a database | `download_phix`, `download_amr_db` (CARD), `cazyme_db_download` (dbCAN), `icescan_models` | `wget`, `tar`, `sha256sum`, `awk`, `grep`, `head`, `cut`, `find` |
| Build a symlink view of a database you already hold | `download_amr_db_local`, `cazyme_db_local`, `secondary_metabolites_db_local`, `genomad_db_local`, `virsorter2_db_local` | `ln`, `mkdir`, `basename` |
| Run one of BacFlux's own helper scripts | `select_contigs`, `build_replicons`, and nine rules of the mobilome module | `python` |
| Filter, copy or stage a FASTA | `filter_contigs`, `finalize_contigs`, `stage_qc_genomes` | `awk`, `grep`, `head`, `cp`, `cat` |

Apart from `wget`, the shell commands above are all standard Unix utilities,
present on any Linux system that can run Snakemake at all. The helper scripts
import nothing outside Python's standard library, so the Python that came with
your Snakemake environment already runs them — there is no `pip install` step.

!!! warning "`wget` is the one that actually bites"

    It is genuinely absent from some minimal conda base environments and from some
    HPC login shells, and all four rules that need it are database downloads.

    Check before you launch, and fix it in the launcher environment if it is
    missing:

    ```bash
    conda activate snakemake
    command -v wget tar sha256sum awk python

    # if wget is not there
    conda install -n snakemake -c conda-forge wget
    ```

    Three of those four downloads can be avoided entirely by pointing BacFlux at a
    copy you already hold: `directories.card_db`, `directories.dbcan_db` and
    `mobilome.icescan.dir`. See [Reference databases](databases.md).

The helper scripts' own tests run against that same Python, so they double as a
check that the launcher environment is sound. `pytest` is not in the environment
created above, so add it first:

```bash
conda install -n snakemake -c conda-forge pytest
pytest workflow/scripts -q        # 457 passed, 2 skipped
```

## 5. Reference databases

Five databases must be downloaded by hand before a run: Bakta, NCBI `core_nt`,
eggNOG, GTDB and Platon. A further six may be pointed at a copy you already hold,
and the rest BacFlux fetches itself on the first run that needs them. Sizes,
download recipes and the config keys that hold the paths are all on
[Reference databases](databases.md).

Delete one of the five required keys from your config and the run stops at parse
time, in seconds, rather than an hour in. Leaving it as the empty string the
template ships is not the same thing: the key is there, parsing carries on, and
the rule that needs the path fails later. The six optional keys get the harder
check — BacFlux confirms the directory exists *and* holds a file only that
database has, so a path aimed at the wrong copy is caught before anything runs.

## Next

- [Quick start](quick-start.md) — clone to first result, in one page.
- [Choosing a mode](choosing-a-mode.md) — which of the four front ends your data belongs in.
- [Reference databases](databases.md) — what to download, and where to put it.
- [Running BacFlux](../reference/running.md) — the command line, CPU budget and restarting.
