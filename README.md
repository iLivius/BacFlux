# BacFlux
A workflow for bacterial whole-genome assembly, QC, annotation, AMR and mobile genetic elements — from Illumina reads, Nanopore reads, both together, or pre-assembled contigs.

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.10.1-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11143917.svg)](https://doi.org/10.5281/zenodo.11143917)
---
```bash
________             _______________              
___  __ )_____ _________  ____/__  /___  _____  __
__  __  |  __ `/  ___/_  /_   __  /_  / / /_  |/_/
_  /_/ // /_/ // /__ _  __/   _  / / /_/ /__>  <  
/_____/ \__,_/ \___/ /_/      /_/  \__,_/ /_/|_|  

BacFlux v2.0.0
Unified genomic analysis of bacterial genomes
Illumina · Nanopore · Hybrid · Pre-assembled contigs

Livio Antonielli, 2026
```

## Synopsis
`BacFlux` is a comprehensive and automated bioinformatics workflow for the processing and analysis of bacterial genomic data. It integrates several powerful tools, each performing a specific task, into a seamless workflow managed by Snakemake.

**v2.0.0 is one workflow with four entry points.** Whatever you sequenced, the data converges on the same downstream analysis, so results are comparable across sequencing strategies:

| `mode:` | Input | Front end |
|---|---|---|
| `illumina` | Illumina paired-end reads | PhiX removal, `fastp`, `SPAdes` |
| `nanopore` | Oxford Nanopore reads | `NanoPlot`, `Filtlong`, `Flye`, `Medaka`, `dnaapler` |
| `hybrid` | Illumina + Nanopore | ONT assembly polished with the short reads (`Polypolish`), plus a `Snippy` comparison of the two assemblies |
| `contigs` | Pre-assembled genomes (FASTA) | contig filtering only — this is the former `FastaFlux` |

From there every mode runs the same shared analysis: decontamination, assembly QC, completeness and contamination, taxonomic placement, functional annotation, secondary metabolites, CAZymes, antimicrobial resistance and virulence screening, plasmid detection, prophage detection, and an aggregated report.

`BacFlux` v2.0.0 also adds an **optional mobilome / AMR-mobility module** (stage `08.mobilome`, **off by default**). For every AMR gene it reports whether the gene sits inside a mobile genetic element and how transferable that element is, on a six-rung ladder from *chromosomal, intrinsic candidate* to *predicted self-transmissible*. See [the mobilome module](#the-mobilome-module).

## Table of Contents
- [What changed in v2.0.0](#what-changed-in-v200)
- [Quick Start](#quick-start)
- [Rationale](#rationale)
- [Description](#description)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running BacFlux](#running-bacflux)
- [The mobilome module](#the-mobilome-module)
- [Licensing and commercial use](#licensing-and-commercial-use)
- [Validation](#validation)
- [Output](#output)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [References](#references)

## What changed in v2.0.0
Read this section if you used `BacFlux` v1.x, `FastaFlux`, `BacFluxL` or `BacFluxL+`.

### Four workflows became one
v1 was four sibling workflows that shared roughly two-thirds of their rules and drifted apart with every change. v2 merges them into a single repository with one Snakefile that dispatches on a `mode:` key in the config:

| v1 | v2 |
|---|---|
| `BacFlux` (Illumina short reads) | `mode: illumina` |
| `FastaFlux` (pre-assembled contigs) | `mode: contigs` |
| `BacFluxL` (Nanopore long reads) | `mode: nanopore` |
| `BacFluxL+` (Illumina + Nanopore) | `mode: hybrid` |

With v2.0.0 the `BacFluxL` and `BacFluxL+` repositories are retired in favour of this one. Their existing DOIs stay valid, so analyses already published with them remain citable and reproducible.

### Your old command
```bash
# v1, short reads
snakemake --sdm conda --jobs 4 --cores 12

# v1, pre-assembled contigs
snakemake --sdm conda --snakefile workflow/FastaFlux --jobs 2 --cores 12
```
becomes, in both cases:
```bash
# v2 — the mode is set in the config file, not on the command line
snakemake --sdm conda --snakefile workflow/Snakefile.v2 --configfile config/config_custom.yaml --cores 12
```

> ⚠ **Your existing `config/config_custom.yaml` is a v1 file and will not run v2.** The v1 config has no `mode:` key, and v2 stops at parse time without one:
>
> ```
> KeyError in file ".../workflow/rules/shared/00_common.smk", line 55: 'mode'
> ```
>
> This is expected, not a bug: `config/config_custom.yaml` is git-ignored, so it is *your* file and a `git pull` never touches it. Start again from the v2 example — `cp config/config_v2.yaml config/config_custom.yaml` — and copy your database paths across by hand. The v2 file is laid out differently (`mode`, `input`, `directories`, `links`, `resources`, `parameters`, `phage`, `mobilome`), so it is worth reading rather than pasting into.

**`--jobs` is gone on purpose, and you should not add it back.** For a local run `--jobs`/`-j` is an alias for `--cores`, not an independent "N jobs of M cores each" setting. Combining the two was measured to allow real CPU oversubscription: two rules each declaring 8 threads, launched with `--jobs 2 --cores 8`, both got their full 8 threads and ran at the same time — 16 real threads against a declared budget of 8. Every CPU-bound rule in v2 declares Snakemake's built-in `threads:`, so **`--cores N` alone is now the complete and correct CPU ceiling**.

### The output directory is renumbered
v1 numbered the shared stages differently in each workflow (taxonomy was `03`, `04` or `09` depending on how many front-end stages came before it). v2 groups all technology-specific work under two fixed parents, so every shared stage has the same number in every mode:

| v1 (short-read) | v2 |
|---|---|
| `01.pre-processing` | `01.reads` |
| `02.assembly` | `02.assembly` |
| `03.post-processing` | `02.assembly/{sample}/eval` and `02.assembly/{sample}/contaminants` |
| `04.taxonomy` | `03.taxonomy` |
| `05.annotation` | `04.annotation` |
| `06.AMR` | `05.amr` |
| `07.plasmids` | `06.plasmids` |
| `08.phages` | `07.phages` |
| — | `08.mobilome` (new, opt-in) |
| `09.report` | `09.report` |

There is no in-place upgrade: point `output_dir` at a fresh directory rather than trying to reuse a v1 one.

### Other changes worth knowing before you run
- **Sample names may now contain underscores.** Three of the four v1 workflows forbade them; v2 allows them everywhere. The exact rule is in [Configuration](#configuration).
- **GTDB moved from R226 to R232, and GTDB-Tk from 2.6.1 to 2.7.2.** GTDB-Tk pins itself to one compatible reference release, so an R226 database will not work with the pinned version. You must download R232; there is no in-place upgrade. Verified on two already-classified genomes: the lineage was identical under both releases.
- **The Bakta database must be v6.0** (Bakta 1.12.0). If you hold a v5.x database, download the new one.
- **The virus caller is selectable.** `VirSorter2` remains the default; `geNomad` is opt-in because it is licensed for academic/non-commercial use only. Choosing `geNomad` also turns on a Platon + geNomad plasmid concordance table.
- **Several databases can now be supplied from a copy you already hold** (`checkv_db`, `vs2_db`, `antismash_db`, `dbcan_db`, `card_db`, `genomad_db`), which skips the download entirely. `BacFlux` only ever reads those paths.
- **The mobilome module is new and off by default.** See [the mobilome module](#the-mobilome-module).
- **`workdir:` is gone.** In v1 Snakemake changed into `output_dir` at parse time. It no longer does, so relative input paths now resolve against the directory you launched from, which is what most people expect. `output_dir` may still be given as a relative path; it is resolved to an absolute one automatically.

What did **not** change: the analysis steps themselves, the thresholds, the decontamination policy, and how the output is to be read. The merge was a consolidation, not a redesign — long-read mode was checked against the published `BacFluxL` baseline and, given the same reads, tool versions, read mode and thread count, produced a byte-identical assembly.

[⬆ Back to Table of Contents](#table-of-contents)

## Quick Start
This gets you started with `BacFlux`. Quick guide:

- Download the latest release:
  ```bash
  # git command
  git clone https://github.com/iLivius/BacFlux.git
  ```

- Install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (if not installed already) and activate the environment:
  ```bash
  # optional, if not installed already
  conda create -c conda-forge -c bioconda -n snakemake snakemake
  # activate Snakemake environment
  conda activate snakemake
  ```

- Configure the workflow. Copy the shipped example config and edit your copy — never the shipped file, so that a `git pull` cannot overwrite your settings:
  ```bash
  # config/config_custom.yaml is the convention here and is git-ignored.
  # Choose another name if you already have one you want to keep.
  cp config/config_v2.yaml config/config_custom.yaml
  ```
  In it you must set: `mode` (which pipeline to run), the input directory for that mode, `output_dir`, and the paths to the databases that are **not** downloaded automatically:

    * bakta_db: path to the [Bakta](https://github.com/oschwengers/bakta?tab=readme-ov-file#database) database directory (**v6.0**)
    * blast_db: path to the [NCBI core nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database directory
    * eggnog_db: path to the [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13#user-content-Installation) diamond database directory
    * gtdbtk_db: path to the [GTDB](https://ecogenomics.github.io/GTDBTk/installing/index.html) database directory (**R232**)
    * platon_db: path to the [Platon](https://github.com/oschwengers/platon?tab=readme-ov-file#database) database directory

- Launch the workflow from the repository root:
   ```bash
  snakemake --sdm conda --snakefile workflow/Snakefile.v2 --configfile config/config_custom.yaml \
            --keep-going --ignore-incomplete --keep-incomplete --cores 24
  ```

  This command uses the following options:

    - --sdm: uses conda for dependency management (one environment per rule, built on first use)

    - --snakefile / --configfile: which workflow and which config to run (see the note below)

    - --keep-going: continues execution even if errors occur in some steps

    - --ignore-incomplete: ignores rules with missing outputs

    - --keep-incomplete: keeps incomplete intermediate files

    - --cores 24: the total CPU budget for the run (adjust as needed). **Do not also pass `--jobs`/`-j`.**

*NOTE ON FILE NAMES: while v2 lives on the `release/v2.0.0` branch, its files are staged under non-conflicting names so the v1 pipeline stays runnable side by side — the Snakefile is `workflow/Snakefile.v2` and the example config is `config/config_v2.yaml`. Passing both explicitly, as above, works either way. Once they take over the plain names (`workflow/Snakefile`, `config/config.yaml`) at the tagged release, `--snakefile` can simply be dropped.*

Refer to the [installation](#installation), [configuration](#configuration) and [running BacFlux](#running-bacflux) sections for detailed instructions.

[⬆ Back to Table of Contents](#table-of-contents)

## Rationale
The analysis of bacterial WGS data often involves a complex series of steps using various bioinformatic tools. Manual execution of this process can be time-consuming, error-prone, and difficult to reproduce. `BacFlux` addresses these challenges by providing a comprehensive and automated Snakemake workflow that streamlines bacterial genomic data analysis.

`BacFlux` integrates several best-in-class bioinformatic tools into a cohesive pipeline, automating tasks from quality control and assembly to annotation, taxonomic classification, identification of antimicrobial resistance genes and viral sequences.

Unifying the four v1 workflows serves the same purpose one level up. The four shared most of their rules, so every fix had to be applied three or four times by hand and the copies drifted. Writing a step once means it behaves identically whether your isolate was sequenced on Illumina, on Nanopore, on both, or arrived as a finished assembly — and it is what makes a cross-technology comparison of the results defensible.

By providing a user-friendly and automated solution, `BacFlux` allows researchers to focus on interpreting the biological meaning of their data.

[⬆ Back to Table of Contents](#table-of-contents)

## Description
Here's a breakdown of the `BacFlux` workflow. Steps 01–02 depend on which mode you run; everything from 03 onwards is identical in all four modes.

01. **Reads (`01.reads`)** — *mode-specific; empty in `contigs` mode.*

    * **Illumina** (`illumina`, `hybrid`): checks raw reads for Illumina phiX contamination using [bowtie2](https://github.com/BenLangmead/bowtie2); filters and removes adapters with [fastp](https://github.com/OpenGene/fastp).

    * **Nanopore** (`nanopore`, `hybrid`): read QC before and after filtering with [NanoPlot](https://github.com/wdecoster/NanoPlot); length/quality filtering with [Filtlong](https://github.com/rrwick/Filtlong). In `hybrid` mode the Illumina reads that survived decontamination are used as Filtlong's reference, so the long reads are scored on agreement with the short-read evidence.

02. **Assembly, decontamination and QC (`02.assembly`)**

    * **Assembly**, per mode:
        - `illumina`: [SPAdes](https://github.com/ablab/spades), then contigs filtered on minimum length (≥500 bp) and coverage (≥2x).
        - `nanopore`: [Flye](https://github.com/mikolmogorov/Flye), reoriented to a conventional start position with [dnaapler](https://github.com/gbouras13/dnaapler) and polished with [Medaka](https://github.com/nanoporetech/medaka).
        - `hybrid`: the Nanopore assembly above, corrected with the Illumina reads using [Polypolish](https://github.com/rrwick/Polypolish). A [Snippy](https://github.com/tseemann/snippy) comparison between the ONT and Illumina assemblies is reported alongside.
        - `contigs`: no assembler. The supplied FASTA is filtered — using the same length/coverage cutoffs when the headers are SPAdes-style and carry those numbers, and by length only when they do not, because inventing a coverage filter would silently delete real sequence.

    * **Decontamination** (all modes): local alignments of contigs against the [NCBI core nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database using [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/), reads (or the contigs themselves) mapped back for coverage, and contaminant contigs identified with [BlobTools](https://github.com/DRL/blobtools). Unless configured otherwise (see [configuration](#configuration)), the output is parsed automatically to discard contaminants based on the relative taxonomic composition of the contigs. Every kept-or-discarded decision, with its reason, is written to `contig_taxonomy_decisions.tsv`.

    * **QC** (all modes): mapping evaluation with [QualiMap](http://qualimap.conesalab.org/) (modes with reads); assembly statistics with [Quast](https://github.com/ablab/quast); genome completeness and contamination with [CheckM](https://github.com/Ecogenomics/CheckM) using taxon-specific markers. In `hybrid` mode both the ONT and the Illumina assembly are carried through QC so the two can be compared directly.

03. **Taxonomic analysis (`03.taxonomy`):**
    * Accurate taxonomic placement with [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) using the curated [GTDB](https://gtdb.ecogenomic.org/) reference database.

04. **Annotation (`04.annotation`):**
    * Annotates contigs using [Bakta](https://github.com/oschwengers/bakta) for functional prediction.
    * Provides further functional annotation with [EggNOG](https://github.com/eggnogdb).
    * Annotates Carbohydrate-Active enZYmes (CAZymes) using the standalone version [run_dbCAN](https://github.com/bcb-unl/run_dbcan) of the [dbCAN3](https://bcb.unl.edu/dbCAN2/) annotation tool.
    * Infers secondary metabolites with [antiSMASH](https://github.com/antismash/antismash).

05. **Antimicrobial resistance (`05.amr`):**
    * Maps filtered reads to the [CARD](https://card.mcmaster.ca/) database with [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/). *Short-read modes only (`illumina`, `hybrid`) — it needs reads.*
    * Screens contigs for antimicrobial resistance and virulence genes using [ABRicate](https://github.com/tseemann/abricate) against eight databases.

06. **Plasmids (`06.plasmids`):**
    * Investigates the presence of plasmids with [Platon](https://github.com/oschwengers/platon) and confirms results with a [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) search.
    * If `geNomad` is opted in as the virus caller, its plasmid calls are compared with Platon's in a concordance table.

07. **Prophages (`07.phages`):**
    * Screens contigs for viral sequences with [VirSorter2](https://github.com/jiarong/VirSorter2) (default) or [geNomad](https://github.com/apcamargo/genomad) (opt-in), followed by [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) for refinement.

08. **Mobilome / AMR mobility (`08.mobilome`) — optional, off by default:**
    * For each AMR gene, reports the mobile-element context and a mobility tier. Adds [AMRFinderPlus](https://github.com/ncbi/amr), [ISEScan](https://github.com/xiezhq/ISEScan) and [CONJscan](https://github.com/gem-pasteur/Macsyfinder_models)/[MacSyFinder](https://github.com/gem-pasteur/macsyfinder), and consumes the plasmid and annotation results from the stages above. See [the mobilome module](#the-mobilome-module).

09. **Reporting (`09.report`):**
    * Parses and aggregates results to generate a report using [MultiQC](https://github.com/MultiQC/MultiQC).

[⬆ Back to Table of Contents](#table-of-contents)

## Installation
BacFlux downloads automatically all dependencies and several databases. However, some external databases require manual download before running the workflow.

1. **Download BacFlux:**

    Head over to the Releases section of the repository.
    Download the latest archive file (typically in .zip or .tar.gz format). This archive contains the `BacFlux` Snakefile, the rule modules, the configuration file and the `envs` environment directory.
    Extract the downloaded archive into your desired location. This will create a directory structure with the necessary files and directories.
    Alternatively, download via command line as:
    ```bash
    # git command
    git clone https://github.com/iLivius/BacFlux.git
    ```

2. **Install Snakemake:**

    `BacFlux` relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to manage the workflow execution. Find the official and complete set of instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install Snakemake as a Conda environment:
    ```bash
    # install Snakemake in a new Conda environment (alternatively, use mamba)
    conda create -c conda-forge -c bioconda -n snakemake snakemake
    ```

    *NOTE: a handful of rules run in the environment you launched from rather than in a conda environment of their own — the ones that only download and unpack a database (`wget`, `tar`, `sha256sum`, `awk`) and the small helper scripts of the mobilome module (`python`, standard library only). All of these are present on a normal Linux system and in the Snakemake conda environment, but `wget` in particular is missing from some minimal conda base environments and HPC login shells, where a database download would otherwise fail with a bare "command not found" after the run has already started.*

3. **Databases:**

    While `BacFlux` automates the installation of all software dependencies, some external databases need to be downloaded manually. If you have already installed them, skip this section and go directly to [configuration](#configuration).

    Here are the required databases and instructions for obtaining them.

    * `Bakta` database (**version 6.0** — Bakta 1.12.0 will refuse an older one):
        ```bash
        # Bakta database comes in two flavours. To download the full database, use the following link (recommended):
        wget -c https://zenodo.org/records/14916843/files/db.tar.xz
        tar -xJf db.tar.xz
        rm db.tar.xz

        # alternatively, download a lighter version
        wget -c https://zenodo.org/records/14916843/files/db-light.tar.xz
        tar -xJf db-light.tar.xz
        rm db-light.tar.xz

        # if the AMRFinderPlus db gives an error, update it by activating the Bakta Conda env and running the following command by targeting the Bakta db directory:
        amrfinder_update --force_update --database db/amrfinderplus-db/
        ```
        *NOTE: according to the [source](https://github.com/oschwengers/bakta?tab=readme-ov-file#database) the light version takes 1.3 GB compressed and 3.9 GB decompressed, whereas the full database takes 30 GB zipped and 84 GB unzipped. Note the archives are `.tar.xz` (hence `tar -xJf`), not `.tar.gz` as in earlier releases. The mobilome module reads `{bakta_db}/amrfinderplus-db/` directly, so a working AMRFinderPlus database inside the Bakta database is what makes that module free of extra downloads.*

    * `NCBI core nt` database, adapted from [here](https://gist.github.com/ppflrs/336e49f8ae3843dc06cc3925940f3024):
        ```bash
        # create a list of all core nt links in the directory designated to host the database (recommended)
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_links.list
        
        # alternatively, create a list of nt links for bacteria only 
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt_prok.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_prok_links.list
       
        # download in parallel, without overdoing it
        cat nt*.list | parallel -j4 'rsync -h --progress rsync://{} .'

        # decompress with multiple CPUs
        find . -name '*.gz' | parallel -j4 'echo {}; tar -zxf {}'

        # get NCBI taxdump
        wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        tar -zxvf taxdump.tar.gz

        # get NCBI BLAST taxonomy
        wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'
        tar -zxvf taxdb.tar.gz

        # get NCBI accession2taxid file
        wget -c 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        gunzip nucl_gb.accession2taxid.gz
        ```
        *NOTE: the complete NCBI core nt database and taxonomy-related files should take around 300 GB of hard drive space (September 2025). The `nodes.dmp` and `names.dmp` files from the taxdump must sit in this same directory — BlobTools reads them from there.*

    * `eggNOG diamond` database:
        ```bash
        # the easiest way is to install a Conda environment with eggnog-mapper, first
        conda create -n eggnog-mapper eggnog-mapper=2.1.13

        # activate the environment
        conda activate eggnog-mapper

        # then, create a directory where you want to install the diamond database for eggnog-mapper 
        mkdir /data/eggnog_db
        #replace /data/eggnog_db with your actual PATH

        # finally, download the diamond db in the newly created directory 
        download_eggnog_data.py --data_dir /data/eggnog_db -y
        ```
        *NOTE: the eggNOG database requires ~50 GB of space.*

    * `GTDB` database (**release R232** — GTDB-Tk 2.7.2 hard-pins itself to one compatible reference release and will reject R226 or older):
        ```bash
        #move first inside the directory where you want to place the database, then download and decompress either the full package or the split package version

        # full package
        wget -c https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
        tar xzvf gtdbtk_r232_data.tar.gz
        rm gtdbtk_r232_data.tar.gz

        # split package (alternative)
        base_url="https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/split_package/gtdbtk_r232_data.tar.gz.part_"
        suffixes=(aa ab ac ad ae af ag ah ai aj ak al)
        printf "%s\n" "${suffixes[@]}" | xargs -n 1 -P 12 -I {} wget "${base_url}{}"
        cat gtdbtk_r232_data.tar.gz.part_* > gtdbtk_r232_data.tar.gz
        tar xzvf gtdbtk_r232_data.tar.gz
        rm gtdbtk_r232_data.tar.gz
        ```
        *NOTE: the compressed archive is ~61 GB and the extracted release directory measured 94 GB here. Point `gtdbtk_db` at the `release232` directory itself. GTDB-Tk 2.7.x reads a pre-sketched skani database that GTDB now ships inside the release, so unlike 2.6.1 it no longer builds (and no longer needs space for) its own ~57 GB sketch cache.*

    * `Platon` database:
        ```bash
        # download the database in a directory of your choice 
        wget https://zenodo.org/record/4066768/files/db.tar.gz
        tar -xzf db.tar.gz
        rm db.tar.gz
        ```
        *NOTE: according to the [source](https://github.com/oschwengers/platon?tab=readme-ov-file#database), the zipped version occupies 1.6 GB and 2.8 GB when unzipped.*

    The remaining databases (PhiX, CARD, CheckV, dbCAN, VirSorter2, antiSMASH, and — if opted in — geNomad) are fetched by the workflow itself. If you already hold any of them, point the matching `directories.*_db` key at your copy and nothing is downloaded; see [configuration](#configuration).

[⬆ Back to Table of Contents](#table-of-contents)

## Configuration
Before running `BacFlux`, copy `config/config_v2.yaml` and edit your copy with a text editor. The file is organized in sections: `mode`, `input`, `directories`, `links`, `resources`, `parameters`, `phage` and `mobilome`. Every key carries an inline comment saying which modes read it, so the file itself is the detailed reference; this section covers the choices that need explaining.

- `mode`

    **Required, and the first thing to set.** One of `illumina`, `nanopore`, `hybrid`, `contigs`. It selects the front-end rules and the sample-discovery logic. Keys belonging to other modes are simply ignored.

- `input`

    One directory per technology; a mode reads only what it needs.

    - **illumina_dir** (`illumina`, `hybrid`): paired-end reads, named `{sample}_R1.<ext>` and `{sample}_R2.<ext>`.
    - **nanopore_dir** (`nanopore`, `hybrid`): long reads, named `{sample}_ont.<ext>`.
    - **contigs_dir** (`contigs`): assemblies, named `{sample}.<fasta|fa|fna>`.

    Conditions on the input files:
    1. Read extensions can only be `fastq`, `fq`, `fastq.gz` or `fq.gz`; contig extensions only `fasta`, `fa` or `fna`.
    2. You can provide multiple samples but the extension must be the same for all files of one kind. Don't mix.
    3. In `hybrid` mode a sample must be present in both directories.

    See an example, below:

    ```bash
    # the input dir contains the PE reads of two strains, PE212-1 and PE253-B, respectively
    ahab@pequod:~/data$ ls -lh
    total 1,6G
    -rw-rw-r-- 1 ahab ahab 379M Apr  8 16:48 PE212-1_R1.fastq.gz
    -rw-rw-r-- 1 ahab ahab 385M Apr  8 16:48 PE212-1_R2.fastq.gz
    -rw-rw-r-- 1 ahab ahab 421M Apr  8 16:48 PE253-B_R1.fastq.gz
    -rw-rw-r-- 1 ahab ahab 433M Apr  8 16:48 PE253-B_R2.fastq.gz
    ```

    **Sample names.** Sample names become filenames, wildcards, locus tags and report labels, so a few characters are refused outright — the run stops at parse time, before anything is computed. The rejected characters, anywhere in the name, are:

    `*` `#` `@` `%` `^` `/` `!` (space) `?` `&` `:` `;` `|` `<` `>`

    Everything else is accepted, **including the underscore** (`_`) and the dot (`.`). This is a change from v1, where three of the four workflows forbade underscores. A dotted file name in `contigs` mode such as `my.genome.fasta` is read as sample `my.genome` — the last dot splits off the extension, earlier dots stay part of the name. The check lives in one place and applies identically to all four modes.

- `directories`

    Update paths based on your file system:

    - **output_dir**: This directory will store all output files generated by `BacFlux`, together with every database the workflow downloads for itself (CheckV, CARD, dbCAN, VirSorter2, antiSMASH, and the optional mobilome databases). Reusing the same output directory for subsequent runs avoids downloading them again. It may be given as a relative path and is resolved to an absolute one automatically.

      *NOTE, changed in v2:* the **conda environments** no longer live here. Snakemake creates them under `.snakemake/conda` in the directory you launch from — in v1 that happened to be `output_dir`, because Snakemake changed into it; v2 does not (see [what changed](#what-changed-in-v200)). Launching from the same directory each time is therefore what avoids rebuilding the ~30 environments; pass `--conda-prefix /some/path` if you would rather keep them somewhere specific, for example shared between projects.

    - **bakta_db**: path to either the light or full (recommended) **v6.0** database of `Bakta`.
    - **blast_db**: path to the `NCBI core nt` (recommended) or prokaryotic database only, and related taxonomic dependencies, see [installation](#installation).
    - **eggnog_db**: path to the diamond database for `eggNOG`.
    - **gtdbtk_db**: path to the **R232** release of `GTDB`.
    - **platon_db**: path to the `Platon` database.

    Optional — set any of these to a copy you already hold and the corresponding download is skipped entirely. `BacFlux` only ever *reads* these paths (through a local symlink view; never in place), and never writes to or deletes them:

    - **checkv_db**: the *parent* directory holding the versioned folder, e.g. `/path/to/checkv/` containing `checkv-db-v1.5/`. Worth setting: the official CheckV database lives on `portal.nersc.gov`, which is unavailable often enough to cost real time.
    - **vs2_db**: what `virsorter setup` produced (~10 GB). Read only when `phage.caller: virsorter2`.
    - **antismash_db**: what `download-antismash-databases` produced.
    - **dbcan_db**: a dbCAN database matching the version in `links.dbcan_link`.
    - **card_db**: an extracted CARD database. Read in `illumina` and `hybrid` only.
    - **genomad_db**: the inner `genomad_db/` directory (the one holding `version.txt` and `genomad_db.dbtype`), read only when `phage.caller: genomad`. **Strongly recommended if you opt into geNomad**: its downloader hard-codes `portal.nersc.gov` with no mirror option, so a local copy is the only fallback when that host is down. `version.txt` must be *readable*, not merely present — in a shared database it easily ends up mode 0640 while the big files beside it are world-readable.

- `links`

    Download URLs for the databases the workflow fetches itself. These should work as they are; change them only if a link breaks or to update a database version.

    - [card_link](https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2): the Comprehensive Antibiotic Resistance Database (`CARD`). Read in `illumina` and `hybrid` only.
    - checkv_link: the `CheckV` database. The default points at an unmodified Zenodo mirror of CheckV's own database, because the official NERSC host is frequently unreachable. Ignored when `directories.checkv_db` is set.
    - [dbcan_link](https://zenodo.org/records/18622157/files/dbcan_db_v5.1.2.tar.gz): the `dbCAN` database for `run_dbcan` v5.1.2.
    - [phix_link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015): the PhiX genome reference used by Illumina as a sequencing control. Read in `illumina` and `hybrid` only.
    - genomad_link / genomad_md5: only read when `phage.caller: genomad`. The default is the geNomad authors' own Zenodo copy of their database, for the same hosting reason as CheckV. Change the MD5 whenever you change the link — a mismatch stops the run.

    *NOTE on the dbCAN pin: `dbcan` is deliberately pinned to 5.1.2, with a version-pinned Zenodo copy of its database, for reproducibility. This is a considered choice, not a workaround: run_dbcan's official mirror only ever hosts the current database, and moving to 5.2.x is not a drop-in swap — subcommand arguments and output schema changed, and on the same input 5.2.9 called 205 CAZyme genes against 5.1.2's 315. The ~110 dropped genes were almost all weak single-tool (`dbCAN_sub`-only) hits that 5.2.9 correctly filters through an e-value fix 5.1.2 ignores, so 5.2.9 is arguably more correct — but the numbers differ materially, while the high-confidence (2–3 tool) calls are preserved either way. Power users should read 5.1.2's weak `dbCAN_sub`-only hits with that in mind.*

- `resources`

    In this section you can specify the hardware resources available to the workflow:

  - threads: the machine scale. Each rule caps its own CPU request at a value appropriate to its tool, against this number — nothing needs configuring per rule.
  - ram_gb: max amount of RAM (read by `illumina` and `hybrid`, for SPAdes and the JVM tools).

  These describe your machine; `--cores N` on the command line is what actually enforces the budget at run time. Set `--cores` to the same value as `threads`, and do not pass `--jobs`.

- `parameters`

    1. **Database selection**: `BacFlux` requires specifying the version of the `NCBI nt` database for `BLAST` operations. You can choose between the `core_nt` and `nt_prok` versions. By default the configuration file is set to use the `core_nt` database. For instructions on installing the `BLAST` database, refer to the [installation](#installation).

    2. **eggNOG `--dbmem`** (`parameters.eggnog.dbmem`, default `false`): eggNOG-mapper's annotation phase does random-access lookups into a 39 GB SQLite database and is the slowest tail of a run. Setting this to `true` loads that database wholly into RAM, at a cost of ~42 GB **per concurrent eggNOG job**. `BacFlux` checks at parse time that `resources.ram_gb` can hold at least one such job, and prints the exact `--resources mem_gb=N` flag to add so that Snakemake limits concurrency to what your RAM actually holds.

    3. **Long-read QC** (`parameters.long_read_qc`, `nanopore` and `hybrid`): **the biggest lever on small-plasmid recovery.** Filtlong ranks reads by a combination of length and quality and then discards the bottom of that ranking. A small plasmid cannot produce reads longer than itself, so if length dominates the ranking, its reads sit at the bottom by construction and are deleted first. Measured on a clinical *K. pneumoniae* isolate whose 5,596 bp Col plasmid was missing from the assembly while sitting complete in the Illumina data: of 614 raw ONT reads mapping to that plasmid, 93 survived at `length_weight: 10` (too few for Flye to assemble it) and 602 at `length_weight: 1`. The defaults are therefore `min_length: 1000`, `keep_percent: 95`, `length_weight: 1` (Filtlong's own default). Raising `length_weight` trades small replicons for chromosome contiguity; if you do it, say so in your methods.

    4. **Long-read assembly and polishing** (`parameters.nanopore` / `parameters.hybrid`):

        - `medaka_model`: `auto` infers the model from the basecaller tag in the FASTQ headers. That only works if the reads carry one — Guppy and Dorado stamp it in, but older or re-headered public data often does not, and Medaka then fails with *"Input file did not contain precisely 1 basecaller model reference"*. The remedy is to name the model explicitly. `FALSE` skips Medaka entirely.
        - **A coupling worth knowing:** with `flye_input_mode: auto` (the default) the *assembler* read mode is chosen from the Medaka model name — a model whose name contains **"fast"** switches Flye to `--nano-raw`, anything else gives `--nano-hq`. So setting `medaka_model` changes the assembler, not just the polisher. To pin the model without touching the assembler, set `flye_input_mode` explicitly. (The published `BacFluxL` baselines were produced with `--nano-raw` through exactly this route.)
        - The Medaka model is validated *before* assembly, so a typo fails in seconds with a table of suggestions instead of after the assembler has run for an hour.

    5. **Decontamination controls**: all four modes use a common taxonomy selector for the `BLAST`/`BlobTools` decontamination step. The default behavior is:

        ```yaml
        decontamination:
          mode: auto
          discard_no_hit: true
          include_genera:
          include_genera_by_sample:
          exclude_genera:
          exclude_genera_file:
          sample_overrides:
        ```

        - **mode**: controls how contigs are selected. Use `auto` to keep the most abundant genus inferred from the BlobTools output, `include` to keep only user-defined genera, `exclude` to remove user-defined genera, or `off` to skip genus-based filtering.
        - **discard_no_hit**: when `true`, contigs with no informative BLAST genus assignment are removed. This is the bacterial default.
        - **include_genera** and **exclude_genera**: optional comma- or semicolon-separated genus lists applied to all samples.
        - **include_genera_by_sample**: optional TSV file with columns `sample` and `genus`, useful when different samples in the same run belong to different target genera.
        - **exclude_genera_file**: optional plain-text file listing genera to remove, one per line or as comma-/semicolon-separated values.
        - **sample_overrides**: optional TSV file with columns `sample`, `mode`, `include_genera`, `exclude_genera`, and `discard_no_hit`. Use it to correct only specific samples without changing global settings.

        The selector writes `contig_taxonomy_decisions.tsv` next to `contigs.list`, so each kept or removed contig can be audited.
        In `auto` and `include` modes, the selector also treats selected genus aliases and retained legacy prefixes as equivalent to reduce false removal caused by BLAST/BlobTools assignments across related or recently reclassified genera. For example, `Bacillus` can retain `Paenibacillus`, `Arthrobacter` can retain `Pseudarthrobacter`, `Pseudoarthrobacter`, or `Paenarthrobacter`, and `Burkholderia` can retain `Paraburkholderia`. Alias-based decisions are marked in `contig_taxonomy_decisions.tsv` with reasons such as `auto_genus_alias` or `included_genus_alias`. `exclude` mode remains exact, because broad alias removal could otherwise remove legitimate target contigs. These aliases are heuristic safeguards, not a formal taxonomic reconciliation system.

        > ⚠ **`auto` mode can remove a genuine small plasmid.** It keeps contigs whose assigned genus matches the sample's dominant genus — right for the common case, but a broad-host-range plasmid can be discarded *precisely because* it is mobile and its closest database relatives sit in another genus. A real example: *K. pneumoniae* ATCC BAA-2146 plasmid pMYS (2,014 bp) was removed as *Escherichia*. Its single best BLAST hit was *K. pneumoniae* at 100% identity over the full length — the correct answer — but BlobTools' `bestsum` rule sums bitscores per taxon across all 95 hits, and because *E. coli* is hugely over-represented in `nt`, *Escherichia* summed to 58,709 against *Klebsiella*'s 16,253. The assignment followed database composition, not biology.
        >
        > Check `contig_taxonomy_decisions.tsv` before trusting the *absence* of a plasmid. If small mobile replicons matter to your question, consider `mode: off` or an explicit include list. The same limitation applies in reverse: a genuine contaminant of the same genus but a different species cannot be caught by a genus-level rule at all.
        >
        > **In `hybrid` mode this bites twice.** Only the Illumina reads mapping to the *selected* contigs are passed to Filtlong as its short-read reference, and ONT reads not covered by that reference score as low quality and are discarded before Flye ever sees them. So a contig dropped here can take its long reads with it, and the sequence disappears from the assembly rather than merely from the taxonomy table. This is independent of `long_read_qc` above; both can delete a plasmid, for different reasons, and both should be ruled out.

- `phage`

    `caller: virsorter2` (default) or `genomad`. See [licensing](#licensing-and-commercial-use) before choosing `genomad`: it is licensed for academic/non-commercial use only. Choosing it also turns on the Platon + geNomad plasmid concordance; on the default, geNomad never runs and the plasmid deliverable is Platon's `verified_plasmids.txt`.

- `mobilome`

    Off by default. See [the mobilome module](#the-mobilome-module).

[⬆ Back to Table of Contents](#table-of-contents)

## Running BacFlux
`BacFlux` can be executed as simply as a Snakefile. Please refer to the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for more details.

`config/config_custom.yaml` below is **your own copy of the v2 example**, made with `cp config/config_v2.yaml config/config_custom.yaml` as described in [Quick Start](#quick-start) — it is git-ignored so your edits survive a `git pull`. If you upgraded from v1 you may already have a file of that name; it is a v1 config and will not run v2 (see [what changed](#what-changed-in-v200)). You can equally point `--configfile` at any other copy you have made; just never edit `config/config_v2.yaml` itself.

```bash
# first, activate the Snakemake Conda environment
conda activate snakemake

# navigate inside the directory where the BacFlux archive was downloaded and decompressed

# check what would run, without running it
snakemake --sdm conda --snakefile workflow/Snakefile.v2 --configfile config/config_custom.yaml -n

# launch the workflow
snakemake --sdm conda --snakefile workflow/Snakefile.v2 --configfile config/config_custom.yaml --cores 24
```

The same command runs all four modes — which one you get is decided by `mode:` inside the config file. `BacFlux` prints the mode, the phage caller, the decontamination policy and whether the mobilome module is on as its first lines of output, so a glance at the header confirms you are running what you meant to run.

*NOTE: starting from Snakemake version 8.4.7, the `--use-conda` option has been deprecated. Use `--software-deployment-method conda` or `--sdm conda` instead. While v2 lives on the `release/v2.0.0` branch, `--snakefile workflow/Snakefile.v2` and `--configfile` are both required, because the v1 `workflow/Snakefile` is still present and would otherwise be picked up. At the tagged release the v2 files take the plain names and `--snakefile` becomes unnecessary.*

**CPU budget.** `--cores N` is the one knob, and it is sufficient on its own: every CPU-bound rule declares Snakemake's built-in `threads:`, which `--cores` enforces automatically. Verified after the conversion from the older custom-resource mechanism: with `--cores 24` the per-rule caps are respected (the VirSorter2 database rule stays at 4, ONT QC and Medaka at 8, adapter trimming at 16, the rest 24); with `--cores 8` Snakemake caps every request at 8 by itself, except the rules whose own lower cap correctly wins. **Do not add `--jobs`/`-j`** — for local execution it is an alias for `--cores`, and passing both can silently reintroduce the oversubscription the ceiling exists to prevent.

**Restarting.** Snakemake only redoes what is missing or out of date, so a run interrupted for any reason can simply be relaunched with the same command. Change a config value and only the affected steps are recomputed. Downloaded databases live under `output_dir` and the conda environments under `.snakemake/conda` in the launch directory, so keeping both stable across runs avoids re-downloading and rebuilding them.

[⬆ Back to Table of Contents](#table-of-contents)

## The mobilome module
**Optional. Off by default. Read the [licensing](#licensing-and-commercial-use) section before turning it on.**

### The question it answers
For each AMR gene detected in the genome: **is it embedded in a mobile genetic element, and if so, how transferable is that element?**

That is the distinction regulators actually ask about. EFSA's framing is *intrinsic* versus *acquired* resistance — intrinsic being roughly species-wide, chromosomal and not transferable, acquired being horizontally gained and possibly able to move again. A resistance gene that is part of a species' normal chromosomal repertoire is a very different risk from the same activity sitting on a conjugative plasmid, because only the second can move into another bacterium.

This module produces the **supporting evidence** for that judgement. It does not produce the judgement, it is not "EFSA-compliant" and does not claim to be, and every mobility call is a **prediction** — the confirmatory experiment is a filter or broth mating assay, not software.

### The mobility ladder
Every AMR gene gets one tier, lowest to highest:

| Tier | Context | What it means |
|:--:|---|---|
| **1** | chromosomal, no mobile-element context | intrinsic candidate |
| **2** | an IS sits upstream, pointing at the gene | the IS can supply an outward-reading ("hybrid") promoter that raises expression. **Expression modulation, not mobilisation** — the gene still cannot move |
| **3** | the gene sits between two copies of the same IS | a composite transposon; the whole block can hop within the cell |
| **4** | the gene sits inside a named unit transposon or integron cassette | mobilisable, with a curated architecture we can name |
| **5** | the gene is on a mobilisable plasmid, or inside an IME | it can move to another cell, but only with a helper element |
| **6** | the gene is on a conjugative plasmid, or inside an ICE | **predicted self-transmissible** |

One case is reported separately and never counted as mobilisation: an IS that has landed **inside** the AMR coding sequence, which usually inactivates the gene (`is_inside_amr_cds`).

### Two things to keep in mind about the output
1. **On a fragmented (short-read) assembly, a located IS count is a floor, not a count.** IS elements are the single biggest cause of contig breaks, because multiple identical copies collapse in the assembly graph. The AMR gene and its flanking IS very often land on different contigs — the exact structure you are trying to detect is what destroyed the assembly. Every row therefore carries `dist_to_contig_end`, `is_at_contig_boundary` and `spans_contigs`, and the IS summary reports what fraction of calls sit at a contig end. Read those.
2. **Published IS-detection false-discovery rates are 8–24% even on curated data.** The output is deliberately tiered evidence (`confidence`: high / medium / low) and never a bare count.

### Turning it on
```yaml
mobilome:
  run: true
```
That is the whole switch. It adds two tools — `ISEScan` and `CONJscan`/`MacSyFinder`, both from bioconda — and surfaces `AMRFinderPlus`, which Bakta already installs and which reads its database straight from `{bakta_db}/amrfinderplus-db/`. **No database path needs configuring.** The CONJscan model package is fetched at run time by `macsydata`, so the machine needs internet access on the first run.

`AMRFinderPlus` is run with `--organism` derived from the GTDB-Tk call where a confident mapping exists, which unlocks curated **point mutations** — the intrinsic, chromosomal, non-transferable determinants tier 1 rests on, and something ABRicate structurally cannot see. There is no official GTDB → AMRFinderPlus organism convention (GTDB's own FAQ states there is no direct translation to NCBI taxa), so `BacFlux` uses a hand-checked, project-local table. It is documented and audited, not authoritative: every decision, **including every refusal**, is written to `{sample}_amrfinder_organism_audit.tsv` with its reason, and a genome with no confident mapping simply gets no `--organism` rather than a guess. Re-check the table when you change the GTDB release. `ABRicate` and the CARD read-mapping leg are unaffected and keep running as before — the three AMR legs are complementary, not redundant.

### Optional layers
Each of these is skipped unless you configure it, and each writes a `PROVENANCE.txt` recording the source, fetch date, checksum and sequence count — because none of these endpoints is versioned, and without that there is no way to say later which release a result came from.

| Layer | Config key | What it adds | Licence |
|---|---|---|---|
| **ICEscan models** | `mobilome.icescan.run` + `url`/`dir` | a second MacSyFinder model set run alongside CONJscan, adding the IME and AICE classes and extra integrase profiles. Over 12 curated IMEs, detections went from 3 to 5 — and detections whose called length is within a factor of two of the published length from 1 to 4 | **CC BY-NC-SA 4.0**, non-commercial |
| **TnCentral** | `mobilome.tncentral.url`/`dir` | curated transposon and integron names — this is what makes **tier 4** reachable at all. A curated hit is not an inference and gets architectures right that the pattern rules miss (notably IS*26*, whose copies sit in direct orientation and break the same-orientation rule tier 3 depends on) | "All Rights Reserved"; no terms page exists |
| **ICEberg 3.0** | `mobilome.iceberg.urls`/`dir` | names the ICE/IME candidates CONJscan already found. Changes no gene's tier — it turns "predicted self-transmissible element" into a name you can look up | no licence or terms statement published |
| **ISOSDB copy number** | `mobilome.isosdb.fasta_url` + `family_map_url` | maps reads against an IS database to estimate how many IS copies the assembler collapsed. Turns the "the count is a floor" warning into a number. Changes no gene's tier. *`illumina` and `hybrid` only — it needs reads* | MIT (from the pseudoR repository) |

### What it outputs
Everything lands in `08.mobilome/{sample}/`. The headline deliverable is:

- **`{sample}_amr_mobility.tsv`** — one row per AMR gene, 46 columns: the gene and its AMRFinderPlus evidence (`amrfinder_method`, identity, coverage, `amr_partial_at_contig_end`), its replicon (chromosome or plasmid id), its `mge_context` (`none` / `is_adjacent` / `composite` / `unit_transposon` / `integron` / `ice` / `ime` / `plasmid`), the measured `distance_bp` and orientation, the IS family and flanking counts, the conjugation machinery behind a tier 5/6 call (`relaxase_type`, `mpf_type`, `machinery_intact`), element boundaries where they could be established (`attL`, `attR`, `boundary_method`), the contig-edge honesty flags, and finally `mobility_tier` and `confidence`.
- **`{sample}_amr_mobility_audit.tsv`** — the decision trail: an explicit reason for every gene that got no mobile-element context, every candidate structure rejected, and every confidence cap applied. Every filtering decision in `BacFlux` gets an audit file with a reason column; this is the mobilome's.

Alongside them: `{sample}_amrfinderplus.tsv` and `_amrfinderplus_mutations.tsv`; `{sample}_is_elements.tsv`, `_is_summary.tsv` (with the fraction of IS calls at a contig boundary) and `_is_discarded.tsv`; `{sample}_ice_candidates.tsv` and `_ice_discarded.tsv`; `{sample}_replicon_calls.tsv`; the raw `isescan/` and `conjscan/` output directories; and, when the optional layers are on, the naming tables and audits plus `{sample}_is_copy_number.tsv`.

### Why the mobilome spans three stages
Plasmids and prophages are mobile genetic elements too, so the layout can read as if `08.mobilome/` were "the mobilome" and `06.plasmids/` and `07.phages/` were something else. They are not. The split is by *what question is being asked*, not by whether the element is mobile:

| Stage | Question | Runs |
|---|---|---|
| `06.plasmids`, `07.phages` | *what replicons and prophages are in this genome?* (detection) | always |
| `08.mobilome` | *is each AMR gene in a mobile element, and how transferable?* (interpretation), plus IS / transposon / integron / ICE detection | opt-in |

Stage 08 **consumes** 06 and 07 rather than re-detecting them.

[⬆ Back to Table of Contents](#table-of-contents)

## Licensing and commercial use
`BacFlux` itself is **MIT licensed**, and that does not change. It ships code and URLs; it never vendors a third-party database or model set. When an optional component is licence-encumbered, you download it yourself, under your own agreement with the licensor — the same arrangement `bakta_db`, `gtdbtk_db`, `platon_db` and the CARD link have always had. Share-alike terms attach to distributed source, not to execution, so running such a tool does not affect `BacFlux`'s own licence.

That protects the *workflow's* licence. It does not tell you whether *your* use is permitted, so:

**The default configuration is fully usable commercially.** Every tool on the default path — including VirSorter2, Platon and CheckV — permits commercial use, and the mobilome module is off.

**Turning the mobilome module on means fetching non-commercial models.**
- The **CONJscan** model package (the HMM profiles and system definitions that detect conjugation machinery, from Institut Pasteur / CNRS) is **CC BY-NC-SA 4.0 — non-commercial, share-alike**. It is fetched whenever `mobilome.run: true`, because without it the module can only ever reach tier 4, and a silently degraded module that still looks complete is worse than an honest switch. MacSyFinder itself, the engine, is GPLv3 and unrestricted.
- The optional **ICEscan** models are a fork of CONJScan by the same authors, under the same **CC BY-NC-SA 4.0** terms.
- Commercial users should either leave `mobilome.run: false` or obtain permission from the model authors.

**geNomad is academic/non-commercial only.** Berkeley Lab licenses it for accredited academic institutions; commercial use requires a separate LBNL licence. (Note the bioconda recipe's `BSD-4-Clause` tag is wrong — the raw LICENSE file is what counts.) It is therefore opt-in, with VirSorter2 as the default. Commercial users should stay on the default.

**The optional naming databases carry unclear terms.** TnCentral publishes an "© All Rights Reserved" notice and no terms page at all — commercial use is not merely unresolved, the site does not address it. ICEberg 3.0 publishes no licence, terms or reuse statement anywhere. Both are opt-in and neither is redistributed by `BacFlux`. ISOSDB is the one mobilome database with no redistribution question attached: it ships in the pseudoR repository under MIT.

**ISfinder** requires written authorisation to download and forbids redistribution; the TnCentral endpoints that carry ISfinder content are deliberately not wired into the config.

Finally, no code from EBI's `mobilome-annotation-pipeline` is used here. Three of its scripts are derived from ICEfinder2 and carry CC BY-NC-SA headers, which would be incompatible with MIT. Reading them to understand an approach is fine and is what was done; only conventions were adopted — Sequence Ontology terms, an element ID format, the discard-with-reason pattern, and two numeric thresholds — and conventions are facts, not expression.

[⬆ Back to Table of Contents](#table-of-contents)

## Validation
The unified workflow was checked against the v1 baselines: given the same reads, tool versions, read mode and thread count, long-read mode reproduced the published `BacFluxL` assembly byte for byte. The helper scripts carry unit tests that run outside Snakemake (`pytest workflow/scripts` — 420 tests at the time of writing).

The mobilome module's ICE/IME caller was benchmarked against **ICEberg's curated coordinates**. The numbers below are the current ones. The methods and the reasoning behind the parameters are written up in `docs/` — `methods_att_and_small_plasmids.md` (the *att*-site search), `methods_icescan_union.md` (the second model set) and `mobilome_worked_example.md` (one genome end to end); the per-element benchmark tables are kept with the benchmark run itself, outside this repository.

**ICE pilot — 18 curated ICEs across 14 genera:** 15 detected. Of the 11 with real chromosomal context, 9 were found, recovering a median 0.48 of the curated element, with a median absolute start offset of 5,141 bp and end offset of 14,912 bp. Of the 7 entries where the deposited record *is* the element, 6 were found, recovering a median 0.67.

**The true ceiling is 17/18, not 18/18.** One entry — SXT(HN1) — is a broken benchmark record: a 19 kb deposit with 17 proteins, no relaxase and no T4SS, standing in for a ~100 kb element. No machinery-based tool can find it. Of the two genuine misses, ICE*Psy*10 has no relaxase, no coupling protein and no VirB4 anywhere in its 161 kb (MacSyFinder's own rejected-candidates file records "quorum of mandatory genes required (3) is not reached: 0"), and ICEB2 in *Mycoplasma bovis* has machinery too divergent for the models to cover.

**IME pilot — 12 curated IMEs across 10 genera, with sizes deliberately straddling the module's size floor:** 5 detected, median 0.69 of the curated element recovered. This is the weaker half of the module and is stated as such.

**Negative control — 12 closed genomes with no curated ICEberg entry, 32.6 Mb in total:** 2 calls, 0 confirmed errors. Nine genomes returned nothing, including *E. coli* K-12 MG1655, *P. aeruginosa* PAO1 and *S. aureus* N315; an archaeon returned nothing too, which is a construction check, since the models are bacterial. Of the two calls, one is ICE*Bs1* in *B. subtilis* 168 — a real element ICEberg has not catalogued, called at 20,571 bp against a ~20.5 kb element, with all four anchor classes, a MOBT relaxase, and a boundary resolved by a tRNA-anchored *att* site at a tRNA-Leu, exactly where ICE*Bs1* is known to integrate. The other is an integrase beside a protein Bakta labels only "DNA-binding protein" but which a T4SS_MOBM model hits at i-evalue 1.8e-84 over 98.6% of the profile; whether that locus is a bona fide IME is a domain question, not a software one.

**This is an upper bound of 2 candidates, not a false-positive rate.** "Absent from ICEberg" is not "contains no element" — and that set originally justified including *B. subtilis* 168 on the grounds that ICE*Bs1* is absent from 168 itself, which is **false**: ICE*Bs1* was discovered there. Counting rather than inspecting would have recorded a textbook-correct detection as a false positive. Every call in a negative control has to be examined, never merely counted.

**Head-to-head against the EBI Mobilome Annotation Pipeline**, over the same 40 genomes and scored with the same script. The caveat has to lead: **these are not two independent callers.** Their ICEfinder2-lite runs MacSyFinder with the ICEscan models — the same engine and model set used here for machinery anchors — so agreement at the detection layer is partly tautological. The genuinely independent parts are boundary refinement and the clustering logic. With that said: across the 30 curated elements of the two pilots (18 ICEs + 12 IMEs) both found 17, `BacFlux` found 3 more (ICE*Vfl*Ind1, Tn*4451*, MTn*Pi*10), the EBI pipeline found none that `BacFlux` missed, and 10 were found by neither — 20/30 against 17/30. Of the 17 both found, the class agrees on 14; all three disagreements are ICEberg-curated ICEs that the other pipeline called IMEs, which is precisely the self-transmissible-versus-needs-a-helper distinction the EFSA framing rests on. Call burden is comparable (65 calls / 2.63 Mb here against 64 / 2.88 Mb), so neither buys recall with volume. On the negative set both made 2 calls, and on *B. subtilis* 168 both returned **base-pair-identical coordinates**, the same 60 bp repeat, the same relaxase and the same mating-pair type, by two different boundary algorithms — independent confirmation that ICE*Bs1* is a real, uncatalogued element.

*One number in that paragraph is easy to quote wrongly, so it is spelled out here: the standalone median of 0.67 above is the pilot-wide figure, over all 6 standalone ICEs `BacFlux` detected. The head-to-head figure is 0.94, and it is computed over only the **5** standalone ICEs **both** callers found — the sixth, ICE*Vfl*Ind1 at 0.29, is one the other pipeline missed entirely, so detecting an extra hard element lowers your own pilot-wide median. Quote 0.67 for `BacFlux` alone and 0.94 only alongside the other caller's 0.94 on the same five.*

*The full account is in `docs/methods_ebi_comparison.md`. If you repeat the comparison, note a conversion trap: the EBI pipeline writes the Sequence Ontology term `integron` both for an ICEfinder2-lite IME and for a genuine IntegronFinder integron. Filtering their GFF on the type column alone silently imports unrelated integrons into the ICE/IME score; the class has to be read from the `mobile_element_type` attribute instead.*

**On vmatch**, which reviewers ask about: the *att*-site search here computes exact maximal repeats in plain Python — vmatch's own semantics. A second implementation built on vmatch's actual data structure (prefix-doubling suffix array, Kasai LCP, cross-flank MEM enumeration) was run against it on every real flank window: 35 of 35 windows, 122 repeats, exact set equality. vmatch would find nothing that is not already found. Separately, bioconda's vmatch declares `license: Unknown / OTHER` and vmatch.de is unreachable, so its terms cannot be established — which is weaker ground than a known-restrictive licence. *One caveat on that measurement: unlike every other number in this section, it cannot currently be re-derived from a file on disk. The result is recorded in the commit message of `4a93d89`, but the diagnostic script that produced it was not kept. It should be re-run and archived with the other validation artefacts before it is relied on in a publication.*

### What the benchmark does not show

Two limitations are not visible in the numbers above, and both change how the output should be read.

**A short boundary becomes a wrong tier, not just a wrong coordinate.** The tables above report whether an element was found and how far its edges landed from the curated ones. But the module's actual output is a *tier per AMR gene*, and a boundary that stops short silently demotes every gene beyond it. Measured on a separate set of eight clinical genomes carrying ICEberg-curated elements: of the **53 AMR genes that sit inside a curated ICE interval, only 15 (28%) reached tier 6**, while **30 (57%) came out at tier 1, "chromosomal, intrinsic candidate"** — the opposite conclusion. `bla`IMP-8 on `CP021851.1` is the clearest case: it lies well inside curated Tn*6397*, but 24,007 bp past the right-hand edge our caller drew, so it is reported as an intrinsic candidate. **A tier 1 call on a genome that has any ICE/IME call in it is therefore weaker evidence than a tier 1 call on a genome that has none**, and the distance from a gene to the nearest element call is worth checking by hand before reading "intrinsic candidate" at face value.

**Tiers 4 and 5 are not equally well tested, and tier 4 has never actually been assigned.** Tier 6 and tier 1 carry the benchmark work described above. The middle of the ladder does not:

- **Tier 4** (inside a *named* transposon or integron) depends entirely on the opt-in TnCentral layer, and **no run retained on disk has ever assigned it**. The naming layer has produced exactly one curated hit on real data — `bla`KPC-2 inside Tn*7247* in a clinical *K. pneumoniae* isolate — and that gene scored **tier 6**, because the transposon sits on a conjugative plasmid and the plasmid evidence outranks it. Tier 4 is reached only when a named element is the *strongest* context available, which in practice means a chromosomal one, and no genome tested so far has produced that combination.
- **Tier 5** is exercised, but only by one of its two routes. Every tier 5 row on disk is a gene on a mobilisable plasmid. **No AMR gene in any run has been assigned an IME context**, so the "inside an IME" half of tier 5 is untested end to end — which matters, because the IME pilot is the weaker half of the module to begin with (5 of 12).

Neither gap is a reason to distrust tiers 1 and 6. Both are a reason to treat a tier 4 or tier 5 call as an unvalidated code path rather than a measured one, and to say so if it appears in a dossier.

[⬆ Back to Table of Contents](#table-of-contents)

## Output
The workflow output reflects the steps described in the [description](#description) section. Here's a breakdown of the subdirectories created within the main output folder, along with explanations of their contents. Stage numbers are the same in all four modes; a stage that a mode cannot produce is simply absent.

- `01.reads`: read QC and filtering. `{sample}/illumina/` holds the [fastp](https://github.com/OpenGene/fastp) (v1.0.1) report and trimmed reads; `{sample}/ont/` holds the [NanoPlot](https://github.com/wdecoster/NanoPlot) (v1.46.2) reports before and after [Filtlong](https://github.com/rrwick/Filtlong) (v0.3.1) filtering. Empty in `contigs` mode.

- `02.assembly`: the assembly and everything used to judge it.
    - The delivered genome is always `{sample}/contigs_final.fasta`, whichever mode produced it — every downstream stage reads that one file.
    - Assembler output: `spades/` ([SPAdes](https://github.com/ablab/spades) v4.2.0), or `flye/` ([Flye](https://github.com/mikolmogorov/Flye) v2.9.6) with `fix_start/` ([dnaapler](https://github.com/gbouras13/dnaapler) v1.3.0), `medaka/` ([Medaka](https://github.com/nanoporetech/medaka) v2.2.2) and, in hybrid mode, `polypolish/` ([Polypolish](https://github.com/rrwick/Polypolish) v0.6.1) and `snps/` ([Snippy](https://github.com/tseemann/snippy) v4.6.0, comparing the ONT and Illumina assemblies).
    - **contaminants**: contig selection based on [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) (v2.16.0) search and [BlobTools](https://github.com/DRL/blobtools) (1.1.1) analysis. Check the `composition` text file for a quick overview of the relative composition of your assembly, and `contig_taxonomy_decisions.tsv` for the per-contig decision and its reason.
    - **eval**: [Quast](https://github.com/ablab/quast) (v5.3.0), [CheckM](https://github.com/Ecogenomics/CheckM) (1.2.4) and [QualiMap](http://qualimap.conesalab.org/) (v2.3) output. `{sample}_qc_genomes.tsv` records which genome is which — relevant in hybrid mode, where the ONT and Illumina assemblies are both evaluated.

- `03.taxonomy`: taxonomic placement performed by [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) (v2.7.2) against GTDB R232.

- `04.annotation`: contains the following sub-directories:
    - **bakta**: accurate annotation output by [Bakta](https://github.com/oschwengers/bakta) (v1.12.0).
    - **eggnog**: functional annotation produced by [EggNOG](https://github.com/eggnogdb) mapper (v2.1.13).
    - **antismash**: secondary metabolites inferred by [antiSMASH](https://github.com/antismash/antismash) (v8.0.4).
    - **dbcan**: carbohydrate-active enzyme and substrate annotation by [dbCAN3](https://github.com/bcb-unl/run_dbcan) (v5.1.2).

- `05.amr`: antimicrobial resistance features are investigated with two complementary approaches:
    - **mapping**: reads filtered by [fastp](https://github.com/OpenGene/fastp) (v1.0.1) are mapped to the CARD database (v4.0.1) using [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (v39.33) with minimum identity = 0.99. Mapping results are parsed and features with a covered length of at least 70% are reported in the `AMR legend` file. *Short-read modes only.*
    - **abricate**: the delivered genome is screened for AMR elements and virulence factors using [ABRicate](https://github.com/tseemann/abricate) (v1.2.0) against eight databases. [EFSA thresholds](https://efsa.onlinelibrary.wiley.com/doi/full/10.2903/j.efsa.2023.8323) are applied: hits must show **≥80% identity and ≥70% gene-length coverage** to be considered.

- `06.plasmids`: the genome is screened for plasmid replicons with [Platon](https://github.com/oschwengers/platon) (v1.7) and results verified by BLAST search to avoid false positives. Contigs ascertained as plasmids are reported in `{sample}/platon/verified_plasmids.txt`. If geNomad is opted in, `{sample}/{sample}_plasmid_concordance.tsv` compares the two callers.

- `07.phages`: the genome is screened for viral sequences with [VirSorter2](https://github.com/jiarong/VirSorter2) (v2.2.4) by default — or [geNomad](https://github.com/apcamargo/genomad) (v1.12.0) if opted in — followed by [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) (v1.0.3) for refinement:
    - **virsorter**: following the instructions provided [here](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=3), viral groups (i.e. dsDNA phage, NCLDV, RNA, ssDNA, and lavidaviridae) are detected with a loose cutoff of 0.5 for maximal sensitivity. Original sequences of circular and (near) fully viral contigs are preserved and passed to the next tool.
    - **checkv**: this second step serves to quality control the results of the previous step to avoid the presence of non-viral sequences (false positives) and to trim potential host regions left at the ends of proviruses.

- `08.mobilome`: **only when `mobilome.run: true`.** Per-sample AMR-mobility tables, IS and ICE/IME calls, and their audit files; see [the mobilome module](#the-mobilome-module). Shared database directories fetched for the optional layers (`conjscan_models/`, `icescan_models/`, `tncentral_db/`, `iceberg_db/`, `isosdb_db/`) sit alongside the per-sample directories.

- `09.report`: [MultiQC](https://github.com/MultiQC/MultiQC) (v1.33) is used to parse and aggregate the results of the following tools:
    1. [fastp](https://github.com/OpenGene/fastp) (v1.0.1) — short-read modes
    2. [NanoPlot](https://github.com/wdecoster/NanoPlot) (v1.46.2) — long-read modes
    3. [QualiMap](http://qualimap.conesalab.org/) (v2.3)
    4. [Quast](https://github.com/ablab/quast) (v5.3.0)
    5. [CheckM](https://github.com/Ecogenomics/CheckM) (1.2.4)
    6. [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) (v2.7.2)
    7. [Bakta](https://github.com/oschwengers/bakta) (v1.12.0)

- `logs` and `benchmarks`: one log and one runtime/memory record per rule, sitting alongside the numbered stages.

[⬆ Back to Table of Contents](#table-of-contents)

## Acknowledgements
This work was originally supported by the [Austrian Science Fund (FWF)](https://www.fwf.ac.at/en/) under Project I6030-B.

## Citation
Antonielli, L., Großkinsky, D. K., Koch, H., Trognitz, F., Sanchez Mejia, A., & Nagel, M. (2024). BacFlux: A workflow for bacterial short-read assembly, QC, annotation, and more. Zenodo. https://doi.org/10.5281/zenodo.11143917

Machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff). The DOI above is the *concept* DOI, which always resolves to the latest release; each tagged release also gets its own version DOI on Zenodo, so cite the version you actually ran when that matters.

**Please also cite the tools whose output you actually used.** `BacFlux` is glue: almost all of the science it reports comes from other people's tools and databases, and the mobilome module in particular is an integrator — it decides how to combine other tools' calls, and every underlying detection belongs to somebody else. The full list, with the licensing notes, is in [`CITATIONS.md`](CITATIONS.md).

## References
1. Abby, S. S., Cury, J., Guglielmini, J., Néron, B., Touchon, M., & Rocha, E. P. C. (2016). Identification of protein secretion systems in bacterial genomes. Scientific Reports, 6, 23080. https://doi.org/10.1038/srep23080

2. Alcock, B. P., Huynh, W., Chalil, R., Smith, K. W., Raphenya, A. R., Wlodarski, M. A., Edalatmand, A., Petkau, A., Syed, S. A., Tsang, K. K., Baker, S. J. C., Dave, M., McCarthy, M. C., Mukiri, K. M., Nasir, J. A., Golbon, B., Imtiaz, H., Jiang, X., Kaur, K., … McArthur, A. G. (2023). CARD 2023: Expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51(D1), D690–D699. https://doi.org/10.1093/nar/gkac920

3. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

4. Blin, K., Shaw, S., Vader, L., Szenei, J., Reitz, Z. L., Augustijn, H. E., Cediel-Becerra, J. D. D., de Crécy-Lagard, V., Koetsier, R. A., Williams, S. E., Cruz-Morales, P., Wongwas, S., Segurado Luchsinger, A. E., Biermann, F., Korenskaia, A., Zdouc, M. M., Meijer, D., Terlouw, B. R., van der Hooft, J. J. J., Ziemert, N., Helfrich, E. J. N., Masschelein, J., Corre, C., Chevrette, M. G., van Wezel, G. P., Medema, M. H., & Weber, T. (2025). antiSMASH 8.0: Extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. Nucleic Acids Research, 53(W1), W32–W38. https://doi.org/10.1093/nar/gkaf334

5. Bouras, G., Grigson, S. R., Papudeshi, B., Mallawaarachchi, V., & Roach, M. J. (2024). Dnaapler: A tool to reorient circular microbial genomes. Journal of Open Source Software, 9(93), 5968. https://doi.org/10.21105/joss.05968

6. Bushnell, B. (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. https://escholarship.org/uc/item/1h3515gn

7. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

8. Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. (2024). Identification of mobile genetic elements with geNomad. Nature Biotechnology, 42(8), 1303–1312. https://doi.org/10.1038/s41587-023-01953-y

9. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825–5829. https://doi.org/10.1093/molbev/msab293

10. Challis, R., Richards, E., Rajan, J., Cochrane, G., & Blaxter, M. (2020). BlobToolKit – Interactive Quality Assessment of Genome Assemblies. G3 Genes|Genomes|Genetics, 10(4), 1361–1374. https://doi.org/10.1534/g3.119.400908

11. Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022). GTDB-Tk v2: Memory friendly classification with the genome taxonomy database. Bioinformatics, 38(23), 5315–5316. https://doi.org/10.1093/bioinformatics/btac672

12. Chen, L., Zheng, D., Liu, B., Yang, J., & Jin, Q. (2016). VFDB 2016: Hierarchical and refined dataset for big data analysis--10 years on. Nucleic Acids Research, 44(D1), D694-697. https://doi.org/10.1093/nar/gkv1239

13. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

14. Coluzzi, C., Garcillán-Barcia, M. P., de la Cruz, F., & Rocha, E. P. C. (2022). Evolution of Plasmid Mobility: Origin and Fate of Conjugative and Nonconjugative Plasmids. Molecular Biology and Evolution, 39(6), msac115. https://doi.org/10.1093/molbev/msac115

15. Cury, J., Touchon, M., & Rocha, E. P. C. (2017). Integrative and conjugative elements and their hosts: Composition, distribution and organization. Nucleic Acids Research, 45(15), 8943–8956. https://doi.org/10.1093/nar/gkx607

16. Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

17. De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: Visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149

18. Doster, E., Lakin, S. M., Dean, C. J., Wolfe, C., Young, J. G., Boucher, C., Belk, K. E., Noyes, N. R., & Morley, P. S. (2020). MEGARes 2.0: A database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Research, 48(D1), D561–D569. https://doi.org/10.1093/nar/gkz1010

19. EFSA Panel on Additives and Products or Substances used in Animal Feed (FEEDAP). (2018). Guidance on the characterisation of microorganisms used as feed additives or as production organisms. EFSA Journal, 16(3), 5206. https://doi.org/10.2903/j.efsa.2018.5206

20. EFSA Panel on Biological Hazards (BIOHAZ). (2023). Statement on how to interpret the QPS qualification on 'acquired antimicrobial resistance genes'. EFSA Journal, 21(10), 8323. https://doi.org/10.2903/j.efsa.2023.8323

21. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

22. Feldgarden, M., Brover, V., Haft, D. H., Prasad, A. B., Slotta, D. J., Tolstoy, I., Tyson, G. H., Zhao, S., Hsu, C.-H., McDermott, P. F., Tadesse, D. A., Morales, C., Simmons, M., Tillman, G., Wasilenko, J., Folster, J. P., & Klimke, W. (2019). Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates. Antimicrobial Agents and Chemotherapy, 63(11), e00483-19. https://doi.org/10.1128/AAC.00483-19

23. Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J. G., Haendiges, J., Haft, D. H., Hoffmann, M., Pettengill, J. B., Prasad, A. B., Tillman, G. E., Tyson, G. H., & Klimke, W. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Scientific Reports, 11, 12728. https://doi.org/10.1038/s41598-021-91456-0

24. Guo, J., Bolduc, B., Zayed, A. A., Varsani, A., Dominguez-Huerta, G., Delmont, T. O., Pratama, A. A., Gazitúa, M. C., Vik, D., Sullivan, M. B., & Roux, S. (2021). VirSorter2: A multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome, 9(1), 37. https://doi.org/10.1186/s40168-020-00990-y

25. Gupta, S. K., Padmanabhan, B. R., Diene, S. M., Lopez-Rojas, R., Kempf, M., Landraud, L., & Rolain, J.-M. (2014). ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrobial Agents and Chemotherapy, 58(1), 212–220. https://doi.org/10.1128/AAC.01310-13

26. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

27. Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085

28. Ingle, D. J., Valcanis, M., Kuzevski, A., Tauschek, M., Inouye, M., Stinear, T., Levine, M. M., Robins-Browne, R. M., & Holt, K. E. (2016). In silico serotyping of E. coli from short read data identifies limited novel O-loci but extensive diversity of O:H serotype combinations within and between pathogenic lineages. Microbial Genomics, 2(7), e000064. https://doi.org/10.1099/mgen.0.000064

29. Jia, B., Raphenya, A. R., Alcock, B., Waglechner, N., Guo, P., Tsang, K. K., Lago, B. A., Dave, B. M., Pereira, S., Sharma, A. N., Doshi, S., Courtot, M., Lo, R., Williams, L. E., Frye, J. G., Elsayegh, T., Sardar, D., Westman, E. L., Pawlowski, A. C., … McArthur, A. G. (2017). CARD 2017: Expansion and model-centric curation of the comprehensive antibiotic resistance database. Nucleic Acids Research, 45(D1), D566–D573. https://doi.org/10.1093/nar/gkw1004

30. Kirsch, J. M., Hryckowian, A. J., & Duerkop, B. A. (2024). A metagenomics pipeline reveals insertion sequence-driven evolution of the microbiota. Cell Host & Microbe, 32(5), 739–754.e4. https://doi.org/10.1016/j.chom.2024.03.005

31. Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology, 37(5), 540–546. https://doi.org/10.1038/s41587-019-0072-8

32. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

33. Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

34. Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake (10:33). F1000Research. https://doi.org/10.12688/f1000research.29032.2

35. Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology, 39(5), 578–585. https://doi.org/10.1038/s41587-020-00774-7

36. Néron, B., Denise, R., Coluzzi, C., Touchon, M., Rocha, E. P. C., & Abby, S. S. (2023). MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes. Peer Community Journal, 3, e28. https://doi.org/10.24072/pcjournal.250

37. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292–294. https://doi.org/10.1093/bioinformatics/btv566

38. Oxford Nanopore Technologies. Medaka. https://github.com/nanoporetech/medaka

39. Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., & Hugenholtz, P. (2022). GTDB: An ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Research, 50(D1), D785–D794. https://doi.org/10.1093/nar/gkab776

40. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043–1055. https://doi.org/10.1101/gr.186072.114

41. Ross, K., Varani, A. M., Snesrud, E., Huang, H., Alvarenga, D. O., Zhang, J., Wu, C., McGann, P., & Chandler, M. (2021). TnCentral: A Prokaryotic Transposable Element Database and Web Portal for Transposon Analysis. mBio, 12(5), e02060-21. https://doi.org/10.1128/mBio.02060-21

42. Schwengers, O., Barth, P., Falgenhauer, L., Hain, T., Chakraborty, T., & Goesmann, A. (2020). Platon: Identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. Microbial Genomics, 6(10), mgen000398. https://doi.org/10.1099/mgen.0.000398

43. Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: Rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11), 000685. https://doi.org/10.1099/mgen.0.000685

44. Seemann, T. (2020). ABRicate. https://github.com/tseemann/abricate

45. Seemann, T. Snippy: Rapid haploid variant calling and core genome alignment. https://github.com/tseemann/snippy

46. Wang, M., Liu, G., Liu, M., Tai, C., Deng, Z., Song, J., & Ou, H.-Y. (2024). ICEberg 3.0: Functional categorization and analysis of the integrative and conjugative elements in bacteria. Nucleic Acids Research, 52(D1), D732–D737. https://doi.org/10.1093/nar/gkad935

47. Wick, R. R. Filtlong: Quality filtering tool for long reads. https://github.com/rrwick/Filtlong

48. Wick, R. R., & Holt, K. E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLoS Computational Biology, 18(1), e1009802. https://doi.org/10.1371/journal.pcbi.1009802

49. Xie, Z., & Tang, H. (2017). ISEScan: Automated identification of insertion sequence elements in prokaryotic genomes. Bioinformatics, 33(21), 3340–3347. https://doi.org/10.1093/bioinformatics/btx433

50. Zankari, E., Hasman, H., Cosentino, S., Vestergaard, M., Rasmussen, S., Lund, O., Aarestrup, F. M., & Larsen, M. V. (2012). Identification of acquired antimicrobial resistance genes. The Journal of Antimicrobial Chemotherapy, 67(11), 2640–2644. https://doi.org/10.1093/jac/dks261

51. Zheng, J., Ge, Q., Yan, Y., Zhang, X., Huang, L., & Yin, Y. (2023). dbCAN3: Automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Research, 51(W1), W115–W121. https://doi.org/10.1093/nar/gkad328

[⬆ Back to Table of Contents](#table-of-contents)
