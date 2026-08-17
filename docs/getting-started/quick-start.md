# Quick start

The shortest path from a clone to a finished genome. Five steps: clone, create the
launcher environment, fill in a config, preview the plan, run it.

The example on this page is **Illumina paired-end reads**, the default `mode`. The
other three modes are the same five steps with a different input directory and one
or two different keys — see [choosing a mode](choosing-a-mode.md).

```bash
# 1. Clone
git clone https://github.com/iLivius/BacFlux.git
cd BacFlux

# 2. Launcher environment (Snakemake only; the analysis tools are installed later,
#    per rule, by conda)
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake

# 3. Make your own copy of the config, then edit the copy
cp config/config.yaml config/config_custom.yaml

# 4. Preview the plan without running anything
snakemake --sdm conda --configfile config/config_custom.yaml --cores 24 -n

# 5. Run it
snakemake --sdm conda --configfile config/config_custom.yaml --cores 24
```

Run these from the repository root. `workflow/Snakefile` is found there
automatically, which is why it never appears on the command line. Paths inside the
config are used as you wrote them, so a relative one resolves against the directory
you launched from, not against the config file — use absolute paths. `output_dir`
is the exception: it is resolved to an absolute path for you.

Steps 1 and 2 in more detail, including which few rules run outside a conda
environment: [installation](installation.md).

!!! note "The config must be named on the command line"

    BacFlux ships no default `configfile:`, deliberately: Snakemake *merges* a
    hard-coded default with the file you pass, so a missing key in your copy would
    be quietly filled in from the shipped example — an empty database path, or a
    download link you did not choose — instead of stopping. Launching without
    `--configfile` gives you:

    ```text
    [BacFlux] No configuration supplied. Pass one explicitly, e.g.:
      snakemake --sdm conda --cores N --configfile config/config.yaml
    ```

## Before step 3 — five databases you supply

Five references are too large, or too tied to one release, for the workflow to
fetch for you. You download them once and give BacFlux the path:

| Config key | Database | Version that is required |
|---|---|---|
| `directories.bakta_db` | Bakta | **v6.0** (Bakta 1.12.1 refuses older ones) |
| `directories.blast_db` | NCBI `core_nt` (or `nt_prok`, named in `parameters.nt_version`), plus the taxonomy files | — |
| `directories.eggnog_db` | eggNOG diamond database | — |
| `directories.gtdbtk_db` | GTDB | **R232** (GTDB-Tk 2.7.2 pins itself to one release) |
| `directories.platon_db` | Platon | — |

Everything else — PhiX, CARD, CheckV, dbCAN, VirSorter2, antiSMASH — is downloaded
by the workflow into `output_dir` on the first run. All of those but PhiX also have
an optional `directories.*_db` key: point it at a copy you already hold and nothing
is fetched. Download recipes, sizes and the six optional keys are on
[reference databases](databases.md).

## Step 3 — the keys that decide a run

Your copy of `config/config.yaml` is long, but almost all of it is tuning that can
stay as shipped. These are the keys an Illumina run cannot start without:

```yaml
mode: illumina

input:
  illumina_dir: /data/reads/illumina        # {sample}_R1.<ext> + {sample}_R2.<ext>

directories:
  output_dir: /data/results/run1            # everything is written here
  bakta_db:   /data/db/bakta/db
  blast_db:   /data/db/NCBI_core_nt
  eggnog_db:  /data/db/eggnog
  gtdbtk_db:  /data/db/release232
  platon_db:  /data/db/platon/db

links:
  phix_link:  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
  card_link:  https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2
  dbcan_link: https://zenodo.org/records/18622157/files/dbcan_db_v5.1.2.tar.gz

resources:
  threads: 16
  ram_gb:  64
```

`dbcan_link` already carries that value in the shipped config; `phix_link` and
`card_link` are empty and **required in `illumina` and `hybrid` mode** — they feed
the PhiX removal step and the read-based CARD screen, neither of which exists in
the other two modes. Every key, with its default and what it trades off, is in
[configuration](../reference/configuration.md).

**Input files are found by pattern, not listed by hand.** BacFlux globs
`illumina_dir` for `{sample}_R1.<ext>` files, rebuilds the R2 path per sample, and
the part before `_R1` becomes the sample name used in every output file:

```text
/data/reads/illumina/
├── PE212-1_R1.fastq.gz
├── PE212-1_R2.fastq.gz
├── PE253-B_R1.fastq.gz
└── PE253-B_R2.fastq.gz          →  samples PE212-1 and PE253-B
```

All files of one kind must share one extension (`fastq`, `fq`, `fastq.gz` or
`fq.gz`); mixing them stops the run. Sample names become filenames and locus tags,
so fifteen characters are refused outright — `*` `#` `@` `%` `^` `/` `!` space `?`
`&` `:` `;` `|` `<` `>` — and everything else, underscores and dots included, is
accepted.

!!! tip "Edit a copy, not the shipped file"

    `config/config_custom.yaml` is the one name `.gitignore` already knows, so a
    copy under that name survives a `git pull` and never lands in a commit. Any
    other name runs just as well but is not ignored. Either way, leave
    `config/config.yaml` alone so you always have the commented template to read.

## Step 4 — the dry run

`-n` (`--dry-run`) makes Snakemake work out everything it would do and then stop.
It costs seconds, and it exercises nearly all of BacFlux's own validation, which
runs while the plan is being built. Skipping it costs whatever the run manages
first — and the first thing a real run does is build conda environments and
download databases.

The header tells you what BacFlux understood. Read it before committing cores:

```text
… BacFlux banner …

Mode: illumina — Genomic analysis of bacterial Illumina reads.
Phage caller: virsorter2 (default). Plasmid stage: Platon-only (geNomad off).
Decontamination mode: auto (discard_no_hit=false).
Using CheckV database from link: 'https://zenodo.org/records/21510554/files/checkv-db-v1.5.tar.gz' (db_id='checkv-db-v1.5').
Sample isolate1 will be processed.
```

Each of those lines is a decision you can still change: the mode, the virus caller,
the decontamination policy, and every database that is being downloaded rather than
read from a copy you hold. When the mobilome module is on it announces itself here
too.

One line per sample means discovery worked; a count that does not match your
sequencing run is a naming mismatch, not a lost file. Below the header, Snakemake
lists the jobs. A single isolate through the default Illumina pipeline, into an
empty output directory, plans **42 jobs**:

```text
Job stats:
job                                count
-------------------------------  -------
AMR_summary                            1
all                                    1
amr_contigs                            8
annotation                             1
…
trim_adapters                          1
viral_quality                          1
virsorter2_db                          1
total                                 42
```

The only count above one is `amr_contigs`: ABRicate runs once per database.

What is checked before the first job starts, and what each failure looks like:

```text
[BacFlux] mode=illumina requires 'input.illumina_dir' to be set in the config.
[BacFlux] input.illumina_dir is not an existing directory: '/nope/reads'
[BacFlux] mode=illumina removes PhiX spike-in reads before assembly, which requires 'links.phix_link' to be set in the config.
[BacFlux] mode=illumina runs the read-based CARD AMR leg, which requires 'links.card_link' to be set in the config.
[BacFlux] directories.checkv_db points at '/data/db/checkv', which is not a directory. …
Missing Illumina mate pair for sample 'PE212-1': /data/reads/illumina/PE212-1_R2.fastq.gz
Sample name 'PE 212' contains unsupported characters.
```

A config with no `mode` at all stops with a `KeyError: 'mode'` from
`workflow/rules/shared/00_common.smk`, rather than guessing a pipeline.

!!! warning "A wrong database path is not always caught at parse time"

    The five required `directories.*_db` paths are read but not opened while the
    plan is built, so a typo in one of them survives the header. Some are then
    caught by the dry run anyway, because they are named as rule inputs — a
    `blast_db` without the NCBI taxonomy dump beside it fails here:

    ```text
    MissingInputException in rule blob_json …
        affected files:
            /data/db/NCBI_core_nt/nodes.dmp
            /data/db/NCBI_core_nt/names.dmp
    ```

    The rest surface when the tool itself runs. A dry run that reaches the job
    table means the plan is sound, not that every database is complete.

## Step 5 — the run

```bash
snakemake --sdm conda --configfile config/config_custom.yaml --cores 24
```

`--cores` is the whole CPU budget and the only knob you need: every CPU-bound rule
declares Snakemake's built-in `threads:`, so `--cores` enforces the ceiling by
itself. **Do not also pass `--jobs`/`-j`** — for a local run it is an alias for
`--cores`, and passing both can let rules oversubscribe the machine. The reasoning
and the measurement are on [running BacFlux](../reference/running.md).

A first run spends a long time before it does any science: a default Illumina run
builds 19 of the 31 conda environments in `workflow/envs/` and downloads the
databases BacFlux fetches for itself. Environments are created under
`.snakemake/conda` in the directory you launched from, so launching from the same
place each time is what stops them being rebuilt.

The run is resumable. If it is interrupted, repeat the same command and Snakemake
redoes only what is missing or out of date. Change a config value and only the
steps downstream of it are recomputed.

## What lands in the output directory

After a successful single-sample run, with the mobilome module off (the default):

```text
/data/results/run1/
├── 01.reads/isolate1/illumina/
│   ├── isolate1_fastp.html               # trimming report; the trimmed FASTQs are temporary
│   └── isolate1_fastp.json
├── 02.assembly/isolate1/
│   ├── contigs_final.fasta               # ← the delivered genome; every later stage reads this file
│   ├── contigs_filt.fasta                # after the length/coverage filter, before decontamination
│   ├── spades/contigs.fasta
│   ├── contaminants/
│   │   ├── contig_taxonomy_decisions.tsv # one row per contig: kept or dropped, and why
│   │   ├── isolate1_composition.txt      # genus composition of the assembly
│   │   ├── isolate1_blastout
│   │   ├── contigs.list
│   │   └── bestscore.blob.blobDB.table.txt
│   └── eval/
│       ├── isolate1_qc_genomes.tsv       # which evaluated genome is which
│       └── quast/  checkm/  qualimap/
├── 03.taxonomy/isolate1/                 # GTDB-Tk classify_wf; the call is in classify/gtdbtk.*.summary.tsv
├── 04.annotation/
│   ├── bakta/isolate1/                   # the annotation: .gff3 .tsv .faa .fna .gbff
│   ├── eggnog/isolate1/
│   ├── antismash/isolate1/               # plus antismash/databases/, fetched once
│   └── dbcan/isolate1/                   # plus dbcan/dbcan_db_v5.1.2/, fetched once
├── 05.amr/
│   ├── abricate/isolate1/                # eight databases, one TSV each, plus AMR_summary.txt
│   └── mapping/isolate1/                 # the read-based CARD leg: covstats, AMR_legend, CARD_report
├── 06.plasmids/isolate1/platon/
├── 07.phages/
│   ├── virsorter/isolate1/  checkv/isolate1/
│   └── checkv_db/  vs2_db/               # fetched once and kept
├── 09.report/multiqc_report.html
└── logs/                                 # one log per job, named rule + sample
```

Stage numbers are the same in all four modes; a stage a mode cannot produce is
simply absent (`01.reads` in `contigs` mode, `05.amr/mapping` wherever there are no
short reads). `08.mobilome` appears only when the module is switched on — see
[turning it on](../mobilome/enabling.md).

Intermediates are declared temporary and deleted as soon as the last rule needing
them has finished, which is why `01.reads` ends up holding only reports and the
PhiX genome and its index disappear entirely. CARD's extracted copy under
`05.amr/card_db/` goes the same way and is fetched again on a fresh run; set
`directories.card_db` to a copy you hold to avoid that.

The files most runs end at:

| Question | File |
|---|---|
| The genome | `02.assembly/{sample}/contigs_final.fasta` |
| Why was that contig kept or dropped? | `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` |
| Is the assembly complete and clean? | `02.assembly/{sample}/eval/checkm/{sample}_checkm_stats.tsv` |
| Assembly statistics | `02.assembly/{sample}/eval/quast/` |
| What species is it? | `03.taxonomy/{sample}/classify/gtdbtk.*.summary.tsv` |
| What genes does it carry? | `04.annotation/bakta/{sample}/` |
| Resistance and virulence, from the contigs | `05.amr/abricate/{sample}/AMR_summary.txt` |
| Resistance, from the reads | `05.amr/mapping/{sample}/{sample}_CARD_report.tsv` |
| Which contigs are plasmids? | `06.plasmids/{sample}/platon/verified_plasmids.txt` |
| One QC page for the whole batch | `09.report/multiqc_report.html` |

The full layout, stage by stage, is in [output files](../reference/output.md).

## Next

- [Choosing a mode](choosing-a-mode.md) — which of the four your data belongs in.
- [Configuration](../reference/configuration.md) — every key, with its default.
- [Running BacFlux](../reference/running.md) — the CPU budget, restarting, where environments and databases live.
- [Troubleshooting](../troubleshooting.md) — the failures that actually happen.
- [Coming from v1](from-v1.md) — if you ran BacFlux, FastaFlux, BacFluxL or BacFluxL+ before.
