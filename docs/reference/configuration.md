# Configuration

The YAML configuration selects the input mode, paths, resources and analysis
settings. BacFlux has no command-line options of its own for analysis parameters,
so this file records every analysis setting of a run, unless a value is overridden
with Snakemake's `--config`. Keep it with the results, together with the BacFlux
version, input identifiers and database provenance. A configuration alone does not
preserve mutable reference data or guarantee reproducibility.

The shipped template is `config/config.yaml`. Fill in the input directory for your
mode, `directories.output_dir`, the five required database paths and, for
`illumina` and `hybrid`, `links.phix_link` and `links.card_link`. The other empty
values belong to optional sources and layers and can stay empty when unused. This page follows
the file's own order, and every section is tagged with the modes that read it —
`[all modes]`, `[nanopore|hybrid]` and so on. Keys belonging to a mode that is
not running are ignored entirely.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **YAML** | the config file's plain-text format; indentation is significant and tabs are refused |
    | **EFSA** | European Food Safety Authority, whose reporting thresholds the ABRicate analysis applies |
    | **GTDB** | Genome Taxonomy Database, the reference GTDB-Tk places each genome in |
    | **AMR** | antimicrobial resistance |
    | **IS** | insertion sequence, the smallest kind of mobile element |
    | **ICE** | integrative and conjugative element — sits in the chromosome and can move itself to another cell |
    | **IME** | integrative mobilisable element — the same, but needs a helper element to move |
    | **AICE** | actinomycete integrative and conjugative element, a third class |
    | ***att* site** | the short repeat marking an integrated element's ends |
    | **HMM** | hidden Markov model, a statistical profile used to recognise a protein family |

## How the file is loaded

The Snakefile declares **no** `configfile:` of its own, so a config must be named
on the command line:

```bash
snakemake --sdm conda --cores 16 --configfile config/config_custom.yaml
```

Start it as a copy of the template, which `.gitignore` already covers so your
edits survive a `git pull`:

```bash
cp config/config.yaml config/config_custom.yaml
```

!!! warning "Your file is the whole config — nothing fills in a missing key"

    Snakemake *merges* a hard-coded `configfile:` with whatever `--configfile`
    supplies rather than replacing it. Had the Snakefile named the shipped
    template as a default, every key you deleted from your own copy would be
    quietly refilled from the template — with its empty-string placeholders and
    its default download links — and the run would carry on against a path that
    does not exist.

    So there is no base layer. Delete a required key and the run stops naming it,
    which is the intended behaviour. Copy the whole template and edit it; do not
    write a short config from scratch.

Launched with no config at all, the run stops in seconds:

```text
[BacFlux] No configuration supplied. Pass one explicitly, e.g.:
  snakemake --sdm conda --cores N --configfile config/config.yaml
```

A single top-level key can still be overridden per run without editing the file
— `--config mode=contigs` — but anything nested belongs in the file.

**Paths.** Every path in the config, `directories.output_dir` included,
resolves against **the directory you launched Snakemake from** — not against the
config file. Since you are told to run from the repository root, a relative
`output_dir` puts your results inside your clone:

| `output_dir` | where the results go |
|---|---|
| `/data/results/run1` | exactly there |
| `results/run1` | inside the clone, at `BacFlux/results/run1` |

Neither is wrong, but only one of them is usually what you meant. `output_dir` is
the one path BacFlux expands to an absolute form for you, so it behaves
predictably once set; the rest (`input.*_dir`, `directories.*_db`, the
decontamination file paths) are used exactly as written. Give all of them as
absolute paths unless you are certain.

## SPAdes k-mer sizes

```yaml
parameters:
  spades_kmers: auto           # auto | a comma-separated list, e.g. "21,33,55,77"
```

Applies to `illumina` and `hybrid`, the two modes that run SPAdes. Leave it on `auto`
unless you have a reason not to.

**`auto` passes no `-k` at all**, which is what lets SPAdes size the ladder from the read
length it measures for itself:

| read length | ladder |
|---|---|
| ≥ 250 bp | 21, 33, 55, 77, 99, 127 |
| ≥ 150 bp | 21, 33, 55, 77 |
| shorter | 21, 33, 55 |

Those three presets have been in SPAdes since at least v3.9 (2016) and are unchanged in
v4.3.0, the version BacFlux pins.

**Giving an explicit list disables that selection.** The list is then used whatever the
reads look like, and the cost is easy to miss: k-mer coverage is not read coverage but
`read_cov × (L − k + 1) / L`. With 145 bp reads, k = 127 keeps only about
**13%** of your coverage, and SPAdes prints the k-mer coverage it saw for each k in
its own log — worth reading once on your own data:

```text
Assembling dataset ... with K=77
  Average coverage = ...
Assembling dataset ... with K=127
  Average coverage = ...        <- typically a small fraction of the K=21 figure
```

SPAdes assembles iteratively and the final contigs come from the **largest** k, so an
over-long ladder means the graph you deliver is the least supported one in the run.

!!! note "When to set it explicitly"

    Two good reasons: pinning the ladder so an assembly is reproducible against a future
    SPAdes whose defaults have moved, or overriding SPAdes on data whose read length
    misrepresents its usable length. Values must be odd, ascending and below 128 —
    SPAdes' own limit — and are checked before any job starts.

    Whichever path is taken, the ladder that ran is recorded in the assembly log.

## What is checked before anything starts

BacFlux checks configuration structure and selected values while Snakemake reads the
workflow. A dry run (`snakemake -n`) catches these errors and resolves declared inputs;
it does not execute tools or validate all database contents. The checks, in the order
the file lists the keys:

| Key | Check | What happens |
|---|---|---|
| any block (`input`, `parameters`, `parameters.decontamination`, `phage`, `mobilome`, `mobilome.tncentral`, …) | its settings indented beneath it, not a single value such as `off` | Stops naming the block; for `decontamination`, `phage`, `mobilome` and the four mobilome layers it also says what to write instead |
| `mode` | present, and one of the four | Missing: stops with a `KeyError`. Present but unrecognised, or a bare `on`/`off`: stops, listing the four |
| `input.<mode>_dir` | set, and an existing directory | Stops naming the key |
| input files | at least one, all sharing one extension, extension recognised | Stops listing the extensions it found |
| sample names | none of the refused characters listed under [Input](#input-all-modes) | Stops naming the offending sample |
| Illumina samples | every R1 has its R2 | Stops naming the missing mate |
| hybrid samples | present in both input directories | Stops listing samples missing from either side |
| `directories.output_dir` and the five required databases | present | Stops with a `KeyError` naming the key |
| `directories.checkv_db` | holds exactly one `*/genome_db/checkv_reps.faa` | Stops, and says to give the parent directory |
| `directories.vs2_db`, `antismash_db`, `dbcan_db`, `card_db`, `genomad_db` | a directory, holding a file only that database has | Stops naming the missing probe file |
| `directories.genomad_db` | every file readable by you; metadata table has the current column set | Stops listing unreadable files, or explaining the schema mismatch |
| `links.dbcan_link` | set, and a `.tar.gz` | Stops |
| `links.card_link`, `links.phix_link` | set whenever the mode has short reads | Stops naming the key and the mode |
| `links.checkv_link`, `links.genomad_link` | a `.tar.gz` if set | Stops |
| `resources.threads`, `ram_gb`, `runtime_hours` | a whole number, 1 or more | Stops, naming the key and the value |
| `parameters.nt_version` | a folder name, not a bare `on`/`off` | Stops |
| `parameters.eggnog.dbmem` | `resources.ram_gb` can hold one job | Stops, and names the number to raise |
| `parameters.spades_kmers` | `auto`, or odd ascending integers below 128 | Stops, naming the offending value |
| `parameters.decontamination.mode`, `long_read_mode` | one of `auto`, `include`, `exclude`, `off`; a bare `off` counts as `off` | Stops |
| `parameters.decontamination.auto_vote` | `bp` or `contigs` | Stops |
| `parameters.decontamination.mode` or `long_read_mode` set to `include` | a genus to include: `include_genera`, `include_genera_by_sample` or a `sample_overrides` file | Stops |
| `parameters.<mode>.flye_input_mode` | one of `auto`, `nano-raw`, `nano-hq` | Stops |
| `phage.caller` | `virsorter2` or `genomad`; there is no `off` | Stops, and says phage calling has no off switch |
| `mobilome.coverage_profile` | a fraction in `(0, 1]` | Stops, and says it is a fraction, not a percentage |
| `mobilome.icescan.run` | has a `url` or a `dir` to use | Stops |
| `mobilome.isosdb.fasta_url` | `family_map_url` set alongside it | Stops |
| any boolean key | recognisable as true/false | Stops naming the value |

Each stop ends the run with exit status 1, so a wrapper script or a cluster
submission sees the failure. Most messages begin with `[BacFlux]`; a missing
required key is reported by Python as a `KeyError` naming it.

!!! note "Unquoted `on`, `off`, `yes` and `no`"

    Snakemake reads the config as YAML 1.1, which turns an **unquoted** `on`,
    `yes` or `true` into the boolean *true*, and `off`, `no` or `false` into
    *false*, before BacFlux sees the value. A word in quotes (`"off"`) always
    stays a word. Where a bare on/off ends up:

    | Where it is written | What BacFlux does |
    |---|---|
    | a switch: `mobilome.run`, `mobilome.icescan.run`, `discard_no_hit`, `eggnog.dbmem`, `short_read_guidance`, `medaka_model_fallback_auto`, `require_trna_boundary_for_high` | takes it as meant |
    | `medaka_model` | `false` skips Medaka, `true` infers the model |
    | `parameters.decontamination.mode` or `long_read_mode` | `off` counts as `off`; `on`, `yes` or `true` names no policy, and the run stops |
    | any other key that takes a word: `mode`, `phage.caller`, `parameters.nt_version`, `spades_kmers`, `auto_vote`, `flye_input_mode` | stops, saying how YAML read the value and listing the words it accepts |
    | a block written as one value, e.g. `decontamination: off` or `mobilome: on` | stops. Decontamination is switched off with `mode: off` inside its block and the mobilome module with `run: false` inside its block; phage calling has no off switch |

`bakta_db`, `eggnog_db`, `gtdbtk_db` and `platon_db` are passed to their tools as
written and are not checked before the run: a path that does not exist passes a dry
run, and the error appears when Bakta, eggNOG-mapper, GTDB-Tk or Platon starts. A
wrong `blast_db` is caught by a dry run, because the `names.dmp` and `nodes.dmp`
files inside it are declared rule inputs; the `parameters.nt_version` subfolder is
not checked. A successful dry run therefore does not prove that every required
database is usable. The optional
`directories.*_db` overrides are the exception — each is validated the moment it
is set, which is why `config/config_v2_test_examples.yaml` leaves them out of its
four minimal per-mode dry-run configs.

What the run prints back while parsing is worth reading: the active mode, the
phage caller, the decontamination policy, one line per discovered sample, and —
where they apply — the eggNOG memory flag to add and the mobilome module's
banner.

## The keys that define a run

Most of the file can stay at its defaults. These are the ones that describe what
the analysis actually is:

| Setting | Choices |
|---|---|
| `mode` | `illumina`, `nanopore`, `hybrid`, `contigs` |
| `input.illumina_dir` / `nanopore_dir` / `contigs_dir` | where the raw inputs are, per technology |
| `directories.output_dir` | where results go |
| `directories.bakta_db`, `blast_db`, `eggnog_db`, `gtdbtk_db`, `platon_db` | the five databases you download yourself |
| `links.dbcan_link`, and `phix_link` + `card_link` in the short-read modes | the downloads the workflow performs |

---

## Mode `[all modes]`

```yaml
mode: illumina                 # illumina | nanopore | hybrid | contigs
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `mode` | Selects the front-end rules (`workflow/rules/<mode>/*.smk`) and the sample-discovery logic | `illumina` | Matched case-insensitively. Anything other than the four names stops the run at parse time. |

The mode is never inferred from the data. Everything from taxonomy onwards is one shared
set of rules, while the mode also determines which read-dependent analyses and metadata
are available — see [choosing a mode](../getting-started/choosing-a-mode.md).

---

## Input `[all modes]`

```yaml
input:
  illumina_dir: ""             # {sample}_R1.<ext> + {sample}_R2.<ext>
  nanopore_dir: ""             # {sample}_ont.<ext>
  contigs_dir:  ""             # {sample}.<fasta|fa|fna>
```

| Key | Read by | What it holds | Notes |
|---|---|---|---|
| `illumina_dir` | `illumina`, `hybrid` | Paired-end reads named `{sample}_R1.<ext>` and `{sample}_R2.<ext>` | Discovery keys off R1; the missing mate of any sample stops the run. |
| `nanopore_dir` | `nanopore`, `hybrid` | Oxford Nanopore Technologies (ONT) reads named `{sample}_ont.<ext>` | |
| `contigs_dir` | `contigs` | Finished assemblies named `{sample}.<ext>` | |

A mode reads only the directories it needs; the others can stay empty.

**Extensions.** Reads may be `fastq`, `fq`, `fastq.gz` or `fq.gz`; assemblies may
be `fasta`, `fa` or `fna`. Every file of one kind in one batch must share a
single extension, because `{sample}.{extension}` would otherwise be ambiguous —
a mixture stops the run and prints the extensions it found.

**Sample names.** The name becomes a filename, a Snakemake wildcard, a Bakta
locus tag and a label in the report, so a handful of characters are refused
outright, anywhere in the name:

`*` `#` `@` `%` `^` `/` `!` (space) `?` `&` `:` `;` `|` `<` `>`

Everything else is accepted, underscores included. In `contigs` mode the *last*
dot splits off the extension, so `my.genome.fasta` is read as sample
`my.genome`; earlier dots stay part of the name and flow safely through paths and
locus tags.

---

## Directories `[all modes]`

```yaml
directories:
  output_dir: ""
  bakta_db:   ""
  blast_db:   ""
  eggnog_db:  ""
  gtdbtk_db:  ""
  platon_db:  ""
  # optional local copies, below
```

### Output

| Key | What it does | Notes |
|---|---|---|
| `output_dir` | Root for every result and log, plus every database the workflow downloads for itself | Created if absent. May be relative — it is the one path resolved to an absolute one automatically. |

Reusing one output directory across runs is what stops the self-downloaded
databases being fetched again. Conda environments do **not** live here: Snakemake
puts them under `.snakemake/conda` in the directory you launch from, or in the
folder `--conda-prefix` names, so launching from the same place each time is what
avoids rebuilding them. See [running BacFlux](running.md).

### The five databases you provide

Required in every mode. None of them has a fallback behind it, so a missing key
stops the run rather than quietly resolving to something wrong.

| Key | What it points at |
|---|---|
| `bakta_db` | A Bakta database, light or full. Also supplies the AMRFinderPlus database at `amrfinderplus-db/latest/`, which the mobilome module uses. |
| `blast_db` | The directory holding the NCBI nucleotide database; the subfolder named by `parameters.nt_version` is what the contamination screen blasts against. |
| `eggnog_db` | The eggNOG-mapper database directory. BacFlux pins eggnog-mapper 2.1.15, which annotates against **eggNOG 5.0**, database release **`emapperdb-5.0.2`**. The version is fixed by the tool rather than chosen by you: `download_eggnog_data.py` builds its URL from `__DB_VERSION__` in the installed release, so the downloader requests the release associated with that tool version. |
| `gtdbtk_db` | A GTDB-Tk release directory. |
| `platon_db` | The Platon plasmid database, release **v1.5.0** ([DOI 10.5281/zenodo.4066768](https://doi.org/10.5281/zenodo.4066768), 2020-10-05). It is versioned separately from Platon itself: Platon's documentation lists database v1.5.0 for Platon 1.5.0 and later, including the 1.8 pinned here, so the database number trailing the tool number is expected rather than a sign of drift. |

Versions and download recipes are on
[reference databases](../getting-started/databases.md).

### Optional local copies of the downloaded databases

Six databases the workflow can fetch itself. Point the matching key at a copy you
already hold and **nothing is downloaded** — the key overrides both the `links.*`
entry and, where the tool has one, its own downloader. These paths are only ever
read: BacFlux builds a local symlink view rather than touching your directory.

| Key | Point it at | Read when | Probe file checked at parse time |
|---|---|---|---|
| `checkv_db` | The **parent** directory holding the versioned folder, e.g. `/path/to/checkv/` containing `checkv-db-v1.5/` | all modes | exactly one `*/genome_db/checkv_reps.faa` |
| `vs2_db` | What `virsorter setup` produced | `phage.caller: virsorter2` | `Done_all_setup` |
| `antismash_db` | What `download-antismash-databases` produced | all modes | `clusterblast/` |
| `dbcan_db` | An extracted dbCAN database matching `links.dbcan_link` | all modes | `dbCAN.hmm` |
| `card_db` | An extracted CARD database | `illumina`, `hybrid` | `aro_index.tsv` |
| `genomad_db` | The `genomad_db/` directory `genomad download-database` produced | `phage.caller: genomad` | `version.txt` |

!!! note "Why CheckV gets a local view rather than being used in place"

    CheckV needs a DIAMOND index (`genome_db/checkv_reps.dmnd`) that the official
    archive does not ship — it is built locally after unpacking. Whatever index
    sits in a shared database was built by whichever DIAMOND that site happened
    to have, and DIAMOND refuses to run a format it does not recognise: CheckV
    then dies at *"Running DIAMOND blastp search… DIAMOND task failed"*, after
    the contamination stage has already succeeded, which reads like a CheckV bug.
    A shared database is also usually not writable, so rebuilding in place is not
    an option. BacFlux therefore symlinks the big read-only files and builds the
    index with its own DIAMOND — about 950 MB and a couple of minutes, once per
    output directory.

!!! warning "geNomad's database gets two extra checks"

    A shared copy can be *partly* readable: some small files, `version.txt`
    among them, can end up mode 0640 while the big files beside them are
    world-readable. A file can be listed and still refuse to open, so checking
    that it exists is not enough: the run would start and die minutes later on a
    permission error. Every file is therefore checked for read access, and the
    offenders are named.

    The second check is on the database *schema*. A geNomad database is coupled
    to the geNomad release, and geNomad does not verify this: it parses
    `genomad_marker_metadata.tsv` positionally from the last columns, so an older
    database shifts every field by one and geNomad dies on
    `ValueError: invalid literal for int() with base 10: '1398618at2'` — a marker
    accession being read as a hallmark count. Verified on 2026-07-24 with
    database v1.7 against geNomad 1.12.0. BacFlux tests for the trailing
    `PREVIOUS_MARKER_ACCESSION` column instead of guessing at version numbers.

---

## Links `[all modes]`

```yaml
links:
  phix_link:   ""
  dbcan_link:  "https://zenodo.org/records/18622157/files/dbcan_db_v5.1.2.tar.gz"
  card_link:   ""
  checkv_link: "https://zenodo.org/records/21510554/files/checkv-db-v1.5.tar.gz"
  genomad_link: "https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz"
  genomad_md5:  "67244b528bb8bed464d1ca147136d33e"
```

| Key | Required in | Default | Notes |
|---|---|---|---|
| `phix_link` | `illumina`, `hybrid` | unset | The PhiX control genome, which may be present in Illumina sequencing data. Those reads are real sequence from another organism and are mapped out before assembly. |
| `dbcan_link` | all modes | dbCAN 5.1.2 on Zenodo | Must be a `.tar.gz`: the checksum URL and the extracted folder name are both derived from it, so the folder always matches the release you pointed at. |
| `card_link` | `illumina`, `hybrid` | unset | CARD, for the read-mapping AMR analysis. Read mapping can detect sequences missing from an assembly, but remains subject to coverage and mapping ambiguity. |
| `checkv_link` | *optional*, all modes | a Zenodo mirror of CheckV's own database | Empty means CheckV downloads its own. Ignored when `directories.checkv_db` is set — and BacFlux says so rather than ignoring the key silently. |
| `genomad_link` | *optional*, `phage.caller: genomad` | the geNomad authors' own Zenodo copy | Must be the database archive (`genomad_db_v*.tar.gz`), not the HMM or MSA archive beside it on the same record. Empty falls back to geNomad's own downloader. |
| `genomad_md5` | with `genomad_link` | the MD5 of the default archive | Zenodo publishes an MD5 rather than the `.sha256` sidecar the other mirrors carry, so it is configured rather than derived. **Change it whenever you change the link.** A mismatch stops rule `genomad_db` before the archive is extracted; an empty value downloads without verification and says so in `logs/genomad_db.log`. |

The CheckV and geNomad mirrors exist for one reason: both databases are served
from `portal.nersc.gov` ([NERSC](https://www.nersc.gov/)), which is unreachable
often enough to cost real time (it was down for the whole of 2026-07-22), and
geNomad's downloader has no option to point elsewhere. Several of these
databases carry terms of their own; BacFlux ships none of the data, and the
licence of each is tabulated on [licensing](../about/licensing.md).

---

## Resources `[all modes]`

```yaml
resources:
  threads: 16
  ram_gb:  64
  runtime_hours: 24
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `threads` | How many cores the machine has to give. Every CPU-bound rule builds its own request against this number. | `16` | Nothing needs configuring per rule. On a cluster, at most the CPUs of one node: the assemblers ask for all of it. |
| `ram_gb` | The memory budget handed to the rules that can be told a limit. On a cluster, also the memory each heavy step asks the scheduler for. | `64` | Not enforced by Snakemake on a single computer — see below. On a cluster, at most one node's memory, and at least about 60, because GTDB-Tk alone needs 55 GB. |
| `runtime_hours` | Cluster only: the time limit each long step asks the scheduler for, in hours | `24` | At most your queue's limit. A step that runs over is stopped by the scheduler; raise the value and run again. Ignored on a single computer, and changing it re-runs nothing. |

All three must be whole numbers, 1 or more. Anything else — `2.5`, a word, `true`,
`0` — stops the run at parse time, naming the key, rather than being rounded or
misread. A key left out takes its default.

These describe your machine. What actually enforces the CPU budget at run time is
`--cores N` on the command line; set it to the same value as `threads`, and see
[running BacFlux](running.md) for why `--jobs` must not be added alongside it. On
a cluster, `ram_gb` and `runtime_hours` become each heavy job's memory and time
request — see [On a cluster](running.md#on-a-cluster).

### How `threads` reaches a rule

Rules ask either for all of `threads`, or for `min(threads, N)` where their tool
stops scaling at N. Both go into Snakemake's built-in `threads:`, which is what
makes a plain `--cores N` a real ceiling.

| Request | Rules |
|---|---|
| all of `threads` | `map_phix`, `illumina_assembly`, `ont_assembly`, `map_contigs`, `map_sel_contigs`, `map_final_contigs` |
| `min(threads, 24)` | `blast_contigs`, `blast_final_contigs`, `genome_assembly_evaluation`, `completeness_and_contamination`, `map_evaluation`, `taxonomic_assignment`, `annotation`, `functional_annotation`, `secondary_metabolites_analysis`, `cazyme_gene_cluster`, `map_amr_db`, `plasmid_search`, `viral_identification_virsorter2`, `viral_quality`, `genomad_end_to_end`, `fix_start`, `short_read_correction`, `compare_hybrid_assemblies`, `long_read_consensus` (nanopore) |
| `min(threads, 16)` | `trim_adapters`, `assembly_depth`, `isosdb_map` |
| `min(threads, 8)` | the NanoPlot rules, `long_read_consensus` (hybrid), `checkv_db`, `checkv_db_local`, and the mobilome rules `amrfinderplus`, `isescan`, `conjscan`, `icescan`, `tncentral_blast`, `iceberg_blast` |
| `min(threads, 4)` | `virsorter2_db` |
| none declared | both ABRicate rules, and the small helper steps — BacFlux's own Python scripts, which run in `envs/python.yaml` |

### How `ram_gb` reaches a rule

| Rule | What it gets |
|---|---|
| `illumina_assembly` | SPAdes `-m`, a hard ceiling, uncapped |
| `map_evaluation` | Qualimap's Java heap, the most memory Java may use: `ram_gb` minus 4 GB, at most 64 GB and at least 1 GB. The 4 GB are left for what Java needs beyond its heap |
| `map_amr_db`, `assembly_depth`, `isosdb_map` | BBMap's `-Xmx`, `min(ram_gb, 32)` GB |
| `functional_annotation` | 42 GB, only when `parameters.eggnog.dbmem` is on |
| the heavy rules, on a cluster | `ram_gb` × 1024 MB as the job's memory request (`mem_mb`); Qualimap asks for its heap plus 4 GB, never more than `ram_gb` |

!!! warning "`ram_gb` is a number handed to tools, not a budget Snakemake polices"

    Snakemake enforces the CPU side on its own, because rules declare the
    built-in `threads:`. Memory figures are gigabytes rather than core counts, so
    they sit in `resources:` under three different names, and are scheduled
    against **only** if the launch line passes the matching one — `--resources
    mem_gb=N` for eggNOG, `java_mem=N` for Qualimap, `ram=N` for SPAdes. Without
    that, two large jobs running at once can each take their full share. eggNOG
    with `--dbmem` is the one case where this matters enough that BacFlux prints
    the exact flag to add. `java_mem=N` also lowers each Qualimap heap to N GB,
    the way to give Qualimap less without lowering `ram_gb` for everything else.

    The three BBMap analyses (`map_amr_db`, `assembly_depth`, `isosdb_map`) take their
    `-Xmx` from a `params:` value rather than a resource, so `--resources` does not
    change it.

    The heavy rules, the BBMap analyses among them, also declare `mem_mb` and
    `runtime`. A cluster executor reads both. A single computer ignores `runtime`
    and uses `mem_mb` only when the launch line passes `--resources mem_mb=N`; each
    heavy step then counts as `ram_gb` × 1024 MB against N.

---

## Parameters

### `nt_version` `[all modes]`

```yaml
parameters:
  nt_version: core_nt          # core_nt | nt_prok
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `nt_version` | Which subfolder of `directories.blast_db` the contamination screen blasts against | `core_nt` | `nt_prok` is the prokaryote-only build: smaller and faster, and enough when contaminants are expected to be bacterial. A config omitting the key falls back to `core_nt`. |

### eggNOG `--dbmem` `[all modes]` *— optional*

```yaml
  eggnog:
    dbmem: false
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `dbmem` | Loads the 39 GB `eggnog.db` wholly into RAM instead of reading it from disk | `false` | Costs about 42 GB **per concurrent eggNOG job**. |

eggNOG-mapper's annotation phase does random-access lookups into that SQLite
database once per seed ortholog, which is why `functional_annotation` is
routinely the last rule still running. `--dbmem` makes those lookups in-memory.

Two things happen when you switch it on. BacFlux refuses at parse time if
`resources.ram_gb` cannot hold even one such job, naming the number to raise. And
it prints the exact `--resources mem_gb=N` flag, with your own `ram_gb` filled
in, so Snakemake will cap concurrent eggNOG jobs at what your memory actually
holds. Without that flag the run still works; concurrency is then bounded by
cores rather than by RAM. More on the annotation stage:
[annotation](../analysis/annotation.md).

### Long-read QC `[nanopore|hybrid]`

```yaml
  long_read_qc:
    min_length: 1000
    keep_percent: 95
    length_weight: 1
```

These settings can affect whether small-plasmid reads survive filtering. Filtlong ranks
reads by

```text
(Length^length_weight × MeanQuality^mean_q_weight)^(1/(lw+mqw)) × WindowQuality
```

and `--keep_percent` deletes the bottom of that ranking. Short plasmid-derived reads can
be disadvantaged by length weighting. Concatemers can produce reads longer than a single
plasmid copy, so this is a tendency rather than an absolute length limit.

| Key | What it does | Default | Reaches Filtlong in |
|---|---|---|---|
| `min_length` | Reads shorter than this are discarded outright | `1000` | both long-read modes |
| `keep_percent` | Keep the best N% of **bases** by the score above | `95` | both long-read modes |
| `length_weight` | How much length dominates the ranking | `1` (Filtlong's own default) | `hybrid` only |

Measured on *Klebsiella pneumoniae* TUM24772 (BioProject PRJNA1168299,
closed genome GCA_043950115.1) — a clinical test isolate, not one of the mobilome
positive controls. Its 5,596 bp Col2 plasmid (CP171790.1) was missing from the
assembly while sitting complete in the Illumina data. Of the 603 raw ONT reads
mapping to it (minimap2 `-x map-ont`, alignment block ≥ 1 kb), the number
surviving each setting pair:

| | `length_weight: 10` | `length_weight: 1` |
|---|--:|--:|
| `keep_percent: 90` (the old value) | 93 | 413 |
| `keep_percent: 95` | 500 | 502 |

In this sample `keep_percent` had the larger effect, so it is the first key to change
when a small plasmid's reads are missing. At the default of 95, changing
`length_weight` is worth two reads; at the old 90 it was worth 320. The plasmid
needed both old values together to disappear, and either change on its own brings
most of it back.

`length_weight` still defaults to 1 because the two keys protect different size
classes. Reads per length band, same sample (the raw set holds 2,954 reads of
1–3 kb and 3,418 of 3–6 kb):

| | 1–3 kb | 3–6 kb |
|---|--:|--:|
| `keep_percent: 90`, `lw: 10` | 0 | 318 |
| `keep_percent: 95`, `lw: 10` | 0 | 2,666 |
| `keep_percent: 90`, `lw: 1` | 123 | 2,331 |
| `keep_percent: 95`, `lw: 1` | 605 | 2,935 |

In this sample, `length_weight: 10` removed the 1–3 kb reads at both tested
`keep_percent` values. Raising `keep_percent` retained more 3–6 kb reads, whereas
lowering `length_weight` also retained some 1–3 kb reads. These observations do not
define universal plasmid-size ranges protected by either setting.

!!! warning "Recovering the reads is necessary, not sufficient"

    None of the settings tried assembled that plasmid — raw unfiltered reads included.
    ONT ligation prep under-represents small circular plasmids roughly four-fold,
    and no filter setting undoes that half of the problem. Read this table as
    "these settings stop the reads being thrown away", not as "and then it
    assembled".

There is no spelling for "off". The number goes straight to Filtlong's
`--keep_percent`, so a blank value arrives as the literal `None` and Filtlong
refuses to start. A value of `100` rejects no read by this percentage; `--min_length`
still applies, as does the 500 Mb `--target_bases` cap (`nanopore`, and `hybrid` with
`short_read_guidance: false`) or, with guidance on, trimming and splitting against the
Illumina reads. Raising either key further rarely gains anything — total
bases barely moved in the test (208 → 227 Mb), so what changes is *which* reads
survive, not how many. If you raise `length_weight` you are trading small
replicons for chromosome contiguity; say so in your methods. Full treatment:
[nanopore mode](../modes/nanopore.md).

### Decontamination `[all modes]`

```yaml
  decontamination:
    mode: auto                 # auto | include | exclude | off
    auto_vote: bp              # bp | contigs
    discard_no_hit: true
    long_read_mode: "off"      # hybrid only: off | auto | include | exclude
    include_genera:
    include_genera_by_sample:
    exclude_genera:
    exclude_genera_file:
    sample_overrides:
```

After BlobTools assigns every contig a genus, the selector keeps or discards
contigs by this policy and writes an audit row for every decision. One policy
covers all four modes, and it runs on the draft assembly, before annotation. In
`illumina`, `contigs` and `nanopore` that draft becomes the delivered genome, so a
contaminant it removes never reaches CheckM or Bakta. In `hybrid` it screens the
Illumina draft; the delivered genome is built from the long reads and is only
judged, a second time, with `long_read_mode`, as a check that changes no contig.

!!! note "`mode` here is the filtering policy, not the pipeline mode"

    `parameters.decontamination.mode` and the top-level `mode` share a name and
    nothing else.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `mode` | `auto` keeps the one genus that wins the vote set by `auto_vote`; `include` keeps only the genera you list; `exclude` removes the ones you list; `off` keeps everything | `auto` | |
| `auto_vote` | How `auto` picks its genus: `bp` adds up the length of each genus's contigs and keeps the genus with the most base pairs, the one most of the assembly belongs to; `contigs` keeps the genus with the most contigs | `bp` | The two agree on a clean isolate. They differ when one genus sits in a few long contigs and another in many short ones — the normal shape of a long-read assembly, where a chromosome plus a few plasmids that BLAST assigns to another genus would outvote the chromosome on a contig count. An exact tie goes to the genus first in alphabetical order, so the answer never changes between runs. The genus Bakta writes as the organism name is counted the same way, over the kept contigs. See [how `auto` picks the genus](../analysis/decontamination.md#how-auto-picks-the-genus). |
| `long_read_mode` | `hybrid` only: the policy of the second screen, which judges the delivered genome (`off`, `auto`, `include`, `exclude`) | `off` | A check, not a filter: it sets the verdict written for each delivered contig in `{sample}_final_contig_taxonomy_decisions.tsv`, and `contigs_final.fasta` is never changed. Uses the same genus lists and `auto_vote` as `mode`. See [the second screen in `hybrid`](../analysis/decontamination.md#the-second-screen-in-hybrid). |
| `discard_no_hit` | Removes contigs BLAST could not place | `true` in the shipped template | Omitting the key entirely resolves to `false`, so keep the line. Applied *before* the mode logic. In `auto` it also drops "no-hit" from the vote that picks the dominant genus. |
| `include_genera` | Genera to keep, for every sample | unset | An inline list: `"GenusA;GenusB"` or `[GenusA, GenusB]`. Commas and semicolons both separate. |
| `include_genera_by_sample` | **A file path** to a two-column TSV (`sample`, `genus`) | unset | The compact form when many samples each need their own include list. A sample may appear on several rows. |
| `exclude_genera` | Genera to discard, for every sample | unset | Same shapes as `include_genera`. |
| `exclude_genera_file` | **A file path** to a genus list, one per line | unset | Blank lines and `#` comments are skipped. |
| `sample_overrides` | **A file path** to a per-sample override TSV | unset | Columns `sample` and `mode` are required; `include_genera`, `exclude_genera` and `discard_no_hit` may follow. A row replaces the run-wide setting for that one sample. |

`include_genera_by_sample`, `exclude_genera_file` and `sample_overrides` are
**paths to files**, not mappings written inline. A YAML block
written there would be stringified and the selector would try to open a file
literally named `{}`.

!!! warning "An empty cell in the override TSV means *keep the global value*"

    Not "match nothing". A row must still carry the full set of tabs where fields
    are blank, or the columns shift and the wrong value is read.

!!! warning "Taxonomic selection can discard a genuine plasmid"

    BlobTools' `bestsum` rule gives each contig the genus with the largest sum of
    BLAST bit scores over all its hits, so the call follows database composition as
    well as similarity. Plasmid pMYS of *K. pneumoniae* ATCC BAA-2146 (2,014 bp;
    NZ_CP006660.1) had *K. pneumoniae* as its single best hit, at 100% identity over
    its full length, but its hits summed to 58,709 for *Escherichia* against 16,253
    for *Klebsiella*. These numbers are bit-score sums, not counts of database
    entries. `auto` kept *Klebsiella* and removed the plasmid
    (`not_auto_genus:Klebsiella`). Broad-host-range plasmids are exposed to this,
    because their closest database relatives often belong to another genus. In
    `illumina`, `nanopore` and `contigs` mode a removed contig is absent from the
    delivered genome.

    In hybrid mode, this selection also affects the Illumina reads used by
    Polypolish. It affects Filtlong only when `short_read_guidance: true`.
    Guidance is disabled by default, and the final hybrid screen removes nothing.

    Check `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` for
    removed contigs of plasmid size whose `assigned_genus` differs from the
    sample's genus. To keep them, set `mode: include` naming both genera,
    `mode: off`, or a `sample_overrides` row for that isolate.
    `discard_no_hit: false` cannot rescue such a contig: it only switches off the
    no-hit prefilter, after which a contig BLAST could not place is still removed by
    `auto` and `include` and kept by `exclude`. Long-read filtering
    (`long_read_qc`, above) can lose a plasmid independently; rule out both.

The selector's alias handling, the audit file's reason codes and the two routes
by which a real plasmid is lost are on
[decontamination](../analysis/decontamination.md).

!!! note "A config written before 2.0.0"

    A single top-level `parameters.genus` is still translated into
    `mode: include` with that genus, so older configs keep working. It is not
    part of the current schema and is not documented beyond this note.

### Long-read assembly and polishing `[nanopore|hybrid]`

```yaml
  nanopore:                    # read only when mode == nanopore
    flye_input_mode: auto      # auto | nano-raw | nano-hq
    medaka_model: auto
    medaka_model_fallback_auto: false
  hybrid:                      # read only when mode == hybrid
    flye_input_mode: auto
    medaka_model: auto
    medaka_model_fallback_auto: false
    short_read_guidance: false # hybrid only: optional Illumina guidance for Filtlong
```

The block is namespaced per mode: a run reads only the block named after itself.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `flye_input_mode` | Which read model Flye assembles under | `auto` | `auto` gives `--nano-hq` unless an explicit *fast* Medaka model is named. |
| `medaka_model` | Which Medaka consensus model polishes the assembly | `auto` | `auto`/`true`/empty infer it from the FASTQ headers; `false` skips Medaka entirely *(set `false` to skip this optional step)*; anything else is taken as a model name. |
| `medaka_model_fallback_auto` | What to do when an explicit model is invalid | `false` | `false` stops with model suggestions. `true` attempts automatic inference. |
| `short_read_guidance` | Hybrid only: use selected Illumina pairs to guide Filtlong, enabling trimming and splitting | `false` | With `false`, Filtlong uses long-read length/quality filtering and the target-base cap, and nothing screens the long reads before Flye; the second screen of the final assembly only reports. With `true`, no `--target_bases` cap is passed, and reads without support from the decontaminated Illumina reads are trimmed, split or discarded. Polypolish uses the selected Illumina reads in either case. |

!!! warning "The two keys are coupled, and it changes the assembly"

    With `flye_input_mode: auto`, an explicit `medaka_model` whose name contains
    **"fast"** also switches the *assembler* to `--nano-raw` — fast basecalling
    means noisier reads. So naming a fast model changes how the genome is
    assembled, not only how it is polished. To pin the polisher and leave the
    assembler alone, set `flye_input_mode` explicitly instead of leaving it on
    `auto`.

!!! warning "`auto` needs a basecaller tag in the reads"

    Inference reads the basecaller name out of the FASTQ headers, so it fails
    outright on reads that carry none — common in public or re-headered data.
    Name the model explicitly for such a dataset.

The model is validated right after read filtering, before the assembler runs, so
a typo fails in seconds rather than an hour in. The same check downloads the
model's weights, on the computer that runs Snakemake, so the polishing job itself
never needs internet ([On a cluster](running.md#cluster-medaka)). Details and the
model table: [nanopore mode](../modes/nanopore.md).

---

## Phage caller `[all modes]`

```yaml
phage:
  caller: virsorter2           # virsorter2 | genomad
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `caller` | Which tool finds viruses and prophages in stage `07.phages` | `virsorter2` | CheckV then scores whatever that caller found, so this changes CheckV's input, not the QC itself. Phage calling cannot be switched off: one of the two callers always runs, and `caller: off` stops the run. |

Choosing `genomad` changes the plasmid stage too. One geNomad run calls viruses
and plasmids together, so opting in also switches on the Platon + geNomad
concordance in `06.plasmids`; on the default, that stage is Platon alone and
geNomad never runs. Read [licensing](../about/licensing.md) before choosing, and
[prophages](../analysis/phages.md) for what the choice does to the output.

---

## Mobilome module `[all modes]` *— optional*

Everything in this block is inert while `run` is `false`, which is the shipped
default. Switching it on adds stage `08.mobilome`: for every AMR gene the
pipeline found, the mobile-element context around it and a mobility tier, from
1 (chromosomal, intrinsic candidate) to 6 (inside an ICE or on a conjugative
plasmid — predicted self-transmissible).

### The switch, and what it costs

```yaml
mobilome:
  run: false
  max_composite_span_bp: 20000
  coverage_profile: 0.5
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `run` | Master switch for the whole module | `false` | Adds 22 rules, two tools (ISEScan and CONJscan/MacSyFinder, both from bioconda), and a model download that needs network access on first use. Those models are non-commercial — see [licensing](../about/licensing.md). |
| `max_composite_span_bp` | Two IS copies further apart than this are not called one composite transposon | `20000` | A convention, not biology. The measured distance is always reported next to the call, so a reader can disagree without re-running anything. |
| `coverage_profile` | The fraction of an HMM profile a protein's alignment must cover before the hit is kept | `0.5` | MacSyFinder's own default, and the value at which every validation result in this documentation was measured; the sweep below also measured 0.4 and 0.3. |

**Nothing else has to be set.** The module runs on the databases BacFlux already
has, plus the CONJscan models it fetches itself. See [turning it
on](../mobilome/enabling.md).

!!! note "`coverage_profile` is not identity and not an E-value"

    A protein can score beyond statistical doubt against a profile and still be
    dropped by this rule alone, because it aligned to only part of the profile.
    That happens for two opposite reasons: it really is a fragment or a decayed
    remnant, or it is a full-length member of a divergent family sharing only the
    catalytic core.

Lowering it was measured end to end at 0.5 / 0.4 / 0.3 over three benchmark sets
— 18 curated ICEs, 12 curated IMEs, and 12 genomes with no curated element at all
(32.6 Mb), with the ICEscan layer on:

| `coverage_profile` | curated ICEs | curated IMEs | calls on the negative set |
|---|--:|--:|--:|
| 0.5 | 15 of 18 | 5 of 12 | 5 (2 without ICEscan) |
| 0.4 | 15 of 18 | 6 of 12 | 7 (4 without ICEscan) |
| 0.3 | 15 of 18 | 6 of 12 | 9 (5 without ICEscan) |

Across 30 curated elements, lowering it recovers exactly one more: a 23 kb IME in
*Faecalibacterium duncaniae*, found in full and correctly classed. Everything
added on the negative controls at 0.4 is a `cime_or_island` at `passive`
mobility, `low` confidence, claiming nothing and carrying no tier — the tiering
absorbs it, which is what it is for. At 0.3 that stops being true: a 21.6 kb IME
is called in *Staphylococcus aureus* N315 at medium confidence with a real
mobility claim. Nothing measured so far supports 0.3 for single isolates. If your
isolates sit in the space where it might help — of 395 curated IMEs, 73 carry a
relaxase hit that is clean on E-value and fails only this rule, 60 of them
`T4SS_MOBT` at about 0.32 coverage — 0.4 is worth trying, and the audit TSV
records a reason for every dropped candidate so you can see what changed.

Both machinery searches use this one value, and must: their hit tables are
merged, so different stringencies would make an element's class depend on which
model set was more permissive.

### The four optional layers

Each needs a source of its own before it does anything, and none of them switches
on merely because `run` is true.

| Layer | Switched on by | What it adds | Modes |
|---|---|---|---|
| `icescan` | `icescan.run: true` **and** a `url` or `dir` (the shipped `url` counts) | A second machinery model set, merged with CONJscan's hits | all |
| `tncentral` | a `url` or `dir` | Curated transposon and integron names — the only thing that makes ladder tier 4 reachable | all |
| `iceberg` | a `urls` entry or `dir` | Names the ICE/IME candidates CONJscan already found | all |
| `isosdb` | `fasta_url` + `family_map_url`, or `dir` | Read-derived IS abundance for comparison with located calls | `illumina`, `hybrid` |

Every layer writes a `PROVENANCE.txt` beside what it fetched, recording the
source, the fetch date, the observed checksum and how much arrived — sequences
for the three databases, HMM profiles for the ICEscan models. None of these URLs
pins a release, so without that file there is no way to say later which one a
result came from.

#### `icescan` — a second machinery model set

```yaml
  icescan:
    run: false
    url: "https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/icefinder2lite/icf2_dbs.tar.gz"
    sha256: ""
    dir: ""
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `run` | Runs MacSyFinder a second time with the ICEscan models and merges the two hit tables | `false` | Switched on with no `url` and no `dir` stops the run — silently producing today's results would hide the fact that the layer never existed. |
| `url` | Where to fetch the model package from | ICEfinder2's database bundle (~61 MB) | Only the MacSyFinder models are extracted; the rest is discarded. Use the `https://` form — the `ftp://` spelling of the same path times out from behind many firewalls. |
| `sha256` | Pins which release you analysed with | unset | Optional but recommended: the URL carries no version. Empty skips the check; the observed checksum is recorded either way. |
| `dir` | A MacSyFinder models directory you already hold, containing an `ICEscan/` folder | unset | Takes precedence over `url`; nothing is downloaded. |

It is a fork of CONJscan by the same authors, one minor version behind the
release BacFlux installs, which is why it runs *alongside* rather than instead:
swapping would lose the MOB relaxase models, the decayed-machinery models and the
whole plasmid set. What it adds is the AICE class, extra integrase profiles, and
many more IME rows — 10 without it against 21 with it, over the 28 benchmark
genomes run both ways. What it does not add is boundaries: MacSyFinder reports
gene ordinals, never base pairs, so every coordinate still comes from BacFlux's
own annotation join, clustering and att-site search. The measurements are on
[validation](../mobilome/validation.md).

#### `tncentral` — curated transposon and integron names

```yaml
  tncentral:
    url: ""
    sha256: ""
    dir: ""
    min_identity: 90.0
    min_reference_coverage: 0.8
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `url` | A ZIP holding `tncentral.fa` — about 512 curated transposons, integrons and IS | unset (layer off) | The server rejects curl's default user-agent, so the download sends a browser one; if that ever stops working the rule stops with an error rather than building an empty database. |
| `sha256` | Pins the release | unset | The endpoint is unversioned. Recorded either way in `tncentral_db/PROVENANCE.txt`. |
| `dir` | A directory already holding `tncentral.fa` | unset | Beats `url`; nothing is fetched. |
| `min_identity` | Percent identity a BLAST hit needs before it may confer a curated **name** | `90.0` | Below this the sequence may well be a relative of the element, but the name would claim more than the data shows. |
| `min_reference_coverage` | Fraction of the **reference element** that must be present | `0.8` | Measured against the reference, not your contig: a 7 kb element inside a 300 kb contig covers 2% of the contig and 100% of itself. A fragment of a transposon is not that transposon. |

Tier 3 is inferred from a pattern: two copies of one IS family within
`max_composite_span_bp`, with the gene between them and both copies in the same
orientation. The orientation test is not applied to the IS6 family (IS26), whose
copies are accepted either way round; the mobility audit then records
`is26_orientation_exemption_applied`. A TnCentral hit instead matches a
characterised, named element: it places the gene inside a known transposon or
integron, but does not show that the element is still active.

!!! warning "Tier 4 is strict, and rarer than it looks"

    Tier 4 needs the TnCentral layer (`url` or `dir`), a hit of at least
    `min_identity` covering at least `min_reference_coverage` of the reference
    element, at least 90% of the gene inside that element, and no ICE, IME or
    plasmid context claiming the gene first. A gene claimed first keeps the
    element's name in `named_element`, with the audit reason
    `inside_named_element_but_higher_tier_applies`. On ATCC BAA-2146 one gene
    reaches tier 4, *bla*<sub>CTX-M-15</sub> inside Tn*Ecp1.1* at 87% of the
    reference; the six other named elements sit on plasmids
    ([Worked example](../mobilome/worked-example.md)). On *Enterobacter hormaechei*
    fifteen rows come from one chromosomal transposon, Tn*SMR478*
    ([Hybrid mode on five closed genomes](../methods_reference_genome_run.md)), so a
    tier-4 count reflects elements as much as genes.

    On fragmented clinical assemblies the same Tn*Ecp1.1* was refused at 12% and
    49% reference coverage, below the threshold, so expect tier 4 far more often on
    closed genomes than on drafts.

#### `iceberg` — which ICE is it?

```yaml
  iceberg:
    urls: []
    dir: ""
    min_identity: 80.0
    min_overlap_fraction: 0.5
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `urls` | One or more FASTA URLs, concatenated into one BLAST database | `[]` (layer off) | A single URL written without a leading dash is accepted too. Note the host: the older `bioinfo-mml.sjtu.edu.cn` path serves the previous ICEberg release, not the current one. The larger file is ~100 MB and the server is slow. |
| `dir` | A directory of the downloaded ICEberg `.fas` files, or the `iceberg_db` folder an earlier BacFlux run built | unset | Beats `urls`; nothing is fetched. Checked when rule `iceberg_db` runs, not at start-up. |
| `min_identity` | Identity floor for a name | `80.0` | Below this the elements are related but not the same, and the name would mislead. |
| `min_overlap_fraction` | How much of **our** candidate the curated element must cover before the two are treated as the same thing | `0.5` | Deliberately lenient next to the transposon cascade: ICEs are mosaic and their cargo varies between strains, so demanding near-completeness would refuse to name exactly the divergent elements a name helps with most. |

ICEberg has no configurable `sha256` key. The concatenated FASTA checksum and fetch
date are recorded in `iceberg_db/PROVENANCE.txt`. This records the data used but does
not ensure that a later fetch from a mutable URL will return identical bytes. Keep a
versioned local copy when repeatability requires it.

This layer only labels; it cannot change any gene's tier. Its real value is
showing how far our boundaries fall short: on *K. pneumoniae* ATCC BAA-2146 the
ICE call spans 54,943 bp against ICEberg's 58,048 bp for the same element — 0.946
of it, stopping 3,138 bp inside its far end. Any AMR gene in that last 3 kb is
scored as though it were outside the ICE. That is the *good* case, a closed
genome where the att search found a tRNA-anchored repeat.

#### `isosdb` — estimating IS abundance { #isosdb-how-many-is-copies-did-the-assembly-lose }

```yaml
  isosdb:
    fasta_url: ""
    family_map_url: ""
    dir: ""
    min_covered_percent: 90.0
    min_copy_number: 0.5
```

`[illumina|hybrid]` only — it needs reads.

| Key | What it does | Default | Notes |
|---|---|---|---|
| `fasta_url` | 22,713 IS nucleotide sequences, already dereplicated at 95% identity | unset (layer off) | About 9.9 MB zipped. |
| `family_map_url` | Entry → IS family (IS3, IS5, IS110 …) | unset | **Required whenever `fasta_url` is set.** `fasta_url` without it stops the run at parse time, rather than the download failing an hour in; `family_map_url` on its own is ignored. |
| `dir` | A directory already holding both files | unset | Beats the URLs. |
| `min_covered_percent` | A database entry must be covered this far end to end before its depth is believed | `90.0` | A partially covered entry is usually a conserved domain shared with another family; averaging it in would inflate every estimate. |
| `min_copy_number` | Below this multiple of the genome baseline the element is treated as absent | `0.5` | Deliberately under 1.0: a real single-copy IS sits near 1×, and sampling noise plus mapping loss routinely push it to 0.6–0.8×. |

The located IS count underestimates the copy number when repeated copies are not
resolved in the assembly. The ISOSDB analysis estimates each family's copy number
from read depth, relative to the length-weighted median depth over the sample's own
contigs, and compares it with the located count. The difference estimates the
copies the assembly collapsed rather than counting them: mapping ambiguity, related
entries and coverage variation also affect it. `collapse_delta` is reported only
when ISOSDB detected the family and the estimate is at least the located count;
otherwise it is `NA` and `db_informative` is `FALSE`. A family with no close ISOSDB
match therefore never yields a negative number: the estimate can show that copies
were collapsed, never that ISEScan over-called.

The result does not change an AMR gene's tier or locate additional copies. Family
summaries reduce sensitivity to which related entry receives an ambiguous read,
but do not remove every mapping or database limitation.

### Confidence and boundary flags

```yaml
  require_trna_boundary_for_high: false
  contig_boundary_bp: 100
```

| Key | What it does | Default | Notes |
|---|---|---|---|
| `require_trna_boundary_for_high` | Also demands a tRNA-anchored att pair before an ICE/IME candidate may be labelled `high` confidence | `false` | Sensible on closed long-read assemblies, where an unresolved boundary really is a warning sign. |
| `contig_boundary_bp` | How close to a contig end counts as "at the boundary" when flagging an **IS** call | `100` | Feeds the `at_contig_boundary` and `fraction_at_contig_boundary` columns of the IS summary. A high fraction indicates that many IS calls occur at contig ends, where context or further copies may be missing. |

Left `false`, confidence answers "how sure are we this **is** an ICE?" — from the
number of anchor classes, whether the machinery is intact, and whether it all
sits on one contig. Where the element *ends* is reported separately, in
`boundary_method`. Set `true` on a fragmented assembly and nearly every call caps
at medium however good the machinery evidence is, because the flanking sequence
simply is not in the contig; that conflates two different questions. Either way,
cargo is never assigned from an unresolved boundary: a de novo direct repeat is
reported but never widens an element (about 16% of arbitrary spans on a real
chromosome produce one by chance), so the interval stays the machinery span — a
floor, never an invention.

!!! note "Three contig-end windows, one of them configurable"

    They answer different questions about different objects, which is why they
    are separate:

    | Window | Question | Set by |
    |---|---|---|
    | 100 bp | Is an **IS call** sitting on the edge of its contig? | `contig_boundary_bp` |
    | 1000 bp | Was there enough flanking sequence to have **seen** a neighbouring element at all? Caps confidence in the mobility table. | fixed in `colocalise.py` |
    | 1000 bp | The same question for an ICE or IME candidate rather than an AMR gene. | fixed in `conjscan_to_ice.py` |

    An IS 500 bp from a contig end is not "at the edge", but you still could not
    have seen its partner 2 kb away — absence of evidence is not evidence of
    absence. The first and third confusingly share the flag name `--boundary-bp`
    while meaning different things; only the first is wired to a config key.

What each column of the output means is on
[reading the output](../mobilome/output.md); what to conclude on a fragmented
assembly is on [draft assemblies](../mobilome/draft-assemblies.md).

---

## Numbers that are not config keys

Some thresholds you might reasonably expect to find in the config are set in the
code instead. This is where each one lives.

### The ICE/IME geometry

`conjscan_ice` passes none of these, so the script's own defaults are what every
measured result quoted on this site was produced at. Changing one means editing
the rule's shell block in `workflow/rules/shared/80_mobilome.smk` — and re-running
the benchmark afterwards.

| Default | What it controls |
|--:|---|
| 15,000 bp | How far apart two machinery genes may sit and still be chained into one candidate element |
| 8,000 bp | Floor on a candidate's length; below it the cluster is dropped with a reason |
| 2,000 bp | The same floor for the IME/AICE architecture, which is genuinely smaller. A separate empirical floor for small machinery clusters |
| 500,000 bp | Ceiling, so a runaway chain of anchors cannot cover a whole replicon |
| 50,000 bp | How far from the machinery an integrase may sit and still be taken as the element's own |
| 50,000 bp | How much sequence either side of the machinery the att-site search reads |

!!! note "Different size floors for admission and boundary search"

    The att search uses an 8,000 bp floor, while IME/AICE machinery clusters can
    be admitted at 2,000 bp. Repeat pairs implying elements shorter than 8 kb are
    rejected. A machinery span below 8 kb does not necessarily imply an element
    of the same length, but this difference limits boundary resolution for
    genuinely small elements. An IME whose only *att* pairs imply less than 8 kb keeps
    `boundary_method: none` and is reported at its machinery span. Four detected
    curated IMEs have machinery spans of 3,840, 4,959, 5,942 and 6,252 bp. Inspect
    their `boundary_method` rather than assuming that detection resolved their
    complete extent.

### Elsewhere in the workflow

| Value | Where | What it is |
|---|---|---|
| `minid=0.76`, ≥70% of reference length covered | `CARD_MIN_IDENTITY` / `CARD_MIN_COVERED`, `00_common.smk` | The CARD read-mapping AMR analysis. Read identity is BBMap's default; the reference-length coverage gate limits partial matches. See below. |
| 80% identity, 70% coverage | `rule amr_contigs`, `50_amr.smk` | ABRicate's EFSA reporting thresholds. Kept on that analysis only — a blanket cutoff would fight AMRFinderPlus's curated per-gene ones. |
| eight databases | `DATABASES`, `00_common.smk` | Which ABRicate databases run: argannot, card, ecoh, ecoli_vf, megares, ncbi, resfinder, vfdb. |
| `--target_bases 500000000` | `rule filter_long_reads`, `nanopore/10_reads.smk` and `hybrid/30_ont_reads.smk` (`FILTLONG_TARGET_BASES`, `00_common.smk`) | The nanopore mode's total-output cap. Hybrid mode also applies the cap when `short_read_guidance: false`. With guidance enabled it uses the selected Illumina reads, `--trim` and `--split 1000` instead; `--length_weight` is passed in hybrid mode. |
| coverage ≥ 2×, length ≥ 500 bp | `FASTA_SEL_CMD`, `00_common.smk` | The contig filter, read out of the SPAdes header. Always applied in `illumina` and `hybrid`; in `contigs` mode only when the first header looks like SPAdes' own (it carries `length_…` or `cov_…`), since inventing a filter for an assembly from anywhere else would silently delete real sequence. |
| 80% of the reference | `EXACT_NAME_REFERENCE_COVERAGE`, `name_ice_elements.py` | Below it, an ICEberg name is suffixed `-like`. Distinct from `iceberg.min_overlap_fraction`, which is measured against our candidate. |

!!! warning "The CARD analysis is a sensitive screen, not a high-identity one"

    Until 2026-08-15 this rule passed BBMap `idfilter=0.99`, and both the code and
    the documentation claimed results were screened at 99% identity. They were
    not. `idfilter` does not filter the primary alignment of a properly-paired
    read (BBMap 39.33, `align2/AbstractMapThread.java`: the guard reads
    `if(!r.paired() && …)`), so every result the analysis ever produced was screened at
    BBMap's default of 0.76. Mapping one set of reads both ways shows it. Run as
    pairs, nearly every alignment that `idfilter=0.99` kept was below 99% identity,
    some far below. Run as single-end, almost none survived.

    The rule now passes `minid=0.76` — which changes no output, and stops the code
    claiming a stringency it never applied. It is deliberately not raised: at a
    real 0.99 the same reads keep almost none of their alignments, which would leave
    this read-based screen with almost nothing to report. Specificity comes from requiring
    ≥70% of the reference gene's **length** to be covered, which is why the loose
    identity floor is used for screening. This gate does not establish a calibrated
    specificity or demonstrate resistance. There is no second pass and no divergent tier.
    See [antimicrobial resistance](../analysis/amr.md).

---

## Where to go next

| Page | Read it for |
|---|---|
| [Reference databases](../getting-started/databases.md) | the databases these keys point at, and how to get them |
| [Quick start](../getting-started/quick-start.md) | a fresh clone to a running job |
| [Running BacFlux](running.md) | what `--cores` does, and how to restart a run |
| [Output files](output.md) | what each numbered stage writes |
| [Enabling the module](../mobilome/enabling.md) · [Optional layers](../mobilome/optional-layers.md) · [Tuning](../mobilome/tuning.md) | the mobilome block in full, layer by layer |
| [Coming from v1](../getting-started/from-v1.md) | which keys changed name or moved |
| [Troubleshooting](../troubleshooting.md) | what to do when a run stops |
