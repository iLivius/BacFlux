# Troubleshooting

Each entry below gives the message you get, why it happens, and what to change.

Run `-n` first, every time. BacFlux does most of its checking while the workflow
is being *read*, so a dry run costs seconds and catches the config problems below
before a single job starts — see
[Always dry-run first](reference/running.md#always-dry-run-first).

## `KeyError: 'mode'`

!!! warning "Symptom"

    The run stops as soon as the config is read, before the line naming the mode:

    ```text
    KeyError in file ".../workflow/rules/shared/00_common.smk", line 81:
    'mode'
    ```

**Cause.** The config has no `mode:` key, which is what selects the pipeline. It is
read with bracket access on purpose (`00_common.smk:81`), so a config without it
stops the run rather than letting BacFlux guess a front end. A config written for
an earlier release will not have it, and `config/config_custom.yaml` is git-ignored
— a `git pull` never updates yours.

**Fix.** Start from the shipped example and copy your database paths across by
hand. The block layout changed as well, so read it rather than paste into it:

```bash
cp config/config.yaml config/config_custom.yaml
```

Details of what moved are in [Coming from v1](getting-started/from-v1.md).

## A database path is wrong

Where a bad path is caught depends on which key it is, and the difference matters:
some cost you seconds, one costs you the length of the run up to that rule.

| Config key | Caught | What you see |
|---|---|---|
| `input.illumina_dir`, `nanopore_dir`, `contigs_dir` | parse time | `[BacFlux] input.contigs_dir is not an existing directory: '/no/such/dir'` |
| `directories.vs2_db`, `antismash_db`, `dbcan_db`, `card_db`, `genomad_db` | parse time | `[BacFlux] directories.vs2_db is '/tmp', but it has no 'Done_all_setup'. Point it at the directory 'virsorter setup' produced …` |
| `directories.checkv_db` | parse time | `[BacFlux] directories.checkv_db is '…', but that directory holds 0 CheckV database(s) (looking for */genome_db/checkv_reps.faa).` |
| `links.dbcan_link`, plus `phix_link` and `card_link` in the modes that need them | parse time | the missing key is named |
| `directories.blast_db` | DAG build (so, a dry run) | `MissingInputException in rule blob_json … affected files: <blast_db>/names.dmp, <blast_db>/nodes.dmp` |
| `directories.bakta_db`, `eggnog_db`, `gtdbtk_db`, `platon_db` | **not checked** | the tool fails when its own rule runs |

Five of the six optional database keys share one check: a file that can only exist
inside a real copy of that particular database — `Done_all_setup` for VirSorter2,
`aro_index.tsv` for CARD, `version.txt` for geNomad, `dbCAN.hmm` for dbCAN,
`clusterblast/` for antiSMASH (`00_common.smk:231-245`). A directory that exists
but holds the wrong thing is therefore still caught. `checkv_db` is the exception,
because it is the one key that points at the *parent* of the versioned database
folder; it is checked for exactly one `*/genome_db/checkv_reps.faa` underneath, so
aiming a level too high or too low is caught as well (`00_common.smk:1427-1446`).

!!! warning "Four of the five required databases are read, not verified"

    `bakta_db`, `eggnog_db`, `gtdbtk_db` and `platon_db` are taken from the config
    as given (`00_common.smk:165-169`). Set all four to paths that do not exist and
    a dry run still builds the whole job list without complaint; the failure comes
    later, inside Bakta, eggNOG-mapper, GTDB-Tk or Platon, which for the first three
    is an hour or more into the run. `blast_db` escapes this only because BlobTools
    reads `names.dmp` and `nodes.dmp` out of it, and those are declared rule inputs.

    Check all five paths by hand before a long batch. What each one should contain
    is on [Reference databases](getting-started/databases.md).

## `wget: command not found`

!!! warning "Symptom"

    A rule that downloads something fails with a bare `command not found`.

**Cause.** Twenty-three of BacFlux's rules declare no conda environment. They run
in whatever environment you launched Snakemake from, on whatever is on its `PATH`.
Four of those rules download with `wget`, which is not part of coreutils and is
genuinely absent from some minimal conda base environments and HPC login shells.

The four sit at different depths, so this can bite at any point: `download_phix`
runs in the first seconds of an `illumina` or `hybrid` run, `download_amr_db`
(CARD) and `cazyme_db_download` (dbCAN) shortly before their own analysis stages,
and `icescan_models` only when the mobilome module is on.

**Fix.** Install it into the launcher environment and relaunch; everything already
finished stays finished.

```bash
conda activate snakemake
command -v wget tar sha256sum awk python
conda install -n snakemake -c conda-forge wget
```

The full list of what the launch environment must provide, and which rules need
it, is in
[Installation](getting-started/installation.md#4-the-tools-bacflux-expects-to-find-and-does-not-install).

## The run stopped, but the exit status was 0

!!! warning "Symptom"

    Snakemake prints a message about your input files and nothing runs — but a
    wrapper script treats the run as successful.

    ```text
    More than one contigs file extension detected:
        fa
        fasta
    ```

**Cause.** Four sample-discovery checks stop the run with `sys.exit(0)` rather than
a non-zero status: an unusable set of file extensions
(`00_common.smk:1358-1369`), a sample name containing one of
`* # @ % ^ / ! ? & : ; | < >` or a space (`00_common.smk:1676`), an Illumina R1
with no R2 beside it (`00_common.smk:1690`), and, in hybrid mode, a sample with
reads of only one kind (`00_common.smk:1734`).

**Fix.** Every sample in one batch must share a single extension, so split mixed
batches into separate runs or rename the odd ones out. If you wrap BacFlux in a
script, do not rely on the exit status alone — check that the expected output
files exist.

## Out of memory

BacFlux's two heaviest steps are the two that place genomes on a reference tree
with `pplacer`: CheckM (`completeness_and_contamination`) and GTDB-Tk
(`taxonomic_assignment`). pplacer's memory scales with its thread count and
reaches tens of gigabytes.

Only GTDB-Tk hands pplacer BacFlux's thread count, through
`--cpus {threads} --pplacer_cpus {threads}` (`30_taxonomy.smk:121-122`). CheckM
gets `-t {threads}` and no `--pplacer_threads` at all, so its placement runs at
CheckM's own default of one thread, which is much the cheaper of the two.

`{threads}` in both rules is `capped_cpus(24)`, which is
`min(resources.threads, 24)` read from your **config**, and Snakemake lowers it
again at run time if `--cores` is smaller (`00_common.smk:421-429`). So if a run
is still being killed there the lever is `resources.threads`, with `--cores` as a
ceiling on top of it — see
[How `threads` reaches a rule](reference/configuration.md#how-threads-reaches-a-rule).
Thread count changes speed and memory only, never the classification.

The two rules are also kept apart on purpose: `taxonomic_assignment` declares
CheckM's stats file as an input purely as a scheduling edge, so they never run on
the same sample at once (`30_taxonomy.smk:91-100`).

!!! warning "Never add `--jobs` alongside `--cores`"

    For a local run `--jobs`/`-j` is an alias for `--cores`, not a second
    "N jobs of M cores each" setting. Passing both switches the core budget off:
    two rules each declaring 8 threads, launched with `--jobs 2 --cores 8`, both
    got their full 8 threads and ran together. If your machine is far more loaded
    than `--cores` should allow, this is why. See
    [`--cores` is the CPU budget](reference/running.md#-cores-is-the-cpu-budget).

eggNOG's `--dbmem` is the other memory lever, and it is off by default. Turning it
on costs about 42 GB per concurrent job, and setting it with `resources.ram_gb`
below that is refused at parse time rather than allowed to fail on the first job
(`00_common.smk:1281-1291`). See
[eggNOG `--dbmem`](reference/running.md#eggnog-dbmem).

## Databases on shared or network storage

Pointing BacFlux at a copy of a database that somebody else installed is supported
and saves a great deal of downloading — every `directories.*_db` path is only ever
read. Three things go wrong with it, and all three have happened here.

**The annotation crawls.** eggNOG-mapper's annotation phase does random-access
lookups into the 39 GB `eggnog.db`, once per seed ortholog. That phase is the slow
tail of a run on any disk — it is why `functional_annotation` is routinely the last
rule still finishing (`00_common.smk:1262-1266`) — and random access is the access
pattern a network filesystem handles worst, so a shared copy is where it hurts
most. Either set `parameters.eggnog.dbmem: true` and pay the RAM, or hold that one
database on local disk. The trade-off is set out in
[Annotation](analysis/annotation.md#-dbmem-the-ram-trade-off).

**CheckV's DIAMOND index belongs to whoever built it.** CheckV's official archive
ships no DIAMOND index; it is built locally after unpacking, and DIAMOND's format
is versioned. Hand a shared copy straight to CheckV and it fails at
`[3/8] Running DIAMOND blastp search... DIAMOND task failed` — *after* the
contamination stage succeeded, which reads as a CheckV bug rather than an index
mismatch. `directories.checkv_db` sidesteps this: BacFlux builds a local view —
symlinks to the big read-only files, plus an index built by its own DIAMOND
(`00_common.smk:186-203`). It costs ~950 MB once per output directory, and nothing
is written to your copy.

**Some files are readable by their owner and not by you.** Seen on a shared geNomad
database: 11 of 27 files were mode 0640, `version.txt` among them, while the big
data files beside them were world-readable. A file you may stat but not read passes
an existence check and kills the run minutes later, so BacFlux tests every file in
that directory for real read access at parse time, then names the first six
offenders and counts the rest (`00_common.smk:300-314`):

```text
[BacFlux] directories.genomad_db is '/shared/db/genomad_db', but 11 file(s) in it are
not readable by you: <first six names> (and 5 more). geNomad needs all of them and
would fail partway through the run. Ask whoever owns that directory to make it
readable (chmod -R a+r), or point at a copy you own.
```

## `medaka_model: auto` cannot infer a model

!!! warning "Symptom"

    A long-read run stops at `check_medaka_model`, seconds after read filtering:

    ```text
    ValueError: Input file did not contain precisely 1 basecaller model reference.
    ```

**Cause.** `auto` runs `medaka tools resolve_model --auto_model consensus_bacteria`
over the filtered reads, which reads the basecaller model out of the FASTQ headers.
Guppy and Dorado stamp it there; older data, re-headered data and much public data
do not carry it at all. Confirmed on a real ONT test set: no basecaller tag in
either the raw FASTQ or the filtlong output, so `auto` cannot work for that dataset
at all.

**Fix.** Name the model explicitly, under the block named after your mode
(`parameters.hybrid` for a hybrid run):

```yaml
parameters:
  nanopore:
    medaka_model: r1041_e82_400bps_sup_v4.2.0
```

!!! warning "An explicit model can change the assembler as well"

    With `flye_input_mode: auto` (the default), the Flye read mode is chosen from
    the Medaka model name: a name containing **`fast`** selects `--nano-raw`,
    anything else `--nano-hq`. So setting `medaka_model` can change how the genome
    is assembled, not just how it is polished. To pin the polisher without touching
    the assembler, set `flye_input_mode` explicitly to `nano-hq` or `nano-raw`.

The failure is deliberately early: `check_medaka_model` runs right after read
filtering and the assembler waits on its output, so a bad model costs seconds
rather than surfacing after Flye has run. Both keys, and what `auto` accepts, are
in
[Configuration](reference/configuration.md#long-read-assembly-and-polishing-nanoporehybrid).

## Medaka reuses an index from a previous run

!!! warning "Symptom"

    You changed a decontamination setting, re-ran a `nanopore` sample, and the
    Medaka log says:

    ```text
    Using the existing fai index file …/assembly_decontam.fasta.fai
    Using the existing mmi index file …/assembly_decontam.fasta.map-ont.mmi
    ```

**Cause.** `medaka_consensus` aligns the reads with `mini_align`, which writes its
`samtools` and `minimap2` indexes **next to the draft** rather than into Medaka's
output directory, and reuses them whenever they already exist (`mini_align`
lines 134-158, in the `medaka=2.2.2` environment). BacFlux does not pass its `-f`
flag, which is what would force a rebuild.

In `nanopore` mode the draft is
`02.assembly/{sample}/contaminants/assembly_decontam.fasta`, a plain file output of
rule `select_contigs` rather than a `directory()`. So re-running the selector
rewrites the FASTA and leaves the two index files beside it, still describing the
contig set you had before, and Medaka aligns against that. `hybrid` mode does not
have the problem: its draft sits inside `fix_start/`, a `directory()` output that
Snakemake wipes before rebuilding, so a change that re-runs `fix_start` takes the
stale index with it.

**Fix.** Delete the two files and re-run. They are rebuilt in seconds.

```bash
rm -f output_dir/02.assembly/*/contaminants/assembly_decontam.fasta.fai \
      output_dir/02.assembly/*/contaminants/assembly_decontam.fasta.map-ont.mmi
```

If the contig set did not change, reusing them is harmless and saves the indexing
step — the message only matters when the draft beside them has been rewritten.

## GTDB-Tk will not accept the reference database

**Cause.** GTDB-Tk hard-pins itself to one compatible reference release in its own
source (`COMPATIBLE_REF_DATA_VERSIONS`). BacFlux pins **GTDB-Tk 2.7.2**, which
requires **GTDB R232**; the 2.6.1 that earlier releases used accepted only R220 and
R226. GTDB-Tk makes the check itself, on startup, so the refusal lands in
`logs/taxonomic_assignment_{sample}.log` rather than at parse time.

**Fix.** Download R232 fresh and point `directories.gtdbtk_db` at the `release232`
directory. There is no in-place upgrade, and taxonomy output cached from an older
release should be treated as stale once that path changes. The download recipe and
sizes are on [Reference databases](getting-started/databases.md); why the two are
pinned together is in
[Taxonomy](analysis/taxonomy.md#version-and-database-are-pinned-together).

## Half the genome left as contamination

!!! warning "Symptom"

    CheckM reports a badly incomplete genome that is **not** contaminated, and
    `logs/select_contigs_{sample}.log` carries the warning the selector prints
    whenever it drops more than 20% of the assembly:

    ```text
    WARNING: decontamination is discarding 52% of the assembly (2,340,112 of
    4,501,776 bp). If the kept assembly then looks incomplete but NOT contaminated,
    the filter has most likely removed genome rather than contamination.
    ```

**Cause.** `auto` mode keeps the genus carried by the most contigs, and BLAST
routinely scatters one organism across several related genus names. When the vote
is close, the wrong half wins. One isolate lost its entire 4.5 Mb chromosome to a
2-2 tie between *Bacillus* and *Peribacillus*, broken alphabetically. The curated
alias table exists so that tie never has to be broken at all
([the genus aliases](analysis/decontamination.md#the-genus-aliases)), but it covers
only the splits it knows about.

**Fix.** Read `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` —
every contig carries the reason it went — and
`{sample}_composition.txt` beside it, which prints each genus's share of the DNA
next to its share of the contig count. When those two columns disagree, the vote
went one way and the sequence the other. Then either set
`decontamination.mode: include` naming both genera, or `mode: off`, and re-run;
Snakemake redoes only the affected steps.

!!! note "If nothing survives at all, the rule fails"

    ```text
    ValueError: No contigs were kept for sample 'X' with mode 'auto'.
    Check taxonomy assignments and decontamination settings.
    ```

    That is deliberate. An empty FASTA would let Snakemake mark the step complete
    and every stage from annotation onwards would run on no sequence. The decisions
    file is written before the check, so it is there to read.

## A plasmid is missing

Two independent steps can delete a real plasmid, and both have to be ruled out.

| Step | Why | Where to look |
|---|---|---|
| Decontamination | `auto` follows BLAST bestsum, which follows database composition rather than biology. A 2,014 bp *K. pneumoniae* plasmid was dropped as *Escherichia* because *E. coli* entries outsummed *Klebsiella* 58,709 to 16,253 | `contaminants/contig_taxonomy_decisions.tsv` |
| Long-read filtering | filtlong ranks reads by quality **and** length, and a plasmid cannot produce reads longer than itself, so its reads sit at the bottom by construction. Of the 603 ONT reads covering a 5,596 bp Col2 plasmid, the hard-coded settings used before these became config keys kept 93; the shipped `keep_percent: 95` and `length_weight: 1` keep 502 | `logs/filter_long_reads_{sample}.log` |

!!! warning "Recovering the reads is necessary, not sufficient"

    On that isolate **no** filtlong setting assembled the plasmid, raw unfiltered
    reads included. ONT ligation prep under-represents small circular plasmids by
    roughly four-fold, and no filter setting undoes that. `keep_percent` is the key
    to reach for; `length_weight` only reaches filtlong in `hybrid` mode, because
    the `nanopore` call caps its output with `--target_bases` instead.

The plasmid is missing from the **assembly**, not merely from the plasmid call, if
it appears in neither. Both routes and their fixes are set out in
[Decontamination](analysis/decontamination.md#5-two-ways-this-step-deletes-a-real-plasmid)
and [Plasmids](analysis/plasmids.md#before-you-conclude-a-plasmid-is-absent).

## The mobilome module cannot fetch its models

!!! warning "Symptom"

    The first run with `mobilome.run: true` fails at rule `conjscan_models`:

    ```text
    You reach the maximum number of request per hour to github.
    Please wait before to try again.
    ```

**Cause.** The CONJscan models are not shipped with BacFlux; they are fetched at
your request by `msf_data install` (`macsydata install` on older MacSyFinder),
which reads `https://api.github.com` over an **unauthenticated** connection and
turns the HTTP 403 you get for exceeding the hourly allowance into that message
(`macsylib/model_package.py:165`, `181-184`). GitHub counts the allowance per IP
address, so a shared machine or an institutional gateway can use it up without you
making a single request yourself.

**Fix.** Wait, then relaunch. The rule fails rather than leaving an empty models
directory behind, so nothing downstream runs on a model set that was never
installed. The models land in `08.mobilome/conjscan_models/`, are shared by every
sample, and are never re-fetched for that output directory once they are there.

This is the only download the module always makes; the naming layers are opt-in,
and they come from ordinary web servers rather than GitHub. Their rules check what
came back, because a redirect to an HTML error page arrives with a 200 and would
otherwise be indexed as an empty database that silently names nothing — so each one
fails loudly and prints the first bytes it received. Turning those layers on is
covered in
[Optional layers](mobilome/optional-layers.md#what-the-download-rules-do-besides-downloading).

## Config (YAML) errors

Snakemake reports any config problem with one generic message:

```text
WorkflowError:
Config file is not valid JSON or YAML. In case of YAML, make sure to not mix
whitespace and tab indentation.
```

That almost always means an indentation or list-syntax mistake, or a stray tab —
not a problem with the values. Validate the file directly to get the exact line:

```bash
python -c "import yaml; yaml.safe_load(open('config/config_custom.yaml'))"
```

The key-by-key reference is in [Configuration](reference/configuration.md).
