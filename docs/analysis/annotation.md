# Annotation

Stage `04.annotation` turns the finished genome into gene calls and functional labels.
Four tools run here, in all four modes, all of them in
`workflow/rules/shared/40_annotation.smk`. Every front end has already converged on one
canonical assembly (`02.assembly/{sample}/contigs_final.fasta`), so nothing in this
stage branches on the sequencing technology.

Bakta goes first and the other three read its output — not the contigs. That is
deliberate: one gene set underlies every functional layer, so a CAZyme cluster, a
biosynthetic cluster and an orthologous group all sit on the same coordinates.

| Tool | Version | Answers | Written to |
|---|---|---|---|
| [Bakta](https://github.com/oschwengers/bakta) | 1.12.0 | where the genes are, and what they are called | `04.annotation/bakta/{sample}/` |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | 2.1.15 | what each protein's orthologous group is known to do | `04.annotation/eggnog/{sample}/` |
| [antiSMASH](https://github.com/antismash/antismash) | 8.0.4 | which biosynthetic gene clusters the genome carries | `04.annotation/antismash/{sample}/` |
| [run_dbcan](https://github.com/bcb-unl/run_dbcan) | 5.1.2 | which sugars and polysaccharides it can process | `04.annotation/dbcan/{sample}/` |

Three of the four are terminal products — you read them, nothing downstream does.
Bakta is the exception: MultiQC summarises it, and the optional mobilome module reads
its proteins and its GFF3.

## Bakta

*Rule `annotation`.*

Bakta calls CDS, tRNA, rRNA, ncRNA and CRISPR arrays on the delivered genome and
writes the whole standard annotation set. The output is declared to Snakemake as a
**directory**, because Bakta names the files inside it; every downstream rule takes
its dependency on the directory and reaches in for the one file it wants.

| | |
|---|---|
| **In** | `02.assembly/{sample}/contigs_final.fasta`; the per-genus composition table from [decontamination](decontamination.md); the replicon table in the long-read modes |
| **Out** | `04.annotation/bakta/{sample}/` — `{sample}.faa`, `.gff3`, `.gbff`, `.tsv` and the rest |
| **Next** | eggNOG-mapper (`.faa`), antiSMASH (`.gbff`), dbCAN (`.faa` + `.gff3`), MultiQC, and the [mobilome module](../mobilome/index.md) (`.faa` + `.gff3`) |

The settings that are worth knowing:

| Flag | Value | What it does |
|---|---|---|
| `--translation-table` | `11` | the bacterial genetic code |
| `--min-contig-length` | `500` | contigs shorter than this are not annotated |
| `--keep-contig-headers` | — | keeps the assembler's contig names, so coordinates in every later stage refer to the same sequence IDs |
| `--locus-tag` | derived | the sample name, with anything outside letters, digits and `_ . -` replaced by `_`, then cut to 24 characters. Bakta rejects some characters in a locus tag; your sample name itself is untouched |
| `--strain` | sample name | what appears in the GenBank header |

The database is `directories.bakta_db` and must be **v6.0** — Bakta 1.12.1 refuses an
older one. See [Reference databases](../getting-started/databases.md).

### The genus hint

Bakta annotates more accurately when it knows roughly what it is looking at, so the
rule derives a genus from the composition table decontamination wrote, and passes it
as `--genus X --species sp.`. It is a hint and nothing more: a wrong genus costs
annotation accuracy, never contigs or reads. Lines for the literal genus `no-hit` are
skipped first, so a poorly-placed sample cannot end up annotated as `--genus no-hit`.
A sample with no usable genus is still annotated, just without the hint.

!!! note "Check which genus was used"

    When the table lists more than one genus, the line picked is **not** the most
    abundant one. Every line carries two labelled figures — see
    [the composition report](decontamination.md#the-composition-report) — so the numeric
    sort that reads them scores each line zero, and the tie then breaks on the whole
    line, reversed: the alphabetically last genus wins. On the four-genus sample
    described there (*Aneurinibacillus* 0.40 of the DNA, *Bacillus* 0.29,
    *Paenibacillus* 0.16, *Brevibacillus* 0.11) Bakta is handed *Paenibacillus*.

    Nothing announces it during the run, so on a mixed assembly read the first line of
    `logs/annotation_{sample}.log`: it records the genus that was used, or says that
    none was found.

### Circular replicons

*`nanopore` and `hybrid` only.*

Bakta treats every sequence as a linear contig unless told otherwise, and that costs
genes: on a sequence declared circular, Pyrodigal may call a gene running across the
origin instead of leaving two partial CDS at the ends. On a closed chromosome that is
typically a handful of genes at position 1, often including *dnaA* itself.

Long-read assemblies are the only ones where the topology is actually measured, so
`build_replicons` writes a small table — topology from Flye's circularity call, type
from dnaapler's marker hit — and the rule adds `--replicons` when that table exists
and is not empty. In `illumina` and `contigs` mode the flag is simply absent. See
[Nanopore mode](../modes/nanopore.md).

## eggNOG-mapper

*Rule `functional_annotation`.*

eggNOG-mapper assigns each protein Bakta predicted to an orthologous group and
attaches what that group is known to do — COG category, GO terms, KEGG KO and pathway,
EC number.

| | |
|---|---|
| **In** | `{sample}.faa` from the Bakta directory |
| **Out** | `04.annotation/eggnog/{sample}/` with the `{sample}.emapper.*` tables |
| **Next** | nobody — you. MultiQC does not read it |

The search runs in DIAMOND mode (`-m diamond`) against `directories.eggnog_db`. A
scratch directory is written alongside the results and deleted by Snakemake when the
rule ends; it is a sibling of the kept output, never inside it, so cleanup can never
reach into a directory you keep.

### `--dbmem`: the RAM trade-off

This is the slow tail of a BacFlux run, and usually the last rule still going. The
reason is the annotation phase: it does random-access lookups into a 39 GB SQLite
database, once per seed ortholog. `--dbmem` loads that database wholly into memory so
the lookups become in-memory ones, and releases it when eggNOG-mapper exits.

```yaml
parameters:
  eggnog:
    dbmem: false        # true = load the 39 GB eggnog.db into RAM
```

| | `dbmem: false` (default) | `dbmem: true` |
|---|---|---|
| Where `eggnog.db` is read from | disk, one random-access lookup per seed ortholog | memory |
| RAM per **concurrent** eggNOG job | eggNOG-mapper's own working set | ~42 GB |
| Worth it when | anything else | the machine has RAM to spare and this rule is what you are waiting for |

Two things happen when you switch it on, both at parse time, before any job starts:

- If `resources.ram_gb` is below 42, the run **stops** with a message naming the key
  to change. One `--dbmem` job that cannot fit would run out of memory the instant it
  started.
- BacFlux prints how many such jobs fit in your declared budget, together with the
  exact `--resources mem_gb=N` flag to add to the launch line.

!!! warning "The flag is not optional if you want the limit enforced"

    `mem_gb` is a gigabyte figure, not a core count, and Snakemake only schedules
    against a named resource when the launch line passes it. Without
    `--resources mem_gb=N`, eggNOG concurrency falls back to the `--cores` bound alone:
    as many eggNOG jobs as your cores allow start together, each claiming 42 GB. See
    [Running BacFlux](../reference/running.md).

    The number comes from `resources.ram_gb`, the budget you declared, not from
    probing free memory. A live probe cannot predict how many jobs will run at once,
    reports the host's RAM inside a container, and on a cluster reads the submit
    node rather than the compute node.

## antiSMASH

*Rules `secondary_metabolites_db` (or `secondary_metabolites_db_local`) and
`secondary_metabolites_analysis`.*

antiSMASH scans the genome for biosynthetic gene clusters — antibiotics, siderophores
and the rest of secondary metabolism — and reports each cluster with its type and its
closest known relative.

| | |
|---|---|
| **In** | `{sample}.gbff` from the Bakta directory; the antiSMASH reference database |
| **Out** | `04.annotation/antismash/{sample}/` |
| **Next** | nobody — you |

It runs with `--taxon bacteria` and `--genefinding-tool none`. That second flag is the
point: **do not re-predict genes, reuse Bakta's calls**, which travel inside the
`.gbff`. It is why this rule consumes the annotated GenBank file rather than the bare
FASTA, and it keeps cluster coordinates on the same gene set as everything else in
this stage.

The reference database is fetched once by `download-antismash-databases` into
`04.annotation/antismash/databases/`, and every sample waits on it. Setting
`directories.antismash_db` to a copy you already hold skips the download: that same
`databases/` directory is then filled with symlinks pointing into your copy.

!!! note "Why a symlink view instead of pointing antiSMASH at your path"

    The database is a `directory()` output, and Snakemake deletes a directory output
    before re-running its rule. If that path were your own shared database, any
    re-run trigger — a changed environment file, a `--forcerun`, an interrupted job —
    would wipe it for everyone using it. Both rules therefore write only inside
    BacFlux's own output tree and read yours. The same pattern is used for dbCAN,
    VirSorter2, CheckV and geNomad.

## dbCAN

*Rules `cazyme_db_download` (or `cazyme_db_local`) and `cazyme_gene_cluster`.*

dbCAN finds the carbohydrate-active enzymes, groups neighbouring ones into CAZyme Gene
Clusters, and predicts the substrate each cluster acts on.

| | |
|---|---|
| **In** | `{sample}.faa` and `{sample}.gff3` from the Bakta directory; the verified dbCAN database |
| **Out** | `04.annotation/dbcan/{sample}/` |
| **Next** | nobody — you |

Four `run_dbcan` sub-commands run in sequence, all writing into that one directory:

1. `CAZyme_annotation` — HMM, DIAMOND and dbCAN-sub over the proteins
2. `gff_process` — put those calls back on the genome via Bakta's GFF3
3. `cgc_finder` — group adjacent CAZymes into clusters
4. `substrate_prediction` — predict what each cluster acts on

Step 2 is passed `--gff_type prodigal`, which is correct for a Bakta GFF3: Bakta calls
its CDS with Pyrodigal.

### The database, and the integrity check

`links.dbcan_link` must point at a `.tar.gz`; the checksum URL and the name of the
extracted directory are both derived from it, so the folder can never disagree with
the link. The download rule fetches the archive and its `.sha256`, compares the two
hashes, and **fails rather than extracting** if they differ. Only on a match does it
extract and write the verified hash into a small marker file, `.verified.sha256`.

The analysis rule waits for that marker file rather than for the database directory
itself, and the difference matters: a directory can be left behind by an interrupted or
never-verified run, whereas the marker is written only after the checksum has passed. So
dbCAN cannot start against a database that failed its check. Setting
`directories.dbcan_db` to a copy you already hold builds a symlink view instead, for the
reason given under antiSMASH above, and writes a marker that says in words `local copy,
not independently checksummed` — the analysis rule tests only that the marker exists, so
the integrity claim holds for the download path alone.

!!! note "Why run_dbcan is pinned to 5.1.2"

    The tool and its database move together, and run_dbcan's own mirror only ever hosts
    the current database — so the shipped link is a version-pinned copy on Zenodo
    ([record 18622157](https://zenodo.org/records/18622157)) instead.

    Upgrading is a project of its own rather than a URL swap: 5.2.x changed both the
    sub-command arguments and the output files, and it changes the answer. On the same
    proteins 5.2.9 called **205** CAZyme genes where 5.1.2 called **315**. The extra
    ones are almost all weak hits found by a single method (`dbCAN_sub`), which 5.2.9
    drops and 5.1.2 keeps. Calls backed by two or three methods are the same either
    way — so treat a single-method hit here with care.

## Resources

| Rule | Threads | Memory |
|---|---:|---|
| `annotation` | capped at 24 | none declared |
| `functional_annotation` | capped at 24 | `mem_gb: 42`, only with `--dbmem` on |
| `secondary_metabolites_analysis` | capped at 24 | none declared |
| `cazyme_gene_cluster` | capped at 24 | none declared |

"Capped at 24" means the rule requests `min(resources.threads, 24)`, so lowering
`threads` lowers all four at once and nothing needs configuring per rule; `--cores N`
caps the total at run time.

The database rules ask for nothing in particular; they fetch or link, once, for the
whole run. `cazyme_db_download` is the one to know about: it runs in the environment
BacFlux was launched from rather than a conda environment of its own, so `wget`, `tar`,
`sha256sum` and `awk` must be available there — see
[Installation](../getting-started/installation.md).

## What to check afterwards

- The Bakta panel in `09.report/multiqc_report.html` — feature counts per sample, side
  by side, which is the quickest way to spot a genome that annotated badly.
- `logs/annotation_{sample}.log` — the genus hint that was used, and the replicon
  table if one was passed.
- `logs/functional_annotation_{sample}.log`, `logs/secondary_metabolites_{sample}.log`,
  `logs/cazyme_{sample}.log` — one per sample per tool.
- `logs/cazyme_db_download.log` — the two downloads. If the rule failed, either a fetch
  did not complete or the tarball's checksum did not match; nothing was extracted
  either way.
