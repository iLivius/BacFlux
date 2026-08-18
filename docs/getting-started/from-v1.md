# Coming from v1

For anyone who ran BacFlux 1.x, FastaFlux, BacFluxL or BacFluxL+. Everything on this
page is about what changed; nothing here is needed for a fresh install.

Version 2.0.0 brings the four workflows into one repository, one mode each. The
analysis steps, the thresholds and the decontamination policy came across as they
were — given the same reads, tool versions, read mode and thread count, each read-based
mode reproduced its predecessor's assembly byte for byte
([Rationale](../about/rationale.md)). What changed is how a run is configured and where
its output lands.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies — the long-read sequencing platform |
    | **DOI** | digital object identifier — the permanent handle a Zenodo record is cited by |
    | **GTDB** | Genome Taxonomy Database, the reference used for taxonomic placement |

| Was | Now |
|---|---|
| `BacFlux` — Illumina reads | `mode: illumina` |
| `FastaFlux` — pre-assembled contigs | `mode: contigs` |
| `BacFluxL` — ONT reads | `mode: nanopore` |
| `BacFluxL+` — Illumina + ONT | `mode: hybrid` |

The `BacFluxL` and `BacFluxL+` repositories are retired in favour of this one. **Their
Zenodo DOIs stay valid**, so an analysis already published with either remains citable
and reproducible — cite the record you actually ran.

## 1. Start your config again

**A v1 config will not run v2.** It has no `mode:` key, and v2 stops at parse time
without one:

```text
KeyError in file ".../workflow/rules/shared/00_common.smk", line 81: 'mode'
```

That is the intended behaviour rather than a bug. `config/config_custom.yaml` is
git-ignored, so it is your file and a `git pull` never touches it — which is also why
an old one can still be sitting there.

```bash
cp config/config.yaml config/config_custom.yaml
```

Then copy your database paths across by hand. The file is organised in sections —
`mode`, `input`, `directories`, `links`, `resources`, `parameters`, `phage`,
`mobilome` — and every key carries an inline comment saying which modes read it, so it
is worth reading rather than pasting into. Every key is also listed in
[Configuration](../reference/configuration.md).

If you cloned fresh, none of this applies: there is no old file to collide with.

## 2. The mode lives in the config, not on the command line

One launch command runs all four modes:

```bash
snakemake --sdm conda --configfile config/config_custom.yaml --cores 24
```

Two things to know about it:

- **`--cores` alone is the CPU ceiling. Do not add `--jobs`/`-j`.** For a local run
  `--jobs` is an alias for `--cores`, not an independent "N jobs of M cores each"
  setting, and passing both lets rules oversubscribe the machine. Every CPU-bound rule
  now declares Snakemake's built-in `threads:`, which is what `--cores` enforces. The
  measurement is on [Running BacFlux](../reference/running.md).
- **`--snakefile` is never needed.** `workflow/Snakefile` is found automatically from
  the repository root, and it dispatches on `mode`.

## 3. Relative paths resolve from where you launch

BacFlux no longer uses Snakemake's `workdir:` directive, so the process never changes
directory. A relative input path now resolves against the directory you launched from,
which is what most people expect. `directories.output_dir` may still be relative; it is
resolved to an absolute path for you.

One consequence is easy to miss: the conda environments have moved. Snakemake
creates them under `.snakemake/conda` in the launch directory. Launching from the same
place each time is what stops the 31 environments being rebuilt; `--conda-prefix
/some/path` puts them somewhere specific, for example shared between projects.

## 4. The output directory is renumbered

All technology-specific work is now grouped under two fixed parents, so every shared
stage has the same number in every mode. Previously each workflow had a different number
of front-end stages, so the shared stages were numbered differently in each.

| Was (short-read v1) | Now |
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

The delivered genome is now always `02.assembly/{sample}/contigs_final.fasta`,
whichever mode produced it, and every stage from `03` onwards reads that one file.

**There is no in-place upgrade.** Point `output_dir` at a fresh directory rather than
trying to reuse an old one. `workflow/scripts/clean_workdir.sh` recognises the four
older layouts as well as the current one, so an old output directory can still be
cleaned down to results and archived.

## 5. Everything else worth knowing before the first v2 run

| Change | What to do |
|---|---|
| **Sample names may contain underscores** (and dots). Fifteen characters are still refused: `* # @ % ^ / ! ? & : ; \| < >` and space | nothing, unless you had renamed files to avoid underscores |
| **GTDB R232 with GTDB-Tk 2.7.2** | download R232; the pinned version rejects older releases |
| **Bakta database v6.0** with Bakta 1.12.1 | download v6.0 |
| **The virus caller is selectable.** VirSorter2 stays the default; geNomad is opt-in, and choosing it also produces a Platon + geNomad plasmid concordance table | nothing, unless you want geNomad |
| **Six databases can be supplied from a copy you already hold** — `checkv_db`, `vs2_db`, `antismash_db`, `dbcan_db`, `card_db`, `genomad_db` | set the ones you have; the download is then skipped |
| **The mobilome module is new and off by default** | see [the mobilome module](../mobilome/index.md) |
| **QUAST now runs in `nanopore` mode**, which never had it, and in `hybrid` it evaluates the delivered genome as well as the Illumina draft | nothing |

One downstream difference is deliberate and appears in the two long-read modes: Bakta
gains or loses a single feature, because it is now handed a `--replicons` table
declaring a closed chromosome circular. That makes its gene caller run in closed mode,
which shifts a couple of marginal start codons. It is 0.02% of the features, and the
better call.

## Next

- [Choosing a mode](choosing-a-mode.md) — which mode your data belongs in.
- [Quick start](quick-start.md) — the first v2 run, start to finish.
- [Output files](../reference/output.md) — what each numbered stage now holds.
- [Citation](../about/citation.md) — which DOI to cite for work already published.
