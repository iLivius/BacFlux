# Reference databases

BacFlux installs its own software but not its own reference data. Three groups:

| Group | Count | What you do |
|---|---|---|
| **Required** — too large or too tied to one release to fetch | 5 | download once, put the path in `directories:` |
| **Fetched for you**, into `output_dir`, on the first run that needs it | 6 | nothing, unless a link needs changing |
| **Fetched, but skippable** — the same six, minus PhiX | 5 (+1) | point `directories.*_db` at a copy you already hold |

Nothing is written back into a database you supply. BacFlux reads your copy through a
local symlink view and never modifies, moves or deletes it.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **GTDB** | Genome Taxonomy Database, the reference used for taxonomic placement |
    | **HPC** | high-performance computing — a shared compute cluster |

## The five you download

| Config key | Database | Version | Size on disk |
|---|---|---|---|
| `directories.bakta_db` | Bakta | **v6.0** | 3.9 GB light, 84 GB full |
| `directories.blast_db` | NCBI `core_nt` (or `nt_prok`) with the taxonomy files | — | ~300 GB |
| `directories.eggnog_db` | eggNOG diamond database | **v5.0.2** | ~50 GB |
| `directories.gtdbtk_db` | GTDB | **R232** | 94 GB extracted |
| `directories.platon_db` | Platon | **v1.5.0** | 2.8 GB |

Two of those versions are requirements. Bakta 1.12.1 refuses a v5.x database, and
GTDB-Tk 2.7.2 pins itself to one GTDB release and rejects R226 or older. There is no
in-place upgrade for either.

Two more are chosen for you rather than by you. The eggNOG release is fixed by the
pinned eggnog-mapper 2.1.15, which builds its own download URL from `__DB_VERSION__`
in the installed release, so the database cannot drift away from the tool. And the
Platon database is versioned separately from Platon itself: v1.5.0 is current for
every Platon from 1.5.0 onwards, including the 1.8 pinned here, so the database
number trailing the tool number is expected.

### Bakta

```bash
# full database (recommended)
wget -c https://zenodo.org/records/14916843/files/db.tar.xz
tar -xJf db.tar.xz && rm db.tar.xz

# or the light version
wget -c https://zenodo.org/records/14916843/files/db-light.tar.xz
tar -xJf db-light.tar.xz && rm db-light.tar.xz
```

Note the archives are `.tar.xz` (hence `-xJf`), not `.tar.gz` as in earlier releases.

!!! tip "If AMRFinderPlus errors on its database, refresh it in place"

    The Bakta archive ships an AMRFinderPlus database alongside the annotation data,
    and the two can fall out of step with the AMRFinderPlus binary in the Bakta conda
    environment. The symptom is a database error from AMRFinderPlus, not a wrong
    result — most often when the mobilome module runs, since that is what reads it.

    Activate the Bakta conda environment and point the update at your Bakta database
    directory:

    ```bash
    amrfinder_update --force_update --database db/amrfinderplus-db/
    ```

    **BacFlux never does this for you.** The workflow only ever reads
    `{bakta_db}/amrfinderplus-db/latest`; it does not check the database's age and
    does not refresh it, so a stale copy stays stale until you run the command above.

!!! note "This is also why the mobilome module needs no database of its own"

    `{bakta_db}/amrfinderplus-db/` is where the module reads AMRFinderPlus from. A
    working Bakta database is the whole prerequisite.

### NCBI core nt

```bash
# NCBI's own downloader. It ships with BLAST+, so it is already in the blast
# environment; --decompress unpacks each volume as it arrives, and every transfer
# is checked against the .md5 NCBI publishes beside it.
update_blastdb.pl --decompress --num_threads 4 core_nt      # or nt_prok

# taxonomy: BlobTools reads nodes.dmp and names.dmp from this same directory
wget -c 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' && tar -zxf taxdump.tar.gz
wget -c 'https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'       && tar -zxf taxdb.tar.gz
wget -c 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
gunzip nucl_gb.accession2taxid.gz
```

!!! warning "rsync no longer works, and ftp:// is a liability"

    NCBI has switched its rsync service off — `rsync://ftp.ncbi.nlm.nih.gov` now
    answers `Connection refused`, so any recipe built on it fails outright.

    `ftp://` still resolves, but it is worth avoiding anyway: sustained FTP transfers
    from an institutional network are a common trigger for rate limiting at NCBI's
    end, and the protocol is blocked outright on many campus networks. Every NCBI URL
    here is `https://`, which is served from the same host and the same paths.

    `update_blastdb.pl` is the route NCBI supports. It resolves the volume list
    itself, so nothing needs updating here when the number of `core_nt` volumes
    changes, and it verifies each volume's md5.

`nodes.dmp` and `names.dmp` **must** sit in the same directory as the BLAST volumes.
`parameters.nt_version` chooses which subfolder is searched: `core_nt` (default) or
`nt_prok`.

### eggNOG

```bash
conda create -n eggnog-mapper eggnog-mapper=2.1.15
conda activate eggnog-mapper
mkdir /data/eggnog_db
download_eggnog_data.py --data_dir /data/eggnog_db -y
```

### GTDB

```bash
wget -c https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
tar xzf gtdbtk_r232_data.tar.gz && rm gtdbtk_r232_data.tar.gz
```

The archive is ~61 GB. Point `gtdbtk_db` at the `release232` directory itself.
GTDB-Tk 2.7.x reads a pre-sketched skani database that GTDB now ships inside the
release, so unlike 2.6.x it neither builds nor needs space for its own ~57 GB sketch
cache.

If the single 61 GB archive is awkward on your connection, GTDB also publishes it
split into 12 parts, which can be pulled in parallel and joined:

```bash
# split package (alternative) — r232 has 12 parts, aa..al
base_url="https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/split_package/gtdbtk_r232_data.tar.gz.part_"
suffixes=(aa ab ac ad ae af ag ah ai aj ak al)
printf "%s\n" "${suffixes[@]}" | xargs -n 1 -P 12 -I {} wget -c "${base_url}{}"

cat gtdbtk_r232_data.tar.gz.part_* > gtdbtk_r232_data.tar.gz
tar xzf gtdbtk_r232_data.tar.gz
rm gtdbtk_r232_data.tar.gz gtdbtk_r232_data.tar.gz.part_*
```

The part count is release-specific — R226 had 14 (`aa`..`an`). `cat` joins whatever
matches the glob without complaint, so a part that failed to download produces a
corrupt archive rather than an error: check you have all 12 before joining.

### Platon

```bash
wget https://zenodo.org/record/4066768/files/db.tar.gz
tar -xzf db.tar.gz && rm db.tar.gz
```

## What the workflow fetches for itself

These land under `output_dir` and are shared by every sample in the run. Reuse the
same output directory across runs of a project and they are downloaded once.

| Database | Where it lands | Link key | Read in |
|---|---|---|---|
| PhiX | `01.reads/phix/` (temporary) | `links.phix_link` | `illumina`, `hybrid` — **required** |
| CARD | `05.amr/card_db/` (temporary) | `links.card_link` | `illumina`, `hybrid` — **required** |
| dbCAN | `04.annotation/dbcan/<version>/` | `links.dbcan_link` | all |
| antiSMASH | `04.annotation/antismash/databases/` | — (its own downloader) | all |
| CheckV | `07.phages/checkv_db/` | `links.checkv_link` | all |
| VirSorter2 | `07.phages/vs2_db/` | — (`virsorter setup`) | unless `phage.caller: genomad` |
| geNomad | `07.phages/genomad_db/` | `links.genomad_link` + `genomad_md5` | only when `phage.caller: genomad` |

`phix_link` and `card_link` ship empty and are required in the two modes that have
short reads. An unset one stops the run at parse time, naming the key.

Two defaults point at Zenodo mirrors rather than the publisher's own host: CheckV's
official database lives on `portal.nersc.gov`, as does geNomad's downloader, and
[NERSC](https://www.nersc.gov/) is unavailable often enough to cost real time. Both
mirrors are unmodified copies. Change `links.genomad_md5` whenever you change
`links.genomad_link` — a mismatch stops the run.

!!! note "The dbCAN pin"

    `dbcan` is pinned to 5.1.2 with a version-pinned Zenodo copy of its database,
    because run_dbcan's official mirror only ever hosts the current one, and moving to
    5.2.x is not a drop-in swap. See [Annotation](../analysis/annotation.md) for what
    the two versions call differently.

## Pointing at a copy you already hold

Set any of these and the corresponding download is skipped entirely.

| Key | Point it at | Worth setting because |
|---|---|---|
| `directories.checkv_db` | the **parent** directory holding the versioned folder, e.g. `/db/checkv/` containing `checkv-db-v1.5/` | the official host is often down |
| `directories.vs2_db` | what `virsorter setup` produced (~10 GB) | only read when `phage.caller: virsorter2` |
| `directories.antismash_db` | what `download-antismash-databases` produced | large |
| `directories.dbcan_db` | a dbCAN database matching the version in `links.dbcan_link` | avoids a repeat download per project |
| `directories.card_db` | an extracted CARD database | CARD is temporary and re-fetched every run otherwise |
| `directories.genomad_db` | the **inner** `genomad_db/` directory — the one holding `version.txt` and `genomad_db.dbtype` | its downloader hard-codes one host with no mirror option, so a local copy is the only fallback |

Each of these gets a harder check than the five required paths: BacFlux confirms the
directory exists **and** holds a file only that database has, so a path aimed at the
wrong copy is caught at parse time rather than by the tool.

!!! warning "`version.txt` must be readable, not merely present"

    In a shared database directory it easily ends up mode 0640 while the big files
    beside it are world-readable. The check reads the file.

## The mobilome layers

The optional mobilome module fetches up to five more data sets — the CONJscan models
whenever the module runs, and ICEscan, TnCentral, ICEberg and ISOSDB only when you
configure a source for each. Each has its own `dir:` key for a copy you already hold,
and each writes a `PROVENANCE.txt` recording what actually arrived. See
[Optional layers](../mobilome/optional-layers.md).

## Practical notes

**Downloads use the launcher environment.** The rules that fetch and unpack a database
declare no conda environment of their own and use `wget`, `tar` and `sha256sum` from
whatever you launched Snakemake in. `wget` is missing from some minimal conda bases
and HPC login shells — check before a long run ([Installation](installation.md)).

**A wrong required path is not always caught at parse time.** The five
`directories.*_db` values are read but not opened while the plan is built. Some are
then caught by the dry run because they are named as rule inputs — a `blast_db` with no
`nodes.dmp` beside it fails there — and the rest surface when the tool runs. A dry run
that reaches the job table means the plan is sound, not that every database is
complete.

**Disk.** A full set is roughly 450 GB, dominated by `core_nt` and GTDB. The light
Bakta database saves 80 GB at some cost in annotation depth. The databases BacFlux
fetches into `output_dir` add about 20 GB.

## Next

- [Quick start](quick-start.md) — the config keys these paths go in.
- [Configuration](../reference/configuration.md) — every key, with its default.
- [Troubleshooting](../troubleshooting.md) — what a missing or half-extracted database
  looks like when it fails.
