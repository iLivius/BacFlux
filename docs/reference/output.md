# Output files

Everything a run writes goes under `directories.output_dir`. BacFlux resolves that
to an absolute path before any job starts, so a relative value in the config still
lands where you expect, and the workflow never changes its working directory.

Inside it the directories are numbered in the order the stages run, so a plain `ls`
reads like the workflow itself. **The numbering is identical in all four modes.**
Everything technology-specific is grouped under `01.reads` and `02.assembly`, and
every shared stage below that gets the same number whichever front end produced the
genome. A stage a mode cannot produce is simply absent — it is never renumbered.

!!! note "The stage numbers changed in 2.0.0"

    They differed per workflow before the four front ends were merged. The old-to-new
    map is on [Coming from v1](../getting-started/from-v1.md); this page describes only
    what a 2.0.0 run writes.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **QC** | quality control — the read, assembly and completeness checks |
    | **AMR** | antimicrobial resistance |
    | **GTDB** | Genome Taxonomy Database, the reference GTDB-Tk places each genome in |
    | **ANI** | average nucleotide identity, a genome-to-genome similarity measure |
    | **EFSA** | European Food Safety Authority, whose reporting thresholds the ABRicate leg applies |
    | **ICE** | integrative and conjugative element — sits in the chromosome and can move itself to another cell |
    | **IME** | integrative mobilisable element — the same, but needs a helper element to move |
    | **COG, GO, KEGG, EC** | the ortholog, gene-ontology, pathway and enzyme classifications eggNOG assigns to a protein |

---

## The files to read first

A full run leaves several thousand files. These are the ones that carry the result.

| File | Answers | Modes |
|---|---|---|
| `02.assembly/{sample}/contigs_final.fasta` | The delivered genome. Every later stage works from this one file. | all |
| `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` | Which contigs were kept, which were dropped, and the reason for each. | all |
| `02.assembly/{sample}/eval/checkm/{sample}_checkm_stats.tsv` | Completeness and contamination of that genome. | all |
| `03.taxonomy/{sample}/classify/gtdbtk.*.summary.tsv` | What the isolate is. | all |
| `04.annotation/bakta/{sample}/{sample}.tsv` | The gene calls, with products. | all |
| `05.amr/abricate/{sample}/AMR_summary.txt` | AMR and virulence hits on the assembly, across eight databases. | all |
| `05.amr/mapping/{sample}/{sample}_CARD_report.tsv` | AMR determinants found in the *reads*, classified by what kind of CARD entry each is. | `[illumina, hybrid]` |
| `06.plasmids/{sample}/platon/verified_plasmids.txt` | Which contigs are plasmids. | all |
| `07.phages/checkv/{sample}/quality_summary.tsv` | Predicted prophages and how complete each one is. | all |
| `08.mobilome/{sample}/{sample}_amr_mobility.tsv` | For each AMR gene: the mobile-element context and a mobility tier. | all, opt-in |
| `09.report/multiqc_report.html` | Read QC, assembly QC, completeness and taxonomy for the whole batch on one page. | all |

!!! warning "One output directory, one batch"

    A sample name is the only thing separating two runs that share an `output_dir`:
    reuse a name and the second run overwrites the first sample's results throughout.
    `09.report/multiqc_report.html` is overwritten unconditionally — the rule passes
    MultiQC `--force`, because without it MultiQC writes `multiqc_report_1.html`
    instead and the rule then fails on a missing output.

    The downloaded databases under `04`, `07` and `08` are *not* per sample and are
    reused, so one output directory per project is cheaper than one per run.

---

## Layout

Bracketed tags follow `config/config.yaml`: `[illumina, hybrid]` means only those
modes write the line. Untagged lines are written in all four.

```text
output_dir/
|-- 01.reads/                                    [illumina, nanopore, hybrid]
|   `-- {sample}/
|       |-- illumina/                            [illumina, hybrid]
|       `-- ont/                                 [nanopore, hybrid]
|-- 02.assembly/
|   `-- {sample}/
|       |-- contigs_final.fasta                  <- the delivered genome
|       |-- contigs_filt.fasta                   [illumina, hybrid, contigs]
|       |-- spades/                              [illumina, hybrid]
|       |-- flye/  fix_start/  medaka/           [nanopore, hybrid]
|       |-- polypolish/  snps/                   [hybrid]
|       |-- {sample}_replicons.tsv  + _audit     [nanopore, hybrid]
|       |-- {sample}_medaka_model.txt            [nanopore, hybrid]
|       |-- contaminants/                        BLAST + BlobTools + the audit
|       `-- eval/                                QUAST, CheckM, Qualimap
|-- 03.taxonomy/{sample}/                        GTDB-Tk
|-- 04.annotation/
|   |-- bakta/{sample}/   eggnog/{sample}/
|   |-- antismash/{sample}/   + antismash/databases/
|   `-- dbcan/{sample}/       + dbcan/<db version>/
|-- 05.amr/
|   |-- abricate/{sample}/                       eight databases + a summary
|   `-- mapping/{sample}/                        [illumina, hybrid]
|-- 06.plasmids/{sample}/platon/
|-- 07.phages/
|   |-- virsorter/{sample}/  + vs2_db/           unless phage.caller: genomad
|   |-- genomad/{sample}/    + genomad_db/       only when phage.caller: genomad
|   `-- checkv/{sample}/     + checkv_db/
|-- 08.mobilome/                                 only when mobilome.run: true
|-- 09.report/                                   multiqc_report.html
`-- logs/                                        one file per rule per sample
```

---

## `01.reads` — read QC and filtering

Absent entirely in `contigs` mode, which has no reads to QC.

| Path | Written by | Contents |
|---|---|---|
| `{sample}/illumina/{sample}_fastp.html` | `trim_adapters` | fastp's adapter- and quality-trimming report, for a human. |
| `{sample}/illumina/{sample}_fastp.json` | `trim_adapters` | The same numbers for MultiQC. |
| `{sample}/illumina/{sample}_sel_R{1,2}.fastq` | `map_sel_contigs` | `[hybrid]` The Illumina pairs that mapped as proper pairs to the decontaminated Illumina assembly. Two rules read them: Filtlong scores the Oxford Nanopore reads against these pairs, so a contig dropped by the contamination screen takes its long reads with it; Polypolish then uses them to correct the final genome. Kept, not temporary — deleting them would force the whole short-read half to re-run. |
| `{sample}/ont/{sample}_filt.fastq` | `filter_long_reads` | The Oxford Nanopore reads Filtlong kept — what Flye assembled and what Medaka polished with. |
| `{sample}/ont/raw_qc/`, `{sample}/ont/filt_qc/` | `raw_long_read_qc`, `filtered_long_read_qc` | NanoPlot read-length and quality profiles, before and after filtering. Read side by side: the difference is what the filter did. |

The trimmed Illumina pairs (`{sample}_trim_R{1,2}.fastq`) and the PhiX reference and
its index under `01.reads/phix/` are all temporary and are gone when the run finishes.
See [Illumina mode](../modes/illumina.md) and [Nanopore mode](../modes/nanopore.md).

---

## `02.assembly` — the genome and everything used to judge it

### The delivered genome

`02.assembly/{sample}/contigs_final.fasta` is the one hand-off between the front end
and everything after it. Every rule in stages `03` to `08` works from this file or
from something built out of it, so nothing downstream has to know which assembler
ran. The one deliberate exception is `hybrid`, where QUAST, CheckM and GTDB-Tk are
also handed the Illumina assembly as a comparator — see `eval/` below.

What it is differs by mode:

| Mode | `contigs_final.fasta` is | Also kept |
|---|---|---|
| `illumina` | the decontaminated SPAdes assembly — decontamination *is* the last assembly step | `contigs_filt.fasta`, the pre-decontamination draft |
| `contigs` | the decontaminated input assembly | `contigs_filt.fasta` |
| `nanopore` | the Medaka consensus of the decontaminated, reoriented Flye assembly | `contaminants/assembly_decontam.fasta` |
| `hybrid` | the Oxford Nanopore assembly after Medaka and Polypolish | `contaminants/contigs_sel.fasta`, the decontaminated Illumina assembly, kept as the QC comparator and as Snippy's reference |

Set `parameters.{nanopore,hybrid}.medaka_model: false` and the consensus step drops
out of the last two rows: `nanopore` then delivers the decontaminated reoriented
assembly as it stands, `hybrid` delivers that assembly after Polypolish alone.

In `contigs` mode, what `contigs_filt.fasta` holds depends on the headers you fed in.
The length and coverage filter reads both numbers out of the header itself, so it can
only run on SPAdes-style headers: those keep their names and lose any contig under
500 bp or 2× coverage. Any other header style has no numbers to read, so nothing is
filtered and the headers are trimmed to their first token instead. Which branch ran
is the first line of `logs/filter_contigs_{sample}.log`.

In `hybrid` the contamination screen runs on the Illumina draft while the delivered
genome comes from the Oxford Nanopore leg, so those two FASTAs have entirely different
contig names. That is not a bug and it is why the plasmid step gets its own BLAST
table (below).

### `contaminants/`

The contamination screen and its audit trail, all under
`02.assembly/{sample}/contaminants/`.

| File | Contents |
|---|---|
| `{sample}_blastout` | The megablast table over the draft assembly: 15 tab-separated columns with the subject title last. |
| `{sample}_final_blastout` | `[nanopore, hybrid]` A second megablast, over `contigs_final.fasta`. It exists because the plasmid step looks contigs up **by ID**, and in the long-read modes the screened contigs and the delivered contigs do not share names. |
| `bestscore.blob.blobDB.table.txt` | The BlobTools per-contig table — taxonomy and coverage. The evidence behind every keep/drop decision, so it is kept. |
| `contig_taxonomy_decisions.tsv` | **The audit file.** Every contig, its assigned genus, and the reason it was kept or dropped. |
| `contigs.list` | The kept contig IDs, one per line. |
| `{sample}_composition.txt` | One line per genus with its share of the DNA and of the contig count, over *every* contig in the BlobTools table — kept or dropped. Also read by the annotation step as Bakta's genus hint. |

`blob.blobDB.json`, the coverage file beside it, and the alignment BAM are all
temporary. Full reasoning, and the two routes by which this step can delete a real
plasmid, on [Decontamination](../analysis/decontamination.md).

### `eval/`

| Path | Contents |
|---|---|
| `{sample}_qc_genomes.tsv` | Which genome is which: bin id, role (`primary` or `comparator`), technology, source path. In `hybrid` this is the file that says which of the two rows is the delivered genome. |
| `quast/` | QUAST assembly metrics — contig count, N50, largest contig, total length, GC. |
| `checkm/{sample}_checkm_stats.tsv` | CheckM completeness and contamination. `lineage.ms` beside it is CheckM's own marker-set file. |
| `qualimap/` | Coverage depth and evenness of the reads on the assembly. Absent in `contigs` mode, where the only alignment is the contigs against themselves and a report on it would say nothing. |

In `hybrid`, QUAST, CheckM and GTDB-Tk each run over **two** genomes per sample —
`{sample}_illumina` and `{sample}_ont` — so one report holds both and the comparison
is the point. The staging directory `eval/genomes/` that carries them there is
temporary. See [Assembly QC](../analysis/assembly-qc.md).

### Front-end working directories

Kept as the assembler left them, for when you need to see what happened.

| Directory | Modes | Contents |
|---|---|---|
| `spades/` | `[illumina, hybrid]` | SPAdes' own output, including `contigs.fasta` before the length and coverage filter. |
| `flye/` | `[nanopore, hybrid]` | `assembly.fasta`, Flye's `assembly_info.txt` (which carries the circularity call), and `ignore_list.txt` — the contigs Flye did not close, which dnaapler must not rotate. |
| `fix_start/` | `[nanopore, hybrid]` | dnaapler's reoriented assembly, the header-trimmed copy, and the reorientation summary. |
| `medaka/` | `[nanopore, hybrid]` | `consensus.fasta` and Medaka's working files. Absent when Medaka is switched off. |
| `polypolish/` | `[hybrid]` | Empty at the end — every file in it is temporary. |
| `snps/` | `[hybrid]` | One numbered sub-directory per Snippy run — `01.flye_snps_dir` through `04.polypolish_snps_dir` — each comparing one Oxford Nanopore stage against the Illumina assembly, plus `SNPs_summary.txt`. With Medaka off, `03` holds a `skipped.txt` instead of a comparison. |
| `{sample}_replicons.tsv` | `[nanopore, hybrid]` | The five-column table handed to Bakta as `--replicons`: topology from Flye, replicon type from dnaapler. `{sample}_replicons_audit.tsv` beside it carries the raw numbers and a reason per contig. |
| `{sample}_medaka_model.txt` | `[nanopore, hybrid]` | The Medaka model that was confirmed or inferred, written early so a bad model fails in seconds rather than after the assembler has run. Absent when Medaka is switched off, along with the rule that writes it. |

---

## `03.taxonomy` — GTDB-Tk

One directory per sample, named by GTDB-Tk itself. The file to read is
`classify/gtdbtk.bac120.summary.tsv` (`gtdbtk.ar53.summary.tsv` for archaea): the
placement, the closest reference and the ANI to it.

MultiQC reads that summary, and — when the mobilome module is on — so does the step
that maps the species to an AMRFinderPlus `--organism`, which is what unlocks
point-mutation detection. See [Taxonomy](../analysis/taxonomy.md).

---

## `04.annotation`

| Path | Contents |
|---|---|
| `bakta/{sample}/` | The full Bakta annotation set — `{sample}.gff3`, `.tsv`, `.faa`, `.fna`, `.gbff` and more. Three later stages reach into this directory: eggNOG takes the `.faa`, antiSMASH the `.gbff`, dbCAN both the `.faa` and the `.gff3`. |
| `eggnog/{sample}/` | `{sample}.emapper.*` — orthologous group per protein, with COG category, GO, KEGG KO and pathway, EC. Terminal: nothing else reads it. |
| `antismash/{sample}/` | Secondary-metabolite clusters, with antiSMASH's own HTML report. Terminal. |
| `antismash/databases/` | The shared antiSMASH reference database: downloaded here, or a symlink view of a copy you already hold. |
| `dbcan/{sample}/` | CAZyme calls, gene clusters and substrate predictions. Terminal. |
| `dbcan/<db version>/` | The shared dbCAN database, likewise. The folder is named from the download URL, so the version is visible on disk. The empty `.verified.sha256` beside it is the marker written once the checksum matched, and its presence is what stops the download repeating. |

eggNOG is the slow tail of a BacFlux run. See [Annotation](../analysis/annotation.md).

---

## `05.amr`

Two legs run independently, and they are meant to disagree at the edges: one sees the
assembly, the other sees the reads.

### `abricate/{sample}/`

One TSV per database — `argannot`, `card`, `ecoh`, `ecoli_vf`, `megares`, `ncbi`,
`resfinder`, `vfdb` — plus `AMR_summary.txt`, which collates a sample's eight tables.
This leg carries the EFSA thresholds: a hit needs ≥80% identity and ≥70% gene-length
coverage.

### `mapping/{sample}/` — `[illumina, hybrid]`

Trimmed reads mapped onto CARD's protein homolog model. Reads are immune to assembly
collapse, so this leg can see determinants that sit on a contig the assembler
destroyed.

| File | Contents |
|---|---|
| `{sample}_CARD_report.tsv` | **The one to read.** One row per CARD sequence that cleared the coverage gate, with CARD's own classification joined on: `aro_accession`, `aro_name`, `covered_percent`, `category`, `resistance_mechanism`, `amr_gene_family`, `drug_class`, `reference_organism`, `note`. `category` is what stops an efflux-pump subunit being counted as an acquired resistance gene. |
| `{sample}_covstats.tsv` | BBMap's raw coverage table, sorted by descending covered percent. The evidence behind the report. |
| `{sample}_AMR_legend.tsv` | Every feature covered ≥70%, with its row from CARD's `aro_index.tsv`. |

!!! warning "What this leg is actually screened at"

    BBMap is passed `minid=0.76` — its own default (`CARD_MIN_IDENTITY`,
    `workflow/rules/shared/00_common.smk:394`). The rule used to pass `idfilter=0.99`
    instead, which BBMap does not apply to the primary alignment of a properly-paired
    read, so 0.76 is what has been in force throughout. Read this as a **sensitive
    screen**: its specificity comes from requiring reads to cover at least **70% of the
    reference gene's length** (`CARD_MIN_COVERED`, same file, line 395), not from
    per-read identity. There is one mapping pass, and no "divergent" tier.

`05.amr/card_db/` holds the CARD database this run used — downloaded there, or a
symlink view of the copy at `directories.card_db`. Either way it is temporary, and
goes once the mapping and the report have both read it. See
[Antimicrobial resistance](../analysis/amr.md).

---

## `06.plasmids`

| Path | Contents |
|---|---|
| `{sample}/platon/verified_plasmids.txt` | The terminal plasmid product on the default path: the contigs Platon called plasmid, each annotated with whether its best nucleotide hit also says "plasmid". That check is **supplementary and never a filter** — a mobile element on a genuinely chromosomal contig can push both signals the same wrong way. |
| `{sample}/platon/contigs_final.tsv` | Platon's per-contig table with its replicon-distribution score. Platon names every file after the input basename, hence `contigs_final.*`. |
| `{sample}/platon/contigs_final.chromosome.fasta` | The contigs Platon called chromosome. |
| `{sample}/{sample}_plasmid_concordance.tsv` | *Optional.* Written only when `phage.caller: genomad`. One row per contig joining Platon's call and geNomad's, with a confidence column: both callers agree → high, one only → medium, a real clash → low. Nothing is dropped; a disagreement is flagged and kept. |

See [Plasmids](../analysis/plasmids.md).

---

## `07.phages`

Which caller runs is set by `phage.caller`, and the choice changes what stage `06`
holds as well as this one.

| Path | Contents |
|---|---|
| `virsorter/{sample}/` | *Default.* VirSorter2's calls at a deliberately loose score cutoff of 0.5; the file CheckV grades is `final-viral-combined.fa`. |
| `genomad/{sample}/` | *Only when opted in.* One geNomad run produces both the virus calls (`contigs_final_summary/contigs_final_virus.fna`) and the plasmid calls (`contigs_final_summary/contigs_final_plasmid_summary.tsv`) that stage `06`'s concordance table needs. |
| `checkv/{sample}/quality_summary.tsv` | **The phage deliverable.** Completeness, contamination and a quality tier per predicted viral sequence. A genome with no viral calls gets a header-only file rather than a failed run. |
| `vs2_db/`, `genomad_db/`, `checkv_db/` | The caller's and CheckV's databases. Each is either downloaded here or built as a symlink view of a copy you already hold — BacFlux only ever reads your copy. |

See [Prophages](../analysis/phages.md).

---

## `08.mobilome` — opt-in

Written only when `mobilome.run: true`. That adds ten rules on the default path; the
module's other twelve belong to four optional layers, each of which stays off until
you configure a source for it. See [Turning it on](../mobilome/enabling.md).

Databases and model sets sit at the top of the stage, fetched once and shared by every
sample: `conjscan_models/`, and — each only when its layer is configured —
`icescan_models/`, `tncentral_db/`, `iceberg_db/`, `isosdb_db/`.

Per sample, `08.mobilome/{sample}/` holds:

| File | Contents |
|---|---|
| `{sample}_amr_mobility.tsv` | **The deliverable.** One row per AMR gene: its mobile-element context, the element's identity where one is known, the contig-edge flags, a mobility tier from 1 to 6 and a confidence. |
| `{sample}_amr_mobility_audit.tsv` | Why every gene that got no context got none. |
| `{sample}_is_summary.tsv` | The honesty metric: how many insertion-sequence calls sit within `mobilome.contig_boundary_bp` of a contig end, and what fraction of the total that is. **Read this before the table above.** |
| `{sample}_is_elements.tsv`, `{sample}_is_discarded.tsv` | One tidy row per insertion sequence, and the dropped rows with a reason. |
| `{sample}_ice_candidates.tsv`, `{sample}_ice_discarded.tsv` | ICE/IME candidates with their anchors, boundaries and class; and the discard trail. |
| `{sample}_amrfinderplus.tsv`, `_amrfinderplus_mutations.tsv` | The AMR calls with coordinates and a method column; point mutations separately. |
| `{sample}_amrfinder_organism.txt`, `_amrfinder_organism_audit.tsv` | Which AMRFinderPlus `--organism` the GTDB-Tk placement mapped to, or none, and why. |
| `{sample}_replicon_calls.tsv` | Chromosome or plasmid, per contig, with the evidence for the call. |
| `{sample}_contig_lengths.tsv` | The yardstick every contig-edge flag is measured against. |
| `isescan/`, `conjscan/` | The tools' own output trees, kept as raw evidence. `icescan/` joins them when that optional second model set is switched on. |
| `{sample}_named_elements.tsv`, `{sample}_named_elements_discarded.tsv`, `{sample}_tncentral_blast.tsv` | *Optional.* Curated transposons and integrons — the layer that makes tier 4 reachable. |
| `{sample}_ice_candidates_named.tsv`, `{sample}_ice_naming.tsv`, `{sample}_iceberg_blast.tsv` | *Optional.* Curated names put on the ICE/IME candidates. Labels only: this layer cannot change any gene's tier. |
| `{sample}_is_copy_number.tsv`, `{sample}_is_copy_number_audit.tsv`, `{sample}_isosdb_covstats.tsv`, `{sample}_assembly_covstats.tsv` | *Optional, `[illumina, hybrid]`.* How many insertion-sequence copies the assembler collapsed, from read depth rather than from the assembly. |

Every column of the mobility table, with real rows, is on
[Reading the output](../mobilome/output.md). What fragmentation does to a call, and
what not to conclude from one, is on
[Draft assemblies](../mobilome/draft-assemblies.md).

---

## `09.report`

| File | Contents |
|---|---|
| `multiqc_report.html` | fastp or NanoPlot, Qualimap, QUAST, CheckM, GTDB-Tk and Bakta for the whole batch, with sample names rewritten into readable panel labels. |
| `multiqc_config.yaml` | The generated config, kept so the renaming rules are inspectable. |
| `multiqc_data/` | MultiQC's parsed numbers behind the report. |

Left alone, MultiQC names a sample after every directory component of the file it came
from. The generated config rewrites those into `<what it is> | <sample>` — so a panel
row reads `assembly QC | S1`, and a hybrid sample's two genomes appear separately as
`completeness Illumina | S1` and `completeness ONT | S1`.

The report directory is deliberately not a declared Snakemake output, so anything you
drop in it survives a re-run.

---

## `logs`

One file per rule invocation, flat, with the sample appended where the rule runs per
sample — `logs/annotation_S1.log`, `logs/multiqc.log`. The stem is usually the rule's
own name, but not always: every mobilome rule logs under a `mobilome_` prefix, so rule
`isescan` writes `logs/mobilome_isescan_S1.log`, and a handful elsewhere are shortened
— `map_amr_db` writes `logs/map_amr_S1.log`. The log is the first place to look when a
rule fails, and the only place some warnings go: the decontamination selector's
"N contigs from the BlobTools table were not found in the FASTA", for instance, is
written there rather than to the console.

There is no `benchmarks/` directory: no rule declares Snakemake's `benchmark:`
directive, so no runtime or memory record is written anywhere.

---

## Temporary files, and directories a re-run wipes

Some outputs are declared temporary: Snakemake deletes them as soon as every rule that
needed them has finished. They are absent from a completed run, which is normal and
not a failure.

| Gone by the end | Why it existed |
|---|---|
| `01.reads/{sample}/illumina/{sample}.{1,2}.fastq` and `{sample}_trim_R{1,2}.fastq` | The PhiX-free pairs and then the fastp-trimmed pairs. The trimmed ones survive until the last rule that reads them is done: the assembler, the coverage alignment, the CARD leg, and in `hybrid` the mapping that picks the Illumina pairs to keep. |
| `01.reads/phix/` | The PhiX reference and its index, needed only while reads are being screened. |
| `02.assembly/{sample}/contaminants/{sample}_map.bam`, the BlobTools JSON and coverage file | Intermediates of the screen; the BlobTools *table* is kept. |
| `02.assembly/{sample}/eval/genomes/` | Staging copies of the assemblies under their bin ids, for QUAST, CheckM and GTDB-Tk. |
| `02.assembly/{sample}/polypolish/` contents | Alignments and the pre-correction draft. |
| `04.annotation/eggnog/{sample}_tmp/` | emapper's scratch space. |
| `05.amr/card_db/`, and `05.amr/card.tar.bz2` on the download path | The CARD database, needed only until the mapping and the report have read it. |
| `08.mobilome/{sample}/assembly_depth_ref/`, `isosdb_ref/` | BBMap indexes for the copy-number leg. |
| `09.report/multiqc_inputs/` | Relabelled CheckM and GTDB-Tk copies staged for MultiQC. |

Nearly every tool directory on this page is a declared Snakemake **directory** output:
the per-sample trees written by Bakta, eggNOG, antiSMASH, dbCAN, QUAST, CheckM,
Qualimap, GTDB-Tk, Platon, the phage caller, CheckV, ISEScan and CONJscan, and every
database directory the workflow fetches. Snakemake deletes a directory output before
re-running the rule that makes it, so anything you put inside one will not survive
that rule re-running. Keep your own files beside them, not in them. `09.report/` is
the deliberate exception, which is why it is safe to drop things there.

---

## See also

- [Configuration](configuration.md) — every key referred to above, and its default.
- [Running BacFlux](running.md) — dry runs, `--cores`, and what a re-run rebuilds.
- [Reference databases](../getting-started/databases.md) — which databases you supply
  and which the workflow fetches into the output directory.
- [Troubleshooting](../troubleshooting.md) — an empty or missing output, and what it
  usually means.
