# Assembly QC

Three questions about the finished genome, one tool each:

| Question | Tool | Answers with |
|---|---|---|
| Is the **assembly** any good? | QUAST 5.3.0 | contig count, N50, largest contig, total length, GC |
| Is the **genome** complete, and is it one organism? | CheckM 1.2.5 | completeness %, contamination %, from lineage-specific single-copy markers |
| Did the **reads** map back sensibly? | Qualimap 2.3 | mean depth, evenness, GC bias, mapping rate, insert size |

The first two need a genome; the third needs the BAM that
[decontamination](decontamination.md#1-coverage-the-reads-mapped-back-onto-the-draft)
already made for BlobTools. All four results — these three plus GTDB-Tk — end up in the
MultiQC report in `09.report`.

Everything here is in `workflow/rules/shared/20_qc.smk` and runs in all four modes,
except Qualimap.

!!! note "Decontamination runs first, and has to"

    CheckM is where a leftover contaminant contig shows up, as contamination *of the
    isolate*. Screening the assembly before QC is what makes the CheckM number mean
    what you think it means.

## Staging: one genome, or two

*Rule `stage_qc_genomes`.*

Three modes deliver one genome per sample. `hybrid` delivers the ONT + Polypolish
genome **and** keeps the decontaminated Illumina draft, so the two can be compared
directly. Rather than duplicate rule bodies or add a wildcard, one small bookkeeping
rule copies this sample's assemblies into a temporary directory under controlled
names, and QUAST, CheckM and GTDB-Tk all read that one directory.

| Mode | Staged as | Role |
|---|---|---|
| `illumina`, `nanopore`, `contigs` | `{sample}.fasta` | primary |
| `hybrid` | `{sample}_illumina.fasta` | comparator (the decontaminated SPAdes draft) |
| `hybrid` | `{sample}_ont.fasta` | primary (the delivered genome) |

Those basenames become the ids CheckM, GTDB-Tk and QUAST report under, so they come
from a single list rather than being typed out per tool. The staging directory itself
is temporary and Snakemake removes it once the last consumer is finished.

Beside it sits a file that is kept: `02.assembly/{sample}/eval/{sample}_qc_genomes.tsv`,
with columns `bin_id`, `role`, `technology`, `source_path`. It says in writing which of
the two hybrid rows is the delivered genome, so nobody has to infer it from a filename
suffix. No rule reads it; it exists to be read by a person.

!!! warning "A hybrid sample does twice the work inside one job"

    CheckM and GTDB-Tk each process both genomes in a single run, under one thread
    reservation. Expect roughly double the wall-clock of a single-genome mode for
    those two rules.

## Assembly metrics (QUAST)

*Rule `genome_assembly_evaluation`.*

```bash
quast eval/genomes/*.fasta -o eval/quast --no-icarus -t {threads}
```

QUAST reports the structural quality of the assembly. It says nothing about biological
completeness — that is CheckM's job — but it says how fragmented the assembly is, and
fragmentation limits everything downstream that depends on gene context. The mobilome
module in particular is built around that constraint; see
[Draft assemblies](../mobilome/draft-assemblies.md).

The shell globs `*.fasta` inside the staged directory, so a hybrid sample gets both
assemblies **side by side in one report** — exactly the comparison a hybrid run is for.
`--no-icarus` skips the interactive contig browser, which nothing here consumes.

Output: `02.assembly/{sample}/eval/quast/` (QUAST names the files inside it).

QUAST labels an assembly by the input file's basename, so the labels are `{sample}`, or
`{sample}_illumina` and `{sample}_ont` in hybrid. The MultiQC report renames them to
`assembly QC | {sample}`, or `assembly Illumina QC` / `assembly ONT QC`.

## Completeness and contamination (CheckM)

*Rule `completeness_and_contamination`.*

CheckM places the genome in its reference tree, picks the set of single-copy marker
genes for that lineage, and counts them.

- Markers **missing** → an incomplete assembly.
- Markers present in **more copies than expected** → either contamination (more than
  one organism in the bin) or a genuine duplication.

This is the standard pass/fail gate for an isolate genome.

```bash
checkm lineage_wf -t {threads} -x fasta eval/genomes/ eval/checkm/
checkm qa -o 2 -t {threads} --tab_table -f {sample}_checkm_stats.tsv lineage.ms eval/checkm/
```

`lineage_wf` places the genome and counts markers; `qa -o 2 --tab_table` writes the
extended per-bin statistics as a parsable TSV. CheckM takes several "bins" in one run
and reports a row for each, which is how the two hybrid genomes are handled.

| Output | What it is |
|---|---|
| `eval/checkm/` | CheckM's own working and result tree |
| `eval/checkm/{sample}_checkm_stats.tsv` | the table MultiQC reads |
| `eval/checkm/lineage.ms` | the marker-set file `lineage_wf` writes and `qa` reads back — the name is CheckM's |

Reading the two numbers together is what makes them useful:

| Completeness | Contamination | Usual meaning |
|---|---|---|
| high | low | a good isolate genome |
| high | high | two organisms in one assembly — a mixed culture that decontamination did not resolve |
| low | low | the decontamination filter probably removed genome, not contamination |
| low | high | poor assembly *and* contamination; start from the read QC |

The third row is the one to watch for after a large-removal warning from the
[selector](decontamination.md#the-two-warnings).

## Mapping quality (Qualimap)

*Rule `map_evaluation` — modes with reads only.*

Qualimap summarises how a sample's reads sit on its own assembly: mean depth and how
even it is, GC bias, mapping rate, and insert size for paired reads. Uneven or
unexpectedly low coverage is an early warning of a mixed culture, a mis-assembly, or
simply not enough data.

```bash
qualimap bamqc -bam {sample}_map.bam --java-mem-size={N}G -nt {threads} \
  -outdir eval/qualimap -outformat html
```

Input is the decontamination BAM — Bowtie2 in the short-read modes, minimap2 `map-ont`
in `nanopore`. That BAM is temporary, so this rule's dependency on it is what keeps it
alive long enough to be used twice.

Output: `02.assembly/{sample}/eval/qualimap/`, an HTML report directory.

The rule is gated on the mode having reads at all, so **`contigs` mode has no Qualimap
report**. That mode does have a BAM, but it is the contigs aligned against themselves
purely to satisfy BlobTools; a mapping-quality report on a self-alignment would say
nothing. It is a mode-specific absence, not a gap.

!!! warning "In hybrid, this panel is Illumina coverage of the Illumina draft"

    The BAM comes from the short-read branch of decontamination: Illumina reads mapped
    to the pre-decontamination SPAdes draft. The ONT reads are never mapped for QC. So
    this is **not** coverage of the delivered ONT genome, and MultiQC labels the panel
    `mapping Illumina QC` to say so. The QUAST, CheckM and GTDB-Tk panels are
    technology-tagged for the same reason.

## Where the two hybrid assemblies get compared

Three places, answering different questions:

| Where | Comparison |
|---|---|
| `eval/quast/` | both assemblies in one QUAST report — contiguity, length, GC |
| `eval/checkm/` and `03.taxonomy/` | one row per genome — completeness, contamination, lineage |
| `02.assembly/{sample}/snps/SNPs_summary.txt` | Snippy variant counts for each ONT polishing stage against the decontaminated Illumina **assembly** as reference; see [Hybrid mode](../modes/hybrid.md) |

## Resources

| Rule | Threads | Memory |
|---|---:|---|
| `stage_qc_genomes` | 1 | negligible — `mkdir` and `cp` |
| `genome_assembly_evaluation` | capped at 24 | modest |
| `completeness_and_contamination` | capped at 24 | tens of GB; pplacer is the driver |
| `map_evaluation` | capped at 24 | `min(resources.ram_gb, 64)` GB handed to the JVM |

"Capped at 24" means the rule asks for `min(resources.threads, 24)`, and `--cores`
caps what Snakemake then grants it. Qualimap's memory figure is
a JVM heap size in gigabytes, not a core count, which is why it is a resource rather
than a thread count; Snakemake only schedules against it if the launch line passes
`--resources java_mem=N`. See [Running BacFlux](../reference/running.md).

CheckM and GTDB-Tk both run pplacer and each can use tens of GB, so
[the taxonomy rule](taxonomy.md) deliberately waits for CheckM rather than running
alongside it on the same sample.

## Output files

All under `02.assembly/{sample}/eval/`.

| Path | Kept? | Read by |
|---|---|---|
| `{sample}_qc_genomes.tsv` | yes | you — which genome is which |
| `quast/` | yes | MultiQC |
| `checkm/{sample}_checkm_stats.tsv` | yes | MultiQC |
| `checkm/` (tree, `lineage.ms`) | yes | `checkm qa` |
| `qualimap/` | yes | MultiQC; absent in `contigs` mode |
| `genomes/` | no | QUAST, CheckM, GTDB-Tk — deleted once all three are done |

## What to check afterwards

- **CheckM completeness and contamination first.** Everything downstream inherits the
  quality of the genome, and the two numbers together tell you which way a problem
  points.
- **QUAST's contig count and N50.** A fragmented assembly is not wrong, but it limits
  what the gene-context stages can say — the mobilome module reports contig-edge flags
  on every row precisely because of this.
- **Qualimap's mean depth and evenness.** A bimodal depth distribution on a
  single-isolate assembly usually means more than one organism, or a plasmid at high
  copy number.
- `09.report/multiqc_report.html` — all of the above for the whole batch on one page.
