# Taxonomy

Stage `03.taxonomy` answers one question — what is this isolate? — with one tool,
GTDB-Tk 2.7.2 against GTDB **R232**.

`classify_wf` finds the ~120 bacterial marker proteins in the genome, aligns them and
places them on the Genome Taxonomy Database reference tree, then refines the call with
skani ANI comparisons against the closest reference genomes. The answer is a full
lineage from domain to species with the ANI evidence behind it.

That matters for a food- and feed-safety workflow: whether an AMR gene found later
reads as intrinsic or acquired depends on having the species right; see
[the mobility ladder](../mobilome/mobility-ladder.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ANI** | average nucleotide identity, the percentage identity between two genomes over the sequence they share |
    | **AMR** | antimicrobial resistance |

## The rule

*Rule `taxonomic_assignment`, `workflow/rules/shared/30_taxonomy.smk`.*

```bash
GTDBTK_DATA_PATH=/path/to/release232 \
gtdbtk classify_wf -x fasta \
  --genome_dir 02.assembly/{sample}/eval/genomes \
  --out_dir 03.taxonomy/{sample} \
  --cpus {threads} --pplacer_cpus {threads}
```

**In:** the staged genome directory that
[`stage_qc_genomes`](assembly-qc.md#staging-one-genome-or-two) filled — this sample's
one genome, or two in `hybrid` mode, named after their bin ids.

**Out:** `03.taxonomy/{sample}/`, a directory, because GTDB-Tk names the files inside
it. The one to read is the summary TSV, `classify/gtdbtk.bac120.summary.tsv`. GTDB-Tk
has moved that file between releases, so the report stage looks inside `classify/`
first and then at the top level.

**Consumed by:** MultiQC, which reads the summary; and — only when the mobilome module
is on — the rule that translates the species call into an AMRFinderPlus `--organism`
name, which is what unlocks point-mutation detection. GTDB and NCBI do not use the same
species names, so that translation goes through a curated lookup table generated per
GTDB release; an unmapped species degrades silently to no `--organism` at all, which is
the common case for environmental isolates. How the table is built and used is on
[the AMR page](amr.md#gtdb-to-amrfinderplus-organism).

The rule also takes CheckM's statistics file as an input. That is a **scheduling edge,
not a data dependency**: GTDB-Tk reads only the genome directory. Both tools run
pplacer and each can use tens of GB, so on a large machine dropping the edge would let
them run concurrently on the same sample and risk an out-of-memory kill. Keeping it
costs nothing and preserves the ordering.

## Hybrid: one run, two rows

A hybrid sample stages two genomes, and one `classify_wf` run classifies both. The
summary comes back with two rows keyed by the staged basenames, `{sample}_illumina` and
`{sample}_ont`. They need no separate output files: the report rewrites the keys into
`taxonomy Illumina | {sample}` and `taxonomy ONT | {sample}`, and
`02.assembly/{sample}/eval/{sample}_qc_genomes.tsv` records which of the two is the
delivered genome.

Two rows that disagree are worth a look. They are the same DNA sequenced two ways, so a
disagreement points at one of the assemblies rather than at the taxonomy.

## GTDB names are not NCBI names

GTDB splits genera and gives them a letter suffix (`Pseudomonas_E`), and it uses
placeholder species names for lineages with no published name (`sp024807945`). That is
a property of GTDB, not a fault — it is what makes the taxonomy consistent with genome
phylogeny.

It does mean the name you get here will not always match a name in NCBI-derived
resources. Where BacFlux needs to cross that gap it does so explicitly, with a curated
lookup, rather than assuming the two vocabularies agree.

## Version and database are pinned together

GTDB-Tk hard-pins itself to **one** compatible reference release, in its own
`COMPATIBLE_REF_DATA_VERSIONS`. Version 2.6.1 accepted only r220 and r226, so GTDB R232
requires 2.7.0 or later; the environment pins **2.7.2**. There is no in-place upgrade
from an older release directory — download R232 and point `directories.gtdbtk_db` at it.
Sizes and the download recipe are in
[Reference databases](../getting-started/databases.md).

The version bump also changed how the skani reference is laid out. Under R226 the
release shipped a raw skani reference and GTDB-Tk built its own sketch cache on first
use. From R232 the release ships that sketch **pre-built**, which is what roughly
halves the space GTDB-Tk needs — the R232 package extracts to 94 GB and no longer
builds a ~57 GB sketch cache beside it. `classify_wf` reads it
straight out of `{gtdbtk_db}/skani/`, and the `--skani_sketch_dir` flag no longer
exists in the command line at all, so there is nothing left for the rule to build or
point at.

## Resources

`--cpus` and `--pplacer_cpus` are both set to `min(resources.threads, 24)` — one knob,
and `--cores` caps what Snakemake grants on top of it. pplacer's
**memory** scales with its thread count and can reach tens of GB, which is why the cap
exists. Thread count changes speed and memory only, never the classification.

This is the single most expensive step in a BacFlux run, which is also why a hybrid
sample is classified in one invocation rather than two: `classify_wf` loads the
reference tree and the skani sketches on every call.

## What to check afterwards

- `03.taxonomy/{sample}/classify/gtdbtk.bac120.summary.tsv` — the lineage, the closest
  reference genome and its ANI. A species-level call needs the ANI to clear GTDB's
  threshold; below it, the row stops at genus.
- The genus here against `{sample}_composition.txt` from
  [decontamination](decontamination.md#the-composition-report). A mismatch means the
  contig filter and the classifier disagree about what the isolate is, and the
  decontamination audit file is where to settle it.
- `09.report/multiqc_report.html` — the classification for the whole batch on one page.
