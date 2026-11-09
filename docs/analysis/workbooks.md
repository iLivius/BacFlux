# Per-strain workbooks

BacFlux writes one Excel workbook per strain with the antimicrobial-resistance (AMR)
genes found in its genome and their genomic context and, with the
[trait profile](traits.md) on, a second one with the strain's predicted traits. Both
land in `10.workbooks/`:

| Workbook | Written | What it holds |
|---|---|---|
| `{sample}_amr.xlsx` | in every run | the resistance genes ABRicate found in CARD, ResFinder and ARG-ANNOT, one row per gene, each with three answers: is it on a plasmid, inside an ICE or IME, inside a phage region? |
| `{sample}_traits.xlsx` | with `traits.run: true` | every trait of the catalogue with its call, the genes behind it, its caveat, a suggested laboratory assay and its references |

*Rules `amr_workbook` and `traits_workbook`, both running
`workflow/scripts/90_report/final_reports.py`.*

| | |
|---|---|
| **In** | the final contigs, GTDB-Tk, the ABRicate tables, Platon, geNomad's plasmid concordance (with `phage.caller: genomad`), the virus caller and CheckV; with `mobilome.run: true`, the mobilome module's replicon, ICE and mobility tables; with `traits.run: true`, the trait tables |
| **Out** | `10.workbooks/` |
| **Next** | user interpretation |

The workbooks make no new call. They read the tables written by the stages before
them, group the database hits of one gene into one row, and compare coordinates to
place each gene in its context. Every row they start from is kept on its own sheet,
so each answer can be traced to the tool that gave it. Every cell has a type:
coordinates and percentages are numbers that sort and filter as such, and a value
such as ABRicate's `GAPS` entry `6/4` stays text instead of being turned into a date.

!!! warning "A prediction, not a test result"

    A resistance gene is not a susceptibility test result, a mobile context does not
    show that the gene transfers, and a trait call does not show that the strain
    expresses the trait. The workbooks keep the evidence, the caveats and the
    suggested assays beside each answer for that reason.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AICE** | actinomycete integrative and conjugative element |
    | **ICE** | integrative and conjugative element: moves itself from cell to cell |
    | **IME** | integrative and mobilisable element: moved by a conjugative element in the same cell |
    | **locus** | one gene in the genome, however many database hits it has |

## Four answers, not two

Each context question is answered with one of four words.

| Answer | Meaning |
|---|---|
| `yes` | the evidence was found |
| `no` | the analysis ran and found none |
| `uncertain` | the tools disagree, or the evidence stops short of a yes |
| `not_assessed` | the analysis that could answer did not run: the mobilome module was off, or Platon failed for this strain |

An analysis that did not run is therefore never written as `no`. The rule behind each
answer is given below, and in the `Data_dictionary` sheet of every workbook.

## The AMR workbook

| Sheet | One row per | Read from |
|---|---|---|
| `Summary` | item | the counts of the sheets below, and the GTDB-Tk classification |
| `AMR_loci` | resistance gene | all the sources below |
| `Database_hits` | ABRicate hit | `05.amr/abricate/{sample}/{card,resfinder,argannot}.tsv` |
| `ICE_IME_elements` | mobile element | `08.mobilome/{sample}/{sample}_ice_candidates.tsv`, or its `_named` copy when ICEberg naming is on |
| `Viral_regions` | virus or prophage region | the virus caller's region table, and CheckV's `quality_summary.tsv` and `contamination.tsv` |
| `Mobility_evidence` | AMRFinderPlus gene | `08.mobilome/{sample}/{sample}_amr_mobility.tsv`, copied unchanged |
| `Mobility_audit` | audit row | `08.mobilome/{sample}/{sample}_amr_mobility_audit.tsv`, copied unchanged |
| `Run_QC` | input or setting | every input, as a path relative to the output folder, and what was read from it |
| `Data_dictionary` | column | the meaning of every column |

Without the mobilome module, the three mobilome sheets hold a single line saying so.

### One gene, one row

The three databases often give one gene different names (`aadA2`, `ANT(3'')-Ia`,
`(AGly)aadA2`), so hits are grouped by their position, never by name. The longest hit
of a gene defines the locus, and a shorter hit joins it when at least 80 % of the
shorter hit lies inside it. Hits are taken longest first and compared with these
longest hits only, so a short hit lying where two genes overlap cannot chain the two
genes into one row. The `locus_id` is built from the sample, the contig and the
coordinates of the longest hit, as in `SAMPLE1|contig_3:3176892-3177704`; a database
version that changes those coordinates changes the ID. Every hit stays on
`Database_hits`, with its identity, coverage and the line of the table it came from.

`reference_accessions` holds the GenBank accessions of the database's reference
sequences, with their coordinates. CARD's ARO numbers are not in ABRicate's output.

`in_resfinder` says whether ResFinder found the gene. ResFinder collects acquired
resistance genes ([Zankari et al. 2012](https://doi.org/10.1093/jac/dks261)), while CARD
and ARG-ANNOT also list genes a species carries anyway, such as efflux pumps, porins and
their regulators; a gene found only by CARD or ARG-ANNOT is therefore worth a look
before it is read as acquired. The `Summary` sheet counts the loci ResFinder found, and
gives ABRicate's version and the date ABRicate lists for its copy of each database
(also on `Run_QC`, from `05.amr/abricate/abricate_databases.tsv`), so each result is
reported with the database it came from.

### On a plasmid?

The answer reads the replicon call of the gene's contig: the mobilome module's
`{sample}_replicon_calls.tsv` when the module ran, otherwise the same call, made by the
same code (`platon_replicons.py`).

| Replicon call of the contig | `plasmid` |
|---|---|
| plasmid, by Platon, by geNomad or by both | `yes` |
| chromosome, by Platon, or inferred from a length of 2 Mb or more | `no` |
| chromosome by Platon or from a length of 2 Mb or more, plasmid by geNomad (call source `conflict`) | `uncertain` |
| unknown: neither tool classified the contig | `uncertain` |
| any, when Platon failed for the strain | `not_assessed`, unless geNomad calls the contig a plasmid |

`plasmid_call_source` says which tools the call rests on: `genomad` when geNomad alone
called a plasmid, `both`, `conflict`, and `platon` for every other call, including a
contig Platon skipped whose call comes from its length (`plasmid_evidence` then says
so). Platon skips every contig of 500 kb or more and every contig under 1 kb, and such a
contig is decided by geNomad's call or by its length, as explained under
[Plasmids](plasmids.md#what-platon-classifies-and-how-bacflux-reads-it). A plasmid is not
transferable by being a plasmid: `plasmid_mobility` gives Platon's reading of its
conjugation and mobilisation genes.

### Inside an ICE or IME?

Decided by coordinates for every gene, from the mobilome module's ICE table.

| Where the gene lies | `ice_ime` |
|---|---|
| inside an ICE, an IME or an AICE | `yes` |
| inside a conjugative region | `uncertain` |
| within 5,000 bp of an ICE, IME or AICE, without lying inside it | `uncertain` |
| anywhere else, a genomic island included | `no` |
| the mobilome module did not run | `not_assessed` |

"Inside" means that at least 90 % of the gene lies within the element's interval, and
"near" means within 5,000 bp: the mobilome module's own rules for its `mge_context`
column. Unless the element's attachment sites were found, its interval is the span of
its machinery genes, not its real ends, so a gene just outside it may still be carried
by it; that is why near is `uncertain` rather than `no`. A conjugative region is
conjugation machinery whose element could not be resolved: no integrase was found
beside it, or it lies on a plasmid, which makes it a conjugative plasmid rather than an
ICE. A genomic island has an integrase but no relaxase, so it cannot be moved by
conjugation on its own; it is named in `ice_ime_evidence`, as is a conjugative region
near the gene.

### Inside a phage region?

The virus caller (VirSorter2, or geNomad with `phage.caller: genomad`) delivers
regions, and [CheckV](phages.md#checkv) splits each prophage region into its viral part
and the host DNA around it.

| Where the gene lies | `phage` |
|---|---|
| on a part CheckV calls viral; a region without a host part is viral as a whole | `yes` |
| in a region of the caller, but only on the part CheckV calls host DNA | `uncertain` |
| in a region of the caller that CheckV did not assess | `uncertain` |
| in no region of the caller | `no` |

The caller's regions are generous on purpose (VirSorter2 runs at a loose cutoff and
leaves the grading to CheckV), so a region can carry tens of kilobases of host DNA at
its ends. `phage_evidence` gives the contig bases CheckV calls viral, and
`Viral_regions` lists both parts of every region.

### The mobility tier

With the mobilome module on, `mobility_tier` and the columns beside it come from the
module's mobility table, which rates the genes AMRFinderPlus found. Those genes are
matched to the database loci by position: same contig, and at least 80 % of the
shorter gene inside the longer one. Resistance genes (`AMR`) and biocide and metal
resistance genes (`STRESS`) are matched; the quaternary-ammonium efflux gene
`qacEdelta1`, for instance, is a `STRESS` gene for AMRFinderPlus and a resistance gene
for the databases. Point mutations are not matched: AMRFinderPlus reports one with the
coordinates of the whole gene it lies in, which a database hit on that gene would
otherwise take for its own. When several rows match one locus, the highest tier is
shown. The tiers are explained under [the mobility ladder](../mobilome/mobility-ladder.md).

### What the AMR workbook leaves out

- The other five ABRicate databases (EcoH, E. coli VF, MEGARes, NCBI and VFDB): their
  tables stay in `05.amr/abricate/{sample}/`.
- The reads mapped to CARD in the short-read modes: see `05.amr/mapping/`.
- AMRFinderPlus's point mutations (for example in `gyrA`), virulence genes, and any
  other AMRFinderPlus gene no database hit matches. With the mobilome module on, they
  are on `Mobility_evidence`, with `locus_ids` `NA`.

## The traits workbook

| Sheet | One row per | Read from |
|---|---|---|
| `Summary` | item | the counts of the calls |
| `Trait_predictions` | trait of the catalogue | `{sample}_trait_profile.tsv`, with the suggested assay and the references of `trait_rules_v1.tsv` |
| `Gene_evidence` | protein and gene with trait evidence | `{sample}_trait_genes.tsv`, with the plasmid and phage answers of each gene |
| `Prediction_audit` | trait | `{sample}_trait_profile_audit.tsv`, copied unchanged |
| `Trait_rules` | trait | `trait_rules_v1.tsv`: the gene sets of each rule |
| `Run_QC`, `Data_dictionary` | | as in the AMR workbook |

The calls and their columns are explained under [Trait profile](traits.md). Each gene's
plasmid and phage answers follow the rules above; the ICE/IME question is not asked for
trait genes.

Before quoting the catalogue, the workbook checks `{sample}_trait_provenance.json`,
written with the trait tables. It holds the SHA-256 fingerprint (a code computed from
every byte of a file) of the two catalogue tables, of EcoFoldDB's descriptions and of
the three trait tables, as they were when the calls were made. If any of these files has
changed since, the workbook stops and asks for rule `trait_profile` to be run again, so
a caveat or a reference is never shown beside a call made with another version of the
catalogue. Only contents are compared, so an output folder that has been moved still
passes. The `Summary` sheet states that laboratory assays were not compared, since the
workbook takes no laboratory results, and that safety is not assessed there: a trait
call says nothing about whether the strain is safe to use.

## Working with the files

A workbook is made again whenever one of its inputs changes, and the new file replaces
the old one, so keep notes in a copy elsewhere, or in a separate table keyed by
`locus_id` or `trait_id`. The paths in a workbook are relative to the output folder, so
it can be passed on without revealing where the analysis ran. A project that needs
other file names can copy the workbooks under those names after the run.
