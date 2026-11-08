# Trait profile

The trait profile reads each genome for traits that matter in plant-growth promotion,
biocontrol, stress tolerance and bioremediation: nitrogen fixation, phosphate
solubilisation, siderophores, plant hormones, osmoprotectants, the breakdown of
aromatic compounds and the steps of denitrification. It serves the purpose BacFlux was
developed for: the genomic characterisation of strains considered for biocontrol,
biostimulant and biofertiliser applications. Each trait is called from a defined set of
genes: for most traits several genes that must all be found, often side by side as in
an operon; for a few, one marker gene.

The layer is off by default; `traits.run: true` switches it on, in all modes. It writes
to `04.annotation/traits/{sample}/`, from rules in
`workflow/rules/shared/40_annotation.smk`, and the summary across samples comes from
`workflow/rules/shared/90_report.smk`. It adds jobs and a section to the MultiQC
report, and changes no other stage's result.

| File | What it is |
|---|---|
| `{sample}_trait_profile.tsv` | **the result** — one row per trait: called, partial or not called, with the genes behind it |
| `{sample}_trait_profile_audit.tsv` | **why** — one row per trait, with a fixed reason for its decision |
| `{sample}_trait_genes.tsv` | **the evidence** — one row per protein and gene: what matched, and whether it counted |
| `{sample}_trait_provenance.json` | **the record** — the fingerprints of the catalogue and of the three tables above, with the calling parameters |
| `09.report/traits_summary.tsv` | every sample's calls side by side, one column per trait |

Read them in this order: the summary, to see the batch at a glance; then a sample's
profile; then its audit and gene table for any call you intend to quote.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **3Di** | Foldseek's structural alphabet: one letter per residue that describes its 3D neighbourhood |
    | **KO** | a KEGG Orthology number, KEGG's identifier for a gene family (`K02588`) |
    | **EC** | an Enzyme Commission number, the identifier of an enzyme reaction (`3.5.99.7`) |
    | **e-value** | how many hits at least this good the search expects by chance; smaller is better |
    | **coverage** | the fraction of the protein (`qcov`) or of the catalogue structure (`tcov`) that the alignment spans |

## How a gene is recognised

*Rules `trait_db` (or `trait_db_local`), `trait_structure_query`,
`trait_structure_search` and `trait_profile`. Off by default; `traits.run: true`
switches them on.*

| | |
|---|---|
| **In** | `{sample}.faa` and `{sample}.tsv` from the Bakta directory; EcoFoldDB; the ProstT5 model from Baktfold's database; `{sample}.json` from antiSMASH; Baktfold's `{sample}.tsv` when `baktfold.run: true` |
| **Out** | `04.annotation/traits/{sample}/` |
| **Next** | `09.report/traits_summary.tsv`, the MultiQC report and `10.workbooks/{sample}_traits.xlsx`; user interpretation |

A gene is recognised in three ways, and one protein can carry all three.

**Structure (S).** Every protein Bakta predicted, up to 4,000 amino acids, is turned
into 3Di by the ProstT5 protein language model, on the CPU, without predicting a full
structure. Foldseek then searches it against EcoFoldDB v2.0 (Ghaly et al. 2025) with
`-s 8`, as EcoFoldDB's own method does. EcoFoldDB holds 205,334 predicted structures,
842 of which form its catalogue (`EcoFoldDB_descriptions.txt`): proteins of described
function; the rest are other Swiss-Prot proteins. A hit is kept only when it passes
the authors' filter:

- the protein's best hit among all 205,334 structures (lowest e-value; the first listed
  wins a tie) is a catalogue entry;
- its e-value is below 1e-10;
- the alignment covers more than 0.8 of the protein and of the catalogue structure.

The first condition keeps look-alikes out: an enzyme that shares a fold with a
catalogue entry, but not its function, usually finds its best hit outside the
catalogue. A kept hit gives `S` to every gene whose `ecofolddb` entries include that
catalogue entry. A longer protein is not searched (`not_searched_too_long`), but its
names are still checked.

**Bakta name (Q).** A protein matches a gene by name when Bakta's gene symbol is
exactly one of the gene's symbols, when one of its KO or EC numbers is listed for the
gene, or when the gene's product pattern is found in Bakta's product. A symbol that
unrelated genes share can be set to count only together with the product pattern or a
KO (`symbol_needs_product`).

**Baktfold name (F).** With `baktfold.run: true`, a protein that Bakta called
"hypothetical protein" (a prefix such as `(pseudo)` aside) is tested the same way,
with the name Baktfold gave it. Baktfold finds its names with the same ProstT5 model
and Foldseek, so `S+F` is two structure searches agreeing, not a sequence-based name
confirming a structure hit. See [Baktfold](annotation.md#baktfold).

**antiSMASH regions** are supporting evidence only
([How a trait is called](#how-a-trait-is-called)).

### Borderline hits

A kept structure hit is flagged borderline when it lies close to a cut-off or to a
rival: its e-value is 1e-11 or above, within ten times the cut-off; its query or target
coverage is 0.85 or less; or the best hit to another catalogue gene, or to a structure
outside the catalogue, has an e-value at most ten times the best one. The flag never
changes a decision. It marks a hit that a small change in the input, such as another
assembly of the same strain, could push across a cut-off or hand to its rival.

## When a gene counts

A protein counts for a gene when a name matches (`Q` or `F`), or when its kept
structure hit (`S`) points to a gene whose `s_alone` is `yes`. Four points qualify
this; the gene table gives each protein's reason.

1. **A veto wins.** A gene can list products that rule it out (`product_exclude`):
   relatives acting on another substrate, unrelated genes with the same symbol, and
   regulators, receptors or binding proteins named after a pathway. When the pattern
   matches Bakta's or Baktfold's product, the protein does not count for that gene,
   whatever else matched, structure included.
2. **Most genes need a name.** For most genes `s_alone` is `no`, so a structure hit
   counts only when a name agrees (reason `structure_only_needs_name`). Often the
   reason is a fold shared by enzymes of other pathways, as for ring-hydroxylating
   (Rieske) oxygenases, non-ribosomal peptide synthetases and HAD-family phosphatases;
   each gene's `notes` give its reason. A pathway name alone, such as "alcaligin
   biosynthesis protein", does not count for a required gene: a cluster's regulators
   and transporters carry it too.
3. **A pseudogene never counts.** Bakta writes `(pseudo)` before the product of a gene
   it found broken, by a frameshift, an internal stop codon or a lost start codon.
   Where its evidence would count, its row reads `pseudogene`, and a set that lacks
   only such genes is `partial`.
4. **One protein, one gene per rule.** A protein counts for two genes of one rule only
   when one source names both, as a fused PqqC/D, DfoC or UreAB does; Bakta's product,
   symbol, KO and EC count as one source. A fusion named after one part only shows the
   other as missing. Otherwise the protein counts once, for the first gene that
   applies:
    - the gene of its structure hit, when a structure hit counts alone for that gene
      (`s_alone` is `yes`);
    - the gene of Bakta's name;
    - the gene of Baktfold's name, last: Baktfold names a protein from its predicted
      structure, while Bakta's name comes from sequence similarity, which usually
      tells related enzymes apart.

    A protein whose structure hit is PcaG, the alpha chain of protocatechuate
    3,4-dioxygenase, but which Bakta names as the beta chain, counts for `pcaG` only:
    one protein cannot fill both places of the pair. The other gene's row reads
    `counted_for_other_gene`, and the protein counts for it in no rule. A pseudogene
    is ranked the same way, but both of its rows keep the reason `pseudogene`, and
    their `detail` names the gene the broken protein is taken as.

Bakta also marks a gene that runs off the end of its sequence: `(5' truncated)`,
`(3' truncated)` or `(partial)`. BacFlux runs Bakta without `--partial`, so on a
linear contig Bakta keeps only genes whose start and stop codons both lie inside it.
These marks therefore appear only on the circular replicons of the `nanopore` and
`hybrid` modes, most often on the replication gene, such as *dnaA*, that dnaapler put
at the start. A circular replicon has no ends, so such a gene counts as usual, with a
note in `detail`. Names are matched without these prefixes or `(pseudo)`, in Bakta's
product and in Baktfold's.

## How a trait is called

Each rule names one or more alternative sets of required genes: `+` joins the genes of
one set and `|` separates the alternatives, as in `ureA+ureC | ureB+ureC`. Sets are
tried in the order listed: the first set that meets the `called` condition gives the
call; failing that, the first one that qualifies as partial, and is not withheld (key
genes, below), is reported.

| Decision | When |
|---|---|
| `called` | every gene of the set counts, and the genes that must lie together do |
| `partial` | the set is complete but its genes are not together where the rule asks; or it is complete except for pseudogenes, which holds for a one-gene set too; or, with two or more genes, at least half of them count (rounded up) |
| `not called` | everything else, including a partial withheld by `key_genes` |

**Co-localisation.** A rule can ask that genes lie together: `all` means the whole set,
otherwise groups separated by `;`, such as
`nahAa+nahAb+nahAc+nahAd; nahB+nahC+nahD+nahE`. A group holds when its genes form a
chain on one contig, each at most `traits.max_gene_gap` positions from the previous
one (default 10: up to nine other genes between). Positions count protein-coding genes
in Bakta's order. A group applies only to a set that holds all its genes, so one rule
can hold routes whose genes are arranged differently. Every contig is treated as
linear: a chain across the start of a circular contig is not joined.

**Key genes.** A rule can name `key_genes`: the genes that set its pathway apart from
look-alikes, such as `desA` for desferrioxamine, since an aerobactin cluster also
matches `desB`, `desC` and `desD` but has no `desA`. A set that lacks only pseudogenes
is always reported. Any other partial is reported only when a key gene is in the
genome, counted or as a pseudogene; otherwise it is withheld, and if no other partial
stands, the trait is `not called` with the reason `partial_withheld`. The idea comes
from MetaCyc's key reactions, which tell related pathway variants apart (Karp et al.
2011).

**Supporting genes** are reported when they count, and no set requires them. Like any
gene, a supporting gene can keep its protein from counting for a required gene of the
same rule ([When a gene counts](#when-a-gene-counts), point 4). An entry
`antiSMASH:<type>`, such as `antiSMASH:NRP-metallophore`, is reported when an
antiSMASH region of exactly that type overlaps, on the same contig, a protein counted
for a required gene.

## The two tables

The rules live in two tables beside the script that applies them, each headed by a
one-line `#` comment that says what it is (skip it when you load a table):

- **`workflow/scripts/40_annotation/trait_rules_v1.tsv`**: 49 rules, one per trait, in
  the order of the seven groups below.
- **`workflow/scripts/40_annotation/trait_genes_v1.tsv`**: the gene dictionary, 274
  entries, one for each gene the rules name, saying how it is recognised.

Each rule follows the primary papers in its `reference` column: the smallest set of
genes that defines the pathway, its alternative routes, and which genes lie together.
Each gene entry takes its structure evidence from EcoFoldDB's curated catalogue, its
KO and EC numbers from KEGG and ExPASy ENZYME only where they belong to that gene
alone, and its names from UniProt and RefSeq, written to match how Bakta names the
protein. Vetoes cover known look-alikes, and `s_alone` is `yes` only where the fold
belongs to that function alone; the `notes` give the reason for each entry.

Editing either table makes the calls again, with the summary and the MultiQC report,
but does not repeat the 3Di step ([Cost and databases](#cost-and-databases)). A fresh
clone of BacFlux does the same in an output folder made before, because it gives the
tables new file dates. Before
calling anything, the script checks both tables: the columns, unique `trait_id`s,
exactly one entry for every gene a rule names and no unused entry, patterns that
compile, `yes`/`no` values, and EcoFoldDB entries present in the catalogue. A failed
check stops the job with a message naming the row.

| Group | Rules | `trait_id` |
|---|--:|---|
| Nitrogen | 3 | `n2_fixation`, `nodulation`, `urea_hydrolysis` |
| Phosphorus | 4 | `mineral_phosphate_solubilisation`, `phytate_degradation`, `phosphate_ester_mineralisation`, `phosphonate_cp_lyase` |
| Iron: siderophores | 11 | `siderophore_pyoverdine`, `siderophore_enterobactin`, `siderophore_bacillibactin`, `siderophore_desferrioxamine`, `siderophore_pyochelin`, `siderophore_achromobactin`, `siderophore_alcaligin`, `siderophore_putrebactin`, `siderophore_rhizobactin_1021`, `siderophore_mycobactin`, `siderophore_pseudomonine` |
| Plant hormones and signals | 9 | `acc_deaminase`, `iaa_iam`, `iaa_ipya`, `iaa_nitrilase`, `gibberellin`, `cytokinin`, `gaba`, `putrescine`, `spermidine` |
| Stress tolerance | 3 | `ectoine`, `glycine_betaine`, `trehalose` |
| Aromatic compounds | 15 | `benzoate_to_catechol`, `catechol_ortho_cleavage`, `catechol_meta_cleavage`, `protocatechuate_ortho_cleavage`, `protocatechuate_meta_cleavage`, `phenol_to_catechol`, `naphthalene_upper_pathway`, `biphenyl_upper_pathway`, `toluene_monooxygenase_aerobic`, `toluene_anaerobic`, `xylene_methyl_monooxygenase`, `phenylacetate_degradation`, `phthalate_degradation`, `terephthalate_degradation`, `benzoyl_coa_anaerobic` |
| Nitrogen cycle (context) | 4 | `nitrate_reduction`, `nitrite_reduction_to_no`, `no_reduction`, `n2o_reduction` |

The rules table, `trait_rules_v1.tsv`:

| Column | Contents |
|---|---|
| `trait_id`, `trait` | the identifier every output uses, and the trait in words |
| `application` | one or more of `biofertiliser`, `biocontrol`, `biostimulant`, `stress tolerance` and `bioremediation`, separated by commas; `context` marks a call that describes the strain's metabolism rather than supporting a product claim ([What a call does not mean](#what-a-call-does-not-mean)) |
| `required` | the alternative gene sets |
| `supporting` | genes and antiSMASH region types that are reported but not required, separated by `|` |
| `colocalised` | empty, `all`, or groups separated by `;` |
| `key_genes` | optional: the genes that set the pathway apart from look-alikes, separated by `|` ([Key genes](#how-a-trait-is-called)) |
| `caveat` | what a call does not show; copied into the profile |
| `assay` | the laboratory test that would confirm the trait |
| `reference` | the primary papers behind the rule |

The gene dictionary, `trait_genes_v1.tsv`:

| Column | Contents |
|---|---|
| `gene`, `protein` | the name the rules use, and the protein in words |
| `ecofolddb` | the catalogue entries that give structure evidence, named by the catalogue's `Gene` and `Pathway/Activity` joined by `@`; one name covers every entry with that `Gene` and `Pathway/Activity`, such as `nifH@N2 fixation`, which covers three structures; empty for a gene recognised by name only |
| `ko`, `ec` | the KO and EC numbers that count as a name match |
| `bakta_gene` | the gene symbols that count, matched exactly |
| `symbol_needs_product` | `yes`: a symbol counts only together with the product pattern or a KO |
| `product_include` | a pattern searched for in the product, ignoring case |
| `product_exclude` | a pattern that vetoes the protein for this gene |
| `s_alone` | `yes`: a structure hit counts without a name |
| `notes` | the biology behind the entry |

Lists inside a cell are separated by `;`, and the patterns are Python regular
expressions.

## Reading the output

`NA` marks a missing value in all three tables.

### `{sample}_trait_profile.tsv`

One row per rule, in table order.

| Column | Contents |
|---|---|
| `sample`, `trait_id`, `trait`, `application` | which genome, which trait |
| `decision` | `called`, `partial` or `not called` |
| `genes` | the counted genes of the called or partial set, each as `gene:locus_tag(evidence)`, with several proteins of one gene separated by commas and genes by semicolons; `NA` for a trait that is not called, and for a partial in which no gene counts |
| `missing` | for a partial, the genes of that set that did not count, separated by a comma and a space: each by its name, or, when only broken copies were found, as `gene:locus_tag(pseudogene)`, several copies joined by a comma alone, as in `nifD:SAMPLE1_01234(pseudogene),SAMPLE1_01240(pseudogene), nifE`. A gene whose only broken copy went to another gene ([When a gene counts](#when-a-gene-counts), point 4) is written by its bare name. `NA` otherwise |
| `supporting` | the supporting genes that counted, and the antiSMASH region types found with their position |
| `borderline` | `yes` when a protein counted for the call carries a borderline structure hit for that gene, `no` when none does; `NA` for a trait that is not called |
| `structure_only` | the genes of the call counted from structure alone, with no name agreeing; `NA` when there are none, and for a trait that is not called |
| `caveat` | what the call does not show, from the rules table |

### `{sample}_trait_profile_audit.tsv`

One row per rule, with the columns `sample`, `trait_id`, `decision`, `reason` and
`detail`. `reason` is a fixed word you can filter on; `detail` names the set and the
genes behind it. For a partial, and for a trait that is not called, `detail` also lists
the required genes of the rule's other sets that were found only as pseudogenes. A gene
whose broken protein went to another gene ([When a gene counts](#when-a-gene-counts),
point 4) is treated as not found, unless another protein counts for it or another
broken protein is taken as it, and is written by its name; its row in the gene table
shows the broken protein.

```bash
cut -f3,4 {sample}_trait_profile_audit.tsv | sort | uniq -c | sort -rn
```

| `reason` | `decision` | Means |
|---|---|---|
| `set_complete` | `called` | one alternative set is complete, and its co-localisation holds |
| `set_complete_not_colocalised` | `partial` | a set is complete, but its genes do not lie together where the rule asks |
| `set_complete_but_pseudogene` | `partial` | each gene of a set counts or is found as a pseudogene, and at least one only as a pseudogene; a one-gene set too. `detail` names each such gene with its locus tag |
| `half_of_set_found` | `partial` | at least half of a set of two or more genes counts |
| `partial_withheld` | `not called` | a set qualifies as partial (other than one lacking only pseudogenes), but no key gene is in the genome, counted or as a pseudogene. `detail` reads `no key gene (desA) found (...)`, with the withheld partial in brackets |
| `too_few_genes` | `not called` | some genes count, but fewer than half of the closest set; `detail` says which are missing |
| `no_gene_found` | `not called` | no required gene counts; `detail` names any required gene found only as a pseudogene |

### `{sample}_trait_genes.tsv`

One row per protein and gene with any evidence, sorted by contig and position: the
place to check a call.

| Column | Contents |
|---|---|
| `sample`, `locus_tag`, `contig`, `start`, `stop`, `strand` | the protein, from Bakta |
| `cds_index` | its position in the contig's gene order, the number co-localisation is measured in |
| `gene`, `traits` | the gene, and the `trait_id` of each rule that uses it |
| `counted` | `yes` or `no` |
| `evidence` | `S`, `Q` and `F` joined by `+`: structure, Bakta name, Baktfold name |
| `bakta_gene`, `bakta_product`, `ko`, `ec` | Bakta's annotation of the protein |
| `baktfold_product` | Baktfold's name for a protein Bakta called hypothetical |
| `ecofolddb_target`, `ecofolddb_gene`, `ecofolddb_pathway` | the best structure hit and its catalogue entry |
| `evalue`, `qcov`, `tcov` | that hit's e-value, as Foldseek printed it, and its two coverages |
| `borderline` | on a row whose evidence includes `S`: why that hit is borderline, or `no`; `NA` on other rows |
| `decision`, `reason`, `detail` | `counted`, `not counted` or `not searched`; a fixed reason; free text |

| `reason` | Means |
|---|---|
| `name_and_structure_agree` | a name match and a kept structure hit to one of the gene's `ecofolddb` entries; counted |
| `name_only` | a name match and no structure evidence for this gene, such as no hit, a failed hit or a kept hit to another entry (`detail` says which); counted |
| `structure_counts_alone` | a kept structure hit and no name, for a gene whose `s_alone` is `yes`; counted |
| `pseudogene` | evidence that would count by one of the three reasons above, on a protein Bakta marks as a pseudogene; not counted. When a different source names the protein for another gene of the same rule, `detail` says which of the two it is taken as ([When a gene counts](#when-a-gene-counts), point 4) |
| `counted_for_other_gene` | evidence that would count, on a protein that counts instead for another gene of the same rule ([When a gene counts](#when-a-gene-counts), point 4); not counted, and `detail` names both genes and why |
| `structure_only_needs_name` | a kept structure hit and no name, for a gene whose `s_alone` is `no`; not counted |
| `vetoed_by_name` | the product matched the gene's `product_exclude`; not counted by any evidence |
| `catalogue_gene_unused` | a kept structure hit to a catalogue entry that no gene of the gene dictionary uses; listed without a gene |
| `not_searched_too_long` | a protein longer than 4,000 amino acids, not searched by structure; listed without a gene, and its names are still checked. `detail` gives its length |

A protein that matched by name but whose best structure hit is a catalogue entry that
failed the filter keeps its name row: the failed hit fills the structure columns, and
`detail` says why it failed.

### `09.report/traits_summary.tsv`

One row per sample and one column per `trait_id`, in table order, each cell `called`,
`partial` or `not called`.

The MultiQC report shows the same calls in its Trait profile section, after Bakta,
from the plot files in `09.report/traits_summary_mqc/`:

- **Call counts**: for each strain, how many traits are called, partial and not called.
  The counts help to find strains in a batch; they do not rank strains by benefit.
- **One heatmap per application label**, strains in rows and traits in columns. Each
  trait appears once, under its exact label, so `biofertiliser, biocontrol` has its own
  panel. The axis shows each trait's `trait_id`, and the panel lists every ID with the
  trait's name.
- **References and assays**: for each trait, the papers behind its rule, a suggested
  laboratory assay and its caveat, from the catalogue.

Dark blue is **called**, light blue **partial** and grey **not called**. The heatmaps
draw them from the numbers 2, 1 and 0, which stand for the three categories, not for a
score of confidence or activity. A separate **References** section of the report lists
the papers of the tools and of the trait methods.

### `{sample}_trait_provenance.json`

Written last by `trait_profile.py`, after the three tables. It records the SHA-256
fingerprint (a code computed from every byte of a file) of the two catalogue tables, of
EcoFoldDB's descriptions and of the three tables, each beside its file name, with the
calling parameters. The MultiQC report and the traits workbook quote the catalogue's
caveats, assays and references beside the calls, so rule `traits_summary` and rule
`traits_workbook` compare these fingerprints with the files they read first: a catalogue
edited after the calls were made, or an edited table, stops them, with rule
`trait_profile` to be run again. Only contents are compared, so an output folder that has
been moved still passes.

## What a call does not mean

!!! warning "A call is genomic evidence, not a phenotype"

    A call says that the genes for a trait are present. It does not show that they
    are expressed, that the proteins are active, or that the strain helps a plant.
    The `caveat` column says what each call leaves open, and the `assay` column of the
    rules table names the test that would confirm it.

**`not called` is not the same as absent.** A gene can be missed: its protein may carry
no name the patterns recognise while the gene needs a name to count, the protein may be
too long to be searched by structure, Bakta's gene finder may not predict a short gene
at all, or the assembly may break the gene. A set that
the rule wants together, split across contigs of a draft assembly, is reported as
`partial` with the reason `set_complete_not_colocalised`, not as called. Read the audit
before reporting a trait as missing.

**`partial` is a lead, not half a trait.** Check `missing` in the profile, and the
gene table for the proteins behind it.

**A call can rest on a single gene.** Some rules accept one marker gene as a set, such
as `acdS` for ACC deaminase or `nosZ` for nitrous oxide reduction. When that gene is
listed in `structure_only`, the call rests on a structure hit alone, with no name
confirming it; only a gene whose `s_alone` is `yes`, such as `nosZ`, can count this
way, so `acdS` always needs a name.

**`context` describes the strain, not a product claim.** The `application` column
gives each trait the use it points to. `context` marks a trait whose call describes
the strain's metabolism but does not, on its own, support a biofertiliser,
biocontrol, biostimulant or bioremediation claim. Three kinds of traits carry it:

- the four denitrification steps, `nitrate_reduction` to `n2o_reduction`, `context`
  only: nitrate reduction is the first step both of denitrification, which loses
  nitrogen from the soil, and of nitrate ammonification, which keeps it as ammonium;
- `siderophore_mycobactin`, `context` only: mycobactin is the siderophore of
  mycobacteria and their relatives, so a call says which genus you are holding
  rather than that the strain supplies a plant with iron;
- the central aromatic routes, catechol and protocatechuate ring cleavage and
  phenylacetate degradation, `bioremediation, context`: many aromatic pollutants end
  up in them, but they are also the ordinary exit for plant phenolics in many soil
  bacteria. A pollutant claim rests on an upstream route, such as
  `naphthalene_upper_pathway` or `toluene_monooxygenase_aerobic`.

The label never changes a decision.

**`borderline: yes`** marks a call in which a counted protein carries a structure hit
close to a cut-off or to a rival. Check that hit in the gene table before quoting the
call.

## Cost and databases

The 3Di step dominates: about 1 h 50 min per genome on 16 CPU cores, and it grows with
the number and length of the proteins. It turns every protein into a 3Di sequence,
not only the hypothetical ones Baktfold looks at, so the trait profile and Baktfold do
not share this step. The Foldseek search then takes about 2 minutes per genome, and
the calls seconds.

The layer needs two databases:

- **Baktfold's database**, for the ProstT5 model inside it, even when `baktfold.run` is
  false: about 22 GB, from `directories.baktfold_db`, or downloaded whole by
  `baktfold install` on the first run.
- **EcoFoldDB v2.0**: a 0.57 GB archive from its Zenodo record, checked against
  `links.ecofolddb_md5` and unpacked to about 0.8 GB in
  `04.annotation/traits/databases/`, or a copy you already hold, from
  `directories.ecofolddb_db`.

!!! warning "Use the Zenodo copy of EcoFoldDB"

    The copy of EcoFoldDB in the authors' GitHub repository names some catalogue
    entries differently from the Zenodo record, and the gene dictionary needs the
    Zenodo names. A copy given in `directories.ecofolddb_db` is checked before any job
    starts: when an entry the gene dictionary uses is missing from its
    `EcoFoldDB_descriptions.txt`, the run stops and names it.

The structure files in `04.annotation/traits/{sample}/structure/` are kept, so the 3Di
step runs once per genome. The Foldseek hit table, about 170 MB per genome, is
temporary: making the calls again, after a table edit for example, repeats only the
2-minute search. Fetching EcoFoldDB again repeats the search and the calls; fetching
or linking Baktfold's database again repeats the 3Di step too. The archive clean-up
(`workflow/scripts/clean_workdir.sh`) removes the structure files and both databases,
and a later run in a cleaned folder repeats the 3Di step. Switching `baktfold.run` on
in a folder that already holds trait profiles makes the calls again, with Baktfold's
names.

## Resources

| Rule | Threads | Memory |
|---|---:|---|
| `trait_structure_query` | capped at 24 | 16 GB requested from a cluster scheduler, or `ram_gb` if lower |
| `trait_structure_search` | capped at 24 | none declared |
| `trait_profile` | 1 | none declared |
| `traits_summary` | 1 | none declared |

"Capped at 24" has the meaning given under
[Annotation](annotation.md#resources). The 3Di step runs ProstT5 on the CPU; like
Baktfold's own rule, it is held to its threads and asks no internet service for the
model, so it runs on a compute node without internet. `trait_db` is a local rule: on a
cluster it runs on the computer that runs Snakemake, which is the one with internet. It
takes `wget`, `tar` and the `gzip` that `tar` calls from `workflow/envs/download.yaml`;
`md5sum`, `awk`, `mkdir`, `basename`, `rm` and `ls` come from the computer, as part of
every Linux base system. The 3Di and search rules run in `workflow/envs/baktfold.yaml`,
which carries ProstT5 and Foldseek, and the two scripts in `workflow/envs/python.yaml`.

## What to check afterwards

- `09.report/traits_summary.tsv`, or the trait heatmaps in
  `09.report/multiqc_report.html`: every sample at a glance.
- `{sample}_trait_profile_audit.tsv`: the reason behind each decision.
- `10.workbooks/{sample}_traits.xlsx`: the calls with their genes, caveats, suggested
  assays and references, and the plasmid and phage context of every gene, in one file
  per strain ([Per-strain workbooks](workbooks.md#the-traits-workbook)).
- `logs/trait_db.log`: the download, and whether the MD5 matched. Nothing is
  extracted when it did not. With `directories.ecofolddb_db` set,
  `logs/trait_db_local.log` takes its place.
- `logs/trait_structure_query_{sample}.log`, `logs/trait_structure_search_{sample}.log`
  and `logs/trait_profile_{sample}.log`: one per sample per step;
  `logs/traits_summary.log` for the summary. A `WARNING` in the trait_profile log
  that the hit table holds no hits means that every trait of that sample is decided
  from names alone, with no structure evidence.

## Licensing and citation

The EcoFoldDB database is published on Zenodo under CC BY 4.0. The authors' own
annotation script, on GitHub, states no licence, and BacFlux does not use it: the
filter on this page is BacFlux's own code, written from the published method. See
[Licensing](../about/licensing.md). Cite EcoFoldDB together with Baktfold, ProstT5 and
Foldseek, which run in this layer ([Citation](../about/citation.md)); the sources of
each rule are in the `reference` column of the rules table.

## References

- Ghaly, T. M., Rajabal, V., Russell, D., Colombi, E., & Tetu, S. G. (2025).
  EcoFoldDB: protein structure-guided functional profiling of ecologically relevant
  microbial traits at the metagenome scale. *Environmental Microbiology*
  27(9):e70178. <https://doi.org/10.1111/1462-2920.70178> — source of the database and
  of the method this layer follows: the structure filter, the 4,000-residue limit and
  the `-s 8` search, as in the authors' own annotation pipeline. The database:
  EcoFoldDB v2.0, <https://doi.org/10.5281/zenodo.16908168>.
- Heinzinger, M. *et al.* (2024). Bilingual language model for protein sequence and
  structure. *NAR Genomics and Bioinformatics* 6(4):lqae150.
  <https://doi.org/10.1093/nargab/lqae150> — ProstT5, which predicts the 3Di
  sequences.
- van Kempen, M. *et al.* (2024). Fast and accurate protein structure search with
  Foldseek. *Nature Biotechnology* 42(2):243–246.
  <https://doi.org/10.1038/s41587-023-01773-0>
- Bouras, G. *et al.* (2026). Baktfold: sensitive protein functional annotation across
  the microbial tree of life using structural information. *bioRxiv* (preprint).
  <https://doi.org/10.64898/2026.03.31.715528> — the tool that runs the ProstT5 step;
  its database supplies the model.
- Karp, P. D., Latendresse, M., & Caspi, R. (2011). The Pathway Tools pathway
  prediction algorithm. *Standards in Genomic Sciences* 5(3):424–429.
  <https://doi.org/10.4056/sigs.1794338> — MetaCyc's key reactions, which tell
  related pathway variants apart: the idea behind `key_genes`.

The co-localisation distance, the borderline flag and each rule's key genes are
conventions of this layer, not values taken from these sources: `key_genes` borrows
only the idea of MetaCyc's key reactions.
