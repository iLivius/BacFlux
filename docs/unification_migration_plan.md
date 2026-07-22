# BacFlux v2.0.0 — Unification Migration Plan

**Status:** design agreed — D1, D3, D4 signed off 2026-07-20; D2, D5–D7 accepted as
recommended. No code written yet.
**Goal:** merge four sibling workflows (BacFlux short-read, FastaFlux contigs, BacFluxL
long-read, BacFluxL+ hybrid) into a **single BacFlux repository** with one Snakefile
that dispatches to mode-specific and shared rule modules — the pattern already proven
in MetaFlux. Released as **BacFlux v2.0.0** on the existing repo/DOI lineage.

This document is the reference for the migration. It records the reasoning so it does
not have to be re-derived. Read §4 (decisions needing sign-off) before starting.

---

## 1. Why (one paragraph)

The four workflows share roughly two-thirds of their rules, an almost-identical conda
env set (same version pins), byte-identical helper scripts, and structurally identical
READMEs. Every change today lands 3–4 times by hand and the copies drift. The upcoming
mobilome/AMR-mobility module (see `mobilome_module_SPEC.md`) would otherwise be copied
into four repos and maintained in parallel forever. Unifying first means the mobilome
module — and every future feature — is written once.

---

## 2. Current state (verified from the code, 2026-07)

| Workflow | Mode name | Input | Tech-specific front end | Snakefile lines | Rules |
|---|---|---|---|---|---|
| BacFlux | `illumina` | Illumina PE reads | phix removal (bowtie2), fastp, SPAdes | 1204 | 31 |
| FastaFlux | `contigs` | pre-assembled FASTA | *none* (starts at contig filtering) | 916 | 22 |
| BacFluxL | `nanopore` | ONT reads | filtlong, NanoPlot, Flye, Medaka, dnaapler | 1265 | 27 |
| BacFluxL+ | `hybrid` | Illumina + ONT | all of the above + Polypolish + Snippy | 1845 | 40 |

**Shared downstream core — ~19 rules present in all four, judged identical or
path-only-different:** `blast_contigs`, `blob_json`, `blob_table`, `select_contigs`,
`completeness_and_contamination`, `taxonomic_assignment`, `annotation`
(`accurate_annotation` in two repos), `functional_annotation`, the two antiSMASH rules,
the two dbCAN rules, `amr_contigs`, `AMR_summary`, `plasmid_search`, `viral_db`,
`viral_identification`, `multiqc`.

**Envs:** version pins are already near-identical across repos (abricate=1.2.0,
antismash=8.0.4, bakta=1.12.0, gtdbtk=2.6.1, …). Only difference found: fastp 1.0.1
(BacFlux) vs 1.1.0 (BacFluxL+). Consolidate to one set, one pin each.

**Scripts:** `select_contigs_by_taxonomy.py` and `clean_workdir.sh` are already shared.

---

## 3. Target architecture

```
config/
    config.yaml                 # single config; `mode:` selects the pipeline
    config_custom.yaml          # keep the dry-run convention
workflow/
    Snakefile                   # thin (~50 lines): banner, configfile, include common,
                                #   glob-include the active mode's front end + all shared,
                                #   `rule all: input: all_targets()`
    rules/
        shared/
            00_common.smk       # mode dispatch, path/stage map, sample discovery,
                                #   all inline helpers (deduplicated), the final-assembly
                                #   abstraction, DATABASES list, awk command constants
            10_decontam.smk     # blast_contigs, blob_json, blob_table, select_contigs
            20_qc.smk           # quast, checkm  (+ qualimap wrapper)
            30_taxonomy.smk     # gtdbtk
            40_annotation.smk   # bakta, eggnog, antismash, dbcan
            50_amr.smk          # abricate (+ CARD read-mapping, gated to read-based modes)
            60_plasmid.smk      # platon (+ geNomad plasmid calls, see D9)
            70_phage.smk        # geNomad (default) or virsorter2 (selectable) + checkv
            80_mobilome.smk     # NEW — added last, per mobilome_module_SPEC.md
            90_report.smk       # multiqc (mode-aware inputs)
        illumina/               # phix, fastp, SPAdes, contig mapping
        nanopore/               # filtlong, NanoPlot, Flye, Medaka, dnaapler, finalize
        hybrid/                 # illumina + nanopore front ends + Polypolish + Snippy
        contigs/                # FastaFlux front end (header-aware contig filter, self-map)
    envs/                       # ONE copy of each env
    scripts/                    # unchanged (already shared) + scripts/mobilome/ later
README.md                       # ONE, with a section per mode
CITATION.cff                    # NEW — machine-readable citation, carries the concept DOI
docs/
    mobilome_module_SPEC.md
    unification_migration_plan.md   # this file
```

The Snakefile mirrors MetaFlux almost exactly: resolve `mode` in `00_common.smk`, then
`for f in sorted(glob(rules/{MODE}/*.smk)): include: f`, then include every `shared/*.smk`,
with `90_report.smk` last so it sees all upstream outputs.

---

## 4. Design decisions — SIGN-OFF NEEDED

These are the choices that shape everything. My recommendation is given for each. The
three that needed the author's call (D1, D3, D4) were **confirmed on 2026-07-20**, all on
the recommended option; the rest are accepted as recommended unless revisited.

### D1. Unified directory numbering — collapse the front end ✓ **DECIDED**

**Problem:** today the numeric prefixes encode stage order, but each mode has a different
number of front-end stages, so the *shared* stages land at different numbers — taxonomy is
`04.taxonomy` in short/long but `09.taxonomy` in hybrid, `03.taxonomy` in contigs. This is
the single biggest reason the "shared" rules differ (only by path).

**Fix:** group *all* tech-specific front-end work under two fixed parents, so the shared
downstream always lands at the same number in every mode:

```
01.reads/{sample}/          # all read QC + filtering (sub-dirs: illumina/, ont/)  — empty in contigs mode
02.assembly/{sample}/       # all assembly, polishing, reorientation, decontamination, assembly QC
                            #   sub-dirs: spades/ flye/ medaka/ fix_start/ polypolish/ snps/
                            #             contaminants/ eval/  → terminates in contigs_final.fasta
03.taxonomy/
04.annotation/              # bakta/ eggnog/ antismash/ dbcan/
05.amr/                     # abricate/ mapping/
06.plasmids/
07.phages/
08.mobilome/                # only when config.mobilome.run
09.report/
```

**Consequence:** this changes the output layout of **all four** current tools. That is a
breaking change — but a *major version* (v2.0.0) is exactly where a layout change belongs,
and it must be called out in the changelog and migration notes. Contigs mode simply has an
empty `01.reads/` and a light `02.assembly/` (filter + decontam only).

Rationale: it makes the ~15 clean-shared rules literally path-identical across modes, which
is what lets them live in `shared/` unmodified.

### D2. Single final-assembly hand-off **[recommended]**

Every mode's front end ends by producing one canonical file:
`02.assembly/{sample}/contigs_final.fasta`. Everything in `03.*`–`08.*` consumes **only**
that file, via a single variable defined in `00_common.smk`. This already matches how
BacFluxL works (`assembly_final.fasta`); we standardize the name and path across modes.
Today the hand-off is named differently per mode (`contigs_sel.fasta` short/contigs,
`assembly_final.fasta` long, `fixed_contigs.fasta` hybrid) — unify to `contigs_final.fasta`.

### D3. Decontamination position ✓ **DECIDED: keep current per-mode order (option A)**

The decontam *rules* (blast/blob/select) are identical and go in `shared/10_decontam.smk`.
But their *position in the DAG* differs and cannot be silently unified:

- `illumina` / `contigs`: decontam is the **last** assembly step (`contigs_sel` = final).
- `nanopore`: assemble → reorient → **decontam** → Medaka polish → finalize.
- `hybrid`: decontam runs on the **Illumina** draft (for the decontam decisions + Snippy
  reference); the *delivered* genome is the ONT+Polypolish assembly.

Two options:
- **(A, recommended for v2.0.0) Keep each mode's current ordering.** The decontam rules are
  shared *definitions*, wired per mode via input functions. Preserves exact current behavior
  → lowest regression risk. Slightly less uniform.
- **(B) Standardize** so decontam always runs on the finished draft. Cleaner, but changes
  hybrid behavior (which currently decontaminates the Illumina assembly, not the ONT final)
  → must be revalidated biologically. Defer to a later minor release.

### D4. VirSorter2 input inconsistency ✓ **DECIDED: decontaminated (`contigs_final`) in all modes**

`viral_identification` consumes different contigs per mode:
- BacFlux (`illumina`): `contigs_filt.fasta` — **pre-decontamination**
- FastaFlux (`contigs`): `contigs_sel.fasta` — decontaminated
- BacFluxL / L+ : `assembly_final` / `fixed_contigs` — final

The short-read version scanning the *pre-decontam* contigs is almost certainly an oversight.
Unifying forces one choice. **Recommendation:** run phage detection on `contigs_final.fasta`
(decontaminated) in all modes, matching three of the four. Flag in the changelog that
`illumina` phage output may change slightly vs v1.x. Confirm you agree this is the intended
behavior and not a deliberate choice to keep phages from discarded contigs.

### D5. Rule-name vocabulary **[recommended]**

Standardize the names that currently differ for the same step:
- `annotation` vs `accurate_annotation` → **`annotation`** (drop the FastaFlux/L+ variant name).
- `genome_assembly` / `assembly` / `illumina_assembly` / `ONT_assembly` → per-mode assembly
  rules can keep descriptive names inside their front-end module; the *output* is always
  `02.assembly/{sample}/contigs_final.fasta`, so downstream never sees the difference.

### D6. Config schema — `mode` + namespaced params **[recommended, see §5]**

### D7. Small cleanups to fold in during migration **[recommended]**

- The 8-database AMR list is hard-coded in **both** `amr_contigs` and `AMR_summary` in every
  repo. Define `DATABASES` once in `00_common.smk`; drive both with `expand()`.
- Standardize the inconsistent resource keys (`cpus`, `cpus_p`, `threads`, `java_mem`,
  ad-hoc `min(CPUS,N)` caps) into one convention.
- The CARD read-mapping leg (`download_amr_db`, `map_amr_db`) exists in `illumina` and
  `hybrid` only (it needs reads); `nanopore` and `contigs` skip it. Gate it on
  "mode produces short reads", not on mode name, so it's explicit.

### D8. Phage caller: VirSorter2 **default**, geNomad **opt-in** ✓ **REVISED 2026-07-22 (licensing)**

**Correction (2026-07-22):** the 2026-07-21 form of this decision made geNomad the DEFAULT on
a false licensing premise. geNomad is NOT BSD-4-Clause — verified from its raw LICENSE it is a
Berkeley Lab **ACADEMIC / NON-COMMERCIAL-USE-ONLY** licence ("User must be an accredited
academic institution"; the bioconda `BSD-4-Clause` recipe tag is wrong). BacFlux is MIT and
(mobilome spec §11) must not force a non-commercial restriction on downstream users, so geNomad
CANNOT be the default/mandatory caller. **Revised (Option A, user-approved):** VirSorter2
(GPLv2, commercial-use-OK) is the DEFAULT; geNomad is OPT-IN (`config.phage.caller: genomad`),
treated exactly like a licence-encumbered database. Verified the default path is fully
permissive: VirSorter2 GPLv2 + Platon GPLv3 + CheckV (LBNL **permissive** BSD, not the
non-commercial variant — checked directly, since CheckV is also Berkeley Lab). geNomad's
technical merits still hold (actively maintained, 97.3% vs 94.7% precision, viruses+plasmids in
one run); the licence, not the science, sets its role. CheckV runs unconditionally on whichever
caller ran.

- `shared/70_phage.smk`: default `config.phage.caller: virsorter2`. geNomad's rules are gated
  behind `if PHAGE_CALLER == "genomad":` — on the default run geNomad's env is never built and
  it never runs (no non-commercial tool is forced on a commercial user). Opting in also enables
  the D9 concordance (geNomad does viruses+plasmids in one run). VirSorter2 (with the
  `virsorter_deps_env` packaging fix baked in) stays as the permissive default. CheckV runs
  downstream of either.
- geNomad's `end-to-end` run also emits plasmid calls in the same invocation — see D9,
  this is the same tool run once, consumed by both `70_phage.smk` (virus role) and
  `60_plasmid.smk` (plasmid second-opinion role).

### D9. Plasmid calls: add Platon + geNomad concordance, keep the BLAST-text check as one more column ✓ **DECIDED 2026-07-21**

BacFlux's current `plasmid_search` rule (Snakefile:985) appends a "verification" pass on
top of Platon's calls: for each Platon-flagged plasmid contig, it greps the first line
matching that contig's ID in the general-purpose contamination-screening `blastout` file
(from `blast_contigs`, BLASTed against nt for a completely different purpose) and checks
whether the subject title contains the literal substring "plasmid".

**Decision: keep it — add to it, don't replace it.** It's non-destructive (annotates,
never discards a Platon call), the user has found it practically useful, and it has a real
(if partial) protective mechanism: it's a reasonable check against *borderline/noisy RDS
calls* — cases where Platon's score is marginal and the contig's overall composition and
best BLAST hit both genuinely look chromosomal, so the check correctly flags it as
suspect.

What it can't catch, and is worth being explicit about: **when a Platon false positive is**
**caused by a mobile genetic element** (IS/transposon/prophage remnant) that a genuinely
chromosomal contig carries — the same element that biases Platon's protein-family-based
RDS score toward "plasmid" can *independently* cause the BLAST top hit to land on a
plasmid-titled entry, because that same element also sits on real deposited plasmids. In
that specific case the two signals aren't actually independent (both are driven by the
same confound), so their agreement is weaker evidence than it looks — and this is exactly
the failure mode most relevant to BacFlux, since IS/Tn-rich contigs are precisely what the
mobilome module spends its effort on. Also worth knowing: Platon's own "Plasmid database
hits" column already does a more rigorous version of the same idea internally (BLAST
against a curated RefSeq-plasmid DB, not generic nt) as part of its accuracy-mode
heuristic — so treat the two as overlapping, not fully independent, signals.

geNomad is the better complement precisely *because* its classifier is architecturally
different (gene-content/marker-based classification) rather than another nt-BLAST
text-match, so its errors are less likely to share Platon's specific mobile-element
confound. Genuine Platon/geNomad agreement is stronger evidence than Platon agreeing with
a re-purposed screening BLAST.

**Design for `shared/60_plasmid.smk`:** a concordance table per contig — `platon_call,
platon_rds, platon_blast_hit (the existing check, kept, clearly labeled as a supplementary
signal, not authoritative), geNomad_call, geNomad_score, geNomad_fdr, agreement,
confidence`. Confidence tiers driven by **Platon/geNomad agreement** (both agree → high;
one calls, other silent → medium; disagreement → flagged, reported not discarded); the
existing BLAST-text flag rides along as one more visible column for manual review, exactly
as the user already uses it today, but no longer the sole basis for "verified." Audit file
follows the existing `contig_taxonomy_decisions.tsv` convention. Open exit condition: drop
the BLAST-text column if/when Platon+geNomad concordance alone proves sufficient in
practice, or if a stronger method replaces it — not before.

---

## 5. Mode-dispatch config schema (proposed)

```yaml
# Which pipeline to run. Selects rules/{mode}/*.smk and the sample-discovery rule.
mode: illumina            # illumina | nanopore | hybrid | contigs

input:
  illumina_dir: ""        # required for illumina, hybrid  (expects {sample}_R1/_R2)
  nanopore_dir: ""        # required for nanopore, hybrid  (expects {sample}_ont)
  contigs_dir:  ""        # required for contigs           (expects {sample}.fasta)

directories:              # databases — unchanged from v1
  output_dir: ""
  bakta_db: ""
  blast_db: ""
  eggnog_db: ""
  gtdbtk_db: ""
  platon_db: ""

links:                    # download URLs — unchanged
  phix_link: ""
  dbcan_link: ""
  card_link: ""
  checkv_link: ""

resources:
  threads: 16
  ram_gb: 64

parameters:
  nt_version: ""
  decontamination:        # shared across all modes — unchanged structure
    mode: ""              # auto | include | exclude | off
    # include_genera / exclude_genera are INLINE genus lists (null | "A;B" | [A, B]).
    # The other three are FILE PATHS the selector opens on disk (null | /path),
    # NOT mappings — a literal {} would stringify to "{}" and be opened as a file.
    # (Corrected 2026-07-21 from an earlier {} form; see plan D9 note / 00_common
    # _as_decontam_text, which deliberately has no dict branch.)
    include_genera:            # null | "GenusA;GenusB" | [GenusA, GenusB]
    include_genera_by_sample:  # null | PATH to a (sample, genus) TSV
    exclude_genera:            # null | inline genus list to discard
    exclude_genera_file:       # null | PATH to a one-genus-per-line file
    sample_overrides:          # null | PATH to a per-sample override TSV
    discard_no_hit: false
  nanopore:               # only read in nanopore/hybrid
    flye_input_mode: ""
    medaka_model: ""
  hybrid:
    flye_input_mode: ""
    medaka_model: ""

phage:                    # decision D8 — reserved in the schema, read by the Stage-3 phage module
  caller: genomad         # genomad (default) | virsorter2

mobilome:                 # NEW module — default OFF (see mobilome_module_SPEC.md §5.1)
  run: false
  # … (populated when WP-A of the mobilome spec begins)
```

`00_common.smk` validates `mode`, checks the required `input.*_dir` for that mode exists,
and runs only the matching `glob_wildcards` (the mode-specific `_R1` / `_ont` / `.fasta`
discovery that today runs unconditionally at parse time in each repo).

---

## 6. Rule → module mapping

**Sharing tier legend:**
`S` = clean shared, identical rule, no mode awareness ·
`P` = shared rule, mode-parameterized (input function or a small branch) ·
`F` = mode-specific front end.

| Rule (unified name) | Target module | Tier | Notes / reconciliation |
|---|---|---|---|
| sample discovery, helpers, `DATABASES`, awk consts, `contigs_final` var | `shared/00_common.smk` | — | dedup the inline helpers from all four repos |
| `filter_contigs` | `contigs/`, `illumina/`, `hybrid/` | P/F | contigs mode has the SPAdes-vs-NCBI header branch; short/hybrid use the `<500bp/<2x` awk; nanopore doesn't use it |
| `blast_contigs` | `shared/10_decontam.smk` | S | |
| `blob_json`, `blob_table` | `shared/10_decontam.smk` | S | |
| `select_contigs` | `shared/10_decontam.smk` | P | rule identical; DAG position per D3 |
| `genome_assembly_evaluation` (quast) | `shared/20_qc.smk` | S | |
| `completeness_and_contamination` (checkm) | `shared/20_qc.smk` | P | **hybrid classifies TWO genomes/sample** (illumina + ont) — needs a mode branch |
| `map_evaluation` / `map_qc` (qualimap) | `shared/20_qc.smk` | P | tool shared; fed by the tech-specific BAM |
| `taxonomic_assignment` (gtdbtk) | `shared/30_taxonomy.smk` | P | same hybrid dual-genome coupling as checkm |
| `annotation` (bakta) | `shared/40_annotation.smk` | S | includes the shared genus-inference awk |
| `functional_annotation` (eggnog) | `shared/40_annotation.smk` | S | |
| `secondary_metabolites_db`, `secondary_metabolites_analysis` (antismash) | `shared/40_annotation.smk` | S | |
| `cazyme_db_download`, `cazyme_gene_cluster` (dbcan) | `shared/40_annotation.smk` | S | sentinel-file pattern preserved |
| `amr_contigs`, `AMR_summary` (abricate) | `shared/50_amr.smk` | S | dedup `DATABASES` (D7) |
| `download_amr_db`, `map_amr_db` (CARD reads) | `shared/50_amr.smk` | P | gate to read-producing modes (illumina, hybrid) |
| `plasmid_search` (platon) + geNomad plasmid concordance | `shared/60_plasmid.smk` | S | retires the ad hoc BLAST-text heuristic — see **D9** (pending confirmation) |
| `viral_db`, `viral_identification` (geNomad default / virsorter2 selectable + checkv) | `shared/70_phage.smk` | P | resolve input per **D4**; caller choice per **D8** |
| `multiqc` | `shared/90_report.smk` | P | mode-aware input set + hybrid illumina/ont relabeling |
| phix rules, `trim_adapters`, SPAdes, `index_contigs`, `map_contigs` | `illumina/` | F | |
| filtlong, NanoPlot ×2, Flye, Medaka (conditional), dnaapler, `finalize_contigs` | `nanopore/` | F | Medaka rule defined under `if USE_MEDAKA` |
| Illumina + Nanopore front ends, Polypolish, Snippy | `hybrid/` | F | largest front end; owns the dual-genome checkm/gtdbtk wiring |
| header-aware `filter_contigs`, self-map `map_contigs` | `contigs/` | F | |

---

## 7. Staged migration sequence (with validation gates)

Each stage ends at a gate. Do not start the next stage until the gate passes.

**Stage 0 — Baseline & safety net (do first, do not skip).**
- Confirm the last standalone commit of each of the four workflows is tagged and archived
  on Zenodo (they are released; just verify tags exist).
- Run each of the four current tools on one small known isolate; capture the output tree
  and checksums of the key deliverables (assembly, GTDB summary, Bakta TSV, ABRicate
  summary, Platon table, VirSorter/CheckV table). **This is the regression oracle** for
  proving the unified tool reproduces v1 behavior. Store under `docs/validation_baselines/`.

**Stage 1 — Skeleton + `00_common.smk`. ✓ DONE 2026-07-21 (gate passed).**
Thin Snakefile, mode dispatch, config schema, deduplicated helpers, `contigs_final`
abstraction. `rule all` only. **Gate:** `snakemake -n` (Snakemake 9.14.6) parses cleanly,
exit 0, for all four modes. Delivered on branch `release/v2.0.0`, staged under
non-conflicting names so the v1 pipeline stays runnable during migration:
`workflow/Snakefile.v2`, `workflow/rules/shared/00_common.smk`, `config/config_v2.yaml`,
`config/config_v2_test_examples.yaml`. Produced by a 13-agent Workflow (extract→reconcile
→draft→5-lens adversarial review); all review findings applied (0 blockers; the 1 major —
output_dir must be writable because `workdir:` enters it at parse time — fixed in the test
configs). Env consolidation deferred to when the first rule module needs an env (Stage 2).
Deliberate v1→v2 behavior changes recorded in code comments: underscore now allowed in
sample names (all modes); `CPUS`/`RAM`/`nt_version` now default instead of KeyError;
illumina now checks R2 mates at parse time; `dbcan_link`/`checkv_link` `.tar.gz`-validated
in all modes; dbCAN DB folder name now derived from the link (was hard-coded).

**Stage 2a — annotation + AMR (Tier S). ✓ DONE 2026-07-21 (gate passed).**
Split from the fuller Stage 2 below because `60_plasmid`/`70_phage` grew new-tool work
(geNomad, D8/D9) that `40_annotation`/`50_amr` do not have. Delivered on `release/v2.0.0`:
`workflow/rules/shared/40_annotation.smk` (bakta, eggnog, antismash db+run, dbcan db+run) and
`workflow/rules/shared/50_amr.smk` (ABRicate per-db + AMR_summary; CARD read-mapping leg
deferred to Stage 3). Faithful v1 ports onto the 00_common API; envs already present.
**Gate:** `snakemake -n` resolves the full 15-job DAG (bakta→eggnog/antismash/dbcan;
amr_contigs×8→AMR_summary) when a stub `contigs_final.fasta` + composition file exist; all
four modes still parse. Built by an 11-agent Workflow (extract→reconcile→draft→5-lens review);
all findings applied (0 blockers). Three 00_common additions this stage required:
`COMPOSITION` constant (single-sources the Bakta genus-file cross-stage edge, like
FINAL_CONTIGS), `ANTISMASH_DB_DIR` constant, and a project-wide `wildcard_constraints:
sample=<discovered names>` (kills the databases/{sample} path-collision class). Deliberate
additive change: AMR_summary now has a log (v1 had none). `cazyme_db_download` kept env-less
per v1 (uses wget from the launch env — see README_notes.md item 3, a deferred decision).

**Stage 2b — plasmid + phage (D8/D9). ✓ DONE 2026-07-22 (gate passed, both callers).**
`workflow/rules/shared/60_plasmid.smk` (Platon always; Platon+geNomad concordance when geNomad
opted in — D9) and `70_phage.smk` (VirSorter2 default / geNomad opt-in + CheckV — D8, revised
for licensing above). New: `envs/genomad.yaml`, `envs/checkv.yaml` (split from virsorter),
`envs/virsorter.yaml` (+`mamba>=1.5` packaging fix), `workflow/scripts/plasmid_concordance.py`
(+18 passing pytest cases), and 00_common additions (`PHAGE_CALLER` default virsorter2,
geNomad/Platon path constants, conditional rule-all plasmid target). **Gate:** `snakemake -n`
resolves the correct DAG for BOTH `phage.caller` values — `virsorter2` (5 jobs, zero geNomad
rules, fully permissive-licence path) and `genomad` (6 jobs, zero VS2 rules, concordance leaf).
Built by a 12-agent Workflow (research→extract→reconcile→draft→5-lens review). The review caught
the geNomad **licensing** issue (see D8 revision) — 0 code blockers otherwise. Concordance-script
majors addressed: loud join-key audit (silent all-medium degradation now warns), ragged-row
logging, honest documented limitation on reverse conflicts. **Verify on first real geNomad run:**
geNomad output filenames/columns, `--restart` need, the reverse-conflict upgrade (read geNomad's
aggregated_classification.tsv), and VS2's `--skip-deps-install`/`--use-conda-off` flags.

### Stage-4 design item (logged 2026-07-22): Bakta `--replicons` for the long-read modes

Deferred to Stage 4 because the inputs come from the assembler front end. Currently BacFlux
passes no `--replicons` table, so Bakta treats even a **closed circular** replicon as an
incomplete linear contig. Per Bakta's docs the two columns differ in kind:

- **`topology` (circular/linear) affects GENE PREDICTION** — *"De novo-prediction via Pyrodigal
  respecting sequences' completeness"*, *"detection & annotation of features spanning sequence
  edges"*. Genes crossing the origin of a closed replicon are currently at risk in
  nanopore/hybrid (Flye+Medaka+dnaapler deliver circularised replicons).
- **`type` (chromosome/plasmid/contig) + `name` are primarily OUTPUT METADATA** (INSDC-shaped
  records; trivial downstream replicon attribution). No evidence found that it changes prediction.

Sources of truth already produced by the pipeline:
- **topology** — Flye `assembly_info.txt` circular flag, already parsed for `ignore_list`
  (`IGNORE_LIST_CMD` in 00_common). A hard fact: populate it always.
- **type** — dnaapler `all` already writes `{sample}_all_reorientation_summary.tsv` with a
  per-contig `Gene_Reoriented` column (**dnaA→chromosome, repA→plasmid, terL→phage**) plus
  identity/coverage. Verified on the real BacFluxL+ test output; BacFlux currently discards it.
  An orthogonal signal to Platon (RDS) and geNomad (gene content) — but caveated: circular
  contigs only (non-circular are `--ignore`d), long-read modes only, absence of a marker is not
  evidence of absence (diverse plasmid replicons lack a recognisable repA), and dnaapler is
  built for reorientation, not classification (off-label use).

**Recommendation:** always set `topology` from Flye; set `type` only where signals concur (e.g.
dnaapler dnaA + Platon chromosome), else the neutral `contig` — so a probabilistic call is never
baked into the annotation artifact. Separately, surface dnaapler's `Gene_Reoriented` as an extra
column in the long-read plasmid concordance, and use `terL` hits as a free cross-check against
the phage stage.

**Stage 2 — Clean-shared tail (Tier S). [original combined scoping, superseded by 2a/2b above]**
`40_annotation`, `50_amr` (abricate leg), `60_plasmid`, `70_phage`. These consume
`contigs_final.fasta` and are the safest (identical rules). **Gate:** with a stub
`contigs_final.fasta`, each rule runs and matches the baseline output for that step.

**Stage 3 — Mode-parameterized shared (Tier P). ✓ DONE 2026-07-22 (gate passed).**
Delivered `workflow/rules/shared/{10_decontam,20_qc,30_taxonomy,90_report}.smk` + the CARD
read-mapping leg appended to `50_amr.smk`. **Contract closure achieved:** `select_contigs` now
produces `COMPOSITION` and `blast_contigs` produces `BLASTOUT`, the two cross-stage inputs
Stages 2a/2b already consumed as stubs. **Gate:** all four modes parse, and from a single
front-end stub the full shared tail builds as ONE 16-job DAG (map_contigs -> blast -> blob ->
select_contigs -> annotation / plasmid / QC / taxonomy / phage / multiqc).

Built by an 11-agent Workflow; the review found 11 majors, all addressed:
- **Hybrid plasmid regression (blocker-class):** the decontam BLAST screens the Illumina draft
  (SPAdes names) while Platon runs on the delivered ONT genome (Flye names), so the contig-ID
  lookup could never match and every hybrid plasmid would silently read "not verified". Fixed
  with `PLASMID_BLASTOUT` + a `blast_final_contigs` rule — and the guard was widened from
  hybrid-only to **both long-read modes** (`NEEDS_FINAL_BLAST`), since nanopore has the same
  mismatch across Medaka.
- **plasmid_search hardening:** restored the v1 long-read form — `grep -F` (fixed string),
  `grep -qi` (NCBI titles capitalise "Plasmid", so the case-sensitive grep silently missed
  them), first-token contig IDs via awk (whole-header IDs can never match), `: >` truncation,
  and Platon exit-code capture.
- **Bakta genus:** skip the literal `no-hit` genus, and run Bakta exactly once with the hint
  applied conditionally — v1's `for`-loop form would skip annotation entirely if no genus
  survived.
- **GTDB-Tk ordering:** restored its scheduling edge on CheckM (v1 guaranteed serialisation;
  without it two pplacer runs can collide and OOM on a large machine).
- **Config merge leak:** Snakemake *merges* a hard-coded `configfile:` with `--configfile`, so
  the v1 config was silently supplying values (verified). The default is removed; the config is
  now required, with an actionable message. This immediately exposed a latent `card_link`
  dependency, now validated by name in 00_common.
- **`mobilome.run: true`** aborted with a MissingInputException naming a directory; now exits
  with a message naming the config key.
- **Report wiring:** QUAST renames are generated from `QC_GENOMES` (they were hand-typed and
  branched on mode, so adding a QC genome updated CheckM/GTDB-Tk but silently not QUAST), and
  staged CheckM/GTDB-Tk rows get prefix-stripping rules — without them MultiQC's `-d` prepended
  the staging path and the hybrid Illumina-vs-ONT relabelling never reached the report.

**Stage 3 — Mode-parameterized shared (Tier P). [original scoping]**
`10_decontam`, `20_qc`, `30_taxonomy`, `50_amr` (CARD leg), `90_report`. Implement the
hybrid dual-genome branch and the D3/D4 decisions here. **Gate:** dry run of all modes still
parses; the dual-genome path resolves correct targets for hybrid.

**Stage 4 — Front ends. ✓ DONE 2026-07-22 (gate passed, all four modes).**
Delivered `workflow/rules/{illumina,nanopore,hybrid,contigs}/*.smk`, plus
`shared/15_replicons.smk` + `workflow/scripts/build_bakta_replicons.py` (+16 passing tests)
for the Bakta `--replicons` work. **This is the stage that flips `all_targets()` from empty to
the full pipeline:** `snakemake -n` now plans a complete run for every mode — illumina 41 jobs,
nanopore 41, hybrid 53, contigs 32. `build_replicons` fires only in the long-read modes (0 jobs
in illumina/contigs); hybrid stages BOTH genomes for dual CheckM/GTDB-Tk.

Review found 7 majors, all addressed. The most serious were data-loss and silent-failure bugs:
- **DRAFT_CONTIGS was nested inside the assembler's `directory()` output** (`SPADES_DIR` was
  derived as `dirname(DRAFT_CONTIGS)`). Snakemake wipes a `directory()` output before re-running
  its rule, so any re-run of the assembler would have silently deleted the hand-off the entire
  contamination screen keys on — reproduced by the reviewer on a minimal workflow. Fixed by
  moving `contigs_filt.fasta` up beside the assembler directory (where contigs mode already put
  it) and spelling `SPADES_DIR` out independently.
- **Medaka auto-model resolution** was switched to the RAW ONT file; v1 used the filtlong output.
  Restored — the raw input may be gzipped while the filtlong output never is, so the draft
  introduced a failure mode v1 could not have.
- **The replicons script treated a total join failure as a warning**, emitting a well-formed
  all-linear table and exiting 0 — Bakta would then annotate as if no table existed, with the
  only trace in an unread log. Zero overlap is now fatal (partial overlap stays a warning, since
  dnaapler only reports contigs it could reorient).
- **Hybrid's Qualimap panel was labelled generically** "mapping QC", though it is Illumina reads
  on the PRE-decontamination SPAdes draft — a reader would take it as coverage of the delivered
  ONT genome. Now tagged "mapping Illumina QC" in hybrid, matching the other tagged panels.

Still open (logged, not blocking): Bakta's honouring of `type=contig` + `topology=circular` is
pin-sensitive and nothing detects a future Bakta change silently reverting it to linear; and
`FASTA_HEAD_CMD` strips the descriptions Bakta uses to infer topology in **contigs** mode, so a
closed Unicycler/NCBI input genome is annotated as linear contigs.

**Stage 4.5 — Real end-to-end validation, illumina mode. ✓ PASSED 2026-07-22.**
Full run on strain CDRTa11 at `/media/data/antonielli_dir/BacFlux_v2_validation/illumina/`
(see its `RUN_NOTES.md` for the complete table and the caveats). 10/10 steps, zero errors,
every stage compared against the v1.3.1 baseline in `BacFlux_test/output_dir`:

- **Both the draft and the decontaminated assembly are byte-identical to v1.3.1**
  (`33b1582f…`, `d8c836d9…`). The v2 restructuring changed the plumbing, not the biology.
- Identical downstream: CheckM (35 contigs / N50 746047), GTDB-Tk (*Arthrobacter*), Bakta
  (4748 features), ABRicate (8 DBs, tables byte-identical), CARD mapping (6052 rows),
  eggNOG (4360), antiSMASH (7 regions), dbCAN (313 CAZymes, 5.1.2 from Zenodo).
- **The phage stage completed for the first time on BacFlux short-read.** v1.3.1 never
  produced VirSorter2 output at all. Zero viral contigs is the correct answer here:
  BacFluxL (same strain) also calls zero; only the more contiguous BacFluxL+ hybrid
  assembly pushes one marginal prophage past the 0.5 threshold.

Two defects found and fixed by this run:
- **`envs/virsorter.yaml` was not self-sufficient.** v2 carried `--use-conda-off` +
  `--skip-deps-install` (introduced in `dbf47b8`; present in NO v1 workflow) without the
  dependency list those flags require, so VirSorter2 died on `No module named 'screed'`.
  Fixed in `1f40303` by installing VirSorter2's own packaged `envs/vs2.yaml` list.
- **A false PASS in the comparison method itself.** The first ABRicate check read a v1 path
  that does not exist, so both sides reported "0 hits, same". Baseline comparisons must
  assert the reference file exists before comparing counts.

Not covered by this run, still to do: nanopore and hybrid validations; the `checkv_db`
DOWNLOAD rule (the database was staged by hand because portal.nersc.gov was unreachable).

**Stage 4 — Front ends, one mode at a time, each a gate. [original scoping]**
Do them in increasing complexity: `illumina` → `contigs` → `nanopore` → `hybrid`.
After each, run the full mode end-to-end on the Stage-0 isolate. **Gate per mode:** output
matches the baseline oracle (allowing for the intentional D1 layout change and the D4 phage
change, which are documented, not regressions).

**Stage 5 — Docs.**
Single README with per-mode sections; `CITATION.cff` with the concept DOI; changelog
documenting the merge, the layout change (D1), and the phage-input change (D4).

**Stage 6 — Release v2.0.0.**
Tag on the BacFlux repo → Zenodo mints a new version DOI under the existing concept DOI
`10.5281/zenodo.11143917`; set the updated title/authors on that record. Archive the
BacFluxL and BacFluxL+ repos (GitHub → read-only) with a banner pointing here; their DOIs
stay valid for reproducibility.

**Stage 7 — Mobilome module.**
Now add it **once** in `shared/80_mobilome.smk` + `scripts/mobilome/`, gated by
`config.mobilome.run`, following `mobilome_module_SPEC.md`.

---

## 8. Risks & must-verify-during-implementation

- **Parse-time side effects.** Each current Snakefile runs banner prints, `glob_wildcards`,
  `sys.exit` guards, and (BacFlux/L) a `global id` mutation *unconditionally* at import.
  All of this must move behind the mode dispatch in `00_common.smk` or it will fire for the
  wrong mode. (MetaFlux already demonstrates the gated pattern.)
- **Conditional/anonymous Medaka rules** (nanopore, hybrid) defined inside `if USE_MEDAKA`
  must be reconciled into the `nanopore/` and `hybrid/` modules without changing DAG identity.
- **Item-level behavior derived from agent summaries, not yet line-verified:** the exact
  hybrid dual-genome checkm/gtdbtk wiring, and the Snippy 4-stage comparison. Read these
  rule bodies in full before porting.
- **Regression discipline:** the Stage-0 oracle is the whole safety story. Every mode gate
  compares against it. Only two intentional differences are allowed and must be documented:
  D1 (directory layout) and D4 (phage input).
- **Readability:** keep each `.smk` flat and heavily commented in the house style. The only
  new abstraction the reader must learn is the mode-dispatch layer in `00_common.smk` — keep
  that one file exemplary.

---

## 9. What this unlocks

One repo, one README, one env set, one DOI lineage. The mobilome module — and every feature
after it — is written and maintained **once** across short-read, long-read, hybrid, and
pre-assembled inputs, instead of four times.
