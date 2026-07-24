# Mobilome WP-A — ground truth (verified empirically 2026-07-24)

Established BEFORE writing any rule, per `mobilome_module_SPEC.md` §1.5 / §12.
Every line here was checked against the live env and real v2 outputs, not inferred.
These findings CORRECT several assumptions in the spec (which pre-dates v2's paths).

## AMRFinderPlus availability — CONFIRMED

- **In the Bakta conda env already** — no new env, no new install.
  `amrfinder --version` → **4.2.7**, alongside `bakta` in the same built env
  (`.snakemake/conda/18ff8047534fdb974565f8170428ceb6_`).

## Spec §12 Q1: does `--annotation_format bakta` work? — YES

- amrfinder 4.2.7 `-a/--annotation_format` accepts: bakta, genbank, microscope,
  patric, pgap, prodigal, prokka, pseudomonasdb, rast, standard.
- So the spec's fallback ("if not accepted, drop -g and run -p + -n only") is NOT
  needed. Use `--annotation_format bakta` with `-g`.

## DB PATH — spec's example rule is WRONG, corrected here

- Spec §4 example: `db = os.path.join(config["bakta_db"], "amrfinderplus-db")`.
- Reality: that parent dir holds VERSIONED subfolders
  (`2024-12-18.1/`, `2026-01-21.1/`) + a `latest -> 2026-01-21.1` symlink.
  The built BLAST index (`AMRProt.fa.phr/.pin/.psq`) lives INSIDE the versioned
  folder, not the parent. Pointing `-d` at the parent fails with
  "The BLAST database for AMRProt.fa ... was not found".
- **CORRECT WP-A path: `{bakta_db}/amrfinderplus-db/latest`** (symlink resolves to
  the current versioned DB, index present). Verified: `--list_organisms` succeeds
  against `.../latest`, fails against the parent.
- No index-building rule needed — Bakta ships the built index in the versioned
  folder. (The dir is also writable, so a rebuild would be possible if ever
  needed, but it is not.)

## Spec §12 Q2: GTDB→AMRFinderPlus organism mapping — 31 curated organisms

`amrfinder --list_organisms -d .../latest` returns exactly these `--organism`
values (build the GTDB-species → this-name lookup against this list; emit
`--organism X` only on a match, else empty string, silently):

  Acinetobacter_baumannii, Bordetella_pertussis, Burkholderia_cepacia,
  Burkholderia_mallei, Burkholderia_pseudomallei, Campylobacter,
  Citrobacter_freundii, Clostridioides_difficile, Corynebacterium_diphtheriae,
  Enterobacter_asburiae, Enterobacter_cloacae, Enterococcus_faecalis,
  Enterococcus_faecium, Escherichia, Haemophilus_influenzae, Helicobacter_pylori,
  Klebsiella_oxytoca, Klebsiella_pneumoniae, Neisseria_gonorrhoeae,
  Neisseria_meningitidis, Pseudomonas_aeruginosa, Salmonella, Serratia_marcescens,
  Staphylococcus_aureus, Staphylococcus_pseudintermedius, Streptococcus_agalactiae,
  Streptococcus_pneumoniae, Streptococcus_pyogenes, Vibrio_cholerae,
  Vibrio_parahaemolyticus, Vibrio_vulnificus

Note: NONE of the 6 screening isolates so far (Paenibacillus, Arthrobacter,
Pseudomonas_E) map to a curated `--organism` — Pseudomonas_E is not P. aeruginosa,
and the others have no entry. So the silent-fallback (no --organism) path is the
COMMON case for this collection and must be well-tested, not an afterthought.

## Extra 4.2.7 capabilities the mobilome can exploit (not in the spec)

- `--mutation_all MUT_ALL_FILE` — reports point mutations even below the report
  threshold. This is the intrinsic/chromosomal-resistance signal ABRicate cannot
  see (spec §4 item 5). Only meaningful with `--organism`, so gated on the mapping.
- `--nucleotide_flank5_output` + `--nucleotide_flank5_size` — emits each hit's 5'
  flanking sequence. Directly useful for the AMR×MGE co-localisation step (WP-D):
  the flank is where an upstream IS/promoter would sit.

## v2 path corrections (spec §1.2 was explicitly UNVERIFIED)

Spec assumed v1 paths; v2 (post-D1 renumbering) reality, from a real run:
- Genome:  `FINAL_CONTIGS = 02.assembly/{sample}/contigs_final.fasta`
           (spec said `.../selected_contigs.fasta` — WRONG)
- Bakta:   `04.annotation/bakta/{sample}/{sample}.{faa,fna,gff3,gbff,tsv}`
           (spec said `05.annotation/...` — renumbered)
- AMR stage: `05.amr/` (spec said `06.AMR/` — renumbered)
- Mobilome stage: `08.mobilome/` (already reserved in config_v2.yaml, run:false)

## Net effect on WP-A

The rule is buildable exactly as the spec envisioned, with three concrete fixes:
1. `-d {bakta_db}/amrfinderplus-db/latest` (not the parent).
2. All paths on the v2 (renumbered) layout, FINAL_CONTIGS not selected_contigs.
3. The GTDB→organism map is a real 31-entry table; the no-match path is the
   common case here and must be graceful + silent.
No new conda env, no new database download, no index-build rule. Cheapest win, as
the spec promised — now with the landmines removed before the first line of rule.

---

# New-tool bioconda availability (verified 2026-07-24) — hard constraint HOLDS

CLAUDE.md rule: every tool must be bioconda-installable. Checked the mobilome's
three new tools:

- **ISEScan** — bioconda, latest **1.7.3** (WP-C, IS detection on contigs). ✓
- **MacSyFinder** — bioconda, latest **2.1.6**. CONJscan runs as a MacSyFinder
  model package (`macsydata install`, a runtime model fetch — engine is bioconda,
  matches spec §3.1). ✓
- **bedtools** — not re-searched (conda search is slow under load); a certainty
  on bioconda, one of the most standard genomics tools. WP-D co-localisation. ✓

Net: the whole A→D + CONJscan plan stays inside bioconda. No tool forces a
departure from the conda-per-rule model. The only non-bioconda fetches are
DATABASES (ISOSDB, TnCentral) — handled by the same user-download + config-URL
pattern BacFlux already uses for bakta_db/blast_db/etc. and just extended to
VS2/antiSMASH/dbCAN/CARD earlier today.

STILL OPEN (spec §12 Q3, not yet checked — needs the pseudoR repo): does ISOSDB
ship a protein FASTA, or must prodigal be run once over the nucleotide FASTA?
Deferred to WP-B/C implementation (not on the WP-A critical path).
