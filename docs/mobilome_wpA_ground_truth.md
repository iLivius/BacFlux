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

---

# AMRFinderPlus 4.2.7 OUTPUT CONTRACT (verified by a real run, 2026-07-24)

Ran, successfully (exit 0, 62 s), on real v2 output — a validation isolate (Pseudomonas_E),
proving spec §12 Q1 empirically:

```
amrfinder -p <bakta>/006.faa -n <asm>/contigs_final.fasta -g <bakta>/006.gff3 \
  --annotation_format bakta --plus \
  -d {bakta_db}/amrfinderplus-db/latest --threads 8 -o out.tsv
```
`--annotation_format bakta` WITH `-g` is accepted. No fallback needed.

## VERBATIM column order (22 columns) — the parser contract

 1 Protein id                     12 Subclass
 2 Contig id                      13 Method
 3 Start                          14 Target length
 4 Stop                           15 Reference sequence length
 5 Strand                         16 % Coverage of reference
 6 Element symbol                 17 % Identity to reference
 7 Element name                   18 Alignment length
 8 Scope                          19 Closest reference accession
 9 Type                           20 Closest reference name
10 Subtype                        21 HMM accession
11 Class                          22 HMM description

**NOTE the spec (and older AMRFinderPlus docs) say "Gene symbol" — in 4.2.7 it is
`Element symbol` / `Element name`. A parser keying on "Gene symbol" silently
finds nothing. Column names above are copied from the real header.**

- Coordinates: `Start`/`Stop`, 1-based inclusive; `Strand` is +/-.
- `Type` = AMR | STRESS | VIRULENCE (the --plus categories).
- `Method` = the free confidence tier: observed BLASTP and HMM on real data;
  the full vocabulary includes EXACTX/BLASTX/PARTIALX/POINTX etc. A PARTIAL*
  method value is how a truncated hit (typically at a contig end) is signalled —
  that is the honest short-read fragmentation flag the spec asks for.

## Real result on a validation isolate (sanity anchor for tests)

5 AMR rows, all Type=AMR, all on contig_1 (the chromosome): emhC, emhB, mexE
(EFFLUX), aac(6') (AMINOGLYCOSIDE), ampC (BETA-LACTAM); Methods BLASTP x2, HMM x3.
Biologically consistent with the ABRicate result for the same isolate (intrinsic
Pseudomonas efflux/AmpC, chromosomal, nothing acquired) — a good regression anchor
and a good demonstration of the intrinsic-vs-acquired framing the module exists to
support.

---

# PLATON OUTPUT CONTRACT (verified from a real run, 2026-07-24)

`06.plasmids/{sample}/platon/contigs_final.tsv` — verbatim header:

```
ID  Length  Coverage  # ORFs  RDS  Circular  Inc Type(s)  # Replication
# Mobilization  # OriT  # Conjugation  # AMRs  # rRNAs  # Plasmid Hits
```

Plus a hard chromosome/plasmid split as FASTA:
- `contigs_final.chromosome.fasta` — contigs Platon calls chromosomal
- `contigs_final.plasmid.fasta`    — contigs Platon calls plasmid

**IMPORTANT for the mobility ladder (a better route than the spec assumed):**
Platon ALREADY reports `# Mobilization`, `# OriT` and `# Conjugation` per contig.
That distinguishes ladder tier 5 (mobilisable — needs a helper) from tier 6
(conjugative — predicted self-transmissible) FOR PLASMID CONTIGS directly, with
no extra tool. Real example, a validation isolate contig_2 (a confirmed 58 kb plasmid):
`# Mobilization = 1`, `# Conjugation = 1`, `# OriT = 0`, `Inc Type(s) = 0`.

So the replicon-call input for WP-D is built as:
- contig in `*.plasmid.fasta`      -> replicon = plasmid
- contig in `*.chromosome.fasta`   -> replicon = chromosome
- a plasmid contig with `# Conjugation > 0`  -> conjugative      (tier 6 eligible)
- a plasmid contig with `# Mobilization > 0` or `# OriT > 0` but no conjugation
                                    -> mobilisable     (tier 5)
- otherwise                         -> non-mobilisable plasmid

CONJscan remains the route for CHROMOSOMAL conjugation machinery (an ICE — the
case that breaks the naive "chromosomal therefore not transferable" assumption),
which is exactly where Platon says nothing because the contig is not a plasmid.
The two are complementary, not redundant.

---

# ISESCAN 1.7.3 CONTRACT (verified by building the env and RUNNING it, 2026-07-24)

## Real CLI (from `isescan.py --help` in the built env — NOT from docs)

```
isescan.py --seqfile SEQFILE --output OUTPUT [--nthread N]
           [--removeShortIS] [--no-FragGeneScan]
```
- `--seqfile` genome FASTA, `--output` output DIRECTORY, `--nthread` CPUs.
- `--removeShortIS` drops incomplete IS (length < 400, or single-copy without a
  perfect TIR). **BacFlux must NOT pass it** — the spec wants partials kept so we
  can tier them ourselves, which is the honest treatment on fragmented assemblies.
- The entry point is `isescan.py` (a script on PATH), not `isescan`.

## GOTCHA 1 — it must run with the ENV's python (libssw.so)

Invoking the script by absolute path while the env's `bin/` is not first on PATH
picks up the SYSTEM python and dies with:
`OSError: libssw.so: cannot open shared object file`.
`libssw.so` IS shipped, at `<env>/lib/libssw.so`; the failure is purely that the
wrong interpreter (and so the wrong library path) was used. Snakemake's `conda:`
directive activates the env properly, so the rule is fine as long as it calls
`isescan.py` by NAME (letting PATH resolve it) and never by absolute path.

## GOTCHA 2 — output is nested under a path derived from the INPUT's location

With `--output <outdir>`, ISEScan does NOT write results flat into `<outdir>`.
It creates `<outdir>/proteome/<name-of-input's-parent-dir>/…`,
`<outdir>/hmm/<same>/…` etc. Running on
`…/02.assembly/386/contigs_final.fasta` produced
`out/proteome/386/contigs_final.fasta.{faa,gff,out,ffn}` — the `386` component
comes from the input file's PARENT DIRECTORY name, not from anything we passed.

Consequence for the rule: declare the output as a DIRECTORY and have the parser
FIND the results file inside it (walk for the expected suffix) rather than
hard-coding a path. In BacFlux the parent dir happens to be `{sample}`, so the
nesting is predictable today — but relying on that would be a silent trap the
moment the input layout changes.

## Runtime

FragGeneScan then hmmsearch; several minutes on a 4.2 Mb genome with 8 threads.
Not a trivial rule — it deserves a real thread allocation.

## ISEScan RESULTS SCHEMA (from the real run — 31 IS on a validation isolate)

Output lands in `<outdir>/<input-parent-dir-name>/<input-basename>.{tsv,csv,gff,sum,is.fna,orf.faa,orf.fna,out,raw}`
(plus `proteome/` and `hmm/` working dirs). The main table is the **`.tsv`**.

VERBATIM 24-column header:
```
 1 seqID       7 ncopy4is   13 irId      19 orfLen
 2 family      8 start1     14 irLen     20 E-value
 3 cluster     9 end1       15 nGaps     21 E-value4copy
 4 isBegin    10 start2     16 orfBegin  22 type
 5 isEnd      11 end2       17 orfEnd    23 ov
 6 isLen      12 score      18 strand    24 tir
```
Key columns for BacFlux:
- `seqID` contig, `isBegin`/`isEnd` 1-based inclusive IS coordinates, `isLen`.
- `family` (IS21, IS3, IS481, IS110, IS256, ISNCY…) and `cluster` (family_NNN).
- **`type` = `c` (complete) or `p` (partial)** — this is the complete/partial flag
  the spec wants tiered rather than discarded (so never pass --removeShortIS).
- **`ncopy4is`** = copy number of that IS in the assembly — directly relevant to
  the collapse problem (multi-copy IS are what break contigs).
- `tir` = the two terminal inverted repeats, colon-separated; `irId`/`irLen` score them.
- `strand`, `orfBegin`/`orfEnd` = the transposase ORF.

There is ALSO a `.gff` with SO-typed features, already using the ontology terms
the spec's §9 output schema adopts:
```
contig_1 ISEScan insertion_sequence        36 1337 . + . ID=contig_1_IS_1;family=IS21;cluster=IS21_259
contig_1 ISEScan terminal_inverted_repeat  36   53 . + . ID=contig_1_IS_1_TIR;parent=contig_1_IS_1
```
And a `.sum` per-family summary (nIS, %Genome, bps4IS, dnaLen) — useful for the QC
line but BacFlux computes its own contig-boundary statistics, which `.sum` lacks.

REAL RESULT, a validation isolate (Arthrobacter, 4.23 Mb, 2 contigs incl. a 58 kb plasmid):
31 IS elements — 18 complete, 13 partial; families IS481 (9), IS3 (6), ISNCY (5),
IS21 (4), IS110 (4), IS256 (3); 0.82% of the genome. **All 31 on contig_1 (the
chromosome); none on the plasmid contig.** Saved as test fixtures:
`workflow/scripts/80_mobilome/testdata/isescan_386_real.{tsv,sum}`.

---

# FIRST END-TO-END RUN OF THE MODULE ON REAL DATA (2026-07-24)

Ran the three scripts against genuine pipeline outputs, no mocks:

**Sample 386** (Arthrobacter, 4.23 Mb, confirmed 58 kb conjugative plasmid):
- `isescan_to_table.py` parsed the REAL ISEScan run: **31 IS, 18 complete /
  13 partial** — matching an independent count of the raw file exactly. QC:
  1 IS within 100 bp of a contig end (fraction 0.032). It auto-discovered the
  results file inside ISEScan's nested output directory.
- Replicon table built from REAL Platon output: contig_1 chromosome,
  contig_2 plasmid + **conjugative** (from Platon's own `# Conjugation` column) —
  i.e. tier-6 eligible, exactly as designed.
- AMRFinderPlus on 386: **0 AMR rows**. The module handled the empty input
  gracefully and wrote a well-formed empty table (exit 0), as required.

**Sample 006** (Pseudomonas_E, 5 real AMR genes) — the decisive test:
all five genes placed at **mobility tier 1 (chromosomal, no MGE context =
intrinsic candidate)**: emhC, emhB, mexE (efflux), aac(6'), ampC. That is the
biologically correct call for Pseudomonas efflux pumps and AmpC, and it is
exactly the intrinsic-vs-acquired evidence the module exists to produce.

**The audit file is the part worth keeping.** Run without an IS table, it did NOT
silently report "no mobile element". It recorded, per gene, that no IS calls were
available, that "absence of context" therefore means "not looked at" rather than
"looked at and found nothing", and it **capped confidence at medium** for that
reason. That distinction is the whole point of the tiered-evidence design.

Deliverable schema: 37 columns, covering everything spec §9 asks for (replicon,
mge_context, distance_bp, orientation, is_family, n_flanking_is,
same_orientation, is_inside_amr_cds, spans_contigs, mobility_tier, confidence)
plus the honest short-read signals (amr_partial_at_contig_end, contig_length,
dist_to_contig_end, is_at_contig_boundary).

---

# CONJSCAN / MACSYFINDER 2.1.6 CONTRACT (built the env and RAN it, 2026-07-24)

## Model package — installable, no licence problem

`msf_data available` lists **CONJScan (2.1.0)**; installed cleanly with:
```
msf_data install --target <models-dir> CONJScan
```
**GOTCHA: `macsydata` is DEPRECATED in 2.1.6 and prints a rename warning — the
current command is `msf_data`.** A rule copied from older docs would emit noise
today and break when the alias is dropped.

The package ships TWO model sets:
- `CONJScan/Chromosome/…` — conjugation machinery on a CHROMOSOME, i.e. the ICE
  case. This is the one that matters most: it breaks the naive assumption that a
  chromosomal gene cannot transfer (spec §2.4).
- `CONJScan/Plasmids/…`   — the plasmid case (Platon already gives us plasmid
  mobility directly, so this is the cross-check, not the primary route).

Models per set: `MOB` (relaxase), `T4SS_type{B,C,F,FA,FATA,G,I,T}` (the MPF
types), and `dCONJ_type*` (decayed/degraded systems — the "overcalling" trap the
spec warns about at §8 Phase 4). Profiles include the relaxase families
`T4SS_MOB{B,C,F,H,M,P1,P2,P3,Q,T,V}` and coupling proteins `T4SS_t4cp1/t4cp2`.

## Real invocation (verified, exit 0)

```
macsyfinder --models CONJScan/Chromosome all \
  --sequence-db <bakta>/{sample}.faa \
  --db-type ordered_replicon \
  --models-dir <models-dir> \
  --out-dir <out> --worker N
```
`ordered_replicon` is correct for a single genome's proteome, and Bakta's `.faa`
is emitted in genome order, so the ordering assumption holds.

## Output — `best_solution.tsv`, verbatim 22 columns

```
 1 replicon      7 sys_loci        13 hit_status       19 hit_begin_match
 2 hit_id        8 locus_num       14 hit_seq_len      20 hit_end_match
 3 gene_name     9 sys_wholeness   15 hit_i_eval       21 counterpart
 4 hit_pos      10 sys_score       16 hit_score        22 used_in
 5 model_fqn    11 sys_occ         17 hit_profile_cov
 6 sys_id       12 hit_gene_ref    18 hit_seq_cov
```
- `gene_name` gives the machinery component (e.g. `T4SS_MOBP1`, `T4SS_t4cp2`) —
  this is what identifies relaxase type vs coupling protein vs MPF.
- `model_fqn` gives the system model (`CONJScan/Chromosome/MOB`).
- `sys_id` groups hits into one system; **`sys_wholeness`** says how complete that
  system is — directly usable as the "machinery_intact / degraded" downgrade the
  spec asks for (§8 Phase 4).
- Also written: `all_systems.tsv`, `best_solution_summary.tsv` (a per-model COUNT
  matrix, one row per replicon), `rejected_candidates.tsv`, `hmmer_results/`.

## REAL RESULT — a validation isolate (Arthrobacter), and why it matters

Two MOB systems found ON THE CHROMOSOME, each `sys_wholeness = 0.667` (i.e.
INCOMPLETE):
- `386_MOB_1`: `T4SS_MOBP1` (relaxase) + `T4SS_t4cp2` (coupling protein)
- `386_MOB_2`: `T4SS_MOBF` (relaxase) + `T4SS_t4cp2` (coupling protein)
No T4SS/MPF system was found (all `T4SS_type*` columns are 0).

Biologically this is exactly the nuance the ladder is designed to express:
relaxase + coupling protein but NO mating-pair-formation apparatus = the element
can be MOBILISED by a helper, but cannot self-transmit. Tier 5, not tier 6 — and
`sys_wholeness = 0.667` is the honest flag that even that is a partial system.
Reporting "conjugation machinery present" as a bare fact here would be misleading;
reporting the type, the wholeness and the missing MPF is the useful answer.

Fixture saved: `workflow/scripts/80_mobilome/testdata/conjscan_386_real_best_solution.tsv`.

## Empty-result behaviour (important for graceful degradation)

Still to confirm on a genome with NO system at all: whether `best_solution.tsv`
is written empty/header-only or omitted. The parser must tolerate BOTH — most
environmental isolates carry no conjugative system, so this is the common path.

---

# ⚠ LICENCE FINDING — CONJScan MODELS ARE CC BY-NC-SA 4.0 (verified 2026-07-24)

**This is NEW and it changes the module's design. The spec (§3.1) listed CONJscan
as a plain "bioconda + macsydata install" adopt, with no licence flag.**

Verified by reading the installed package's own metadata, not a web claim:
`<models>/CONJScan/metadata.yml` says verbatim:

```
license: CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
copyright: 2018-2026, Institut Pasteur, CNRS
vers: 2.1.0
```
(a `LICENSE` file ships alongside it.)

## What this does and does NOT mean

- The **engine** (MacSyFinder 2.1.6, bioconda) is not the problem.
- The **models** (HMM profiles + XML definitions) are non-commercial + share-alike.
- They are **data, fetched at runtime, never vendored into this repo** — so under
  the spec §11 database rule, BacFlux's own MIT licence is NOT contaminated. This
  is the same position as `bakta_db`, `blast_db`, TnCentral, etc.
- BUT: making CONJscan a MANDATORY step would force every downstream user —
  including commercial ones — through a non-commercial-only resource. That is
  precisely the argument that made **geNomad opt-in rather than default** (§3.1,
  D8/D9 revision 2026-07-22).

## Consequence for the design — CONJscan must be OPT-IN, not automatic

Treat it exactly like geNomad:
- a separate config switch (e.g. `mobilome.conjscan.run`), **default false**;
- a clear README notice stating the models are CC BY-NC-SA 4.0, academic /
  non-commercial use, and that the user fetches them under their own agreement
  with Institut Pasteur/CNRS;
- the mobility ladder must DEGRADE GRACEFULLY without it. That is workable
  because Platon already supplies plasmid mobility (`# Conjugation`,
  `# Mobilization`, `# OriT` — see the Platon contract above), so tiers 5 and 6
  remain reachable for PLASMID-borne AMR with no CONJscan at all.
  What is lost when CONJscan is off is specifically the **ICE case**: conjugation
  machinery sitting on the CHROMOSOME. Those genes then stay at tier 1 with an
  audit reason saying the ICE check was not run — honest, not silently wrong.

## Second finding from the same research — the ordered_replicon caveat

`--db-type ordered_replicon` means "one complete genome". MacSyFinder takes the
replicon name from the sequence-db FILE BASENAME, so a multi-contig draft is
treated as ONE pseudo-replicon and the clustering can join anchors ACROSS contig
boundaries. On a fragmented short-read assembly that is a real false-positive
route. Mitigation: keep it, but hard-flag any system whose hits span contigs
(the `spans_contigs` column already exists in the deliverable) and cap those at
low confidence — consistent with how the rest of the module treats
contig-boundary evidence.

---

# DATABASE STATUS (verified live, 2026-07-24) — spec §12 Q3 ANSWERED

## ISOSDB — spec §12 Q3 resolved: NUCLEOTIDE ONLY, and it is MIT

- Source: the pseudoR repo itself; no Zenodo/figshare DOI exists.
  `https://github.com/joshuakirsch/pseudoR/raw/main/ISOSDB.V3.fna.zip`
  (~9.9 MB zip -> `ISOSDB.V3.fna`, 33 MB, **22,713 nucleotide sequences**).
- Repo LICENSE = **MIT**. Freely redistributable and auto-fetchable — no
  encumbrance, so this one can be default-on when the module runs.
- **NO protein FASTA anywhere in the repo.** So a protein leg would need prodigal
  run once over the nucleotide file. Per spec §5.4 BacFlux only needs the
  NUCLEOTIDE form (for the read-mapping copy-number leg), so this costs nothing
  today — the question is simply closed.
- Bonus files worth using: `ISOSDB_names.tsv` (name lookup) and
  `IS_fam_annot.txt` (family annotation) for the naming cascade, plus `IRs.fa`
  with a pre-built BLAST v5 index (150 bp IS terminal ends) — directly useful for
  the partial-IS / contig-end signal.

## TnCentral — all endpoints live; the spec listed only half of them

Six endpoints, all HTTP 200, application/zip, no auth, no click-through:

| endpoint    | main file                          | contains ISfinder? |
|-------------|------------------------------------|--------------------|
| `nc/tn`         | tncentral.fa                   | no  |
| `nc/tn_in`      | tncentral_integrall.fa         | no  |
| `nc/tn_is`      | tncentral_isfinder.fa          | YES |
| `nc/tn_in_is`   | tncentral_integrall_isfinder.fa| YES |
| **`prot/tn`**   | tncentral.prot.fa              | **no — the spec missed this one** |
| `prot/tn_is`    | tncentral_isfinder.prot.fa     | YES |

- Format confirmed empirically as **BLAST v4** (`.nhr/.nin/.nsq`, no `.ndb/.ntf/.nto`),
  so spec §5.3's dump-and-rebuild-as-v5 recipe stands.
- Undated URLs as warned, BUT every internal file has mtime **2025-05-16** — a
  usable de facto version stamp to record alongside the sha256.
- Terms: **no licence page, no terms-of-use, no commercial statement** — only the
  footer "© TnCentral 2024 - All Rights Reserved". Spec §5.2's "commercial use
  unverified — email them before claiming anything in the README" REMAINS OPEN.
- **Practical win: `prot/tn` is an ISfinder-FREE protein database**, so the
  blastp/diamond family-assignment leg can be default-on rather than opt-in.
  Spec §5.1 should gain `tncentral_prot_url` beside `tncentral_nt_url`.

## ISfinder — terms unchanged, still cannot be automated

Verbatim, still current: *"It is not permitted to download the ISfinder database
without written authorization. Moreover, it is also not permitted at any time to
distribute the database to third parties..."* — so ISfinder content stays behind a
config URL the user must supply themselves, exactly as the spec required.

## What WP-B may automate

- **Auto-fetch when the module runs:** ISOSDB (MIT), TnCentral `nc/tn`,
  `nc/tn_in`, `prot/tn` — no ISfinder content, no stated restriction beyond
  copyright. Record URL + sha256 + fetch date to a provenance TSV.
- **Config URL, default null, user pastes it:** `nc/tn_is`, `nc/tn_in_is`,
  `prot/tn_is` (ISfinder terms apply).
