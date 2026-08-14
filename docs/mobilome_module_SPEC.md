# BacFlux — Mobilome & AMR-Mobility Module: Implementation Spec

**Status:** **IMPLEMENTED** in v2.0.0 (as of 2026-07-31). This is now a *design
record*, not a work order. See the box below before acting on anything in it.
**Author of design notes:** Livio Antonielli (iLivius), with an AI brainstorming session
**Target repo:** `https://github.com/iLivius/BacFlux` (MIT) — written when BacFlux had long-read siblings BacFluxL / BacFluxLplus; v2.0.0 merged all of them into this one repository.
**Document purpose:** give a coding agent everything needed to reconstruct context and start implementing without re-deriving the reasoning.

---

> ## ⚠ How to read this document now
>
> This spec was written **before** the module existed, and it still reads like a
> plan. The module is built: **22 rules in
> `workflow/rules/shared/80_mobilome.smk`**, with the helper package in
> `workflow/scripts/80_mobilome/` and a passing test suite
> (`pytest workflow/scripts` — 420 tests). Some sections were updated in place as
> they were implemented (§5.3 and §5.5 carry dated correction boxes); most were
> not.
>
> **Three things in here are stale by construction, and are flagged where they
> appear:**
>
> 1. **The v1 output paths in §1.2** (`05.annotation`, `06.AMR`, `07.plasmids`,
>    `08.phages`, `04.taxonomy`) were renumbered in v2. The real paths were
>    established in [`mobilome_wpA_ground_truth.md`](mobilome_wpA_ground_truth.md)
>    — read that, not §1.2.
> 2. **The rule-file path in §10** (`workflow/rules/mobilome.smk`) is not where
>    the code lives.
> 3. **"BacFluxL only"**, in §2.3, §6, §8 and §13, meant "the long-read sibling
>    repository". There is no such repository any more. Those parts run in **all
>    four v2 modes**; the *caveat* behind the restriction — that they are much less
>    reliable on a fragmented short-read assembly — is still entirely real and is
>    the reason the module reports contig-edge flags on every row.
>
> For what the module actually does and how well it works, read the README's
> *Validation* section, [`methods_att_and_small_plasmids.md`](methods_att_and_small_plasmids.md),
> [`methods_icescan_union.md`](methods_icescan_union.md),
> [`methods_ebi_comparison.md`](methods_ebi_comparison.md) and
> [`mobilome_worked_example.md`](mobilome_worked_example.md).

---

## 0. TL;DR for the agent

Build a new optional Snakemake module that answers one question per sample:

> **For each AMR gene detected, is it embedded in a mobile genetic element, and if so, how transferable is that element?**

It does this by integrating outputs from tools BacFlux already runs (Bakta, ABRicate, Platon, GTDB-Tk) plus three new ones (AMRFinderPlus surfaced properly, ISEScan, CONJscan), then intersecting coordinates and emitting a tiered evidence table.

Nothing here requires new heavy dependencies. The single biggest new piece of original code is an att-site/direct-repeat search (~200 lines Python).

---

## 1. Existing BacFlux context the agent must know

### 1.1 Repo layout
```
config/          config.yaml
workflow/        Snakefile (main, short-read), FastaFlux (pre-assembled contigs entry point), envs/
miscellaneous/
```
Run with `snakemake --sdm conda`. Conda-per-rule is the dependency model. Everything must stay installable via bioconda.

### 1.2 Relevant existing outputs

> ⚠️ **SUPERSEDED — DO NOT USE THESE PATHS.** They were guessed from the v1 README
> and never corrected here. The ground truth was established by reading the
> Snakefile and is written up in
> [`mobilome_wpA_ground_truth.md`](mobilome_wpA_ground_truth.md); **that is the
> document to use.**
>
> Two things changed. First, the original warning below ("run the greps in §1.5
> and paste the results in") was never actioned, so the table is still guesswork.
> Second, and more importantly, **v2 renumbered every stage**, so even a corrected
> v1 table would now be wrong:
>
> | this table says (v1) | v2 actual |
> |---|---|
> | `04.taxonomy` | `03.taxonomy` |
> | `05.annotation` | `04.annotation` |
> | `06.AMR` | `05.amr` |
> | `07.plasmids` | `06.plasmids` |
> | `08.phages` | `07.phages` |
> | — | `08.mobilome` |
>
> The tool versions below are also v1: v2 runs **GTDB-Tk 2.7.2 against GTDB R232**,
> not 2.6.1 / R226. The table is kept only because §1.5 and the surrounding
> reasoning refer to it.

**Known facts (from v1.3.1 release notes, verified):**
- ABRicate runs as a **named wildcard rule `amr_contigs`** across multiple databases (anonymous per-database rules were replaced in v1.3.1). So there are *several* ABRicate outputs per sample — the mobilome module must choose one (or merge) rather than assume a single file.
- `workflow/scripts/` already exists (e.g. `select_contigs_by_taxonomy.py`) — the proposed `workflow/scripts/80_mobilome/` package fits the existing convention.
- `config/config_custom.yaml` exists alongside `config/config.yaml`; dry-runs are done against the custom one.
- **Path layouts differ across the family.** BacFluxLplus numbers directories differently (`09.taxonomy`, `10.annotation`, `12.plasmids`) and adds dnaapler / Polypolish / Snippy stages; BacFluxL uses older tool versions (ABRicate v1.0.1, antiSMASH 7.1.0) and runs both Prokka and Bakta. **Never hard-code a path across repos** — resolve via config or a per-repo path map.


| Output | Path (under `output_dir`) | Produced by |
|---|---|---|
| Selected (decontaminated) contigs | `02.assembly/{sample}/…selected…fasta` | SPAdes + BlobTools selector |
| Contig taxonomy audit | `contig_taxonomy_decisions.tsv` | decontamination selector |
| Bakta annotation | `05.annotation/bakta/{sample}/{sample}.{gff3,tsv,faa,fna,gbff}` | Bakta v1.12.0 |
| ABRicate AMR/virulence | `06.AMR/…` | ABRicate v1.2.0, EFSA thresholds ≥80% id / ≥70% cov |
| CARD read mapping | `06.AMR/AMR_mapping/…` | BBMap v39.33 vs CARD v4.0.1, minid 0.99, ≥70% covered length |
| Plasmid calls | `07.plasmids/…` incl. `verified plasmids` file | Platon v1.7 + BLAST verification |
| Prophage | `08.phages/…` | VirSorter2 v2.2.4 + CheckV v1.0.3 |
| Taxonomy | `04.taxonomy/…` | GTDB-Tk v2.6.1 (GTDB R226) |

### 1.3 Existing config keys already available
`bakta_db`, `blast_db`, `eggnog_db`, `gtdbtk_db`, `platon_db`, plus `resources.threads`, `resources.ram_gb`.

### 1.4 Conventions already in use in BacFlux (follow them)
- **Audit every filtering decision to a TSV** (see `contig_taxonomy_decisions.tsv`). The mobilome module must do the same.
- Databases that are large or license-encumbered are **user-downloaded**, path given in `config.yaml`, documented in README.
- Modules degrade gracefully rather than hard-failing.

### 1.5 Commands to establish ground truth (run these first)

```bash
cd /path/to/BacFlux
grep -n "^rule " workflow/Snakefile workflow/FastaFlux            # rule inventory
sed -n '/^rule bakta/,/^rule /p'      workflow/Snakefile          # Bakta outputs (.faa .gff3 .tsv)
sed -n '/^rule amr_contigs/,/^rule /p' workflow/Snakefile         # ABRicate: dbs, wildcards, thresholds
sed -n '/^rule platon/,/^rule /p'     workflow/Snakefile          # plasmid outputs
sed -n '/^rule gtdbtk/,/^rule /p'     workflow/Snakefile          # taxonomy summary path
grep -n "selected" workflow/Snakefile | head -40                  # exact selected-contigs filename
ls workflow/envs/                                                 # env naming convention
sed -n '1,60p' config/config.yaml                                 # key names + nesting style
```
Paste the results into §1.2 and delete the warning above.

*(Done, but not here: this was carried out against the real Snakefile and written
up in [`mobilome_wpA_ground_truth.md`](mobilome_wpA_ground_truth.md) instead of
being pasted back into §1.2. Note also that these commands are themselves v1 —
in v2 the rules live in `workflow/rules/{shared,illumina,nanopore,hybrid,contigs}/*.smk`
rather than in one `workflow/Snakefile`. At the time this was written the config to
read was `config/config_v2.yaml`; the v1/v2 cutover has since happened and it is
simply `config/config.yaml` now — see `unification_migration_plan.md` Stage 1.)*

---

## 2. Core conceptual model

This section is the reasoning that took the longest to establish. Do not re-litigate it.

### 2.1 Two different meanings of "reference" — do not conflate
- **IS reference database** (ISfinder / ISOSDB / TnCentral): a catalogue of element sequences. Used to *name* things.
- **Reference genome**: a second genome acting as a coordinate system representing the state *without* an insertion.

BacFlux takes arbitrary user isolates with **no comparator genome**. Therefore all tools that require an external reference genome are **out of scope**: IS_mapper, MGEfinder, panISa, ISCompare, ISdetector. They are excluded because they need a *comparison*, not because they need a database.

### 2.2 The three-layer detection philosophy
Sensitivity increases, specificity/nameability decreases, going down:
1. **Nucleotide homology** to named elements → precise, gives a citable name, blind to novelty.
2. **Protein / profile HMM** of catalytic domains → generalises to divergent elements, gives family not name.
3. **Structural signature** (TIRs, direct repeats, copy number ≥2, att sites) → most sensitive, no identity.

Best practice: **HMM finds candidates → structure validates → nucleotide names.**

Empirical anchor: in the ISOSDB paper, **97.5%** of ISOSDB transposases had protein homologs in ISfinder but only **37.9%** had nucleotide homologs. Protein search is ~2.5× more sensitive at recovering relationship to known elements.

Reality check: published benchmarks put IS-detection FDR at **8–24%** even on curated datasets. Report tiered confidence, never a bare count.

### 2.3 Short-read fragmentation is the central constraint
IS elements are the primary cause of contig breaks. On SPAdes assemblies:
- Multi-copy IS collapse → located IS counts are a **floor, not a count**.
- An AMR gene and its flanking IS usually land on **different contigs** — i.e. the exact structure being detected is what destroyed the assembly.

Consequences (hard rules):
- Composite-transposon and ICE-boundary calling → **BacFluxL only**.
- BacFlux (short-read) reports the weaker, honest signals: distance to contig end, IS-at-contig-boundary flag, read-depth-derived copy number.

> **Revised in v2.0.0.** The first bullet is no longer how this is enforced.
> There is no separate long-read repository to restrict things to — v2 is one
> workflow with four entry points, and composite-transposon and ICE-boundary
> calling run in **all four**, including `illumina`. Gating a whole capability on
> the sequencing technology turned out to be the wrong instrument: what actually
> degrades the call is *assembly fragmentation*, which is correlated with the
> technology but not determined by it (a good short-read assembly of a low-repeat
> genome can beat a poor long-read one).
>
> What v2 does instead is report the fragmentation directly, per call, so the
> reader can judge: `dist_to_contig_end`, `at_contig_boundary`, `spans_contigs`,
> and a hard rule that **anything spanning contigs is capped at `low` confidence**
> regardless of how good its machinery looks. The second bullet — report the
> weaker, honest signals — is therefore the part that survived, and it now applies
> everywhere rather than only in short-read mode.

### 2.4 IS vs transposon vs ICE — the definitions drive the tool choice
- **IS** = encodes only what it needs to move. By definition carries **no** passenger genes. ISfinder covers this.
- **Transposon** (composite or unit/Tn3-family) = carries passenger genes incl. AMR. **TnCentral** covers this, not ISfinder.
- **ICE/IME** = integrates into the chromosome **and** encodes conjugation machinery → self-transmissible. Breaks the naive "chromosomal ⇒ non-transferable" assumption. ICEberg/CONJscan cover this.

There is no hidden ISfinder feature that tells you an IS carries AMR. The analysis must be built.

### 2.5 The deliverable: a mobility ladder per AMR gene
1. Chromosomal, no MGE context → intrinsic candidate
2. IS adjacent/upstream, oriented at the gene → expression modulation (hybrid promoter), not mobilisation
3. Inside IS-flanked composite → mobilisable within the cell
4. Inside unit transposon / integron cassette → mobilisable, named architecture
5. On a mobilizable plasmid → transferable with helper
6. **Inside an ICE, or on a conjugative plasmid → predicted self-transmissible**

Special case: IS inserted *inside* an AMR CDS → likely inactivation. Report separately; do not pollute the mobilisation count.

### 2.6 Regulatory framing (EFSA)
WGS is mandatory for strain characterisation (FEEDAP 2018; extended to all food-chain microorganisms by the 2024 EFSA statement). The regulatory hook is the **intrinsic vs acquired** distinction — intrinsic ≈ species-wide, chromosomal, non-transferable; acquired ≈ horizontally gained via plasmid/transposon.

EFSA does not literally require "report the replicon and flanking IS". It requires evidencing whether a determinant is acquired and transferable. This module provides that evidence.

**Language discipline:** always write "**predicted** self-transmissible", never "transmissible". Never label output "EFSA-compliant"; frame as "supporting evidence for intrinsic/acquired classification". Confirmatory experiment is a filter/broth mating assay — say so in the docs.

---

## 3. Tool decisions

### 3.1 ADOPT

| Tool | Role | Install | Notes |
|---|---|---|---|
| **AMRFinderPlus** | primary AMR input for mobility module | **already present** in Bakta conda env; DB at `{bakta_db}/amrfinderplus-db/` | see §4 — this is the cheapest win |
| **ISEScan** | IS detection on contigs | bioconda | self-contained pHMMs, no external DB, works on drafts, flags complete vs partial |
| **CONJscan** (MacSyFinder module) | conjugation machinery → ICE/IME mobility class | bioconda + `macsydata install` | conservative, model-based, no giant DB. **Added 2026-07-30:** an optional second model set, **ICEscan** (from ICEfinder2, CC BY-NC-SA, default off), now runs alongside it and is UNIONed — it supplies the IME and AICE classes CONJScan 2.1.0 has no model for. It is a *fork* of CONJScan, not an independent tool, and four of its integrase profiles are deliberately distrusted. Full rationale, measurements and limitations: `docs/methods_icescan_union.md`. |
| **ISOSDB** | redistributable IS nucleotide DB | from pseudoR repo | read-mapping copy number + sensitive detection |
| **TnCentral (± ISfinder ± Integrall)** | naming + transposon/integron layer | direct download, see §5 | |
| **geNomad** | phage caller + 2nd-opinion plasmid caller alongside Platon — **role under review, see licensing note** | bioconda; **⚠ LICENSE: Berkeley Lab ACADEMIC / NON-COMMERCIAL USE ONLY** (verified from the raw LICENSE 2026-07-22 — NOT BSD-4-Clause; the bioconda recipe's `BSD-4-Clause` tag is wrong). Clause 4: "NON-COMMERCIAL USE, purposes ONLY. User must be an accredited academic institution." Commercial use requires a separate LBNL license. | **Licensing correction 2026-07-22 (supersedes the 2026-07-21 "default/mandatory geNomad" framing of D8/D9).** BacFlux is MIT and the project's hard rule (§11) is to impose no non-commercial restriction on downstream users — so geNomad CANNOT be the default/mandatory caller: that would force every user (incl. commercial) into a non-commercial tool. BacFlux invokes geNomad (does not vendor its code), so BacFlux's own MIT is not contaminated, but a mandatory geNomad step is not commercially usable. **Design decision pending (see unification plan D8/D9 revision):** make geNomad OPT-IN (default off) with a clear non-commercial README notice — same pattern BacFlux already uses for licence-encumbered databases (§11) — with VirSorter2 (permissive) as the default phage caller and the D9 plasmid concordance degrading to Platon-only when geNomad is off. Technical merits still hold (v1.12.0 actively maintained, 97.3% precision vs VirSorter2 94.7%, one run does viruses+plasmids); the licence, not the science, constrains the role. |

### 3.2 CONSIDER / OPTIONAL
| Tool | When | Notes |
|---|---|---|
| **digIS** | BacFluxL, opt-in | novelty discovery; heavier deps, ~2021 vintage |
| **MobileElementFinder** | BacFluxL, cross-check only | native composite-transposon flagging; reference-bound, DTU CGE. **NOT the same tool as MGEfinder** |
| **MOB-suite (`mob_recon`/`mob_typer`)** | if plasmid mobility typing wanted | gives conjugative/mobilizable/non-mobilizable — directly maps to ladder tiers 5–6 |
| **hAMRonization** | if ABRicate + AMRFinderPlus outputs need reconciling | PHA4GE standard |

### 3.3 REJECT (and why — do not revisit)
| Tool | Reason |
|---|---|
| IS_mapper, MGEfinder, panISa, ISCompare, ISdetector | require an external reference genome / comparative design |
| MGEfinder specifically | only detects insertions *absent* from the reference → F1 = 0.01 on *S. sonnei* benchmark |
| pseudoR | built for fragmented metagenome assemblies; on clonal isolates it measures within-population heterogeneity, not inventory. Heavy R+Bowtie2+blastn+mosdepth stack. **Still needed as the ISOSDB download source** |
| mobileOG-db | protein-level MGE context only; duplicates Bakta; manual DB download |
| ICEfinder 2.0 / icefinder-opt | no installation section, no conda recipe, hard-coded `config.ini` paths, CentOS 7. **vmatch** dependency has an UNVERIFIABLE licence → cannot be redistributed with an MIT workflow. ⚠ **Corrected 2026-07-28:** the earlier claim that vmatch "is not on bioconda" was WRONG — `bioconda/vmatch 2.3.1` exists (linux-64, osx-64). The rejection stands on the licence alone: the recipe declares `license: Unknown / OTHER` and vmatch.de is unreachable, so its terms cannot be established, which §11 makes a blocker. See `docs/methods_att_and_small_plasmids.md`. NOTE ICEfinder2 is CC BY-NC-SA 4.0 — its source may be READ to establish an algorithm (§11 permits this) but never copied. |
| EBI mobilome-annotation-pipeline as a dependency | Nextflow, metagenome-oriented. **Use as escape hatch / design reference only.** ⚠ **Corrected 2026-07-30:** the claim that it "currently does not run AMRFinderPlus" was WRONG — see the correction at §11. The rejection stands on the remaining grounds (Nextflow runtime, metagenome/MAG orientation), not on missing AMR calling. |

---

## 4. Work package A — surface AMRFinderPlus (do this first, ~half a day)

**Rationale:** Bakta already calls AMRFinderPlus internally but only surfaces gene name/product. The full report adds five things the mobility module needs:
1. `Method` column (EXACTX / BLASTX / PARTIALX / HMM / POINTX…) → a free confidence tier
2. Partial-at-contig-end flags → directly addresses the fragmentation problem
3. Element type / subtype (AMR / STRESS / VIRULENCE)
4. Class / subclass
5. Point mutations (with `--organism`) → the *intrinsic, non-transferable* category ABRicate structurally cannot see

**No new conda env, no new database download.**

```python
rule amrfinderplus:
    input:
        faa = f"{OUT}/05.annotation/bakta/{{sample}}/{{sample}}.faa",
        gff = f"{OUT}/05.annotation/bakta/{{sample}}/{{sample}}.gff3",
        fna = f"{OUT}/02.assembly/{{sample}}/selected_contigs.fasta",
    output:
        tsv = f"{OUT}/06.AMR/amrfinderplus/{{sample}}_amrfinder.tsv",
    params:
        db = lambda w: os.path.join(config["bakta_db"], "amrfinderplus-db"),
        organism = lambda w: amrfinder_organism(w.sample),   # see below; "" if unmapped
    conda: "../envs/bakta.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        amrfinder -p {input.faa} -n {input.fna} -g {input.gff} \
          --annotation_format bakta \
          -d {params.db} --plus --threads {threads} \
          {params.organism} -o {output.tsv}
        """
```

**Verify:** that `--annotation_format bakta` is accepted by the amrfinder version pinned in the Bakta env. If not, drop `-g` and run `-p` + `-n` only (still yields contig coordinates).

**`--organism` from GTDB-Tk:** BacFlux already runs GTDB-Tk. Build a small YAML/TSV lookup mapping GTDB species → AMRFinderPlus curated organism names (~30 organisms, NCBI nomenclature, which does not always match GTDB). Emit `--organism X` when mapped, empty string otherwise. This unlocks point-mutation detection. Fallback must be silent and graceful.

**Do not remove ABRicate.** Three complementary legs:
- BBMap→CARD (reads): assembly-collapse-immune
- ABRicate (contigs): multi-DB breadth + existing EFSA-threshold reporting
- AMRFinderPlus (contigs): structured input to the mobility module

Keep EFSA thresholds on the ABRicate leg only — applying a blanket 80/70 cutoff would fight AMRFinderPlus's curated per-gene cutoffs.

---

## 5. Work package B — databases

### 5.1 Config additions
```yaml
mobilome:
  run: false                       # master switch, default OFF
  isescan:
    remove_short_is: false         # keep partials; tier them ourselves

  # Layer 1 — redistributable, always on when mobilome.run
  isosdb_url: "<from pseudoR repo: https://github.com/joshuakirsch/pseudoR>"

  # Layer 2 — no ISfinder terms attached, safe default
  tncentral_nt_url: "https://tncentral.ncc.unesp.br/api/download_blast/nc/tn"

  # Layer 3 — opt-in, ISfinder terms apply, user must set
  tncentral_isfinder_nt_url:   null   # .../api/download_blast/nc/tn_in_is
  tncentral_isfinder_prot_url: null   # .../api/download_blast/prot/tn_is

  ice:
    run: false
    window_bp: 15000               # anchor clustering gap
    min_element_bp: 8000
    max_element_bp: 500000
  composite:
    max_span_bp: 20000             # configurable, not biology
```

### 5.2 Licensing rules (important — BacFlux is MIT)
- **Never vendor** TnCentral/ISfinder files in the repo. Fetch at install time via a `download_db` rule, document terms in README.
- ISfinder's stated terms require written authorisation to download and forbid redistribution; TnCentral offers direct download links but carries an "All Rights Reserved" notice. **Open question:** commercial use is unverified — email TnCentral/ISfinder contacts before making claims in the README.
- ISOSDB is openly licensed → safe to bundle or auto-fetch.
- The `thanhleviet/ISfinder-sequences` GitHub scrape must **not** be used — unofficial, ~2018 vintage, unauthorised redistribution.

### 5.3 BLAST DB handling

> ✅ **IMPLEMENTED 2026-07-28** as rules `tncentral_db` and `iceberg_db`
> (`80_mobilome.smk`). The paragraph below was written from the TnCentral web
> page and is **wrong about the format** — corrected inline. Everything here was
> re-established by actually fetching the endpoints; see §5.5.

TnCentral DBs are built with makeblastdb 2.6.0+ (**v4 format**) and the endpoints are **undated**. Always:
```bash
blastdbcmd -db <db> -entry all -out <db>.fasta
makeblastdb -in <db>.fasta -dbtype nucl -out <db>_v5 -blastdb_version 5
sha256sum <db>.fasta > <db>.sha256      # record checksum + fetch date in provenance
```

**Correction:** the `nc/tn` endpoint does not serve a bare BLAST database. It
serves a **ZIP** containing `tncentral.fa` (a plain FASTA) *plus* a v4 index
alongside it. So no `blastdbcmd` round-trip is needed — unzip, discard the
shipped `.n*` index files, and `makeblastdb` the FASTA directly as v5. The ZIP
members are dated (2025-05-16 in the release fetched), which is a better
provenance anchor than the spec assumed existed.

### 5.5 Endpoint ground truth (established 2026-07-28 by fetching them)

The spec's URLs were **not blocked by licensing** — they were stale or
mislabelled, which is why this layer sat unbuilt for so long.

| Claim in this spec | Reality |
|---|---|
| TnCentral endpoint returns 403 | **Bot block on curl's default user-agent.** With a browser UA it serves 5.47 MB. The download rules set one. |
| TnCentral is a BLAST v4 database | **ZIP holding `tncentral.fa`**, 512 sequences, plus a v4 index |
| Three TnCentral endpoints | **Six.** `nc/tn` = TnCentral only (the safe default, and the one wired in); `nc/tn_in` adds Integrall; the `*_is` variants add ISfinder content and its terms |
| ICEberg at `bioinfo-mml.sjtu.edu.cn` | **404.** ICEberg 3.0 lives at `tool2-mml.sjtu.edu.cn/ICEberg3/`; the old host still serves ICEberg **2.0**, which is how the stale URL kept looking plausible |

**⚠ Both nucleotide FASTAs are MALFORMED and must be repaired before indexing.**
Some deflines are glued onto the *end* of a sequence line instead of starting
their own:

```
...gtgcagccgtcttctgaaaacgaca>In1223-KX784502
```

| | records at line start | `>` characters | broken |
|---|---|---|---|
| TnCentral `nc/tn` (2025-05-16) | 512 | 533 | **21** |
| ICEberg `ICE_seq_all` + `IME_seq_all` (2023-06-01) | 1773 | 1774 | **1** |
| ISOSDB V3 | 22713 | 22713 | clean |

The damage runs **both ways**, and the second direction is the dangerous one:

- the affected elements are **invisible** to `makeblastdb` — for TnCentral that
  includes Tn*7* itself and eleven integrons, the very class tier 4 exists to name;
- and the records they were glued to become **chimeric**, absorbing the defline
  text plus the next element's sequence. `In_Tn6162` measured **41,492 bp instead
  of 8,911** — 4.7× its true length. Because the naming cascade tests coverage as
  `alignment / slen`, an inflated `slen` makes those elements almost impossible
  to name, and it fails *silently*.

The shipped v4 index has the same 21 records missing, since it was built from the
unrepaired file — an independent reason to rebuild rather than reuse it. Both
`download_db` rules now split embedded deflines onto their own lines and then
**assert** that every `>` begins a line, failing rather than indexing a file whose
shape they do not understand. Repairing TnCentral recovered `In781_p` on the
positive control (99.8% identity, 82% coverage) — a real integron that the
unrepaired database could not see. (The note originally named that control
"KPNIH1"; this project has used both *K. pneumoniae* KPNIH1 (`CP008827.1`) and
ATCC BAA-2146 (`CP006659.2`) as positive controls and repeatedly conflated them,
so the strain is left unstated rather than guessed.)

(The `download_full/fa` endpoint ships 513 well-formed per-element files and is a
possible alternative source, but it is a year older — 2024-05-16 versus
2025-05-16 — and carries 20 fewer elements, so repairing the newer file is
preferred.)

**Deflines** — the parsers are written against these, so they are recorded verbatim:

```
TnCentral   >Tn4401b-JX560992              <NAME>-<ACCESSION>, split on the LAST hyphen
            >IS1133_Tn10_IS903B-CP000602.1  (names contain hyphens and underscores)
            >Tn7246-                        (accession may be empty)

ICEberg     >ICEberg|1174|ICEKpnATCCBAA-2146-1|GenBank|CP006659.2|4603840..4661887 ...
            pipe-delimited; field 2 = element name, field 4 = source accession
```

**Licensing, re-checked verbatim (still: never vendor, URL only):**
- TnCentral: *"© TnCentral 2024 - All Rights Reserved"*. There is **no terms page at all**, so §5.2's "commercial use is unverified" is not merely unresolved — the site does not address it. The citation request is Ross *et al.* 2021, mBio.
- ICEberg 3.0 (released 2023-06-01): **no licence, terms or reuse statement anywhere.** Every page carries only *"Copyright © 2023 All Rights Reserved by Microbial Bioinformatics Group in MML, SJTU."*
- ISOSDB: the pseudoR repo is **MIT**, confirming §5.2. `ISOSDB.V3.fna.zip` (22,713 sequences) and `IS_fam_annot.txt` (family per element) are fetchable straight from the repo.

Both naming layers are therefore **off unless a URL is configured**, and each
writes a `PROVENANCE.txt` recording source, fetch date, checksum and sequence
count — because neither endpoint is versioned, and without that there is no way
to say later which release a result came from.
The dumped FASTA is also the input for minimap2/diamond legs.

### 5.4 Which DB for which operation
| Operation | Input | DB | Answers |
|---|---|---|---|
| map (`bwa mem`/`minimap2`) | reads | ISOSDB nt | copy number (collapse-immune) |
| blastn | ISEScan element seqs | TnCentral+ISfinder nt | element name, family |
| blastp/diamond | ISEScan transposase ORFs | TnCentral+ISfinder prot | family when nt fails |
| blastn | whole contigs | full TnCentral (± Integrall) | composite/Tn3 transposons, integrons |

**Naming cascade:** blastn ≥90% id over ≥80% length → named element. Unassigned → protein search → family only. Still nothing → structure only.

---

## 6. Work package C — IS detection

### BacFlux (short-read)
1. `ISEScan` on selected contigs. Run **without** `--removeShortIS`; tier partials ourselves.
2. QC metric: fraction of IS hits within *N* bp of a contig end.
3. Read-mapping leg against ISOSDB → per-element copy number = element depth / median assembly depth. **Report the delta** vs ISEScan's located count: this quantifies assembler collapse and makes the located inventory honestly a floor. (Structurally the same rule pattern as the existing BBMap→CARD leg.)
   - Dereplicate DB (ISOSDB already CD-HIT 95%), filter on MAPQ + identity, handle multi-mapping explicitly, report at family **and** element level.

### BacFluxL/Lplus
Same core, plus per-replicon IS burden (chromosome vs each plasmid), composite-transposon detection, IS×CDS intersection for pseudogene candidates. Optional: digIS (novelty), MobileElementFinder (cross-check).

*In v2 this heading means "`nanopore` and `hybrid` mode", and the split is softer
than it looks: per-replicon IS burden, composite-transposon detection and the
IS×CDS intersection all run in every mode. What genuinely differs is how much you
should believe them, which is why every call carries its contig-edge flags. The
one leg that really is mode-restricted is the ISOSDB read-mapping copy-number
step — it needs reads, so it is `illumina` and `hybrid` only. digIS and
MobileElementFinder were not adopted.*

---

## 7. Work package D — AMR × MGE co-localisation

Inputs: AMRFinderPlus TSV, ISEScan GFF, TnCentral contig BLAST, Platon/MOB-suite replicon call, Bakta GFF3.

```bash
bedtools window   -a amr.bed -b is.bed -w {composite.max_span_bp} > pairs.tsv
bedtools closest  -a amr.bed -b is.bed -d -D a                    > nearest.tsv
bedtools intersect -a amr.bed -b is.bed -f 0.5 -wo                > disrupted.tsv
```

Composite call = short parser over `pairs.tsv`: same contig, ≥2 IS hits, **same family**, **same orientation**, AMR gene between them, span ≤ threshold.

**Pitfalls (hard-code these):**
- **IS26 / IS6 family breaks the orientation rule** — forms translocatable units with copies in *direct* orientation. Needs an explicit exception or the single most clinically important architecture is missed.
- A **TnCentral named hit should override** the pattern-based call — the curated entry has the true architecture and gets IS26 right for free.
- Distance thresholds are convention, not biology. Always report the actual distance alongside the tier.

---

## 8. Work package E — ICE-lite (~9–12 working days)

*Was "BacFluxL only". **Built, and it runs in all four v2 modes** — a hybrid dry
run prints "Mobilome module: ON" like any other. Phases 0–6 are implemented in
`workflow/scripts/80_mobilome/` (`att_search.py`, `conjscan_to_ice.py`,
`colocalise.py`); Phase 7 was carried out and is written up in the README's
Validation section and in `methods_ebi_comparison.md`. The single biggest
departure from the plan below is Phase 3 — see the superseded-probe box in it.*

Design principle: **build an evidence integrator, not an ICE finder.** Only genuinely new algorithm is att-site search.

**Phase 0 — data contract (½ d).** One loader per input returning a tidy frame: `contig, start, end, strand, type, label`. Everything downstream uses this shape.

**Phase 1 — anchors (1 d).** Extract into one table:
- relaxase, T4CP, T4SS/MPF ← CONJscan
- integrase ← Bakta product regex (`tyrosine recombinase|Phage_integrase|XerC/D|serine recombinase`)
→ `{sample}_ice_anchors.tsv`

**Phase 2 — candidate seeding (1 d).** Cluster anchors on the same contig within `window_bp`. Keep clusters with ≥ relaxase or T4SS. Provisional interval = min→max anchor, snapped to CDS boundaries. Drop < `min_element_bp`.

**Phase 3 — att-site detection (3–4 d).** The only real algorithm.
- *Mode A, tRNA-anchored (high precision):* ICEs integrate at tRNA 3′ ends. Probe = last 15–25 bp of a nearby tRNA (strand-aware); search for a second copy (≤1 mismatch) on the opposite side within window → attL/attR.
  > ⚠ **SUPERSEDED 2026-07-28 — a FIXED probe length cannot work.** The ICE*Kp*
  > direct repeat is **17 bp** (Lam *et al.* 2018), so a 25 bp probe is
  > structurally incapable of finding it; measured on two clinical *K. pneumoniae*
  > isolates, nothing was found at 25 bp and pairs appeared at 18. Bacterial ICE
  > repeats span ~10–60 bp, so no single length is defensible. Every reference
  > tool uses a variable-length search: ICEfinder2 runs `vmatch -l 15`, DEPhT and
  > DBSCAN-SWA BLAST the two flanks against each other. BacFlux now does the
  > latter (`blastn -task blastn-short -dust no`, ≥15 bp, ranked by bitscore),
  > with tRNA proximity demoted to a scoring bonus. Full evidence, parameters and
  > citations: `docs/methods_att_and_small_plasmids.md`.
- *Mode B, de novo:* k-mer dict of the left flank window (k=25 stepping down to 12), scan right flank for same-orientation matches; rank by length desc, resulting interval within [min,max]_element_bp, symmetry of flank distance.
- Use `pyfaidx`/Biopython; a plain dict k-mer index is fast enough at these window sizes — **this is why vmatch is not needed**.
- **Critical:** mask the interval against ISEScan calls before the de novo scan, or IS terminal repeats flood the candidate pairs. This is the most likely failure mode.
- Output `attL, attR, att_seq, tRNA, boundary_method ∈ {tRNA, denovo, none}`.

**Phase 4 — classification (½ d).** Pure rules:

| Integrase | Relaxase | T4SS | Class | Mobility |
|:-:|:-:|:-:|---|---|
| ✓ | ✓ | ✓ | ICE | predicted self-transmissible |
| ✓ | ✓ | ✗ | IME | mobilisable (needs helper) |
| ✓ | ✗ | ✗ | CIME / island | passive |
| ✗ | ✓ | ✓ | conjugative region, unbounded | report, don't call ICE |

Plus truncation check: relaxase or VirB4 with HMM alignment covering <70% of the profile → `degraded`, downgrade mobility. Decayed ICEs are common; this is where naive tools overcall.

**Phase 5 — cargo & context (1 d).** Intersect final intervals with AMRFinderPlus, ISEScan, Platon/MOB-suite, Bakta CDS.

**Phase 6 — confidence flags (1 d).** Per element: `n_anchor_classes`, `boundary_method`, `machinery_intact`, `spans_contigs`, `dist_to_contig_end`, `confidence`. Rule: `high` = 4 anchor classes + tRNA-anchored + single contig + intact. **Anything spanning contigs is capped at `low` regardless.**

**Phase 7 — validation (3–4 d, DO NOT SKIP).**
- Positive controls: *V. cholerae* SXT/R391; *K. pneumoniae* ICE*Kp*; *E. faecalis* Tn*916* (Gram+, tetM); *E. coli* ICE*Ec2*. Score element recall and boundary bp offset vs published coordinates.
- Negative control: genomes with no reported ICE.
- Cross-check against the **ICEfinder 2.0 web server** — benchmark with zero installation.
- Unit tests on synthetic contigs with planted att repeats at known offsets.

**Sequencing note:** assemble the Phase 7 positive controls *before* starting Phase 3.

**Usable v0 after ~3 days** = Phases 0–2 + 4, no boundaries. Already answers "is there conjugation machinery and what class" — tiers 5–6 of the ladder.

---

## 9. Output schema

Adopt EBI Mobilome Annotation Pipeline conventions (facts/design, not code — see §11):

**ID format:** `contig_id|mge_type-start:end`

**Sequence Ontology terms for the GFF `Type` column:**
| Type | SO term |
|---|---|
| `insertion_sequence` | SO:0000973 |
| `inverted_repeat_element` | SO:0000481 |
| `direct_repeat` | SO:0000314 |
| `conjugative_integron` (= ICE) | SO:0000371 |
| `integron` | SO:0000365 |
| `plasmid` | SO:0000155 |
| `prophage` | SO:0001006 |

**Thresholds borrowed from MAP:** discard MGEs <500 bp; discard predictions with no CDS; assign a CDS to an MGE at ≥0.9 coverage of the boundaries.

**Discard-with-reason file** — mandatory, matches BacFlux's existing `contig_taxonomy_decisions.tsv` philosophy:
`{sample}_discarded_mge.tsv` with an explicit reason per dropped prediction (`mge<500bp`, `no_cds`, `spans_contigs`, `machinery_degraded`, …).

**Main deliverable table** — one row per AMR gene:
```
sample, contig, replicon(chromosome|plasmid_id), amr_gene, amr_class, amr_subclass,
amrfinder_method, amr_partial_at_contig_end,
mge_context(none|is_adjacent|composite|unit_transposon|integron|ice|ime|plasmid),
mge_id, mge_name(if TnCentral/ICEberg match), distance_bp, orientation,
is_family, n_flanking_is, same_orientation,
relaxase_type, mpf_type, machinery_intact,
attL, attR, boundary_method, spans_contigs,
mobility_tier(1-6), confidence(high|medium|low)
```

Example report line (target quality):
> `blaCTX-M-15` — chromosomal, contained within predicted ICE (contig 3, 84.2 kb, single contig); integrase + MOBF relaxase + T4CP + VirB4 all intact; flanking DRs at tRNA-Gly; 94% id to ICEberg ICE*Ec2*. **Confidence: high. Predicted self-transmissible.**

---

## 10. Snakemake structure

> ⚠ **The path below is wrong; the layout was planned before v2's rule directories
> existed.** The module is implemented in **`workflow/rules/shared/80_mobilome.smk`**
> (22 rules), alongside the other shared stages
> (`00_common`, `10_decontam`, `20_qc`, `30_taxonomy`, `40_annotation`, `50_amr`,
> `60_plasmid`, `70_phage`, `90_report`). The helper package landed where planned,
> at `workflow/scripts/80_mobilome/`, though with different file names — the real ones
> are `att_search.py`, `conjscan_to_ice.py`, `colocalise.py`, `isescan_to_table.py`
> and `name_elements.py`, plus their `test_*.py` siblings.

```
workflow/rules/mobilome.smk        # ACTUAL: workflow/rules/shared/80_mobilome.smk
    isescan
    isescan_name_tncentral       # blastn + blastp cascade
    isosdb_read_mapping          # copy number
    amrfinderplus                # WP-A
    conjscan
    ice_anchors                  # runs in all four v2 modes; see the §2.3 revision
    ice_candidates
    ice_boundaries               # script: att_search.py
    ice_classify
    amr_mge_colocalisation
    mobilome_report

workflow/scripts/80_mobilome/       # importable package + thin CLI, pytest-able outside Snakemake
    __init__.py  loaders.py  att_search.py  classify.py  colocalise.py  report.py
```
Gate everything behind `config["mobilome"]["run"]` (default `false`) and `config["mobilome"]["ice"]["run"]`.

---

## 11. Attribution & licensing (read before copying anything)

**HARD CONSTRAINT:** EBI's `mobilome-annotation-pipeline` names three scripts as ICEfinder2-derived and licenses those algorithms under **CC BY-NC-SA 4.0** (their modifications: Apache-2.0):
- `bin/ice_boundary_refinement.py`
- `bin/map_tools/icefinder_process.py`
- `bin/prescan_to_fasta.py`

CC BY-NC-SA is **non-commercial + share-alike → incompatible with MIT**. Copying or adapting that code contaminates BacFlux's licence and imposes a non-commercial restriction on downstream users.

**Two separate issues — do not conflate them:**

*Databases (TnCentral/ISfinder):* solved by not shipping the data. The workflow distributes **code and URLs only**; the user downloads under their own agreement with the licensor. Default-off + README notice is good practice and matches how BacFlux already handles `bakta_db`, `blast_db`, `eggnog_db`, `gtdbtk_db`, `platon_db` and the CARD link.

*Code (CC BY-NC-SA):* **default-off does NOT cure this.** Share-alike attaches to distributed source, not to execution. Code copied into the repo is licence-contaminated whether or not a config flag ever runs it. The only remedies are: don't copy it, or relicense the module. This project takes the first route.

**Permitted:** adopting conventions (SO terms, ID format, discard-reason pattern, 500 bp / 0.9 coverage thresholds) — these are facts and design decisions, not expression. Reading the scripts to understand the algorithm and then implementing independently is also fine, but keep the separation real.

**CITATIONS.md structure:**
- *Tools & databases used:* ISEScan (Xie & Tang 2017), CONJscan/MacSyFinder, AMRFinderPlus (Feldgarden et al. 2021), ABRicate, TnCentral (Ross et al. 2021), ISfinder (Siguier et al. 2006), ISOSDB/pseudoR (Kirsch et al. 2024, Cell Host Microbe), ICEberg 3.0 (Wang et al. 2024, NAR), ICEfinder 2.0
- *Design influence & further analysis:* EBI Mobilome Annotation Pipeline (Apache-2.0, with its own ICEfinder2 attribution) — frame as **complementary scope**, not "simplified alternative": MAP goes deep on the mobilome for metagenomes/MAGs; BacFlux covers assembly→annotation→AMR→plasmids→mobility for single isolates.

> ⚠ **CORRECTED 2026-07-30 — the paragraph below is FALSE. Do not put it in the README.**
> It claimed MAP "deliberately does not run AMRFinderPlus (gene-level AMR/virulence
> subworkflow in development)". Established by reading the pipeline at revision
> `3aa408d` (2026-07-23): **MAP does run AMRFinderPlus, on by default.**
> `subworkflows/ebi-metagenomics/amr_annotation/main.nf` includes the nf-core
> `AMRFINDERPLUS_RUN` module (alongside DeepARG, RGI/CARD, hAMRonization and an
> `AMRINTEGRATOR` that joins AMR calls to the mobilome). It is included at
> `workflows/mobilomeannotation.nf:39`, invoked at line 385, has publishDir rules in
> `conf/modules.config`, ships test fixtures, and is opt-**out** via
> `--skip_amrfinderplus`. AMRFinderPlus is also cited in their `CITATIONS.md`.
>
> **What this changes.** The "we are not duplicating their roadmap" argument is dead —
> MAP already does AMR × mobilome integration. The *complementary scope* framing in the
> bullet above still stands, but it now rests solely on the real distinction: MAP goes
> deep on the mobilome for **metagenomes and MAGs**; BacFlux covers
> assembly → annotation → AMR → plasmids → mobility for **single isolates**. Frame it
> that way and nothing else. See `CITATIONS.md`.
>
> The §11 HARD CONSTRAINT above is **unaffected** — re-verified at the same revision:
> `bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py` and
> `bin/prescan_to_fasta.py` are all still present and still carry the header
> *"Original ICEfinder2 work licensed under CC BY-NC-SA 4.0"*. Still never copy them.

~~Note in the README that MAP's current release deliberately does not run AMRFinderPlus (gene-level AMR/virulence subworkflow in development), so this module is not duplicating their roadmap.~~

**Escape hatch to document for users:** "for full ICE delimitation and deeper mobilome analysis, run the EBI Mobilome Annotation Pipeline; BacFluxL writes `{sample}_contigs.fna` and `{sample}.gbk` ready for it." Consuming its `mobilome.gff.gz` output carries no licence contamination.

---

## 12. Open questions / decisions not yet made
1. **`--annotation_format bakta`** support in the pinned amrfinder version — verify empirically.
2. **GTDB → AMRFinderPlus organism mapping table** — needs building; ~30 curated organisms; GTDB names ≠ NCBI names in places.
3. **Does ISOSDB ship a protein FASTA?** Verify in the pseudoR repo. If not, run prodigal over the nucleotide FASTA once (nucleotide is the form actually needed for BacFlux's use cases anyway).
4. ~~**geNomad adoption**~~ — **RESOLVED 2026-07-21, REVISED 2026-07-22 (licensing)**, see §3.1. geNomad is an **OPT-IN** phage caller and 2nd-opinion plasmid caller; **VirSorter2 is the default**, because geNomad is licensed academic/non-commercial-only while BacFlux is MIT (§11). Both roles live in the always-on 06.plasmids / 07.phages stages, not inside the mobilome module — the module *consumes* their calls. Details and sourcing: `docs/unification_migration_plan.md` §4 D8/D9.
5. **PLSDB** — if reintroduced, use for *annotation* (which known plasmid group / replicon type), **not** classification. Interpret Mash **distance** (≈1−ANI) and **shared-hash fraction**, not the p-value (which is the probability of shared hashes by chance and collapses to ~0 for any real match). `mash dist` for one-contig-vs-DB with per-sequence processing; `mash screen` only for mixtures/read sets. Beware: plasmids and chromosomes share IS/Tn/AMR content, so a PLSDB hit on an IS-rich contig is expected and uninformative.

---

## 13. Implementation order (recommended)

1. **WP-A** AMRFinderPlus surfacing — ½ day, immediate value, zero new deps
2. **WP-B** database rules + checksums/provenance — 1 day
3. **WP-C** ISEScan + naming cascade + read-depth copy number — 2–3 days
4. **WP-D** AMR × IS co-localisation + mobility tiers 1–4 — 2–3 days
5. **CONJscan** rule (tiers 5–6, no boundaries) — ½ day
6. **WP-E** ICE-lite — 9–12 days, only if boundaries are actually wanted (planned as "BacFluxL only"; built, and runs in all four v2 modes)

Steps 1–5 deliver a complete, defensible mobility ladder. Step 6 is the expensive refinement.

*All six were completed for v2.0.0. Two notes on how the finished ladder differs
from the plan: step 4's "tiers 1–4" in practice reaches tier 3 without the
optional TnCentral layer, because **tier 4 requires a curated name** and that
layer is opt-in; and **tier 4 has never actually been assigned** in any retained
run, because on the one occasion the naming layer did match a curated transposon
the gene was on a conjugative plasmid and scored tier 6 instead. See the README's
"What the benchmark does not show".*
