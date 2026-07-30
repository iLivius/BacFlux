# Citations

BacFlux is glue. Almost all of the science it reports comes from other people's
tools and databases, and this file is where they are credited.

**If you publish results from BacFlux, cite the tools whose output you actually
used** — not just BacFlux. The mobilome module in particular is an integrator:
it decides how to combine other tools' calls, and every underlying detection
belongs to somebody else.

Two sections, following the structure set out in
`docs/mobilome_module_SPEC.md` §11:

1. **Tools and databases used** — code that runs, or data that is read.
2. **Design influence** — work we learned from and adopted decisions from,
   without using their code.

---

## 1. Tools and databases used

### Assembly, annotation and taxonomy

- **SPAdes** — Prjibelski A, Antipov D, Meleshko D, Lapidus A, Korobeynikov A. (2020)
  *Using SPAdes de novo assembler.* Current Protocols in Bioinformatics 70:e102.

- **Bakta** — Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. (2021)
  *Bakta: rapid and standardized annotation of bacterial genomes via alphanumerical identifiers.*
  Microbial Genomics 7(11):000685. <https://doi.org/10.1099/mgen.0.000685>

- **CheckM** — Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. (2015)
  *CheckM: assessing the quality of microbial genomes recovered from isolates,
  single cells, and metagenomes.* Genome Research 25:1043–1055.

- **GTDB-Tk** — Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH. (2022)
  *GTDB-Tk v2: memory friendly classification with the genome taxonomy database.*
  Bioinformatics 38(23):5315–5316.

- **GTDB** — Parks DH, Chuvochina M, Rinke C, Mussig AJ, Chaumeil P-A, Hugenholtz P. (2022)
  *GTDB: an ongoing census of bacterial and archaeal diversity through a
  phylogenetically consistent, rank normalized and complete genome-based taxonomy.*
  Nucleic Acids Research 50(D1):D785–D794.

- **BlobTools** — Laetsch DR, Blaxter ML. (2017)
  *BlobTools: Interrogation of genome assemblies.* F1000Research 6:1287.

### Antimicrobial resistance

- **AMRFinderPlus** — Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG,
  Haendiges J, Haft DH, et al. (2021) *AMRFinderPlus and the Reference Gene Catalog
  facilitate examination of the genomic links among antimicrobial resistance,
  stress response, and virulence.* Scientific Reports 11:12728.
  <https://doi.org/10.1038/s41598-021-91456-0>

- **ABRicate** — Seemann T. *ABRicate: mass screening of contigs for antimicrobial
  and virulence genes.* <https://github.com/tseemann/abricate>

- **CARD** — Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, et al. (2023)
  *CARD 2023: expanded curation, support for machine learning, and resistome
  prediction at the Comprehensive Antibiotic Resistance Database.*
  Nucleic Acids Research 51(D1):D690–D699.

- **BBMap / BBTools** — Bushnell B. *BBMap: A Fast, Accurate, Splice-Aware Aligner.*
  <https://sourceforge.net/projects/bbmap/>

### Mobile genetic elements

- **ISEScan** — Xie Z, Tang H. (2017) *ISEScan: automated identification of
  insertion sequence elements in prokaryotic genomes.* Bioinformatics 33(21):3340–3347.
  <https://doi.org/10.1093/bioinformatics/btx433>

- **MacSyFinder v2** — Néron B, Denise R, Coluzzi C, Touchon M, Rocha EPC, Abby SS. (2023)
  *MacSyFinder v2: Improved modelling and search engine to identify molecular
  systems in genomes.* Peer Community Journal 3:e28.

- **CONJScan models** — cite all three; they are the models BacFlux actually runs:
  > Coluzzi C, Garcillán-Barcia MP, de la Cruz F, Rocha EPC. (2022)
  > *Evolution of plasmid mobility: origin and fate of conjugative and
  > non-conjugative plasmids.* Molecular Biology and Evolution 39(6):msac115.

  > Cury J, Touchon M, Rocha EPC. (2017) *Integrative and conjugative elements
  > and their hosts: composition, distribution and organization.*
  > Nucleic Acids Research 45(15):8943–8956. <https://doi.org/10.1093/nar/gkx607>

  > Abby SS, Cury J, Guglielmini J, Néron B, Touchon M, Rocha EPC. (2016)
  > *Identification of protein secretion systems in bacterial genomes.*
  > Scientific Reports 6:23080. <http://dx.doi.org/10.1038/srep23080>

  *Licence: CC BY-NC-SA 4.0 (Institut Pasteur / CNRS) — academic, non-commercial
  use only. BacFlux does not ship these models; they are fetched at run time.*

- **ICEscan models / ICEfinder2** — the optional second machinery model set
  (`mobilome.icescan.run`, default off), which supplies the IME and AICE classes:
  > Wang M, Goh Y-X, Tai C, Wang H, Deng Z, Ou H-Y. (2024)
  > *ICEberg 3.0: functional categorization and analysis of the integrative and
  > conjugative elements in bacteria.* Nucleic Acids Research 52(D1):D732–D737.
  > <https://doi.org/10.1093/nar/gkad935>

  ICEscan is a **fork of CONJScan** by the same Institut Pasteur authors — its own
  `metadata.yml` still identifies it as CONJScan 2.0.1 — so when it is enabled,
  **cite the CONJScan references above as well**. Same CC BY-NC-SA 4.0 terms; also
  not shipped. Rationale, measurements and limitations:
  `docs/methods_icescan_union.md`.

- **ICEberg 3.0** — element sequences and the curated ICE/IME coordinates used as
  benchmark ground truth. Wang *et al.* 2024, as above.

- **ISfinder** — Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M. (2006)
  *ISfinder: the reference centre for bacterial insertion sequences.*
  Nucleic Acids Research 34(Database issue):D32–D36.
  *Terms require written authorisation to download and forbid redistribution.
  Opt-in only; BacFlux ships a URL, never the data.*

- **ISOSDB / pseudoR** — Kirsch JM, Hryckowian AJ, Duerkop BA. (2024)
  *A metagenomics pipeline reveals insertion sequence-driven evolution of the
  microbiota.* Cell Host & Microbe 32(5):739–754.
  *ISOSDB is openly licensed (pseudoR repo is MIT).*

- **TnCentral** — Ross K, Varani AM, Snesrud E, Huang H, Alvarenga DO, Zhang J, et al. (2021)
  *TnCentral: a prokaryotic transposable element database and web portal for
  transposon analysis.* mBio 12(5):e0206021.
  *"© TnCentral — All Rights Reserved"; no terms page exists. Opt-in only.*

- **Platon** — Schwengers O, Barth P, Falgenhauer L, Hain T, Chakraborty T, Goesmann A. (2020)
  *Platon: identification and characterization of bacterial plasmid contigs in
  short-read draft assemblies exploiting protein sequence-based replicon
  distribution scores.* Microbial Genomics 6(10):mgen000398.

- **VirSorter2** — Guo J, Bolduc B, Zayed AA, Varsani A, Dominguez-Huerta G,
  Delmont TO, et al. (2021) *VirSorter2: a multi-classifier, expert-guided
  approach to detect diverse DNA and RNA viruses.* Microbiome 9:37.

- **CheckV** — Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. (2021)
  *CheckV assesses the quality and completeness of metagenome-assembled viral
  genomes.* Nature Biotechnology 39:578–585.

- **geNomad** — Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, et al. (2024)
  *Identification of mobile genetic elements with geNomad.* Nature Biotechnology 42:1303–1312.
  *⚠ Licence: Berkeley Lab **academic / non-commercial use only** — not BSD, despite
  the bioconda recipe's tag. Opt-in, default off; VirSorter2 is the default phage
  caller. See spec §3.1.*

### Supporting tools

- **HMMER** — Eddy SR. (2011) *Accelerated profile HMM searches.*
  PLoS Computational Biology 7:e1002195.

- **BLAST+** — Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K,
  Madden TL. (2009) *BLAST+: architecture and applications.* BMC Bioinformatics 10:421.

- **bedtools** — Quinlan AR, Hall IM. (2010) *BEDTools: a flexible suite of
  utilities for comparing genomic features.* Bioinformatics 26(6):841–842.

- **Snakemake** — Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH,
  Sochat V, et al. (2021) *Sustainable data analysis with Snakemake.*
  F1000Research 10:33.

- **Bioconda** — Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J,
  Tomkins-Tinch CH, et al. (2018) *Bioconda: sustainable and comprehensive
  software distribution for the life sciences.* Nature Methods 15:475–476.

---

## 2. Design influence

Work that shaped how the mobilome module is built, without contributing code to it.

### EBI Mobilome Annotation Pipeline

> EBI-Metagenomics `mobilome-annotation-pipeline` (Apache-2.0)
> <https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline>
> — which carries its own attribution to ICEfinder2.

**Complementary scope, not a simplified alternative.** These two pipelines answer
different questions and it would be wrong to present BacFlux as a lighter version
of theirs:

| | EBI Mobilome Annotation Pipeline | BacFlux mobilome module |
|---|---|---|
| **Input** | metagenomes and MAGs | single bacterial isolates |
| **Depth on the mobilome** | deep — full ICE delimitation, integron and phage annotation, mobilome GFF as the product | one focused question: for each AMR gene, is it in a mobile element, and how transferable is that element |
| **Scope of the workflow** | mobilome annotation | reads → assembly → decontamination → annotation → taxonomy → AMR → plasmids → phage → AMR mobility |
| **Runtime** | Nextflow | Snakemake, conda-per-rule |

What we took from them — all of it facts and design decisions, which spec §11
expressly permits:

- that **ICEscan** exists, ships the IME and AICE models, and is worth running;
- Sequence Ontology terms for the mobilome GFF `Type` column;
- the `contig_id|mge_type-start:end` element ID format;
- the discard-with-reason file pattern;
- the 500 bp minimum element size and 0.9 CDS-coverage assignment thresholds.

**No code was copied.** Three of their scripts —
`bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py` and
`bin/prescan_to_fasta.py` — are derived from ICEfinder2 and carry
**CC BY-NC-SA 4.0** headers. That licence is non-commercial and share-alike, so
copying or adapting any of it would contaminate BacFlux's MIT licence and impose
a non-commercial restriction on everyone downstream. Reading it to understand an
algorithm is fine and is what was done; the separation is real and deliberate
(spec §11).

**Escape hatch for users.** For full ICE delimitation and deeper mobilome
analysis, run their pipeline — BacFlux writes `{sample}_contigs.fna` and
`{sample}.gbk` ready for it, and consuming its `mobilome.gff.gz` output carries
no licence consequences.

> **Correction, 2026-07-30.** Earlier BacFlux documentation said this pipeline
> "deliberately does not run AMRFinderPlus". That is **false** — it runs
> AMRFinderPlus by default, alongside DeepARG, RGI/CARD and hAMRonization, and
> integrates the results with the mobilome. The complementary-scope framing above
> rests on the metagenome/MAG versus single-isolate distinction, and on nothing else.

### ICEfinder2 and the att-site literature

ICEfinder2's approach to boundary detection (variable-length repeat search rather
than a fixed-length probe) was established by reading its source and the
surrounding literature, then implemented independently. The evidence, the
parameters and the negative results are recorded in
`docs/methods_att_and_small_plasmids.md`.

### Regulatory framing

- **EFSA FEEDAP Panel** (2018) *Guidance on the characterisation of microorganisms
  used as feed additives or as production organisms.* EFSA Journal 16(3):5206.
- **EFSA** (2024) *Statement on how to interpret the QPS qualification on
  'absence of acquired antimicrobial resistance genes'.* EFSA Journal.

BacFlux provides **supporting evidence** for the intrinsic-versus-acquired
distinction these documents require. It is not "EFSA-compliant" and does not
claim to be, and every mobility call is a **prediction** — confirmation requires
a filter or broth mating assay.
