# Citations

Most of the science BacFlux reports is produced by other people's tools and
databases, and this file is where they are credited.

**If you publish results from BacFlux, cite the tools whose output you actually
used** — not just BacFlux. The mobilome module in particular is an integrator:
it decides how to combine other tools' calls, and every underlying detection
belongs to somebody else.

Two sections, because using another project's code or data is not the same as
learning from its design:

1. **Tools and databases used** — code that runs, or data that is read.
2. **Design influence** — work we learned from and adopted decisions from,
   without using their code.

---

## 1. Tools and databases used

Every program BacFlux runs and every database it reads is listed below, with the
papers its authors ask you to cite. Where they ask for more than one paper, or for
the programs and databases their tool uses in turn, the entry gives those as well;
for Bakta, whose list is long, it points to that list. Generic utilities (wget,
curl, tar, gzip, bzip2, unzip) and the Python library that writes the workbooks
(XlsxWriter) are not listed.

### Reads

- **fastp** (adapter and quality trimming of the Illumina reads in `illumina` and
  `hybrid` mode) — Chen S. (2025) *fastp 1.0: An ultra-fast all-round tool for FASTQ
  data quality control and preprocessing.* iMeta 4(5):e70078.
  <https://doi.org/10.1002/imt2.70078>
  Its README lists two earlier papers as well: Chen S, Zhou Y, Chen Y, Gu J. (2018)
  *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics
  34(17):i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>; Chen S. (2023)
  *Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication
  using fastp.* iMeta 2(2):e107. <https://doi.org/10.1002/imt2.107>

- **Bowtie 2** (in `illumina` and `hybrid` mode: removes the PhiX reads, and maps the
  reads back onto the draft assembly for the BlobTools coverage track; in `hybrid`
  mode it also maps them onto the decontaminated assembly to pull out the clean read
  pairs) — Langmead B, Salzberg SL. (2012) *Fast gapped-read alignment with Bowtie 2.*
  Nature Methods 9(4):357–359. <https://doi.org/10.1038/nmeth.1923>
  Its manual asks for a second paper as well: Langmead B, Wilks C, Antonescu V,
  Charles R. (2019) *Scaling read aligners to hundreds of threads on general-purpose
  processors.* Bioinformatics 35(3):421–432.
  <https://doi.org/10.1093/bioinformatics/bty648>

- **PhiX control genome** (the spike-in removed from the Illumina reads before
  assembly in `illumina` and `hybrid` mode; `links.phix_link`) — Escherichia phage
  phiX174, NCBI RefSeq assembly GCF_000819615.1 (sequence NC_001422.1), the file the
  quick start names.

- **NanoPlot** (length and quality profiles of the raw ONT reads and of the reads
  Filtlong keeps, in `nanopore` and `hybrid` mode) — De Coster W, Rademakers R.
  (2023) *NanoPack2: population-scale evaluation of long-read sequencing data.*
  Bioinformatics 39(5):btad311. <https://doi.org/10.1093/bioinformatics/btad311>

- **fastplong** (ONT adapter and barcode trimming in `nanopore` and `hybrid` mode;
  `parameters.<mode>.ont_adapter: off` skips it) — OpenGene.
  <https://github.com/OpenGene/fastplong>. It has no publication of its own; its README
  lists the two fastp papers under "citations": Chen S. (2023) *Ultrafast one-pass
  FASTQ data preprocessing, quality control, and deduplication using fastp.* iMeta
  2(2):e107. <https://doi.org/10.1002/imt2.107>; Chen S, Zhou Y, Chen Y, Gu J. (2018)
  *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics
  34(17):i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>
  *MIT licence.*

- **Filtlong** (selects the ONT reads that go to assembly by length and quality, in
  `nanopore` and `hybrid` mode; with `parameters.hybrid.short_read_guidance: true` it
  judges quality by agreement with the decontaminated short reads instead of the Phred
  scores, and trims and splits reads where they disagree) — Wick RR. *Filtlong:
  quality filtering tool for long reads.* <https://github.com/rrwick/Filtlong>. It has
  no publication of its own.

### Assembly and polishing

- **SPAdes** (short-read assembly in `illumina` and `hybrid` mode, with `--isolate`) —
  Prjibelski A, Antipov D, Meleshko D, Lapidus A, Korobeynikov A. (2020) *Using SPAdes
  De Novo Assembler.* Current Protocols in Bioinformatics 70(1):e102.
  <https://doi.org/10.1002/cpbi.102>
  SPAdes asks for this, its latest paper, and allows its older papers as well, among
  them Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, et al.
  (2012) *SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell
  Sequencing.* Journal of Computational Biology 19(5):455–477.
  <https://doi.org/10.1089/cmb.2012.0021>

- **Flye** (long-read assembly in `nanopore` and `hybrid` mode) — Kolmogorov M,
  Yuan J, Lin Y, Pevzner PA. (2019) *Assembly of long, error-prone reads using repeat
  graphs.* Nature Biotechnology 37(5):540–546.
  <https://doi.org/10.1038/s41587-019-0072-8>

- **dnaapler** (`dnaapler all` rotates each circular contig of the ONT assembly to
  start at dnaA, repA, terL or COG1474, in `nanopore` and `hybrid` mode) — Bouras G,
  Grigson SR, Papudeshi B, Mallawaarachchi V, Roach MJ. (2024) *Dnaapler: A tool to
  reorient circular microbial genomes.* Journal of Open Source Software 9(93):5968.
  <https://doi.org/10.21105/joss.05968>
  Its README also asks you to consider citing its dependencies where relevant. It
  always searches with **MMseqs2**: Steinegger M, Söding J. (2017) *MMseqs2 enables
  sensitive protein sequence searching for the analysis of massive data sets.* Nature
  Biotechnology 35(11):1026–1028. <https://doi.org/10.1038/nbt.3988>
  When a contig's best hit does not begin with a valid start codon, it also calls
  **Pyrodigal** to find the gene, and Pyrodigal asks for its own paper and Prodigal's:
  Larralde M. (2022) *Pyrodigal: Python bindings and interface to Prodigal, an
  efficient method for gene prediction in prokaryotes.* Journal of Open Source
  Software 7(72):4296. <https://doi.org/10.21105/joss.04296>; Hyatt D, Chen G-L,
  LoCascio PF, Land ML, Larimer FW, Hauser LJ. (2010) *Prodigal: prokaryotic gene
  recognition and translation initiation site identification.* BMC Bioinformatics
  11:119. <https://doi.org/10.1186/1471-2105-11-119>

- **Medaka** (long-read consensus polishing of the ONT assembly in `nanopore` and
  `hybrid` mode; `parameters.<mode>.medaka_model: false` skips it) — Oxford Nanopore
  Technologies. <https://github.com/nanoporetech/medaka>. It has no publication of its
  own, and its README asks for no citation.

- **BWA** (BWA-MEM: `bwa mem -a` aligns each Illumina mate to the ONT assembly for
  Polypolish in `hybrid` mode; Snippy runs it too) — Li H. (2013) *Aligning sequence
  reads, clone sequences and assembly contigs with BWA-MEM.* arXiv:1303.3997v2
  [q-bio.GN]. <https://doi.org/10.48550/arXiv.1303.3997>

- **Polypolish** (corrects the ONT assembly with the Illumina reads in `hybrid` mode:
  `polypolish filter`, then `polypolish polish`) — Wick RR, Holt KE. (2022)
  *Polypolish: Short-read polishing of long-read bacterial genome assemblies.* PLOS
  Computational Biology 18(1):e1009802. <https://doi.org/10.1371/journal.pcbi.1009802>
  Its README also lists the follow-up paper that describes Polypolish v0.6.0, the
  release that introduced the `filter` and `polish` subcommands BacFlux runs:
  Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. (2024) *How low can
  you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies.*
  Microbial Genomics 10(6):001254. <https://doi.org/10.1099/mgen.0.001254>

- **Snippy** (compares each stage of the ONT assembly with the Illumina assembly in
  `hybrid` mode, run with `--ctgs`) — Seemann T. (2015) *snippy: fast bacterial
  variant calling from NGS reads.* <https://github.com/tseemann/snippy>

### Quality control and decontamination

- **QUAST** (assembly statistics of the finished genome, in every mode; both
  assemblies in `hybrid` mode) — Gurevich A, Saveliev V, Vyahhi N, Tesler G. (2013)
  *QUAST: quality assessment tool for genome assemblies.* Bioinformatics
  29(8):1072–1075. <https://doi.org/10.1093/bioinformatics/btt086>
  QUAST's manual names this paper for QUAST in general, and for its latest version,
  5.3.0, the one BacFlux runs, it asks for: Mikheenko A, Saveliev V, Hirsch P,
  Gurevich A. (2023) *WebQUAST: online evaluation of genome assemblies.* Nucleic Acids
  Research 51(W1):W601–W606. <https://doi.org/10.1093/nar/gkad406>

- **CheckM** (completeness and contamination of the finished genome, in every mode;
  both assemblies in `hybrid` mode) — Parks DH, Imelfort M, Skennerton CT,
  Hugenholtz P, Tyson GW. (2015) *CheckM: assessing the quality of microbial genomes
  recovered from isolates, single cells, and metagenomes.* Genome Research
  25(7):1043–1055. <https://doi.org/10.1101/gr.186072.114>
  CheckM runs **pplacer**, **Prodigal** and **HMMER**, and its authors recommend
  citing them too: Matsen FA, Kodner RB, Armbrust EV. (2010) *pplacer: linear time
  maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed
  reference tree.* BMC Bioinformatics 11:538.
  <https://doi.org/10.1186/1471-2105-11-538>; for Prodigal, the paper CheckM names:
  Hyatt D, LoCascio PF, Hauser LJ, Uberbacher EC. (2012) *Gene and translation
  initiation site prediction in metagenomic sequences.* Bioinformatics
  28(17):2223–2230. <https://doi.org/10.1093/bioinformatics/bts429>; HMMER as listed
  under Supporting tools.

- **Qualimap** (mapping-quality report of the reads on the draft assembly, in
  `illumina`, `nanopore` and `hybrid` mode; in `hybrid`, the Illumina reads on the
  Illumina draft) — Okonechnikov K, Conesa A, García-Alcalde F. (2016) *Qualimap 2:
  advanced multi-sample quality control for high-throughput sequencing data.*
  Bioinformatics 32(2):292–294. <https://doi.org/10.1093/bioinformatics/btv566>

- **BlobTools** (joins the BLAST hits and the coverage track into one taxon call per
  contig for the decontamination screen, in every mode, and for the screen of the
  delivered genome in `hybrid` mode) — Laetsch DR, Blaxter ML. (2017) *BlobTools:
  Interrogation of genome assemblies.* F1000Research 6:1287.
  <https://doi.org/10.12688/f1000research.12232.1>

- **minimap2** (read coverage for BlobTools and Qualimap in `nanopore` mode, and for
  BlobTools on the delivered genome in `hybrid` mode; the contigs mapped to themselves
  in `contigs` mode) — Li H. (2018) *Minimap2: pairwise alignment for nucleotide
  sequences.* Bioinformatics 34(18):3094–3100.
  <https://doi.org/10.1093/bioinformatics/bty191>
  Its README asks for this paper and/or a second one: Li H. (2021) *New strategies to
  improve minimap2 alignment accuracy.* Bioinformatics 37(23):4572–4574.
  <https://doi.org/10.1093/bioinformatics/btab705>

- **SAMtools** (converts, sorts and indexes the alignments behind the BlobTools
  coverage track in every mode; in `hybrid` mode it also writes out the read pairs
  that map as proper pairs to the decontaminated assembly) — Danecek P, Bonfield JK,
  Liddle J, Marshall J, Ohan V, Pollard MO, et al. (2021) *Twelve years of SAMtools
  and BCFtools.* GigaScience 10(2):giab008.
  <https://doi.org/10.1093/gigascience/giab008>

- **NCBI core_nt** (the BLAST database searched for the decontamination screens and
  the plasmid check, in every mode; `parameters.nt_version` can name another NCBI BLAST
  database, such as `nt_prok`) — National Center for Biotechnology Information.
  *BLAST databases.* <https://ftp.ncbi.nlm.nih.gov/blast/db/>
  NCBI makes core_nt from nt by leaving out nuclear eukaryotic chromosome assemblies,
  and its BLAST pages name no paper for either. Its citation guidance names its yearly
  article on its databases, the newest being Sayers EW, Bolton EE, Fine AM, Kelly C,
  Kim S, Landrum M, et al. (2026) *Database resources of the National Center for
  Biotechnology Information in 2026.* Nucleic Acids Research 54(D1):D20–D27.
  <https://doi.org/10.1093/nar/gkaf1060>; and the GenBank release notes ask users of
  GenBank data for the GenBank paper: Sayers EW, Cavanaugh M, Frisse L,
  Pruitt KD, Schneider VA, Underwood BA, et al. (2025) *GenBank 2025 update.* Nucleic
  Acids Research 53(D1):D56–D61. <https://doi.org/10.1093/nar/gkae1114>

- **NCBI Taxonomy** (`nodes.dmp` and `names.dmp`, which BlobTools reads from
  `directories.blast_db` to name the taxon of every BLAST hit, in every mode, and the
  `taxdb` files beside them, which give the BLAST tables their taxon names) —
  Schoch CL, Ciufo S, Domrachev M, Hotton CL, Kannan S, Khovanskaya R, et al. (2020)
  *NCBI Taxonomy: a comprehensive update on curation, resources and tools.* Database
  2020:baaa062. <https://doi.org/10.1093/database/baaa062>
  NCBI's citation page lists a second paper as well: Cox E, Tsuchiya MTN, Ciufo S,
  Torcivia J, Falk R, Anderson WR, et al. (2025) *NCBI Taxonomy: enhanced access via
  NCBI Datasets.* Nucleic Acids Research 53(D1):D1711–D1715.
  <https://doi.org/10.1093/nar/gkae967>

- **Adapter sequences** (the ONT adapter scan in `nanopore` and `hybrid` mode, which
  `parameters.<mode>.ont_adapter: off` skips, and the adapter screen of the finished
  genome in every mode) — Oxford Nanopore Technologies, *Chemistry Technical
  Document*, version CHTD_500_v1_revAV_24July2026, sections 14 and 15
  (<https://nanoporetech.com/document/chemistry-technical-document>): the adapter and
  the 96 native and 96 rapid barcodes shipped in
  `workflow/data/adapters/ont_chtd.fasta`. National Center for Biotechnology
  Information, *The UniVec Database*
  (<https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/>; files at
  <https://ftp.ncbi.nlm.nih.gov/pub/UniVec/>): the PacBio and MGI records shipped in
  `workflow/data/adapters/univec_pacbio_mgi.fasta`, and the Illumina records of build
  10.1, which the workflow downloads.
  The 50 bp after every rapid barcode are the right end of phage Mu, identical to the
  last 50 bp of RefSeq NC_000929.1 (derived from GenBank AF083977.1). A genome with a
  Mu prophage holds them too, so the screen lists them alone inside a contig without
  counting them: Morgan GJ, Hatfull GF, Casjens S, Hendrix RW. (2002) *Bacteriophage
  Mu genome sequence: analysis and comparison with Mu-like prophages in Haemophilus,
  Neisseria and Deinococcus.* Journal of Molecular Biology 317(3):337–359.
  <https://doi.org/10.1006/jmbi.2002.5437>
  *ONT states no licence for its document's sequences; they are reproduced as ONT
  publishes them, with the source in every record. UniVec is a work of the US
  Government; Illumina's oligonucleotide sequences are copyright Illumina, Inc., so
  BacFlux downloads them on your computer and does not redistribute them.*

### Annotation and taxonomy

- **Bakta** (gene calling and annotation of the finished genome, every mode; BacFlux
  expects database v6.0, <https://doi.org/10.5281/zenodo.14916843>) — Schwengers O,
  Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A. (2021) *Bakta: rapid and
  standardized annotation of bacterial genomes via alignment-free sequence
  identification.* Microbial Genomics 7(11):000685.
  <https://doi.org/10.1099/mgen.0.000685>
  Bakta runs tRNAscan-SE, Aragorn, Infernal, PILER-CR, Pyrodigal, DIAMOND, BLAST+,
  PyHMMER, AMRFinderPlus and, for its genome plot, pyCirclize; its database draws on
  sources such as UniProt,
  RefSeq, Rfam, Pfam, KEGG and ISfinder; its README asks you to cite these primary
  sources as well where you rely on them, and lists the papers.

- **Baktfold** (names, from predicted structure, for proteins Bakta leaves
  hypothetical; optional, `baktfold.run: true`; its ProstT5 step also runs for the
  trait profile, `traits.run: true`) — Bouras G, Lim SW, Durr L, Vreugde S,
  Goesmann A, Edwards RA, et al. (2026) *Baktfold: sensitive protein functional
  annotation across the microbial tree of life using structural information.* bioRxiv
  2026.03.31.715528. <https://doi.org/10.64898/2026.03.31.715528>
  Baktfold runs **ProstT5** and **Foldseek**, cite both with it:
  Heinzinger M, Weissenow K, Sanchez JG, Henkel A, Mirdita M, Steinegger M, et al.
  (2024) *Bilingual language model for protein sequence and structure.* NAR Genomics
  and Bioinformatics 6(4):lqae150. <https://doi.org/10.1093/nargab/lqae150>;
  van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, et al. (2024)
  *Fast and accurate protein structure search with Foldseek.* Nature Biotechnology
  42(2):243–246. <https://doi.org/10.1038/s41587-023-01773-0>
  Its database, v0.0.1, the one Baktfold 0.3.0 installs
  (<https://doi.org/10.5281/zenodo.17347516>), is built from SwissProt, the AlphaFold
  Protein Structure Database, PDB and CATH; if you report which of them a name came
  from, cite that database as well (Baktfold's README lists the papers).

- **EcoFoldDB** (the reference structures the trait profile searches; optional,
  `traits.run: true`) — Ghaly TM, Rajabal V, Russell D, Colombi E, Tetu SG. (2025)
  *EcoFoldDB: Protein Structure-Guided Functional Profiling of Ecologically Relevant
  Microbial Traits at the Metagenome Scale.* Environmental Microbiology
  27(9):e70178. <https://doi.org/10.1111/1462-2920.70178>
  The trait profile searches its database, EcoFoldDB v2.0
  (<https://doi.org/10.5281/zenodo.16908168>), with **Foldseek**, after Baktfold's
  **ProstT5** step has turned each protein of up to 4,000 amino acids into a string of
  Foldseek's 3Di structural alphabet; cite Baktfold, ProstT5 and Foldseek as listed
  under Baktfold above. The trait profile also reads the KEGG Orthology numbers in
  Bakta's annotation, which Bakta's database takes from KEGG's KOfam profiles; for
  those, Bakta's README gives: Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K,
  Kanehisa M, Goto S, et al. (2020) *KofamKOALA: KEGG Ortholog assignment based on
  profile HMM and adaptive score threshold.* Bioinformatics 36(7):2251–2252.
  <https://doi.org/10.1093/bioinformatics/btz859>
  The primary papers behind each trait rule are in the `reference` column of
  `workflow/scripts/40_annotation/trait_rules_v1.tsv`.
  *The database is CC BY 4.0. The EcoFoldDB repository's code carries no licence, and
  BacFlux does not use it.*

- **eggNOG-mapper** (orthology-based functional labels for Bakta's proteins: COG, GO,
  KEGG and EC; every mode) — Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P,
  Huerta-Cepas J. (2021) *eggNOG-mapper v2: Functional Annotation, Orthology
  Assignments, and Domain Prediction at the Metagenomic Scale.* Molecular Biology and
  Evolution 38(12):5825–5829. <https://doi.org/10.1093/molbev/msab293>
  Its authors ask you to cite the database as well; BacFlux reads **eggNOG 5.0**
  (release 5.0.2): Huerta-Cepas J, Szklarczyk D, Heller D, Hernández-Plaza A,
  Forslund SK, Cook H, et al. (2019) *eggNOG 5.0: a hierarchical, functionally and
  phylogenetically annotated orthology resource based on 5090 organisms and 2502
  viruses.* Nucleic Acids Research 47(D1):D309–D314.
  <https://doi.org/10.1093/nar/gky1085>
  They also ask you to cite the program behind the search step: BacFlux runs the
  DIAMOND mode (`-m diamond`), so cite **DIAMOND** too (under Supporting tools).

- **antiSMASH** (biosynthetic gene clusters on Bakta's annotated genome; every mode) —
  Blin K, Shaw S, Vader L, Szenei J, Reitz ZL, Augustijn HE, et al. (2025)
  *antiSMASH 8.0: extended gene cluster detection capabilities and analyses of
  chemistry, enzymology, and regulation.* Nucleic Acids Research 53(W1):W32–W38.
  <https://doi.org/10.1093/nar/gkaf334>
  BacFlux fetches antiSMASH's own database bundle with its `download-antismash-databases`
  command.

- **run_dbCAN / dbCAN3** (CAZymes, CAZyme gene clusters and their predicted
  substrates; every mode; run_dbCAN and PyHMMER are installed from PyPI with pip, not
  from a conda channel) — Zheng J, Ge Q, Yan Y, Zhang X, Huang L, Yin Y. (2023)
  *dbCAN3: automated carbohydrate-active enzyme and substrate annotation.* Nucleic
  Acids Research 51(W1):W115–W121. <https://doi.org/10.1093/nar/gkad328>
  The run_dbCAN README asks for two earlier dbCAN papers as well: Zhang H, Yohe T,
  Huang L, Entwistle S, Wu P, Yang Z, et al. (2018) *dbCAN2: a meta server for
  automated carbohydrate-active enzyme annotation.* Nucleic Acids Research
  46(W1):W95–W101. <https://doi.org/10.1093/nar/gky418>; Huang L, Zhang H, Wu P,
  Entwistle S, Li X, Yohe T, et al. (2018) *dbCAN-seq: a database of
  carbohydrate-active enzyme (CAZyme) sequence and annotation.* Nucleic Acids
  Research 46(D1):D516–D521. <https://doi.org/10.1093/nar/gkx894>

- **GTDB-Tk** (taxonomic placement, `classify_wf`; every mode) — Chaumeil P-A,
  Mussig AJ, Hugenholtz P, Parks DH. (2022) *GTDB-Tk v2: memory friendly
  classification with the genome taxonomy database.* Bioinformatics 38(23):5315–5316.
  <https://doi.org/10.1093/bioinformatics/btac672>; Chaumeil P-A, Mussig AJ,
  Hugenholtz P, Parks DH. (2020) *GTDB-Tk: a toolkit to classify genomes with the
  Genome Taxonomy Database.* Bioinformatics 36(6):1925–1927.
  <https://doi.org/10.1093/bioinformatics/btz848>
  Its authors strongly encourage citing the programs it runs as well. `classify_wf`
  first compares every genome with the GTDB representative genomes using **skani**; a
  genome that this screen does not assign to a GTDB species is then placed in the
  reference tree with **Prodigal**, **HMMER** (under Supporting tools) and
  **pplacer**: Shaw J, Yu YW. (2023) *Fast and robust metagenomic sequence comparison
  through sparse chaining with skani.* Nature Methods 20(11):1661–1665.
  <https://doi.org/10.1038/s41592-023-02018-3>; for Prodigal, the paper GTDB-Tk
  names: Hyatt D, Chen G-L, LoCascio PF, Land ML, Larimer FW, Hauser LJ. (2010)
  *Prodigal: prokaryotic gene recognition and translation initiation site
  identification.* BMC Bioinformatics 11:119.
  <https://doi.org/10.1186/1471-2105-11-119>; Matsen FA, Kodner RB, Armbrust EV.
  (2010) *pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement
  of sequences onto a fixed reference tree.* BMC Bioinformatics 11:538.
  <https://doi.org/10.1186/1471-2105-11-538>

- **GTDB** (the reference taxonomy GTDB-Tk places each genome in; BacFlux expects
  release R232) — Parks DH, Chaumeil P-A, Mussig AJ, Rinke C, Chuvochina M,
  Hugenholtz P. (2026) *GTDB release 10: a complete and systematic taxonomy for
  715 230 bacterial and 17 245 archaeal genomes.* Nucleic Acids Research
  54(D1):D743–D754. <https://doi.org/10.1093/nar/gkaf1040>; Parks DH, Chuvochina M,
  Rinke C, Mussig AJ, Chaumeil P-A, Hugenholtz P. (2022) *GTDB: an ongoing census of
  bacterial and archaeal diversity through a phylogenetically consistent, rank
  normalized and complete genome-based taxonomy.* Nucleic Acids Research
  50(D1):D785–D794. <https://doi.org/10.1093/nar/gkab776>
  GTDB's own list asks for four more papers as well: Chuvochina M, Mussig AJ,
  Chaumeil P-A, Skarshewski A, Rinke C, Parks DH, et al. (2023) *Proposal of names
  for 329 higher rank taxa defined in the Genome Taxonomy Database under two
  prokaryotic codes.* FEMS Microbiology Letters 370:fnad071.
  <https://doi.org/10.1093/femsle/fnad071>; Rinke C, Chuvochina M, Mussig AJ,
  Chaumeil P-A, Davín AA, Waite DW, et al. (2021) *A standardized archaeal taxonomy
  for the Genome Taxonomy Database.* Nature Microbiology 6(7):946–959.
  <https://doi.org/10.1038/s41564-021-00918-8>; Parks DH, Chuvochina M, Chaumeil P-A,
  Rinke C, Mussig AJ, Hugenholtz P. (2020) *A complete domain-to-species taxonomy for
  Bacteria and Archaea.* Nature Biotechnology 38(9):1079–1086.
  <https://doi.org/10.1038/s41587-020-0501-8>; Parks DH, Chuvochina M, Waite DW,
  Rinke C, Skarshewski A, Chaumeil P-A, et al. (2018) *A standardized bacterial
  taxonomy based on genome phylogeny substantially revises the tree of life.* Nature
  Biotechnology 36(10):996–1004. <https://doi.org/10.1038/nbt.4229>
  The newest paper covers releases up to R226, so name the release you used (R232) in
  your methods. Two small tables BacFlux ships, which map GTDB species to
  AMRFinderPlus organisms for the mobilome module
  (`workflow/scripts/80_mobilome/gtdb_organism_equivalences.tsv` and
  `gtdb_organism_genus_rules.tsv`), were computed from GTDB R232's metadata.

### Antimicrobial resistance and virulence

- **ABRicate** (screens the finished contigs for resistance and virulence genes in
  every mode, once against each of its own copies of CARD, NCBI, ResFinder and VFDB) —
  Seemann T. *ABRicate: mass screening of contigs for antimicrobial and virulence
  genes.* <https://github.com/tseemann/abricate>
  ABRicate asks you to also cite the database behind each table you report: CARD,
  ResFinder and VFDB, listed in this section, and for `ncbi` the AMRFinderPlus paper.
  The ABRicate version and the date of each database copy are in
  `05.amr/abricate/abricate_databases.tsv`.

- **AMRFinderPlus and the NCBI Reference Gene Catalog** (Bakta runs AMRFinderPlus in
  every mode for its AMR gene annotations; with `mobilome.run: true` BacFlux also runs
  it itself, from the Bakta environment against the database inside
  `directories.bakta_db`; ABRicate's `ncbi` database, used in every mode, is built from
  the catalogue's core resistance genes) — Feldgarden M, Brover V,
  Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, et al. (2021) *AMRFinderPlus and
  the Reference Gene Catalog facilitate examination of the genomic links among
  antimicrobial resistance, stress response, and virulence.* Scientific Reports
  11(1):12728. <https://doi.org/10.1038/s41598-021-91456-0>
  NCBI asks you to state the software and database versions you used. The
  AMRFinderPlus wiki also lists the paper on its predecessor, AMRFinder: Feldgarden M,
  Brover V, Haft DH, Prasad AB, Slotta DJ, Tolstoy I, et al. (2019) *Validating the
  AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance
  Genotype-Phenotype Correlations in a Collection of Isolates.* Antimicrobial Agents
  and Chemotherapy 63(11):e00483-19. <https://doi.org/10.1128/AAC.00483-19>

- **CARD** (resistance genes: ABRicate's `card` database in every mode; in `illumina`
  and `hybrid` mode also the reference the trimmed reads are mapped onto with BBMap,
  from `links.card_link` or your own copy in `directories.card_db`) — Alcock BP,
  Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, et al. (2023) *CARD 2023:
  expanded curation, support for machine learning, and resistome prediction at the
  Comprehensive Antibiotic Resistance Database.* Nucleic Acids Research
  51(D1):D690–D699. <https://doi.org/10.1093/nar/gkac920>
  CARD's own list of papers also gives earlier ones, among them Jia B, Raphenya AR,
  Alcock B, Waglechner N, Guo P, Tsang KK, et al. (2017) *CARD 2017: expansion and
  model-centric curation of the comprehensive antibiotic resistance database.*
  Nucleic Acids Research 45(D1):D566–D573. <https://doi.org/10.1093/nar/gkw1004>

- **ResFinder** (acquired resistance genes; ABRicate's `resfinder` database, every
  mode) — Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, et al.
  (2012) *Identification of acquired antimicrobial resistance genes.* Journal of
  Antimicrobial Chemotherapy 67(11):2640–2644. <https://doi.org/10.1093/jac/dks261>
  That is the paper the ResFinder database asks for. ABRicate's README gives
  ResFinder 4.0 for the same database, so cite it too: Bortolaia V, Kaas RS, Ruppe E,
  Roberts MC, Schwarz S, Cattoir V, et al. (2020) *ResFinder 4.0 for predictions of
  phenotypes from genotypes.* Journal of Antimicrobial Chemotherapy
  75(12):3491–3500. <https://doi.org/10.1093/jac/dkaa345>

- **VFDB** (virulence factors; ABRicate's `vfdb` database, VFDB's core set of
  experimentally verified factors, every mode) — Chen L, Zheng D, Liu B, Yang J,
  Jin Q. (2016) *VFDB 2016: hierarchical and refined dataset for big data
  analysis—10 years on.* Nucleic Acids Research 44(D1):D694–D697.
  <https://doi.org/10.1093/nar/gkv1239>
  That is the paper ABRicate's README gives. VFDB's website lists its papers without
  asking for one; the newest is Zhou S, Liu B, Zheng D, Chen L, Yang J. (2025)
  *VFDB 2025: an integrated resource for exploring anti-virulence compounds.* Nucleic
  Acids Research 53(D1):D871–D877. <https://doi.org/10.1093/nar/gkae968>

- **BBMap** (BBTools; in `illumina` and `hybrid` mode, maps the trimmed reads onto
  CARD, and, with `mobilome.run: true` and ISOSDB configured, onto ISOSDB and the
  genome to estimate IS copy number) — Bushnell B. (2014) *BBMap: A Fast, Accurate,
  Splice-Aware Aligner.* 9th Annual Genomics of Energy & Environment Meeting, Walnut
  Creek, CA (report LBNL-7065E). <https://www.osti.gov/biblio/1241166>; software:
  <https://bbmap.org>, <https://github.com/bbushnell/BBTools>

### Mobile genetic elements, plasmids and phages

- **ISEScan** (finds insertion sequences on the final contigs; mobilome module,
  `mobilome.run: true`) — Xie Z, Tang H. (2017) *ISEScan: automated identification
  of insertion sequence elements in prokaryotic genomes.* Bioinformatics
  33(21):3340–3347. <https://doi.org/10.1093/bioinformatics/btx433>

- **MacSyFinder v2** (runs the CONJScan models on Bakta's proteins when
  `mobilome.run: true`, and the ICEscan models as well when
  `mobilome.icescan.run: true`) — Néron B, Denise R, Coluzzi C, Touchon M,
  Rocha EPC, Abby SS. (2023) *MacSyFinder v2: Improved modelling and search engine
  to identify molecular systems in genomes.* Peer Community Journal 3:e28.
  <https://doi.org/10.24072/pcjournal.250>

- **CONJScan models** (the conjugation-machinery models MacSyFinder runs, version
  2.1.0, fetched when `mobilome.run: true`) — cite all three; the package's own
  `metadata.yml` lists them, and they are the models BacFlux actually runs:
  > Coluzzi C, Garcillán-Barcia MP, de la Cruz F, Rocha EPC. (2022)
  > *Evolution of plasmid mobility: origin and fate of conjugative and
  > nonconjugative plasmids.* Molecular Biology and Evolution 39(6):msac115.
  > <https://doi.org/10.1093/molbev/msac115>

  > Cury J, Touchon M, Rocha EPC. (2017) *Integrative and conjugative elements
  > and their hosts: composition, distribution and organization.*
  > Nucleic Acids Research 45(15):8943–8956. <https://doi.org/10.1093/nar/gkx607>

  > Abby SS, Cury J, Guglielmini J, Néron B, Touchon M, Rocha EPC. (2016)
  > *Identification of protein secretion systems in bacterial genomes.*
  > Scientific Reports 6:23080. <https://doi.org/10.1038/srep23080>

  *Licence: CC BY-NC-SA 4.0 (Institut Pasteur / CNRS) — academic, non-commercial
  use only. BacFlux does not ship these models; they are fetched at run time.*

- **ICEscan models / ICEfinder2** (the optional second model set, which supplies the
  IME and AICE classes; `mobilome.icescan.run: true`, default off) — ICEscan is a
  fork of CONJScan 2.0.1 (Institut Pasteur / CNRS) distributed inside the database
  bundle of ICEfinder2, the ICE finder of the ICEberg group (Shanghai Jiao Tong
  University); its own `metadata.yml` lists the three CONJScan papers. When it is
  on, cite the three CONJScan references above and the ICEberg 3.0 paper, which
  describes ICEfinder 2.0 and which the ICEberg site asks users to cite:
  > Wang M, Liu G, Liu M, Tai C, Deng Z, Song J, et al. (2024)
  > *ICEberg 3.0: functional categorization and analysis of the integrative and
  > conjugative elements in bacteria.* Nucleic Acids Research 52(D1):D732–D737.
  > <https://doi.org/10.1093/nar/gkad935>

  Licence: **CC BY-NC-SA 4.0**, the same terms as the CONJScan models it forks;
  not shipped by BacFlux. Rationale, measurements and limitations:
  `docs/methods_icescan_union.md`.

- **ICEberg 3.0** (names the ICE and IME candidates, `mobilome.iceberg.urls` or
  `.dir`; its curated ICE/IME coordinates are also the benchmark ground truth) —
  Wang *et al.* 2024, as above.

- **ISfinder** (BacFlux's configuration downloads nothing from it, but Bakta's
  database integrates annotations from the ISfinder database, so cite it when you
  report Bakta's insertion-sequence annotations; cite it too if you point
  `mobilome.tncentral.url` at a TnCentral endpoint that bundles ISfinder content) —
  Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M. (2006) *ISfinder: the
  reference centre for bacterial insertion sequences.* Nucleic Acids Research
  34(Database issue):D32–D36. <https://doi.org/10.1093/nar/gkj014>
  ISfinder also asks for its database URL: <https://www-is.biotoul.fr>
  *Terms require written authorisation to download and forbid redistribution.
  BacFlux's configuration offers no download URL for it and never ships the data; see
  [Licensing](https://ilivius.github.io/BacFlux/about/licensing/#isfinder-is-not-wired-in).*

- **ISOSDB / pseudoR** (the IS sequences reads are mapped to for IS copy number;
  `mobilome.isosdb.fasta_url` + `family_map_url`, or `.dir`; `illumina` and `hybrid`
  mode) — Kirsch JM, Hryckowian AJ, Duerkop BA. (2024) *A metagenomics pipeline
  reveals insertion sequence-driven evolution of the microbiota.*
  Cell Host & Microbe 32(5):739–754.e4. <https://doi.org/10.1016/j.chom.2024.03.005>
  *ISOSDB is MIT licensed, in the pseudoR repository.*

- **TnCentral** (names curated transposons and integrons, which makes mobility tier 4
  reachable; `mobilome.tncentral.url` or `.dir`) — Ross K, Varani AM, Snesrud E,
  Huang H, Alvarenga DO, Zhang J, et al. (2021) *TnCentral: a prokaryotic
  transposable element database and web portal for transposon analysis.*
  mBio 12(5):e02060-21. <https://doi.org/10.1128/mBio.02060-21>
  *TnCentral's pages carry an "All Rights Reserved" notice; no terms page exists.
  Opt-in only.*

- **Platon** (sorts every contig of the finished genome into plasmid or chromosome;
  always on, every mode; BacFlux's database page points to its database,
  <https://doi.org/10.5281/zenodo.4066768>) — Schwengers O, Barth P, Falgenhauer L,
  Hain T, Chakraborty T, Goesmann A. (2020) *Platon: identification and
  characterization of bacterial plasmid contigs in short-read draft assemblies
  exploiting protein sequence-based replicon distribution scores.* Microbial Genomics
  6(10):mgen000398. <https://doi.org/10.1099/mgen.0.000398>
  Platon's database carries PlasmidFinder's incompatibility-group probes, MOBscan's
  relaxase HMMs and MOB-suite's oriT sequences, and Platon asks you to cite those too:
  Carattoli A, Zankari E, García-Fernández A, Voldby Larsen M, Lund O, Villa L, et al.
  (2014) *In Silico Detection and Typing of Plasmids using PlasmidFinder and Plasmid
  Multilocus Sequence Typing.* Antimicrobial Agents and Chemotherapy 58(7):3895–3903.
  <https://doi.org/10.1128/AAC.02412-14>; Garcillán-Barcia MP, Redondo-Salvo S,
  Vielva L, de la Cruz F. (2020) *MOBscan: Automated Annotation of MOB Relaxases.* In:
  Horizontal Gene Transfer. Methods in Molecular Biology 2075:295–308.
  <https://doi.org/10.1007/978-1-4939-9877-7_21>; Robertson J, Nash JHE. (2018)
  *MOB-suite: software tools for clustering, reconstruction and typing of plasmids
  from draft assemblies.* Microbial Genomics 4(8):e000206.
  <https://doi.org/10.1099/mgen.0.000206>

- **VirSorter2** (the default virus and prophage caller, every mode;
  `phage.caller: virsorter2`) — Guo J, Bolduc B, Zayed AA, Varsani A,
  Dominguez-Huerta G, Delmont TO, et al. (2021) *VirSorter2: a multi-classifier,
  expert-guided approach to detect diverse DNA and RNA viruses.* Microbiome 9(1):37.
  <https://doi.org/10.1186/s40168-020-00990-y>
  Its authors also publish a protocol for running it, which BacFlux's phage page
  links: Guo J, Vik D, Adjie Pratama A, Roux S, Sullivan M. (2021) *Viral sequence
  identification SOP with VirSorter2 v3.* protocols.io.
  <https://doi.org/10.17504/protocols.io.bwm5pc86>

- **CheckV** (grades the completeness and contamination of whatever the virus caller
  found, every mode) — Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S,
  Kyrpides NC. (2021) *CheckV assesses the quality and completeness of
  metagenome-assembled viral genomes.* Nature Biotechnology 39(5):578–585.
  <https://doi.org/10.1038/s41587-020-00774-7>

- **geNomad** (opt-in, `phage.caller: genomad`: replaces VirSorter2 as the virus
  caller, and its plasmid calls are compared with Platon's) — Camargo AP, Roux S,
  Schulz F, Babinski M, Xu Y, Hu B, et al. (2024) *Identification of mobile genetic
  elements with geNomad.* Nature Biotechnology 42(8):1303–1312.
  <https://doi.org/10.1038/s41587-023-01953-y>
  *⚠ Licence: Berkeley Lab **academic / non-commercial use only** — not BSD, despite
  the bioconda recipe's tag. Opt-in, default off; VirSorter2 is the default phage
  caller. See [Licensing](https://ilivius.github.io/BacFlux/about/licensing/).*

### Supporting tools

- **HMMER** (version 3.4; the profile-HMM search inside many of the tools BacFlux
  runs: CheckM, antiSMASH, Platon, CheckV and the AMRFinderPlus inside Bakta in every
  mode; GTDB-Tk for a genome its skani screen does not assign; Bakta's own searches and
  run_dbCAN through the PyHMMER library; VirSorter2, the default virus caller; and
  ISEScan, MacSyFinder and AMRFinderPlus in the mobilome module. antiSMASH also runs
  HMMER 2.3.2's `hmmpfam2` for some of its analyses) — its authors ask for the web
  site and the version: <http://hmmer.org/>, HMMER 3.4. Where a journal needs a
  paper: Eddy SR. (2011) *Accelerated profile HMM searches.* PLoS Computational
  Biology 7(10):e1002195. <https://doi.org/10.1371/journal.pcbi.1002195>
  PyHMMER asks to be cited together with HMMER: Larralde M, Zeller G. (2023)
  *PyHMMER: a Python library binding to HMMER for efficient sequence analysis.*
  Bioinformatics 39(5):btad214. <https://doi.org/10.1093/bioinformatics/btad214>

- **BLAST+** (megablast against NCBI core_nt for the decontamination screens and the
  plasmid check, and blastn-short for the adapter screen, in every mode; with the
  mobilome's naming layers on, it also searches TnCentral and ICEberg; Bakta,
  AMRFinderPlus, ABRicate, Platon and ISEScan run it themselves too) — Camacho C,
  Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. (2009) *BLAST+:
  architecture and applications.* BMC Bioinformatics 10:421.
  <https://doi.org/10.1186/1471-2105-10-421>

- **DIAMOND** (protein alignment, every mode: BacFlux builds CheckV's DIAMOND index
  with it, eggNOG-mapper runs in its DIAMOND mode, run_dbCAN uses it, and Bakta,
  Platon, CheckV and antiSMASH (its tailoring-enzyme search) search with it) —
  Buchfink B, Reuter K, Drost H-G. (2021)
  *Sensitive protein alignments at tree-of-life scale using DIAMOND.* Nature Methods
  18(4):366–368. <https://doi.org/10.1038/s41592-021-01101-x>

- **MultiQC** (the aggregated HTML report, in every mode) — Ewels P, Magnusson M,
  Lundin S, Käller M. (2016) *MultiQC: summarize analysis results for multiple tools
  and samples in a single report.* Bioinformatics 32(19):3047–3048.
  <https://doi.org/10.1093/bioinformatics/btw354>

- **Snakemake** (the workflow engine; BacFlux needs 9.14.6 or newer and its
  installation commands install 9.27.1) — Mölder F, Jablonski KP, Letcher B,
  Hall MB, Tomkins-Tinch CH, Sochat V, et al. (2021) *Sustainable data analysis
  with Snakemake.* F1000Research 10:33.
  <https://doi.org/10.12688/f1000research.29032.2>

- **Bioconda** (the channel the Snakemake launcher and almost every per-rule conda
  environment install from) — Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J,
  Tomkins-Tinch CH, et al. (2018) *Bioconda: sustainable and comprehensive software
  distribution for the life sciences.* Nature Methods 15(7):475–476.
  <https://doi.org/10.1038/s41592-018-0046-7>

- **conda-forge** (the channel every per-rule environment and the Snakemake launcher
  install from, beside Bioconda) — conda-forge community. (2015) *The conda-forge
  Project: Community-based Software Distribution Built on the conda Package Format
  and Ecosystem.* Zenodo. <https://doi.org/10.5281/zenodo.4774216>

---

## 2. Design influence

Work that shaped how BacFlux is built, without contributing code to it.

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

What we took from them — all of it facts and design decisions, not code (see
[Licensing](https://ilivius.github.io/BacFlux/about/licensing/#what-this-repository-deliberately-does-not-contain)):

- that **ICEscan** exists, ships the IME and AICE models, and is worth running;
- Sequence Ontology terms for element types (BacFlux's `mge_type` column);
- the `contig_id|mge_type-start:end` element ID format;
- the discard-with-reason file pattern;
- the 0.9 CDS-coverage assignment threshold. Their 500 bp minimum element size is
  deliberately not applied to IS calls, because a short call is usually a partial IS
  cut by a contig break.

**No code was copied.** Three of their scripts —
`bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py` and
`bin/prescan_to_fasta.py` — are derived from ICEfinder2 and carry
**CC BY-NC-SA 4.0** headers. That licence is non-commercial and share-alike, so
copying or adapting any of it would contaminate BacFlux's MIT licence and impose
a non-commercial restriction on everyone downstream. Reading it to understand an
algorithm is fine and is what was done; the separation is real and deliberate.

**Escape hatch for users.** For full ICE delimitation and deeper mobilome
analysis, run their pipeline on the finished genome BacFlux delivers,
`02.assembly/{sample}/contigs_final.fasta`; consuming its `mobilome.gff.gz` output
carries no licence consequences.

Their pipeline also runs AMRFinderPlus by default, alongside DeepARG, RGI/CARD and
hAMRonization, and integrates the results with the mobilome. The difference in scope
set out above rests on metagenomes and MAGs versus single isolates, and on nothing
else.

### ICEfinder2 and the att-site literature

ICEfinder2's approach to boundary detection (variable-length repeat search rather
than a fixed-length probe) was established by reading its source and the
surrounding literature, then implemented independently. The evidence, the
parameters and the negative results are recorded in
`docs/methods_att_and_small_plasmids.md`.

### NCBI FCS-adaptor and the adapter screen

> National Center for Biotechnology Information. *Foreign Contamination Screen (FCS),
> FCS-adaptor.* <https://github.com/ncbi/fcs>

The adapter screen of the finished genome counts an adapter hit as being at a contig
end when it lies within 100 bp of the end of a linear contig: the distance within
which FCS-adaptor extends a trim to the contig end
(<https://github.com/ncbi/fcs/wiki/FCS-adaptor-output>). BacFlux runs no FCS code and
ships none of its data; run FCS-adaptor itself before a GenBank submission.

### IS calling accuracy

- **digIS** (not run by BacFlux) — Puterová J, Martínek T. (2021) *digIS: towards
  detecting distant and putative novel insertion sequence elements in prokaryotic
  genomes.* BMC Bioinformatics 22:258. <https://doi.org/10.1186/s12859-021-04177-6>
  Cited on the mobilome page for the IS false-discovery rates it measured, ISEScan's
  among them, and for ISEScan's partial hits over DnaA, DnaC and DnaI: part of the
  evidence behind the confidence level on every IS call and the removal of false
  IS21 calls over the replication genes.

### Regulatory framing

- **EFSA Scientific Committee** (2025) *Guidance on the characterisation of
  microorganisms in support of the risk assessment of products used in the food
  chain.* EFSA Journal 23(11):e9705. <https://doi.org/10.2903/j.efsa.2025.9705>
  Its Appendix A lists the chapters of the 2018 FEEDAP guidance (EFSA Journal
  16(3):5206, <https://doi.org/10.2903/j.efsa.2018.5206>) that it may supersede,
  in whole or in part.
- **EFSA BIOHAZ Panel** (2023) *Statement on how to interpret the QPS qualification
  on 'acquired antimicrobial resistance genes'.* EFSA Journal 21(10):8323.
  <https://doi.org/10.2903/j.efsa.2023.8323>
- **EFSA** (2024) *EFSA statement on the requirements for whole genome sequence
  analysis of microorganisms intentionally used in the food chain.* EFSA Journal
  22(8):e8912. <https://doi.org/10.2903/j.efsa.2024.8912> — the source of the 80%
  identity and 70% coverage reporting thresholds the ABRicate screen uses, and of
  the requirement to search at least two maintained, curated AMR databases.

BacFlux provides **supporting evidence** for the intrinsic-versus-acquired
distinction these documents require. It is not "EFSA-compliant" and does not
claim to be, and every mobility call is a **prediction** — confirmation requires
a filter or broth mating assay.
