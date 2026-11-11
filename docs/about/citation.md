# Citation and references

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR** | antimicrobial resistance |
    | **WGS** | whole-genome sequencing |
    | **EFSA** | European Food Safety Authority |
    | **QPS** | Qualified Presumption of Safety, EFSA's assessment framework for listed taxonomic units subject to qualifications |
    | **ICE** | integrative and conjugative element — moves itself between cells |
    | **IME** | integrative mobilisable element — needs a helper element to move |
    | **MAG** | metagenome-assembled genome, a genome reconstructed from a mixed community |

## Citing BacFlux

> Antonielli, L. (2026). *BacFlux: a unified Snakemake workflow for bacterial
> whole-genome sequencing.* Zenodo. <https://doi.org/10.5281/zenodo.11143917>

That DOI is the **concept DOI** — it always resolves to the newest release. Every
tagged release also gets its own **version DOI** on Zenodo. Cite the version DOI
when the exact code matters, which for a methods section it usually does.

The first Zenodo record, v1.1.5 (2024,
<https://doi.org/10.5281/zenodo.11143918>), also credits D. K. Großkinsky, H. Koch,
F. Trognitz, A. Sanchez Mejia and M. Nagel; cite it for work done with that version.

Machine-readable metadata is in
[`CITATION.cff`](https://github.com/iLivius/BacFlux/blob/main/CITATION.cff);
GitHub reads that file to build a citation for you.

!!! note "The retired sibling workflows keep their DOIs"
    v2.0.0 folds `BacFluxL` and `BacFluxL+` into this repository. Their Zenodo
    records stay valid, so an analysis already published with one of them remains
    citable and reproducible — cite the record you actually ran. See
    [Coming from v1](../getting-started/from-v1.md).

## Cite the tools as well

BacFlux combines established tools and reference databases with its own processing
and interpretation rules. Cite the components that contributed to the reported
results, as well as the workflow.

Which references you need depends on the mode you ran and on which optional
stages you switched on, so the tables below are grouped that way. The complete
list, with the licensing notes attached, is in
[`CITATIONS.md`](https://github.com/iLivius/BacFlux/blob/main/CITATIONS.md).

### Shared stages { #every-run-whatever-the-mode }

| Step | Tool or database | Reference |
|---|---|---|
| Workflow engine | Snakemake; Bioconda (one conda environment per rule) | 44, 30 |
| Decontamination | BLAST+ against NCBI core nt; BlobTools | 10, 38 |
| Coverage track for BlobTools | Bowtie 2 (short reads), minimap2 (long reads or contigs), SAMtools | 39, 40, 42, 21 |
| Assembly statistics | QUAST 5.3.0 | 32, 43 |
| Completeness and contamination | CheckM 1.2.5 | 52 |
| Mapping QC (modes with reads) | Qualimap 2.3 | 48 |
| Adapter screen of the finished genome | BLAST+ against ONT's adapters and barcodes (ONT Chemistry Technical Document) and the PacBio, MGI and Illumina adapters of NCBI UniVec | 10, 50, 45 |
| Taxonomic placement | GTDB-Tk 2.7.2 against GTDB **R232** | 14, 13, 53, 54 |
| Annotation | Bakta 1.12.1 (database **v6.0**) | 58 |
| Functional annotation | eggNOG-mapper 2.1.15 against eggNOG 5.0.2, in its DIAMOND mode | 12, 35, 8 |
| Secondary metabolites | antiSMASH 8.0.4 | 3 |
| CAZymes | run_dbCAN 5.1.2 (dbCAN3) | 68, 67, 34 |
| AMR and virulence on contigs | ABRicate 1.4.0 — plus the databases below | 60 |
| Plasmids | Platon 1.8, with a supplementary BLAST+ line beside each call | 57, 10 |
| Prophages | VirSorter2 2.2.4 (default) or geNomad 1.12.0 (`phage.caller: genomad`), then CheckV 1.1.1, with the DIAMOND index BacFlux builds for it, to assess the predicted viral sequences | 31, 11, 46, 8 |
| Aggregated report | MultiQC 1.33 | 26 |

**ABRicate's four databases.** The screen runs once per database and the four
tables are kept side by side, so cite the ones you actually read.

| `db` | Database | Reference |
|---|---|---|
| `card` | CARD | 2 |
| `ncbi` | NCBI AMRFinderPlus Reference Gene Catalog | 28 |
| `resfinder` | ResFinder | 66, 4 |
| `vfdb` | VFDB | 15, 69 |

### What your mode adds

| `mode:` | Front end | Reference |
|---|---|---|
| `illumina` | Bowtie 2 (PhiX removal), fastp 1.0.1, SPAdes 4.3.0 | 39, 40, 17, 18, 55 |
| `nanopore` | NanoPlot 1.47.4 (NanoPack2), fastplong 0.7.1 (adapter and barcode trimming; `parameters.{nanopore,hybrid}.ont_adapter: off` skips it), Filtlong 0.3.1, Flye 2.9.6, dnaapler 1.4.0, Medaka 2.2.2 *(optional — `parameters.{nanopore,hybrid}.medaka_model: false` skips it)* | 22, 49, 16, 18, 63, 37, 5, 51 |
| `hybrid` | both front ends above — the short reads are assembled too, as the comparator genome — then BWA + Polypolish 0.7.1 to correct the ONT assembly with them, and Snippy 4.6.0 to compare the two | the two rows above, plus 41, 64, 6, 59 |
| `contigs` | no assembler — contig filtering only, with minimap2 self-mapping for the coverage track | 42 |

**Short-read modes only.** `illumina` and `hybrid` also map the trimmed reads
onto CARD v4.0.2 with BBMap 39.33 (refs 9, 2). Read
[Antimicrobial resistance](../analysis/amr.md) before quoting that analysis: it is a
sensitive screen whose specificity comes from the reference-length coverage gate
rather than from read identity; it does not confirm phenotypic resistance.

### Optional stages and layers

Each of these is off until you turn it on, and each brings its own citations.

| Switch | What runs | Reference |
|---|---|---|
| `phage.caller: genomad` | geNomad 1.12.0 replaces VirSorter2 as the virus caller, and its plasmid calls are compared with Platon's | 11 |
| `baktfold.run: true` | Baktfold 0.3.0, which turns each of Bakta's hypothetical proteins into a string of Foldseek's 3Di structural alphabet with the ProstT5 language model, then searches them with Foldseek against SwissProt, AlphaFold Database clusters, PDB and CATH to suggest names — cite all three; if you report which database a name came from, cite that database as well (Baktfold's README lists the papers) | 7, 33, 61 |
| `traits.run: true` | the trait profile: every protein of up to 4,000 amino acids is turned into a string of Foldseek's 3Di structural alphabet by Baktfold's ProstT5 step, then searched with Foldseek against EcoFoldDB v2.0 — cite EcoFoldDB, Baktfold, ProstT5 and Foldseek; the papers behind each rule are in the `reference` column of its rules table ([Trait profile](../analysis/traits.md)) | 29, 7, 33, 61 |
| `mobilome.run: true` | AMRFinderPlus, run from the Bakta environment against the AMRFinderPlus database that ships inside `bakta_db`; ISEScan 1.7.3 and MacSyFinder 2.1.6, both HMMER searches; and the **CONJscan** models MacSyFinder runs — cite all three CONJscan papers; they are the models that run | 28, 27, 65, 47, 23, 1, 19, 20 |
| `mobilome.icescan.run: true` | the **ICEscan** model set, distributed with ICEfinder2 and itself a fork of CONJscan — so cite the ICEberg paper *and* the three CONJscan papers | 62, 1, 19, 20 |
| `mobilome.tncentral.url` or `.dir` | TnCentral, to name curated transposons and integrons | 56 |
| `mobilome.iceberg.urls` or `.dir` | ICEberg 3.0 sequences, to name the ICE and IME candidates | 62 |
| `mobilome.isosdb.fasta_url` + `family_map_url`, or `.dir` | ISOSDB read mapping with BBMap, for IS copy number | 36, 9 |

The mobility tiers are organised around EFSA's distinction between intrinsic and
acquired resistance genes (refs 24, 25). They provide genomic-context evidence; they
do not establish that a gene is intrinsic and do not demonstrate transfer.
[Antimicrobial resistance](../analysis/amr.md) lists the two further EFSA documents —
the 2024 WGS statement behind the 80%/70% thresholds, and the *Bacillus* catalogue —
and the [intrinsic/acquired methods page](../methods_amr_intrinsic_acquired.md) sets
out what the evidence supports.

!!! note "ISfinder reaches a run only through Bakta's database"
    Bakta's database integrates annotations from the ISfinder database, so cite
    ISfinder when you report Bakta's insertion-sequence annotations. BacFlux itself
    downloads nothing from ISfinder: TnCentral publishes further endpoints whose
    content carries ISfinder's terms, and those are deliberately not wired into
    `config.yaml`. If you point `mobilome.tncentral.url` at one of them yourself,
    cite ISfinder as well:

    > Siguier, P., Perochon, J., Lestrade, L., Mahillon, J., & Chandler, M.
    > (2006). ISfinder: the reference centre for bacterial insertion sequences.
    > *Nucleic Acids Research*, 34(Database issue), D32–D36.
    > <https://doi.org/10.1093/nar/gkj014>

    Terms for every optional layer are listed on [Licensing](licensing.md); what
    each layer adds is on
    [Optional layers](../mobilome/optional-layers.md).

## Acknowledgements

This work was originally supported by the
[Austrian Science Fund (FWF)](https://www.fwf.ac.at/en/) under Project I6030-B.

Much of v2.0.0 — merging the four v1 workflows into one, building the mobilome
module and writing this documentation — was done with
[Claude Code](https://claude.com/claude-code). Anthropic provided six months of
Claude Max through their Open Source programme, and the scale of the v2 rewrite
would not have been realistic without it.

## References

1. Abby, S. S., Cury, J., Guglielmini, J., Néron, B., Touchon, M., & Rocha, E. P. C. (2016). Identification of protein secretion systems in bacterial genomes. *Scientific Reports*, 6, 23080. <https://doi.org/10.1038/srep23080> (CONJscan model set)
2. Alcock, B. P., et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D690–D699. <https://doi.org/10.1093/nar/gkac920>
3. Blin, K., et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. *Nucleic Acids Research*, 53(W1), W32–W38. <https://doi.org/10.1093/nar/gkaf334>
4. Bortolaia, V., et al. (2020). ResFinder 4.0 for predictions of phenotypes from genotypes. *Journal of Antimicrobial Chemotherapy*, 75(12), 3491–3500. <https://doi.org/10.1093/jac/dkaa345> (the ResFinder paper ABRicate's README gives for its `resfinder` database)
5. Bouras, G., Grigson, S. R., Papudeshi, B., Mallawaarachchi, V., & Roach, M. J. (2024). Dnaapler: a tool to reorient circular microbial genomes. *Journal of Open Source Software*, 9(93), 5968. <https://doi.org/10.21105/joss.05968>
6. Bouras, G., Judd, L. M., Edwards, R. A., Vreugde, S., Stinear, T. P., & Wick, R. R. (2024). How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. *Microbial Genomics*, 10(6), 001254. <https://doi.org/10.1099/mgen.0.001254> (Polypolish 0.6 and later, whose `filter` and `polish` steps BacFlux runs)
7. Bouras, G., et al. (2026). Baktfold: sensitive protein functional annotation across the microbial tree of life using structural information. *bioRxiv* (preprint). <https://doi.org/10.64898/2026.03.31.715528>
8. Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods*, 18(4), 366–368. <https://doi.org/10.1038/s41592-021-01101-x>
9. Bushnell, B. (2014). *BBMap: a fast, accurate, splice-aware aligner.* 9th Annual Genomics of Energy & Environment Meeting, Walnut Creek, CA (report LBNL-7065E). <https://www.osti.gov/biblio/1241166>
10. Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. <https://doi.org/10.1186/1471-2105-10-421>
11. Camargo, A. P., et al. (2024). Identification of mobile genetic elements with geNomad. *Nature Biotechnology*, 42(8), 1303–1312. <https://doi.org/10.1038/s41587-023-01953-y>
12. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution*, 38(12), 5825–5829. <https://doi.org/10.1093/molbev/msab293>
13. Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics*, 36(6), 1925–1927. <https://doi.org/10.1093/bioinformatics/btz848>
14. Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022). GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics*, 38(23), 5315–5316. <https://doi.org/10.1093/bioinformatics/btac672>
15. Chen, L., Zheng, D., Liu, B., Yang, J., & Jin, Q. (2016). VFDB 2016: hierarchical and refined dataset for big data analysis — 10 years on. *Nucleic Acids Research*, 44(D1), D694–D697. <https://doi.org/10.1093/nar/gkv1239>
16. Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. *iMeta*, 2(2), e107. <https://doi.org/10.1002/imt2.107> (cited for fastplong, as its authors ask)
17. Chen, S. (2025). fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta*, 4(5), e70078. <https://doi.org/10.1002/imt2.70078> (the paper fastp's README asks for)
18. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>
19. Coluzzi, C., Garcillán-Barcia, M. P., de la Cruz, F., & Rocha, E. P. C. (2022). Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative plasmids. *Molecular Biology and Evolution*, 39(6), msac115. <https://doi.org/10.1093/molbev/msac115> (CONJscan model set)
20. Cury, J., Touchon, M., & Rocha, E. P. C. (2017). Integrative and conjugative elements and their hosts: composition, distribution and organization. *Nucleic Acids Research*, 45(15), 8943–8956. <https://doi.org/10.1093/nar/gkx607> (CONJscan model set)
21. Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. <https://doi.org/10.1093/gigascience/giab008>
22. De Coster, W., & Rademakers, R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics*, 39(5), btad311. <https://doi.org/10.1093/bioinformatics/btad311> (NanoPlot)
23. Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology*, 7(10), e1002195. <https://doi.org/10.1371/journal.pcbi.1002195> (the profile-HMM engine inside CheckM, GTDB-Tk, antiSMASH, Platon, CheckV, AMRFinderPlus, VirSorter2, ISEScan and MacSyFinder, and, through PyHMMER, inside Bakta and run_dbCAN)
24. EFSA Panel on Biological Hazards (BIOHAZ). (2023). Statement on how to interpret the QPS qualification on 'acquired antimicrobial resistance genes'. *EFSA Journal*, 21(10), 8323. <https://doi.org/10.2903/j.efsa.2023.8323> (the method for telling intrinsic from acquired)
25. EFSA Scientific Committee. (2025). Guidance on the characterisation of microorganisms in support of the risk assessment of products used in the food chain. *EFSA Journal*, 23(11), e9705. <https://doi.org/10.2903/j.efsa.2025.9705> (the cited 2025 edition; Appendix A describes its relationship to the 2018 FEEDAP guidance)
26. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. <https://doi.org/10.1093/bioinformatics/btw354>
27. Feldgarden, M., et al. (2019). Validating the AMRFinder tool and resistance gene database by using antimicrobial resistance genotype–phenotype correlations in a collection of isolates. *Antimicrobial Agents and Chemotherapy*, 63(11), e00483-19. <https://doi.org/10.1128/AAC.00483-19>
28. Feldgarden, M., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728. <https://doi.org/10.1038/s41598-021-91456-0>
29. Ghaly, T. M., Rajabal, V., Russell, D., Colombi, E., & Tetu, S. G. (2025). EcoFoldDB: protein structure-guided functional profiling of ecologically relevant microbial traits at the metagenome scale. *Environmental Microbiology*, 27(9), e70178. <https://doi.org/10.1111/1462-2920.70178> (database: EcoFoldDB v2.0, <https://doi.org/10.5281/zenodo.16908168>, CC BY 4.0)
30. Grüning, B., et al. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nature Methods*, 15, 475–476. <https://doi.org/10.1038/s41592-018-0046-7>
31. Guo, J., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome*, 9(1), 37. <https://doi.org/10.1186/s40168-020-00990-y>
32. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. <https://doi.org/10.1093/bioinformatics/btt086>
33. Heinzinger, M., et al. (2024). Bilingual language model for protein sequence and structure. *NAR Genomics and Bioinformatics*, 6(4), lqae150. <https://doi.org/10.1093/nargab/lqae150> (ProstT5)
34. Huang, L., et al. (2018). dbCAN-seq: a database of carbohydrate-active enzyme (CAZyme) sequence and annotation. *Nucleic Acids Research*, 46(D1), D516–D521. <https://doi.org/10.1093/nar/gkx894>
35. Huerta-Cepas, J., et al. (2019). eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. *Nucleic Acids Research*, 47(D1), D309–D314. <https://doi.org/10.1093/nar/gky1085>
36. Kirsch, J. M., Hryckowian, A. J., & Duerkop, B. A. (2024). A metagenomics pipeline reveals insertion sequence-driven evolution of the microbiota. *Cell Host & Microbe*, 32(5), 739–754.e4. <https://doi.org/10.1016/j.chom.2024.03.005> (pseudoR, the source of ISOSDB)
37. Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. <https://doi.org/10.1038/s41587-019-0072-8> (Flye)
38. Laetsch, D. R., & Blaxter, M. L. (2017). BlobTools: interrogation of genome assemblies. *F1000Research*, 6, 1287. <https://doi.org/10.12688/f1000research.12232.1>
39. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359. <https://doi.org/10.1038/nmeth.1923>
40. Langmead, B., Wilks, C., Antonescu, V., & Charles, R. (2019). Scaling read aligners to hundreds of threads on general-purpose processors. *Bioinformatics*, 35(3), 421–432. <https://doi.org/10.1093/bioinformatics/bty648>
41. Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, 1303.3997v2 [q-bio.GN]. <https://doi.org/10.48550/arXiv.1303.3997> (BWA-MEM, the BWA algorithm that aligns the reads Polypolish reads)
42. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. <https://doi.org/10.1093/bioinformatics/bty191>
43. Mikheenko, A., Saveliev, V., Hirsch, P., & Gurevich, A. (2023). WebQUAST: online evaluation of genome assemblies. *Nucleic Acids Research*, 51(W1), W601–W606. <https://doi.org/10.1093/nar/gkad406> (the paper QUAST asks for with its latest version)
44. Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. <https://doi.org/10.12688/f1000research.29032.2>
45. National Center for Biotechnology Information. *The UniVec Database.* <https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/> (files at <https://ftp.ncbi.nlm.nih.gov/pub/UniVec/>; the Illumina records of build 10.1, downloaded by the adapter screen, and the PacBio and MGI records shipped with BacFlux)
46. Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology*, 39(5), 578–585. <https://doi.org/10.1038/s41587-020-00774-7>
47. Néron, B., Denise, R., Coluzzi, C., Touchon, M., Rocha, E. P. C., & Abby, S. S. (2023). MacSyFinder v2: improved modelling and search engine to identify molecular systems in genomes. *Peer Community Journal*, 3, e28. <https://doi.org/10.24072/pcjournal.250>
48. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. *Bioinformatics*, 32(2), 292–294. <https://doi.org/10.1093/bioinformatics/btv566>
49. OpenGene. *fastplong: ultrafast preprocessing and quality control for long reads.* <https://github.com/OpenGene/fastplong> (cite it through the two fastp papers, as its README asks)
50. Oxford Nanopore Technologies. *Chemistry Technical Document*, version CHTD_500_v1_revAV_24July2026, sections 14 and 15. <https://nanoporetech.com/document/chemistry-technical-document> (the ONT adapter and barcode sequences shipped with BacFlux)
51. Oxford Nanopore Technologies. *Medaka.* <https://github.com/nanoporetech/medaka>
52. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043–1055. <https://doi.org/10.1101/gr.186072.114>
53. Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., & Hugenholtz, P. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*, 50(D1), D785–D794. <https://doi.org/10.1093/nar/gkab776> (cite the release you ran — BacFlux expects **R232**)
54. Parks, D. H., Chaumeil, P.-A., Mussig, A. J., Rinke, C., Chuvochina, M., & Hugenholtz, P. (2026). GTDB release 10: a complete and systematic taxonomy for 715 230 bacterial and 17 245 archaeal genomes. *Nucleic Acids Research*, 54(D1), D743–D754. <https://doi.org/10.1093/nar/gkaf1040> (GTDB's newest paper; name the release you ran, R232, in your methods)
55. Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. *Current Protocols in Bioinformatics*, 70(1), e102. <https://doi.org/10.1002/cpbi.102> (the paper SPAdes asks for)
56. Ross, K., et al. (2021). TnCentral: a prokaryotic transposable element database and web portal for transposon analysis. *mBio*, 12(5), e02060-21. <https://doi.org/10.1128/mBio.02060-21>
57. Schwengers, O., Barth, P., Falgenhauer, L., Hain, T., Chakraborty, T., & Goesmann, A. (2020). Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. *Microbial Genomics*, 6(10), mgen000398. <https://doi.org/10.1099/mgen.0.000398>
58. Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*, 7(11), 000685. <https://doi.org/10.1099/mgen.0.000685>
59. Seemann, T. (2015). *snippy: fast bacterial variant calling from NGS reads.* <https://github.com/tseemann/snippy>
60. Seemann, T. (2020). *ABRicate: mass screening of contigs for antimicrobial and virulence genes.* <https://github.com/tseemann/abricate>
61. van Kempen, M., et al. (2024). Fast and accurate protein structure search with Foldseek. *Nature Biotechnology*, 42(2), 243–246. <https://doi.org/10.1038/s41587-023-01773-0>
62. Wang, M., et al. (2024). ICEberg 3.0: functional categorization and analysis of the integrative and conjugative elements in bacteria. *Nucleic Acids Research*, 52(D1), D732–D737. <https://doi.org/10.1093/nar/gkad935> (also the paper to cite for ICEfinder2 / ICEscan)
63. Wick, R. R. *Filtlong: quality filtering tool for long reads.* <https://github.com/rrwick/Filtlong>
64. Wick, R. R., & Holt, K. E. (2022). Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Computational Biology*, 18(1), e1009802. <https://doi.org/10.1371/journal.pcbi.1009802>
65. Xie, Z., & Tang, H. (2017). ISEScan: automated identification of insertion sequence elements in prokaryotic genomes. *Bioinformatics*, 33(21), 3340–3347. <https://doi.org/10.1093/bioinformatics/btx433>
66. Zankari, E., et al. (2012). Identification of acquired antimicrobial resistance genes. *Journal of Antimicrobial Chemotherapy*, 67(11), 2640–2644. <https://doi.org/10.1093/jac/dks261> (ResFinder)
67. Zhang, H., et al. (2018). dbCAN2: a meta server for automated carbohydrate-active enzyme annotation. *Nucleic Acids Research*, 46(W1), W95–W101. <https://doi.org/10.1093/nar/gky418>
68. Zheng, J., Ge, Q., Yan, Y., Zhang, X., Huang, L., & Yin, Y. (2023). dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research*, 51(W1), W115–W121. <https://doi.org/10.1093/nar/gkad328>
69. Zhou, S., Liu, B., Zheng, D., Chen, L., & Yang, J. (2025). VFDB 2025: an integrated resource for exploring anti-virulence compounds. *Nucleic Acids Research*, 53(D1), D871–D877. <https://doi.org/10.1093/nar/gkae968> (VFDB's newest paper)

## Design influence

Work that shaped BacFlux without contributing code to it; the distinction is
spelled out on [Licensing](licensing.md).

- **EBI Mobilome Annotation Pipeline** (Apache-2.0).
  <https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline> — from which
  BacFlux adopted decisions rather than code: that ICEscan exists and is worth
  running, Sequence Ontology terms for element types (BacFlux's `mge_type`
  column), the `contig_id|mge_type-start:end` element ID format, the
  discard-with-reason file and the 0.9 CDS-coverage rule. Their 500 bp minimum
  element size is deliberately not applied to IS calls, because a short call is
  usually a partial IS cut by a contig break. The two pipelines answer different
  questions: theirs goes deep on the mobilome of metagenomes and MAGs, BacFlux
  covers reads through to AMR mobility for single isolates.
- **ICEfinder2**, whose variable-length repeat search — rather than a
  fixed-length probe — is why BacFlux's att-site search works the way it does.
  Established by reading it, then implemented independently.
- **NCBI FCS-adaptor** (<https://github.com/ncbi/fcs>), whose rule for a contig end
  the adapter screen follows: an adapter hit within 100 bp of the end of a linear
  contig counts as being at that end. BacFlux runs no FCS code; run FCS-adaptor
  itself before a GenBank submission.

- Puterová, J. & Martínek, T. (2021) digIS: towards detecting distant and putative
  novel insertion sequence elements in prokaryotic genomes. *BMC Bioinformatics*
  22:258. <https://doi.org/10.1186/s12859-021-04177-6> — benchmarked here for the
  IS false-discovery figures quoted on the mobilome pages, not run by BacFlux.

## Licence

BacFlux is released under the
[MIT Licence](https://github.com/iLivius/BacFlux/blob/main/LICENSE). The tools,
model sets and databases it invokes are distributed under their own terms, which
are listed on [Licensing](licensing.md).
