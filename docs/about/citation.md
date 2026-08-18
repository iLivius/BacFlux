# Citation and references

## Citing BacFlux

> Antonielli, L., Großkinsky, D. K., Koch, H., Trognitz, F., Sanchez Mejia, A., &
> Nagel, M. (2024). *BacFlux: A workflow for bacterial short-read assembly, QC,
> annotation, and more.* Zenodo. <https://doi.org/10.5281/zenodo.11143917>

That DOI is the **concept DOI** — it always resolves to the newest release. Every
tagged release also gets its own **version DOI** on Zenodo. Cite the version DOI
when the exact code matters, which for a methods section it usually does.

Machine-readable metadata is in
[`CITATION.cff`](https://github.com/iLivius/BacFlux/blob/main/CITATION.cff);
GitHub and Zenodo read that file to build a citation for you.

!!! note "The retired sibling workflows keep their DOIs"
    v2.0.0 folds `BacFluxL` and `BacFluxL+` into this repository. Their Zenodo
    records stay valid, so an analysis already published with one of them remains
    citable and reproducible — cite the record you actually ran. See
    [Coming from v1](../getting-started/from-v1.md).

## Cite the tools as well

BacFlux is glue. Almost all of the science it reports comes from other people's
tools and databases, and the mobilome module in particular is an **integrator**:
it decides how to combine other tools' calls, and every underlying detection
belongs to somebody else.

Which references you need depends on the mode you ran and on which optional
stages you switched on, so the tables below are grouped that way. The complete
list, with the licensing notes attached, is in
[`CITATIONS.md`](https://github.com/iLivius/BacFlux/blob/main/CITATIONS.md).

### Every run, whatever the mode

| Step | Tool or database | Reference |
|---|---|---|
| Workflow engine | Snakemake; Bioconda (one conda environment per rule) | 37, 24 |
| Decontamination | BLAST+ against NCBI core nt; BlobTools | 7, 33 |
| Coverage track for BlobTools | Bowtie 2 (short reads), minimap2 (long reads or contigs), SAMtools | 34, 35, 15 |
| Assembly statistics | QUAST 5.3.0 | 27 |
| Completeness and contamination | CheckM 1.2.5 | 42 |
| Mapping QC (modes with reads) | QualiMap 2.3 | 40 |
| Taxonomic placement | GTDB-Tk 2.7.2 against GTDB **R232** | 10, 43 |
| Annotation | Bakta 1.12.1 (database **v6.0**) | 46 |
| Functional annotation | eggNOG-mapper 2.1.15 against eggNOG 5.0.2 | 9, 28 |
| Secondary metabolites | antiSMASH 8.0.4 | 4 |
| CAZymes | run_dbCAN 5.1.2 (dbCAN3) | 54 |
| AMR and virulence on contigs | ABRicate 1.2.0 — plus the databases below | 47 |
| Plasmids | Platon 1.8, with a supplementary BLAST+ line beside each call | 45, 7 |
| Prophages | VirSorter2 2.2.4, with CheckV 1.0.3 grading what it found | 25, 38 |
| Aggregated report | MultiQC 1.33 | 21 |

**ABRicate's eight databases.** The screen runs once per database and the eight
tables are kept side by side, so cite the ones you actually read.

| `db` | Database | Reference |
|---|---|---|
| `argannot` | ARG-ANNOT | 26 |
| `card` | CARD | 2, 30 |
| `ecoh` | EcOH — *E. coli* O- and H-antigen loci | 29 |
| `ecoli_vf` | *E. coli* virulence factors | ships with ABRicate — 47 |
| `megares` | MEGARes 2.0 | 17 |
| `ncbi` | NCBI AMRFinderPlus Reference Gene Catalog | 23 |
| `resfinder` | ResFinder | 53 |
| `vfdb` | VFDB | 11 |

### What your mode adds

| `mode:` | Front end | Reference |
|---|---|---|
| `illumina` | Bowtie 2 (PhiX removal), fastp 1.0.1, SPAdes 4.2.0 | 34, 12, 3 |
| `nanopore` | NanoPlot 1.46.2 (NanoPack), Filtlong 0.3.1, Flye 2.9.6, dnaapler 1.4.0, Medaka 2.2.2 *(optional — `parameters.{nanopore,hybrid}.medaka_model: false` skips it)* | 16, 50, 32, 5, 41 |
| `hybrid` | both front ends above — the short reads are assembled too, as the comparator genome — then BWA + Polypolish 0.7.1 to correct the ONT assembly with them, and Snippy 4.6.0 to compare the two | the two rows above, plus 36, 51, 48 |
| `contigs` | no assembler — contig filtering only, with minimap2 self-mapping for the coverage track | 35 |

**Short-read modes only.** `illumina` and `hybrid` also map the trimmed reads
onto CARD v4.0.1 with BBMap 39.33 (refs 6, 2, 30). Read
[Antimicrobial resistance](../analysis/amr.md) before quoting that leg: it is a
sensitive screen, and its specificity comes from the reference-length coverage
gate rather than from read identity.

### Optional stages and layers

Each of these is off until you turn it on, and each brings its own citations.

| Switch | What runs | Reference |
|---|---|---|
| `phage.caller: genomad` | geNomad 1.12.0 replaces VirSorter2 as the virus caller, and its plasmid calls are compared with Platon's | 8 |
| `mobilome.run: true` | AMRFinderPlus, run from the Bakta environment against the AMRFinderPlus database that ships inside `bakta_db`; ISEScan 1.7.3 and MacSyFinder 2.1.6, both HMMER searches; and the **CONJscan** models MacSyFinder runs — cite all three CONJscan papers; they are the models that run | 23, 22, 52, 39, 18, 1, 13, 14 |
| `mobilome.icescan.run: true` | the **ICEscan** model set, distributed with ICEfinder2 and itself a fork of CONJscan — so cite the ICEberg paper *and* the three CONJscan papers | 49, 1, 13, 14 |
| `mobilome.tncentral.url` or `.dir` | TnCentral, to name curated transposons and integrons | 44 |
| `mobilome.iceberg.urls` or `.dir` | ICEberg 3.0 sequences, to name the ICE and IME candidates | 49 |
| `mobilome.isosdb.fasta_url` + `family_map_url`, or `.dir` | ISOSDB read mapping with BBMap, for IS copy number | 31, 6 |

The mobility ladder is built around EFSA's intrinsic-versus-acquired distinction
(refs 19, 20). [Antimicrobial resistance](../analysis/amr.md) lists the two
further EFSA documents — the 2024 WGS statement behind the 80%/70% reporting
thresholds, and the *Bacillus* catalogue — and says what BacFlux may and may not
claim from any of them.

!!! note "ISfinder is not part of a standard run"
    TnCentral publishes further endpoints whose content carries ISfinder's terms.
    Those are deliberately not wired into `config.yaml`, so nothing BacFlux
    fetches by default touches ISfinder data. If you point
    `mobilome.tncentral.url` at one of them yourself, cite ISfinder as well:

    > Siguier, P., Perochon, J., Lestrade, L., Mahillon, J., & Chandler, M.
    > (2006). ISfinder: the reference centre for bacterial insertion sequences.
    > *Nucleic Acids Research*, 34(Database issue), D32–D36.
    > <https://doi.org/10.1093/nar/gkj014>

    Terms for every optional layer are listed on [Licensing](licensing.md); what
    each layer buys you is on
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
3. Bankevich, A., et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology*, 19(5), 455–477. <https://doi.org/10.1089/cmb.2012.0021>
4. Blin, K., et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. *Nucleic Acids Research*, 53(W1), W32–W38. <https://doi.org/10.1093/nar/gkaf334>
5. Bouras, G., Grigson, S. R., Papudeshi, B., Mallawaarachchi, V., & Roach, M. J. (2024). Dnaapler: a tool to reorient circular microbial genomes. *Journal of Open Source Software*, 9(93), 5968. <https://doi.org/10.21105/joss.05968>
6. Bushnell, B. (2014). *BBMap: a fast, accurate, splice-aware aligner.* <https://escholarship.org/uc/item/1h3515gn>
7. Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. <https://doi.org/10.1186/1471-2105-10-421>
8. Camargo, A. P., et al. (2024). Identification of mobile genetic elements with geNomad. *Nature Biotechnology*, 42(8), 1303–1312. <https://doi.org/10.1038/s41587-023-01953-y>
9. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution*, 38(12), 5825–5829. <https://doi.org/10.1093/molbev/msab293>
10. Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022). GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics*, 38(23), 5315–5316. <https://doi.org/10.1093/bioinformatics/btac672>
11. Chen, L., Zheng, D., Liu, B., Yang, J., & Jin, Q. (2016). VFDB 2016: hierarchical and refined dataset for big data analysis — 10 years on. *Nucleic Acids Research*, 44(D1), D694–D697. <https://doi.org/10.1093/nar/gkv1239>
12. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>
13. Coluzzi, C., Garcillán-Barcia, M. P., de la Cruz, F., & Rocha, E. P. C. (2022). Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative plasmids. *Molecular Biology and Evolution*, 39(6), msac115. <https://doi.org/10.1093/molbev/msac115> (CONJscan model set)
14. Cury, J., Touchon, M., & Rocha, E. P. C. (2017). Integrative and conjugative elements and their hosts: composition, distribution and organization. *Nucleic Acids Research*, 45(15), 8943–8956. <https://doi.org/10.1093/nar/gkx607> (CONJscan model set)
15. Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. <https://doi.org/10.1093/gigascience/giab008>
16. De Coster, W., D'Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666–2669. <https://doi.org/10.1093/bioinformatics/bty149> (NanoPlot)
17. Doster, E., et al. (2020). MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. *Nucleic Acids Research*, 48(D1), D561–D569. <https://doi.org/10.1093/nar/gkz1010>
18. Eddy, S. R. (2011). Accelerated profile HMM searches. *PLoS Computational Biology*, 7(10), e1002195. <https://doi.org/10.1371/journal.pcbi.1002195> (the profile-HMM engine inside ISEScan, MacSyFinder, VirSorter2 and run_dbCAN)
19. EFSA Panel on Biological Hazards (BIOHAZ). (2023). Statement on how to interpret the QPS qualification on 'acquired antimicrobial resistance genes'. *EFSA Journal*, 21(10), 8323. <https://doi.org/10.2903/j.efsa.2023.8323> (the method for telling intrinsic from acquired)
20. EFSA Scientific Committee. (2025). Guidance on the characterisation of microorganisms in support of the risk assessment of products used in the food chain. *EFSA Journal*, 23(11), e9705. <https://doi.org/10.2903/j.efsa.2025.9705> (the guidance in force; its Appendix A lists what it supersedes in the 2018 FEEDAP guidance)
21. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. <https://doi.org/10.1093/bioinformatics/btw354>
22. Feldgarden, M., et al. (2019). Validating the AMRFinder tool and resistance gene database by using antimicrobial resistance genotype–phenotype correlations in a collection of isolates. *Antimicrobial Agents and Chemotherapy*, 63(11), e00483-19. <https://doi.org/10.1128/AAC.00483-19>
23. Feldgarden, M., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728. <https://doi.org/10.1038/s41598-021-91456-0>
24. Grüning, B., et al. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nature Methods*, 15, 475–476. <https://doi.org/10.1038/s41592-018-0046-7>
25. Guo, J., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome*, 9(1), 37. <https://doi.org/10.1186/s40168-020-00990-y>
26. Gupta, S. K., et al. (2014). ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. *Antimicrobial Agents and Chemotherapy*, 58(1), 212–220. <https://doi.org/10.1128/AAC.01310-13>
27. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072–1075. <https://doi.org/10.1093/bioinformatics/btt086>
28. Huerta-Cepas, J., et al. (2019). eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. *Nucleic Acids Research*, 47(D1), D309–D314. <https://doi.org/10.1093/nar/gky1085>
29. Ingle, D. J., et al. (2016). In silico serotyping of *E. coli* from short read data identifies limited novel O-loci but extensive diversity of O:H serotype combinations within and between pathogenic lineages. *Microbial Genomics*, 2(7), e000064. <https://doi.org/10.1099/mgen.0.000064> (the `ecoh` database)
30. Jia, B., et al. (2017). CARD 2017: expansion and model-centric curation of the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 45(D1), D566–D573. <https://doi.org/10.1093/nar/gkw1004>
31. Kirsch, J. M., Hryckowian, A. J., & Duerkop, B. A. (2024). A metagenomics pipeline reveals insertion sequence-driven evolution of the microbiota. *Cell Host & Microbe*, 32(5), 739–754.e4. <https://doi.org/10.1016/j.chom.2024.03.005> (pseudoR, the source of ISOSDB)
32. Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. <https://doi.org/10.1038/s41587-019-0072-8> (Flye)
33. Laetsch, D. R., & Blaxter, M. L. (2017). BlobTools: interrogation of genome assemblies. *F1000Research*, 6, 1287. <https://doi.org/10.12688/f1000research.12232.1>
34. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359. <https://doi.org/10.1038/nmeth.1923>
35. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. <https://doi.org/10.1093/bioinformatics/bty191>
36. Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*, 25(14), 1754–1760. <https://doi.org/10.1093/bioinformatics/btp324> (BWA, the aligner Polypolish reads)
37. Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. <https://doi.org/10.12688/f1000research.29032.2>
38. Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology*, 39(5), 578–585. <https://doi.org/10.1038/s41587-020-00774-7>
39. Néron, B., Denise, R., Coluzzi, C., Touchon, M., Rocha, E. P. C., & Abby, S. S. (2023). MacSyFinder v2: improved modelling and search engine to identify molecular systems in genomes. *Peer Community Journal*, 3, e28. <https://doi.org/10.24072/pcjournal.250>
40. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. *Bioinformatics*, 32(2), 292–294. <https://doi.org/10.1093/bioinformatics/btv566>
41. Oxford Nanopore Technologies. *Medaka.* <https://github.com/nanoporetech/medaka>
42. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research*, 25(7), 1043–1055. <https://doi.org/10.1101/gr.186072.114>
43. Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., & Hugenholtz, P. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*, 50(D1), D785–D794. <https://doi.org/10.1093/nar/gkab776> (cite the release you ran — BacFlux 2.0.0 expects **R232**)
44. Ross, K., et al. (2021). TnCentral: a prokaryotic transposable element database and web portal for transposon analysis. *mBio*, 12(5), e02060-21. <https://doi.org/10.1128/mBio.02060-21>
45. Schwengers, O., Barth, P., Falgenhauer, L., Hain, T., Chakraborty, T., & Goesmann, A. (2020). Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. *Microbial Genomics*, 6(10), mgen000398. <https://doi.org/10.1099/mgen.0.000398>
46. Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*, 7(11), 000685. <https://doi.org/10.1099/mgen.0.000685>
47. Seemann, T. (2020). *ABRicate: mass screening of contigs for antimicrobial and virulence genes.* <https://github.com/tseemann/abricate>
48. Seemann, T. *Snippy: rapid haploid variant calling and core genome alignment.* <https://github.com/tseemann/snippy>
49. Wang, M., et al. (2024). ICEberg 3.0: functional categorization and analysis of the integrative and conjugative elements in bacteria. *Nucleic Acids Research*, 52(D1), D732–D737. <https://doi.org/10.1093/nar/gkad935> (also the paper to cite for ICEfinder2 / ICEscan)
50. Wick, R. R. *Filtlong: quality filtering tool for long reads.* <https://github.com/rrwick/Filtlong>
51. Wick, R. R., & Holt, K. E. (2022). Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLoS Computational Biology*, 18(1), e1009802. <https://doi.org/10.1371/journal.pcbi.1009802>
52. Xie, Z., & Tang, H. (2017). ISEScan: automated identification of insertion sequence elements in prokaryotic genomes. *Bioinformatics*, 33(21), 3340–3347. <https://doi.org/10.1093/bioinformatics/btx433>
53. Zankari, E., et al. (2012). Identification of acquired antimicrobial resistance genes. *Journal of Antimicrobial Chemotherapy*, 67(11), 2640–2644. <https://doi.org/10.1093/jac/dks261> (ResFinder)
54. Zheng, J., Ge, Q., Yan, Y., Zhang, X., Huang, L., & Yin, Y. (2023). dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research*, 51(W1), W115–W121. <https://doi.org/10.1093/nar/gkad328>

## Design influence

Work that shaped the mobilome module without contributing code to it — the
distinction matters, and it is spelled out on [Licensing](licensing.md).

- **EBI Mobilome Annotation Pipeline** (Apache-2.0).
  <https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline> — from which
  BacFlux adopted decisions rather than code: that ICEscan exists and is worth
  running, the Sequence Ontology terms for the mobilome GFF, the
  `contig_id|mge_type-start:end` element ID format, the discard-with-reason file,
  the 500 bp minimum element size and the 0.9 CDS-coverage rule. The two
  pipelines answer different questions: theirs goes deep on the mobilome of
  metagenomes and MAGs, BacFlux covers reads through to AMR mobility for single
  isolates.
- **ICEfinder2**, whose variable-length repeat search — rather than a
  fixed-length probe — is why BacFlux's att-site search works the way it does.
  Established by reading it, then implemented independently.

- Puterová, J. & Martínek, T. (2021) digIS: towards detecting distant and putative
  novel insertion sequence elements in prokaryotic genomes. *BMC Bioinformatics*
  22:258. <https://doi.org/10.1186/s12859-021-04177-6> — benchmarked here for the
  IS false-discovery figures quoted on the mobilome pages, not run by BacFlux.

## Licence

BacFlux is released under the
[MIT Licence](https://github.com/iLivius/BacFlux/blob/main/LICENSE). The tools,
model sets and databases it invokes are distributed under their own terms, which
are listed on [Licensing](licensing.md).
