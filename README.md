# BacFlux
A workflow for bacterial short reads assembly, QC, annotation, and more.

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.4.7-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11143917.svg)](https://doi.org/10.5281/zenodo.11143917)

```bash
________             _______________              
___  __ )_____ _________  ____/__  /___  _____  __
__  __  |  __ `/  ___/_  /_   __  /_  / / /_  |/_/
_  /_/ // /_/ // /__ _  __/   _  / / /_/ /__>  <  
/_____/ \__,_/ \___/ /_/      /_/  \__,_/ /_/|_|  
                                                  
BacFlux v1.1.7

May 2024
```

![BacFlux DAG](miscellaneous/BacFlux_v1.1.x_DAG.png)

## Authors and Contributors
[AIT Austrian Institute of Technology, Center for Health & Bioresources](https://www.ait.ac.at/en/research-topics/bioresources)

- Livio Antonielli
- Dominik K. Großkinsky
- Hanna Koch
- Friederike Trognitz

[IPK Leibniz Institute of Plant Genetics and Crop Plant Research, Cryo and Stress Biology](https://www.ipk-gatersleben.de/forschung/genbank/cryo-und-stressbiologie)
- Manuela Nagel
- Alexa Sanchez Mejia

## Synopsis
`BacFlux` is a comprehensive and automated bioinformatics workflow designed specifically for the processing and analysis of bacterial genomic data sequenced with Illumina technology. It integrates several powerful tools, each performing a specific task, into a seamless workflow managed by a Snakemake script.
The pipeline accepts paired-end reads as input and subjects them to a series of analyses including steps for quality control, assessment of genome completeness and contamination, taxonomic placement, annotation, inference of secondary metabolites, screening for antimicrobial resistance and virulence genes, and investigation of plasmid presence.

## Table of Contents
- [Quick Start](#quick-start)
- [Rationale](#rationale)
- [Description](#description)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running BacFlux](#running-bacflux)
- [Output](#output)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [References](#references)

## Quick Start
This guide gets you started with `BacFlux`. Here's a quick guide:

- Download the latest release:
  ```bash
  #git command
  git clone https://github.com/iLivius/BacFlux.git
  ```

- Configure the `config.yaml`: Specify the input directory containing the raw sequencing data (i.e. paired-end FASTQ files: *strain-1_R1.fq.gz*, *strain-1_R2.fq.gz*) and the desired location for the analysis outputs, respectively. `BacFlux` relies on external databases for some analyses. Some of them are not automatically installed and the `config.yaml` must be edited with the path to the following downloaded databases:

    * blast_db: path to the [NCBI nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database directory
    * eggnog_db: path to the [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Installation) diamond database directory
    * gtdbtk_db: path to the [GTDB](https://ecogenomics.github.io/GTDBTk/installing/index.html) database directory
    * bakta_db: path to the [Bakta](https://github.com/oschwengers/bakta?tab=readme-ov-file#database) database directory
    * platon_db: path to the [Platon](https://github.com/oschwengers/platon?tab=readme-ov-file#database) database directory

- Install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (if not installed already) and activate the environment:
  ```bash
  #optional, if not installed already
  conda create -c conda-forge -c bioconda -n snakemake snakemake
  #activate Snakemake environment
  conda activate snakemake
  ```
Now you're all set to run `BacFlux`! Refer to the [installation](#installation), [configuration](#configuration) and [running BacFlux](#running-bacflux) sections for detailed instructions.

- Here is an example. Within the main workflow directory, launch the pipeline as follows:
   ```bash
  snakemake --sdm conda --keep-going --ignore-incomplete --keep-incomplete --cores 50
  ```

  This  command uses the following options:

    --sdm: uses conda for dependency management

    --keep-going: continues execution even if errors occur in some steps

    --ignore-incomplete: ignores rules with missing outputs

    --keep-incomplete: keeps incomplete intermediate files

    --cores 50: cap the amount of local CPUs at this value (adjust as needed).

## Rationale
The analysis of bacterial WGS data often involves a complex series of steps using various bioinformatic tools. Manual execution of this process can be time-consuming, error-prone, and difficult to reproduce. `BacFlux` addresses these challenges by providing a comprehensive and automated Snakemake workflow that streamlines bacterial genomic data analysis.

`BacFlux` integrates several best-in-class bioinformatic tools into a cohesive pipeline, automating tasks from quality control and assembly to annotation, taxonomic classification, identification of resistance genes and viral sequences.

By providing a user-friendly and automated solution, `BacFlux` empowers researchers to focus on interpreting the biological meaning of their data.

## Description
Here's a breakdown of the `BacFlux` workflow:

01. **Preprocessing:**

    * Checks raw reads for Illumina phiX contamination using [bowtie2](https://github.com/BenLangmead/bowtie2).

    * Filters and removes adapters from reads using [fastp](https://github.com/OpenGene/fastp).

02. **Assembly:**
    * Assembles filtered reads into contigs with [SPAdes](https://github.com/ablab/spades).

03. **Quality Control, Contamination and Completeness Assessment:**
    * Filters contigs based on minimum length (at least 500 bp) and coverage (2x).
    * Maps filtered reads back to contigs, using [bowtie2](https://github.com/BenLangmead/bowtie2) and [samtools](https://github.com/samtools/samtools), and analyzes the resulting BAM file with [QualiMap](http://qualimap.conesalab.org/).
    * Performs local alignments of contigs against the [NCBI nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database using [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
    * Checks for contaminant contigs with [BlobTools](https://github.com/DRL/blobtools). Unless otherwise specified (see [configuration](#configuration) section for more details), the output of this step will be parsed automatically to discard contaminants based on the relative taxonomic composition of the contigs.   
    * Evaluates genome assembly quality with [Quast](https://github.com/ablab/quast).
    * Assesses genome completeness and contamination with [CheckM](https://github.com/Ecogenomics/CheckM) using taxon-specific markers.

04. **Taxonomic Analysis:**
    * Performs accurate taxonomic placement with [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) using a curated reference [database](https://gtdb.ecogenomic.org/).

05. **Annotation:**
    * Annotates contigs using [Prokka](https://github.com/tseemann/prokka) and [Bakta](https://github.com/oschwengers/bakta) for functional prediction.
    * Provides further functional annotation with [EggNOG](https://github.com/eggnogdb).
    * Infers secondary metabolites with [antiSMASH](https://github.com/antismash/antismash).

06. **Antimicrobial Resistance (AMR):**
    * Maps filtered reads to the [CARD](https://card.mcmaster.ca/) database with [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/). 
    * Screens contigs for antimicrobial resistance and virulence genes using [ABRicate](https://github.com/tseemann/abricate).

07. **Plasmids:**
    * Investigates the presence of plasmids with [Platon](https://github.com/PlatONnetwork) and confirms results with [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) search.

08. **Prophages:**
    * Screens contigs for viral sequences with [VirSorter2](https://github.com/jiarong/VirSorter2), followed by [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) for refinement.

09. **Reporting:**
    * Parses and aggregates results to generate a report using [MultiQC](https://github.com/MultiQC/MultiQC).

## Installation
BacFlux downloads automatically all dependencies and several databases.  However, some external databases require manual download before running the workflow.

1. **Download BacFlux:**

    Head over to the Releases section of the repository.
    Download the latest archive file (typically in .zip or .tar.gz format). This archive contains the `BacFlux` Snakefile script, the `config.yaml` configuration file, and an `envs` environment directory.
    Extract the downloaded archive into your desired location. This will create a directory structure with the necessary files and directories.
    Alternatively, download via command line as:
    ```bash
    #git command
    git clone https://github.com/iLivius/BacFlux.git
    ```

2. **Install Snakemake:**

    `BacFlux` relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to manage the workflow execution. Find the official and complete set of instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Alternatively, to install Snakemake as a Conda environment:
    ```bash
    #install Snakemake in a new Conda environment (alternatively, use mamba)
    conda create -c conda-forge -c bioconda -n snakemake snakemake
    ```

3. **Databases:**

    While `BacFlux` automates the installation of all software dependencies, some external databases need to be downloaded manually. Unless you have installed them already. In that case, skip this paragraph and jump to the [configuration](#configuration) section. 

    Here are the required databases and instructions for obtaining them.

    * `NCBI nt` database, adapted from [here](https://gist.github.com/ppflrs/336e49f8ae3843dc06cc3925940f3024):
        ```bash
        #create a list of all nt links in the directory designated to host the database (recommended)
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_links.list
        
        #alternatively, create a list of nt links for bacteria only 
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt_prok.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_prok_links.list
       
        #download in parallel, without overdoing it
        cat nt*.list | parallel -j4 'rsync -h --progress rsync://{} .'

        #decompress with multiple CPUs
        find . -name '*.gz' | parallel -j4 'echo {}; tar -zxf {}'

        #get NCBI taxdump
        wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        tar -zxvf taxdump.tar.gz

        #get NCBI BLAST taxonomy
        wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'
        tar -zxvf taxdb.tar.gz

        #get NCBI accession2taxid file
        wget -c 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        gunzip nucl_gb.accession2taxid.gz
        ```
        *NOTE: the complete NCBI nt database and taxonomy-related files should take around 400 GB of hard drive space.*

    * `eggNOG diamond` database:
        ```bash
       #the easiest way is to install a Conda environment with eggnog-mapper, first
       conda create -n eggnog-mapper eggnog-mapper=2.1.12

       #activate the environment
       conda activate eggnog-mapper

       #then, create a directory where you want to install the diamond database for eggnog-mapper 
       mkdir /data/eggnog_db
       #change /data/eggnog_db with your actual PATH

       #finally, download the diamond db in the newly created directory 
       download_eggnog_data.py --data_dir /data/eggnog_db -y
        ```
        *NOTE: the eggNOG database requires ~50 GB of space.*

    * `GTDB` database:
        ```bash
        #move first inside the directory where you want to place the database, then download and decompress the file
        wget -c https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
        tar -xzvf gtdbtk_r214_data.tar.gz
        rm gtdbtk_r214_data.tar.gz
        ```
        *NOTE: the data needed to run GTDB-Tk will occupy around 85 GB of space.*

    * `Bakta` database:
        ```bash
        #Bakta database comes in two flavours. To download the full database, use the following link (recommended):
        wget -c https://zenodo.org/records/10522951/files/db.tar.gz
        tar -xzf db.tar.gz
        rm db.tar.gz

        #alternatively, download a lighter version
        wget https://zenodo.org/record/10522951/files/db-light.tar.gz
        tar -xzf db-light.tar.gz
        rm db-light.tar.gz
        ```
        *NOTE: according to the [source](https://github.com/oschwengers/bakta?tab=readme-ov-file#database) the light version should take 1.4 GB compressed and 3.4 GB decompressed, whereas the full database should get 37 GB zipped and 71 GB unzipped.*

    * `Platon` database:
        ```bash
        #download the database in a directory of your choice 
        wget https://zenodo.org/record/4066768/files/db.tar.gz
        tar -xzf db.tar.gz
        rm db.tar.gz
        ```
        *NOTE: according to the [source](https://github.com/oschwengers/platon?tab=readme-ov-file#database), the zipped version occupies 1.6 GB and 2.8 GB when unzipped.*

## Configuration
Before running `BacFlux`, you must edit the `config.yaml` file with a text editor. The file is organized in different sections: `links`, `directories`, `resources` and `parameters`, respectively.  
 
- `links`

    This section should work fine as it is, therefore it is recommandable to change the `links` only if not working or to update the database versions:

    - [phix_link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz): Path to the PhiX genome reference used by Illumina for sequencing control.
    - [card_link](https://card.mcmaster.ca/download/0/broadstreet-v3.2.9.tar.bz2): Path to the Comprehensive Antibiotic Resistance Database (`CARD`)
    - [checkv_link](https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz): Path to the `CheckV` database for viral genome quality assessment

- `directories`

    Update paths based on your file system: 

    - **fastq_dir**: directory containing the paired-end reads of your sequenced strains, in FASTQ format. You can provide as many as you like but at the following conditions:
        1. Files can only have the following extensions: either `fastq`, `fq`, `fastq.gz` or `fq.gz`.
        2. You can provide multiple samples but the extension should be the same for all files. So, don't mix files with different extensions.
        3. Sample names should be formatted as follows: mystrain_R1.fq and mystrain_R2.fq, strain-1_R1.fq and strain-1_R2.fq, strain2_R1.fq, strain2_R2.fq. In this example, `BacFlux` will interpret the name of each strain as: "mystrain", "strain-1", and "strain2", respectively. Strain names cannot contain underscores. See another example, below:
            ```bash
            #the input dir contains the PE reads of two strains, PE212-1 and PE253-B, respectively
            ahab@pequod:~/data$ ls -lh
            total 1,6G
            -rw-rw-r-- 1 ahab ahab 379M Apr  8 16:48 PE212-1_R1.fastq.gz
            -rw-rw-r-- 1 ahab ahab 385M Apr  8 16:48 PE212-1_R2.fastq.gz
            -rw-rw-r-- 1 ahab ahab 421M Apr  8 16:48 PE253-B_R1.fastq.gz
            -rw-rw-r-- 1 ahab ahab 433M Apr  8 16:48 PE253-B_R2.fastq.gz
            ```
             
    - **out_dir**: This directory will store all output files generated by `BacFlux`. Additionally, by default, `BacFlux` will install required software and databases here, within Conda environments. Reusing this output directory for subsequent runs avoids reinstalling everything from scratch.

    - **blast_db**: path to the whole `NCBI nt` (recommended) or prokaryotic database only, and related taxonomic dependencies, see [installation](#installation).

    - **eggnog_db**: path to the diamond database for `eggNOG`.

    - **gtdbtk_db**: path to the R214 release of `GTDB`.

    - **bakta_db**: path to either the light or full (recommended) database of `Bakta`.

    - **platon_db**: path to the `Platon` database.

- `resources`
    
    In this section you can specify the hardware resources available to the workflow: 

  - threads: max number of CPUs used by each rule
  - ram_gb: max amount of RAM used (SPAdes only).

- `parameters`

    `BacFlux` offers an optional parameter that allows you to specify the bacterial `genus` of contigs you want to keep in the final assembly. Leaving this field blank will work in most cases: `BacFlux` will automatically retain contigs matching the most abundant taxon based on relative composition, determined by BLAST analysis. This approach is often effective but has limitations. For instance, it might not provide the best resolution at species level due to the reliance on the sum of best scores of BLAST hits. Additionally, this method can be misleading if the contaminating organism belongs to the same genus as your target organism or if you have co-cultured closely related species or strains. If this parameter causes more troubles than benefits in your specific case, just remove the `genus` option from the `config.yaml` file. 

    - **Using** the `genus` parameter: 
  
        If a contaminant is ascertained to be more abundant than your target organism, you can re-run the workflow after reviewing the assembly [output](#output). Specify the `genus` of the desired bacterial taxon you want to keep in during the re-run. 
    
    - **Disabling** the `genus` filtering:

        If either the automatic inference of contaminant contigs or the manual selection of the desired taxon are still not working for you, simply delete the `genus` option from the `parameters`. In this case, only contigs tagged as "no-hit" after BLAST search will be filtered out.

## Running BacFlux
`BacFlux` can be executed as simply as a Snakefile. Please refer to the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for more details.
```bash
# first, activate the Snakemake Conda environment
conda activate snakemake

# navigate inside the directory where the BacFlux archive was downloaded and decompressed

# launch the workflow
snakemake --sdm conda --cores 50
```
*NOTE:  Starting from Snakemake version 8.4.7, the --use-conda option has been deprecated. Instead, you should now use --software-deployment-method conda or --sdm conda.*

## Output
The workflow output reflects the steps described in the [description](#description) section. Here's a breakdown of the subdirectories created within the main output folder, along with explanations of their contents:

- `01.pre-processing`: QC and statistics of raw reads and trimmed reads, produced by [fastp](https://github.com/OpenGene/fastp) (v0.23.4).

- `02.assembly`: Content output by [SPAdes](https://github.com/ablab/spades) (v3.15.5). In addiction to the raw contigs, you will find also the filtered contigs (>500bp and at least 2x) and the selected contigs, which are the contigs selected after BLAST search and decontamination (see `parameters` in the [configuration](#configuration) section above). The follow-up applications used during the worflow will either use selected contigs (i.e. for annotation purposes) or raw, filtered and selected contigs (i.e. to evaluate the genome completenness and contamination).

- `03.post-processing`: Contains the following sub-directories:
    - **mapping_evaluation**: [QualiMap](http://qualimap.conesalab.org/) (v2.3) output based on filtered contigs.
    - **contaminants**: Contig selection based on [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) (v2.15.0) search and [BlobTools](https://github.com/DRL/blobtools) (1.1.1) analysis. Check the `composition` text file for a quick overview of the relative composition of your assembly.
    - **assembly_evaluation**: [Quast](https://github.com/ablab/quast) (v5.2.0) output based on selected contigs.
    - **completenness_evaluation**: [CheckM](https://github.com/Ecogenomics/CheckM) (1.2.2) output based on raw, filtered and selected contigs.

- `04.taxonomy`: Taxonomic placement of raw, filtered and selected contigs, performed by [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) (v2.3.2).

- `05.annotation`: Contains the following sub-directories:
    - **prokka**: Legacy annotation performed by [Prokka](https://github.com/tseemann/prokka) (v1.14.5).
    - **bakta**: Accurate annotation outputted by [Bakta](https://github.com/oschwengers/bakta) (v1.9.3).
    - **eggnog**: Functional annotation produced by [EggNOG](https://github.com/eggnogdb) mapper (v2.1.12).
    - **antismash**: Secondary metabolites inferred by [antiSMASH](https://github.com/antismash/antismash) (v7.1.0).

- `06.AMR`: Antimicrobial resistance features are investigated with two complementary approaches:
    - **AMR_mapping**: reads filtered by [fastp](https://github.com/OpenGene/fastp) (v0.23.4) are mapped to the CARD database (v3.2.9.) using [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (v39.06) with minimum identiy = 0.99. Mapping results are parsed and features with a covered length of at least 70% are reported in the `AMR legend` file.
    - **ABRicate**: selected contigs are screened for the presence of AMR elements and virulence factors, using [ABRicate](https://github.com/tseemann/abricate) (v1.0.1).

- `07.plasmids`: Selected contigs are screened for the presence of plasmid replicons with [Platon](https://github.com/oschwengers/platon) and results verified by BLAST search to avoid false positive. Contigs ascertained as plasmids are reported in the `verified plasmids` file.

- `08.phages`: Filtered contigs are screened for the presence of viral sequences using [VirSorter2](https://github.com/jiarong/VirSorter2) (v2.2.4), followed by [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) (v1.0.1) for refinement:
    - **virsorter**: Following the instructions provided [here](https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=3), viral groups (i.e. dsDNA phage, NCLDV, RNA, ssDNA, and lavidaviridae) are detected with a loose cutoff of 0.5 for maximal sensitivity. Original sequences of circular and (near) fully viral contigs are preserved and passed to the next tool.
    - **checkv**: This second step serves to quality control the results of the previous step to avoid the presence of non-viral sequences (false positive) and to trim potential host regions left at the ends of proviruses.

- `09.report`: [MultiQC](https://github.com/MultiQC/MultiQC) (v1.17) is used to parse and aggregate the results of the following tools:
    1. [fastp](https://github.com/OpenGene/fastp) (v0.23.4)
    2. [QualiMap](http://qualimap.conesalab.org/) (v2.3)
    3. [Quast](https://github.com/ablab/quast) (v5.2.0)
    4. [Prokka](https://github.com/tseemann/prokka) (v1.14.5)
    5. [Bakta](https://github.com/oschwengers/bakta) (v1.9.3)

## Acknowledgements
This work was supported by the [Austrian Science Fund (FWF)](https://www.fwf.ac.at/en/) [Project I6030-B].

## Citation
Antonielli, L., Großkinsky, D. K., Koch, H., Trognitz, F., Sanchez Mejia, A., & Nagel, M. (2024). BacFlux: A workflow for bacterial short reads assembly, QC, annotation, and more. Zenodo. https://doi.org/10.5281/zenodo.11143917

## References
1. Alcock, B. P., Huynh, W., Chalil, R., Smith, K. W., Raphenya, A. R., Wlodarski, M. A., Edalatmand, A., Petkau, A., Syed, S. A., Tsang, K. K., Baker, S. J. C., Dave, M., McCarthy, M. C., Mukiri, K. M., Nasir, J. A., Golbon, B., Imtiaz, H., Jiang, X., Kaur, K., … McArthur, A. G. (2023). CARD 2023: Expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51(D1), D690–D699. https://doi.org/10.1093/nar/gkac920

2. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

3. Blin, K., Shaw, S., Augustijn, H. E., Reitz, Z. L., Biermann, F., Alanjary, M., Fetter, A., Terlouw, B. R., Metcalf, W. W., Helfrich, E. J. N., van Wezel, G. P., Medema, M. H., & Weber, T. (2023). antiSMASH 7.0: New and improved predictions for detection, regulation, chemical structures and visualisation. Nucleic Acids Research, 51(W1), W46–W50. https://doi.org/10.1093/nar/gkad344

4. Bushnell, B. (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. https://escholarship.org/uc/item/1h3515gn

5. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

6. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825–5829. https://doi.org/10.1093/molbev/msab293

7. Challis, R., Richards, E., Rajan, J., Cochrane, G., & Blaxter, M. (2020). BlobToolKit – Interactive Quality Assessment of Genome Assemblies. G3 Genes|Genomes|Genetics, 10(4), 1361–1374. https://doi.org/10.1534/g3.119.400908

8. Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2022). GTDB-Tk v2: Memory friendly classification with the genome taxonomy database. Bioinformatics, 38(23), 5315–5316. https://doi.org/10.1093/bioinformatics/btac672

9. Chen, L., Zheng, D., Liu, B., Yang, J., & Jin, Q. (2016). VFDB 2016: Hierarchical and refined dataset for big data analysis--10 years on. Nucleic Acids Research, 44(D1), D694-697. https://doi.org/10.1093/nar/gkv1239

10. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

11. Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

12. Doster, E., Lakin, S. M., Dean, C. J., Wolfe, C., Young, J. G., Boucher, C., Belk, K. E., Noyes, N. R., & Morley, P. S. (2020). MEGARes 2.0: A database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Research, 48(D1), D561–D569. https://doi.org/10.1093/nar/gkz1010

13. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England), 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

14. Feldgarden, M., Brover, V., Haft, D. H., Prasad, A. B., Slotta, D. J., Tolstoy, I., Tyson, G. H., Zhao, S., Hsu, C.-H., McDermott, P. F., Tadesse, D. A., Morales, C., Simmons, M., Tillman, G., Wasilenko, J., Folster, J. P., & Klimke, W. (2019). Validating the AMRFinder Tool and Resistance Gene Database by Using Antimicrobial Resistance Genotype-Phenotype Correlations in a Collection of Isolates. Antimicrobial Agents and Chemotherapy, 63(11), e00483-19. https://doi.org/10.1128/AAC.00483-19

15. Guo, J., Bolduc, B., Zayed, A. A., Varsani, A., Dominguez-Huerta, G., Delmont, T. O., Pratama, A. A., Gazitúa, M. C., Vik, D., Sullivan, M. B., & Roux, S. (2021). VirSorter2: A multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome, 9(1), 37. https://doi.org/10.1186/s40168-020-00990-y

16. Gupta, S. K., Padmanabhan, B. R., Diene, S. M., Lopez-Rojas, R., Kempf, M., Landraud, L., & Rolain, J.-M. (2014). ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes. Antimicrobial Agents and Chemotherapy, 58(1), 212–220. https://doi.org/10.1128/AAC.01310-13

17. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics (Oxford, England), 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

18. Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085

19. Ingle, D. J., Valcanis, M., Kuzevski, A., Tauschek, M., Inouye, M., Stinear, T., Levine, M. M., Robins-Browne, R. M., & Holt, K. E. (2016). In silico serotyping of E. coli from short read data identifies limited novel O-loci but extensive diversity of O:H serotype combinations within and between pathogenic lineages. Microbial Genomics, 2(7), e000064. https://doi.org/10.1099/mgen.0.000064

20. Jia, B., Raphenya, A. R., Alcock, B., Waglechner, N., Guo, P., Tsang, K. K., Lago, B. A., Dave, B. M., Pereira, S., Sharma, A. N., Doshi, S., Courtot, M., Lo, R., Williams, L. E., Frye, J. G., Elsayegh, T., Sardar, D., Westman, E. L., Pawlowski, A. C., … McArthur, A. G. (2017). CARD 2017: Expansion and model-centric curation of the comprehensive antibiotic resistance database. Nucleic Acids Research, 45(D1), D566–D573. https://doi.org/10.1093/nar/gkw1004

21. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), Article 4. https://doi.org/10.1038/nmeth.1923

22. Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., & Köster, J. (2021). Sustainable data analysis with Snakemake (10:33). F1000Research. https://doi.org/10.12688/f1000research.29032.2

23. Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology, 39(5), Article 5. https://doi.org/10.1038/s41587-020-00774-7

24. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292–294. https://doi.org/10.1093/bioinformatics/btv566

25. Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A., & Hugenholtz, P. (2022). GTDB: An ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. Nucleic Acids Research, 50(D1), D785–D794. https://doi.org/10.1093/nar/gkab776

26. Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043–1055. https://doi.org/10.1101/gr.186072.114

27. Schwengers, O., Barth, P., Falgenhauer, L., Hain, T., Chakraborty, T., & Goesmann, A. (2020). Platon: Identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. Microbial Genomics, 6(10), mgen000398. https://doi.org/10.1099/mgen.0.000398

28. Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: Rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11), 000685. https://doi.org/10.1099/mgen.0.000685

29. Seemann, T. (2014). Prokka: Rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068–2069. https://doi.org/10.1093/bioinformatics/btu153

30. Seemann, T. (2020). ABRicate. https://github.com/tseemann/abricate

31. Zankari, E., Hasman, H., Cosentino, S., Vestergaard, M., Rasmussen, S., Lund, O., Aarestrup, F. M., & Larsen, M. V. (2012). Identification of acquired antimicrobial resistance genes. The Journal of Antimicrobial Chemotherapy, 67(11), 2640–2644. https://doi.org/10.1093/jac/dks261
