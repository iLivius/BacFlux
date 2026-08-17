# BacFlux v2.0.0 — illumina front end, assembly.
#
# Turn the cleaned read pairs into contigs, then drop the contigs that cannot be
# real isolate sequence. SPAdes --isolate is tuned for a single high-coverage
# bacterial culture, and its multi-k run resolves both
# low-complexity and repeat-rich regions in one pass.
#
# Stage chain: TRIM_R1/TRIM_R2 → illumina_assembly → SPADES_CONTIGS →
# filter_contigs → DRAFT_CONTIGS → the BLAST/BlobTools screen in
# shared/10_decontam.smk → DECONTAM_CONTIGS.
#
# illumina_assembly  : one SPAdes run per sample. Declares its own working
#                      directory as its only directory() output, never the
#                      sample directory.
# filter_contigs     : two awk passes that keep contigs of at least 500 bp and at
#                      least 2x coverage, reading both numbers off SPAdes' own
#                      FASTA headers.
#
# Where this mode ends. In illumina mode, decontamination IS the last assembly
# step: 00_common sets DECONTAM_CONTIGS = FINAL_CONTIGS (the same Python object),
# so the shared `select_contigs` rule writes the canonical delivered genome and
# this front end must NOT declare FINAL_CONTIGS anywhere. Declaring it here would
# give Snakemake two rules producing one file.
#
# Everything referenced here comes from 00_common.smk: TRIM_R1/TRIM_R2,
# SPADES_DIR, SPADES_CONTIGS, DRAFT_CONTIGS, FASTA_LIN_CMD, FASTA_SEL_CMD, LOGS,
# CPUS, RAM.


# ──────────────────── Short-read assembly (SPAdes) ─────────────
# Takes in: r1, r2 = TRIM_R1 / TRIM_R2, the fastp-trimmed, PhiX-free pairs written
#           by rule trim_adapters (illumina/10_reads.smk).
# Does:     one SPAdes de novo assembly per sample in --isolate mode, over the six
#           the k-mer ladder set by parameters.spades_kmers (default auto, which
#           lets SPAdes size it from the reads). OMP_NUM_THREADS is exported so the
#           OpenMP parts of SPAdes obey the same thread budget as the -t flag;
#           without it they default to every core on the machine, which
#           oversubscribes a multi-sample run. -m is a hard RAM ceiling in GB.
# Produces:
#   dir     = SPADES_DIR, 02.assembly/{sample}/spades/ — SPAdes' own working and
#             output tree. Declared so Snakemake owns it; no other rule reads it.
#   contigs = SPADES_CONTIGS, SPAdes' own spades/contigs.fasta inside that tree.
# Consumed by: filter_contigs, below, which reads SPADES_CONTIGS.
#
# Both the thread count and the RAM ceiling are left uncapped, as in v1, because
# SPAdes is the single most expensive step in this mode and scales with both.
#
# LAYOUT LANDMINE — do not "simplify" the output back to the sample directory.
# v1 declared directory("02.assembly/{sample}"). Under the v2 layout (D1) that
# directory ALSO holds contaminants/, eval/ and (in the long-read modes)
# fix_start/, medaka/ and snps/, written by OTHER rules. Snakemake deletes and
# recreates a directory() output when the rule re-runs, so declaring the sample
# directory would wipe every other stage's results on a re-run. The assembler
# therefore declares only its own sub-directory, SPADES_DIR.
#
# (v1 rule name: `genome_assembly` in BacFlux, `illumina_assembly` in BacFluxL+.
#  D5 unifies on `illumina_assembly`, so the SPAdes rule has ONE name in both the
#  illumina and hybrid front ends. The log filename changes with it.
#  v1 message: "--- SPAdes: genome assembly. ---")
rule illumina_assembly:
    input:
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        dir = directory(SPADES_DIR),
        contigs = SPADES_CONTIGS,
    params:
        # Empty when parameters.spades_kmers is auto, so SPAdes chooses the ladder from
        # the read length it measures. See the note on SPADES_KMER_FLAG in 00_common.smk.
        kmers = SPADES_KMER_FLAG,
    conda:
        "../../envs/spades.yaml"
    threads: CPUS
    resources:
        ram = RAM
    log:
        LOGS + "/illumina_assembly_{sample}.log"
    shell:
        """
        OMP_NUM_THREADS={threads} \
        spades.py {params.kmers} --isolate \
          --pe1-1 {input.r1} \
          --pe1-2 {input.r2} \
          -o {output.dir} \
          -t {threads} \
          -m {resources.ram} > {log} 2>&1
        """


# ──────────────────── Contig length and coverage filter ────────
# A pure-isolate assembly has a long tail of tiny, low-coverage contigs —
# sequencing noise, chimeras, and fragments of whatever else was in the tube. They
# add nothing to the annotation but do inflate the contig count, the CheckM
# contamination estimate and the BLAST screening time. This rule drops them.
#
# Takes in: contigs = SPADES_CONTIGS, the raw output of rule illumina_assembly
#           above.
# Does:     two awk passes, both defined once in 00_common.smk, applying the v1
#           cutoffs — keep a contig only if it is at least 500 bp AND has at least
#           2x coverage, both read straight off SPAdes' own FASTA headers:
#             FASTA_LIN_CMD — put each record's sequence on ONE line, so the
#                             selector below can read a header and its sequence as
#                             a header/getline pair.
#             FASTA_SEL_CMD — keep a record when SPAdes' own header says
#                             length >= 500 and coverage >= 2.0. The header looks
#                             like ">NODE_1_length_12345_cov_67.8", so splitting it
#                             on "_" puts the length in field 4 and the coverage in
#                             field 6 — hence awk -F"_".
# Produces: contigs = DRAFT_CONTIGS, 02.assembly/{sample}/contigs_filt.fasta — the
#           hand-off that the contamination screen in shared/10_decontam.smk reads.
# Consumed by: index_contigs, blast_contigs, blob_json and select_contigs, all in
#              shared/10_decontam.smk. map_contigs reads the Bowtie2 index that
#              index_contigs builds from this file, not the FASTA itself.
#
# The {NAME:q} form lets Snakemake quote each awk program safely for the shell;
# writing the awk inline would need every brace doubled.
#
# conda: NONE — cat and awk only, as in v1.
#
# (v1 message: "Remove short contigs <500 bp and low coverage <2x.")
rule filter_contigs:
    input:
        contigs = SPADES_CONTIGS,
    output:
        contigs = DRAFT_CONTIGS,
    shell:
        """
        cat {input.contigs} | \
        awk {FASTA_LIN_CMD:q} | \
        awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}
        """
