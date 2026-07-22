# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — illumina front end, assembly (rules/illumina/20_assembly.smk)
#
# The biology: turn the cleaned read pairs into contigs, then throw away the
# contigs that cannot be real isolate sequence. SPAdes in --isolate mode is tuned
# for a single high-coverage bacterial culture, and its multi-k approach
# (21..127) resolves both low-complexity and repeat-rich regions in one run.
#
# Data flow through this module:
#
#   TRIM_R1 / TRIM_R2 ──► illumina_assembly (SPAdes) ──► SPADES_CONTIGS
#   (from 10_reads.smk)                                       │
#                                                             ▼
#                                                     filter_contigs
#                                                     (>=500 bp, >=2x cov)
#                                                             │
#                                                             ▼
#                                                       DRAFT_CONTIGS
#                                                             │
#                          ┌──────────────────────────────────┤
#                          ▼                                  ▼
#           shared/10_decontam.smk: BLAST + BlobTools    index/map_contigs
#           -> select_contigs -> DECONTAM_CONTIGS
#
# WHERE THIS MODE ENDS. In illumina mode, decontamination IS the last assembly
# step: 00_common sets DECONTAM_CONTIGS = FINAL_CONTIGS (the same Python object),
# so the shared `select_contigs` rule writes the canonical delivered genome and
# this front end must NOT declare FINAL_CONTIGS anywhere. Declaring it here would
# give Snakemake two rules producing one file.
#
# Everything referenced here comes from 00_common.smk: TRIM_R1/TRIM_R2,
# SPADES_DIR, SPADES_CONTIGS, DRAFT_CONTIGS, FASTA_LIN_CMD, FASTA_SEL_CMD, LOGS,
# CPUS, RAM.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: illumina_assembly — de novo assembly of the isolate (SPAdes) ───────
# Takes in: TRIM_R1 / TRIM_R2, the fastp-trimmed, PhiX-free pairs.
# Does:     one SPAdes run per sample in --isolate mode over six k-mer sizes.
#           OMP_NUM_THREADS is exported so the OpenMP parts of SPAdes obey the
#           same thread budget as the -t flag (without it they default to every
#           core on the machine, which oversubscribes a multi-sample run).
# Produces: SPADES_DIR (the whole SPAdes working directory) and, inside it,
#           contigs.fasta.
# Consumed by: filter_contigs, below.
#
# LAYOUT LANDMINE — do not "simplify" the output back to the sample directory.
# v1 declared directory("02.assembly/{sample}"). Under the v2 layout (D1) that
# directory ALSO holds contaminants/, eval/ and (in the long-read modes)
# fix_start/, medaka/, snps/ — written by OTHER rules. A directory() output is
# deleted and recreated by Snakemake when the rule re-runs, so declaring the
# sample directory would wipe every other stage's results on a re-run. The
# assembler therefore declares only its own sub-directory, SPADES_DIR.
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
    conda:
        "../../envs/spades.yaml"
    resources:
        # Uncapped as in v1: SPAdes is the single most expensive step in this mode
        # and scales with both cores and RAM. RAM is a hard ceiling in GB (-m).
        cpus = CPUS,
        ram = RAM
    log:
        LOGS + "/illumina_assembly_{sample}.log"
    priority: 10
    shell:
        """
        OMP_NUM_THREADS={resources.cpus} \
        spades.py -k 21,33,55,77,99,127 --isolate \
          --pe1-1 {input.r1} \
          --pe1-2 {input.r2} \
          -o {output.dir} \
          -t {resources.cpus} \
          -m {resources.ram} > {log} 2>&1
        """


# ── Rule: filter_contigs — drop short and low-coverage contigs ───────────────
# Biology: a SPAdes assembly of a pure isolate has a long tail of tiny,
# thinly-covered contigs — sequencing noise, chimeras, and fragments of whatever
# else was in the tube. They add nothing to the annotation but do inflate the
# contig count, the CheckM contamination estimate and the BLAST screening time.
# The cutoffs are the v1 ones: at least 500 bp AND at least 2x coverage.
#
# Takes in: SPADES_CONTIGS (raw SPAdes output).
# Does:     two small awk passes, both defined once in 00_common:
#             FASTA_LIN_CMD — put each record's sequence on ONE line, so the
#                             selector below can read a header and its sequence
#                             as a header/getline pair.
#             FASTA_SEL_CMD — keep a record when SPAdes' own header says
#                             length >= 500 and coverage >= 2.0. The header looks
#                             like ">NODE_1_length_12345_cov_67.8", so splitting
#                             it on "_" puts the length in field 4 and the
#                             coverage in field 6 — hence awk -F"_".
#           The {NAME:q} form lets Snakemake quote each awk program safely for the
#           shell; writing the awk inline would need every brace doubled.
# Produces: DRAFT_CONTIGS — the Stage-4 hand-off that shared/10_decontam.smk
#           screens (BLAST + BlobTools + selector).
# Consumed by: index_contigs, map_contigs, blast_contigs, blob_json and
#              select_contigs, all in shared/10_decontam.smk.
#
# conda: NONE — cat and awk only, as in v1.
#
# (v1 message: "Remove short contigs <500 bp and low coverage <2x.")
rule filter_contigs:
    input:
        contigs = SPADES_CONTIGS,
    output:
        contigs = DRAFT_CONTIGS,
    priority: 9
    shell:
        """
        cat {input.contigs} | \
        awk {FASTA_LIN_CMD:q} | \
        awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}
        """
