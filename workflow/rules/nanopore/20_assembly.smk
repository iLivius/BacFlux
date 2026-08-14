# BacFlux v2.0.0 — nanopore front end, assembly.
#
# ONT reads are long enough to span the repeats that break a short-read assembly,
# so Flye usually closes a bacterial chromosome into a single circular contig. A
# circular sequence has no natural start, though, and Flye's arbitrary start point
# makes two assemblies of the same isolate look different. dnaapler rotates each
# replicon so it begins at a conventional landmark — dnaA for a chromosome, repA
# for a plasmid, terL for a phage — which makes assemblies comparable and puts the
# origin of replication where a reader expects it.
#
# Stage chain: FILT_LONG → ont_assembly → FLYE_CONTIGS (plus FLYE_INFO and the
# ignore list derived from it) → fix_start → DNAAPLER_FIXED, which IS
# DRAFT_CONTIGS → the screen in shared/10_decontam.smk → DECONTAM_CONTIGS →
# Medaka.
#
# ont_assembly  : one Flye run plus one awk pass over Flye's own
#                 assembly_info.txt, to list the contigs dnaapler must not
#                 rotate.
# fix_start     : `dnaapler all`, then two awk passes that trim and linearise its
#                 FASTA. Also keeps DNAAPLER_SUMMARY, which tells
#                 shared/15_replicons.smk which contigs are chromosomes and which
#                 are plasmids.
#
# D3 ORDERING, preserved from v1: dnaapler runs BEFORE the contamination screen
# and Medaka runs AFTER it. That is why DRAFT_CONTIGS (what the screen reads) is
# the dnaapler output, and why nanopore/30_polish.smk polishes DECONTAM_CONTIGS.
#
# Everything referenced here comes from 00_common.smk: FILT_LONG, FLYE_DIR,
# FLYE_CONTIGS, FLYE_INFO, FLYE_IGNORE_LIST, FLYE_INPUT_MODE, USE_MEDAKA,
# MEDAKA_MODEL_RESOLVED, DNAAPLER_DIR, DNAAPLER_REORIENTED, DNAAPLER_FIXED
# (== DRAFT_CONTIGS), DNAAPLER_SUMMARY, IGNORE_LIST_CMD, FASTA_HEAD_CMD,
# FASTA_LIN_CMD, LOGS, CPUS, capped_cpus.


# ──────────────────────── Long-read assembly (Flye) ────────────
# Takes in:
#   filt_long = FILT_LONG            the filtlong-selected reads, from rule
#                                    filter_long_reads (nanopore/10_reads.smk).
#   medaka_ok = MEDAKA_MODEL_RESOLVED from check_medaka_model
#                                    (shared/12_medaka_check.smk), or [] when
#                                    Medaka is off. A GATE only — the shell never
#                                    reads it; see the inline note on the input.
# Does:     one Flye de novo assembly, then one awk pass over Flye's own summary
#           table. --iterations 5 is Flye's internal long-read polishing, run five
#           times (the v1 value). {FLYE_INPUT_MODE} is --nano-hq or --nano-raw,
#           resolved once at parse time in 00_common.smk section 8: in auto mode
#           only an explicit *fast* Medaka model selects --nano-raw (fast
#           basecalling means noisier reads), otherwise --nano-hq, and
#           parameters.nanopore.flye_input_mode overrides the choice outright.
# Produces:
#   flye_dir     = FLYE_DIR, 02.assembly/{sample}/flye/ — Flye's own working and
#                  result tree
#   flye_contigs = FLYE_CONTIGS, assembly.fasta — the assembly itself
#   flye_info    = FLYE_INFO, assembly_info.txt — per contig length, coverage,
#                  circularity, repeat status
#   ignore_list  = FLYE_IGNORE_LIST, ignore_list.txt — the contig names Flye did
#                  NOT flag circular. The awk line in the shell writes this one;
#                  the three above are Flye's own files, which we only declare.
# Consumed by: fix_start below, which reads the contigs and the ignore list, and
#              build_replicons (shared/15_replicons.smk), which turns FLYE_INFO's
#              circularity column into the topology column of the Bakta
#              --replicons table. v1 declared FLYE_INFO as an output too, but
#              nothing read it.
#
# In hybrid mode the twin of this rule in hybrid/40_ont_assembly.smk also feeds
# FLYE_CONTIGS to compare_hybrid_assemblies (hybrid/50_polish.smk).
#
# Why the ignore list: dnaapler rotates a sequence, which only makes sense for a
# circular molecule. Rotating a linear contig would move its true ends into the
# middle. IGNORE_LIST_CMD (00_common.smk) reads assembly_info.txt and prints the
# name of every contig whose "circ." column is not "Y"; fix_start then tells
# dnaapler to leave those alone. The file is legitimately EMPTY (0 bytes) when
# every contig is circular — dnaapler handles that, verified in v1 runs.
#
# LAYOUT LANDMINE (same as SPAdes in illumina mode): declare directory(FLYE_DIR),
# never the sample folder. The sample folder also holds contaminants/, eval/,
# fix_start/ and medaka/ written by other rules, and Snakemake wipes and recreates
# a directory() output on a re-run.
#
# (v1 rule name: `assembly` in BacFluxL, `ONT_assembly` in BacFluxL+; D5 unifies
#  on `ont_assembly`. v1 message: "--- Flye: Genome assembly with long reads. ---")
rule ont_assembly:
    input:
        filt_long = FILT_LONG,
        # GATE (not used in the shell): don't spend an hour assembling if the
        # Medaka model is already known bad. check_medaka_model writes this after
        # read filtering; [] when Medaka is off, so there is no gate then.
        medaka_ok = MEDAKA_MODEL_RESOLVED if USE_MEDAKA else [],
    output:
        flye_dir = directory(FLYE_DIR),
        flye_contigs = FLYE_CONTIGS,
        flye_info = FLYE_INFO,
        ignore_list = FLYE_IGNORE_LIST,
    params:
        input_mode = FLYE_INPUT_MODE,
        iterations = 5,
    conda:
        "../../envs/flye.yaml"
    threads: CPUS
    log:
        LOGS + "/ont_assembly_{sample}.log"
    priority: 10
    shell:
        """
        flye \
          {params.input_mode} \
          {input.filt_long} \
          --out-dir {output.flye_dir} \
          --threads {threads} \
          --iterations {params.iterations} > {log} 2>&1

        # Contigs Flye did NOT call circular (column "circ." != "Y") → dnaapler
        # must not rotate them. Positional awk, unchanged from v1.
        awk {IGNORE_LIST_CMD:q} {output.flye_info} > {output.ignore_list}
        """


# ──────────────────── Replicon rotation (dnaapler) ─────────────
# `dnaapler all` searches every contig for the three canonical start genes at once
# — dnaA (chromosome), repA (plasmid), terL (phage terminase) — plus cog1474 for
# archaeal-type origins, and rotates the sequence so the best hit begins at
# position 1 on the forward strand. Two consequences are used later: assemblies of
# the same strain become directly comparable, and the marker that was found is
# itself evidence of what kind of replicon the contig is.
#
# Takes in:
#   flye_contigs = FLYE_CONTIGS      the Flye assembly, from ont_assembly above.
#   ignore_list  = FLYE_IGNORE_LIST  the contigs Flye did not call circular, from
#                                    the same rule — the ones not to rotate.
# Does:     dnaapler at e-value 1e-10 with a fixed seed of 42, so a re-run is
#           reproducible. Two awk passes then run over dnaapler's FASTA:
#             FASTA_HEAD_CMD — cut each header at the first whitespace token. This
#                              matters far beyond tidiness: Bakta's --replicons
#                              table, Platon, geNomad and the BlobTools/BLAST
#                              screen all join on that first token, so trimming
#                              here is what makes every later join work.
#             FASTA_LIN_CMD  — one line per sequence.
# Produces:
#   dnaapler_dir     = DNAAPLER_DIR, 02.assembly/{sample}/fix_start/
#   dnaapler_contigs = DNAAPLER_REORIENTED, {sample}_reoriented.fasta — dnaapler's
#                      own untouched output, kept alongside the trimmed one
#   fixed_contigs    = DNAAPLER_FIXED, {sample}_fixed.fasta — the same sequences
#                      with headers trimmed and one line per sequence. In nanopore
#                      mode this file IS DRAFT_CONTIGS; 00_common.smk asserts it.
#   summary          = DNAAPLER_SUMMARY, {sample}_all_reorientation_summary.tsv —
#                      which start gene was found on each contig, and how well
# Consumed by: the contamination screen in shared/10_decontam.smk, whose
#              map_contigs, blast_contigs, blob_json and select_contigs all read
#              DRAFT_CONTIGS; and build_replicons (shared/15_replicons.smk), which
#              reads DNAAPLER_SUMMARY. No rule reads DNAAPLER_DIR or
#              DNAAPLER_REORIENTED — they are kept for inspection only.
#
# DNAAPLER_SUMMARY is NEW as a declared output. v1 produced
# {sample}_all_reorientation_summary.tsv and never declared or used it, even
# though its Gene_Reoriented / Coverage / Identity_Percentage columns already
# say, per contig, which start gene was found and how convincingly. That is
# exactly the input build_replicons (shared/15_replicons.smk) needs to tell Bakta
# which contigs are chromosomes and which are plasmids, so we now declare it.
#
# (v1 message: "--- dnaapler: Re-orient replicons. ---")
rule fix_start:
    input:
        flye_contigs = FLYE_CONTIGS,
        ignore_list = FLYE_IGNORE_LIST,
    output:
        dnaapler_dir = directory(DNAAPLER_DIR),
        dnaapler_contigs = DNAAPLER_REORIENTED,
        fixed_contigs = DNAAPLER_FIXED,
        summary = DNAAPLER_SUMMARY,
    params:
        evalue = 1e-10,
        seed = 42,
    conda:
        "../../envs/dnaapler.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/fix_start_{sample}.log"
    priority: 9
    shell:
        """
        dnaapler all \
          -i {input.flye_contigs} \
          -p {wildcards.sample} \
          -e {params.evalue} \
          --seed_value {params.seed} \
          -t {threads} \
          -o {output.dnaapler_dir} \
          --ignore {input.ignore_list} \
          --force > {log} 2>&1

        # Trim headers to their first token, then put each sequence on one line.
        awk {FASTA_HEAD_CMD:q} {output.dnaapler_contigs} | \
        awk {FASTA_LIN_CMD:q} > {output.fixed_contigs}
        """
