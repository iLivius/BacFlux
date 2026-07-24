# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — hybrid front end, ONT assembly + consensus
# (rules/hybrid/40_ont_assembly.smk)
#
# The ONT leg proper: assemble the filtlong-selected long reads, rotate each
# circular replicon to a conventional start, and polish the consensus. The
# delivered genome comes out of the NEXT module (50_polish.smk, Polypolish).
#
# NOTE: `ont_assembly` and `fix_start` below are IDENTICAL to
# rules/nanopore/20_assembly.smk — keep them in sync. `long_read_consensus`
# differs from nanopore's in exactly two places, both flagged in its comment.
#
# Data flow through this module:
#
#   FILT_LONG ──► ont_assembly (Flye) ──► FLYE_CONTIGS ──┐
#   (30_ont_reads)      │                                │
#                       ├──► FLYE_INFO ──► FLYE_IGNORE_LIST
#                       │    (circularity)      │        │
#                       │                       ▼        ▼
#                       │                     fix_start (dnaapler)
#                       │           ┌─────────────┼──────────────┐
#                       │           ▼             ▼              ▼
#                       │     DNAAPLER_     DNAAPLER_FIXED  DNAAPLER_SUMMARY
#                       │     REORIENTED          │               │
#                       │                         ▼               │
#                       │            long_read_consensus (Medaka) │
#                       │              [only when USE_MEDAKA]     │
#                       │                         │               │
#                       │                         ▼               │
#                       │                 MEDAKA_CONSENSUS        │
#                       │                         │               │
#                       │                         ▼               │
#                       │        short_read_correction (50_polish.smk)
#                       │                    -> FINAL_CONTIGS     │
#                       └─────────────────────────────────────────┴──►
#                                     shared/15_replicons.smk (build_replicons)
#
# THE ONT LEG IS NEVER DECONTAMINATED HERE. Unlike nanopore mode, Medaka polishes
# the PRE-screen reoriented assembly (DNAAPLER_FIXED), because the ONT reads were
# already filtered against the decontaminated Illumina reads in 30_ont_reads.smk.
# There is no ONT-side screen to wait for, and adding one would be redundant.
#
# Everything referenced here comes from 00_common.smk.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: ont_assembly — long-read de novo assembly (Flye) ───────────────────
# NOTE: identical to rules/nanopore/20_assembly.smk — keep in sync.
#
# Takes in: FILT_LONG (short-read-guided filtlong output, 30_ont_reads.smk).
# Does:     Flye in {FLYE_INPUT_MODE} (--nano-hq or --nano-raw, resolved once in
#           00_common section 8) with 5 internal polishing iterations, then one
#           awk pass over Flye's own summary table.
# Produces: FLYE_DIR plus FLYE_CONTIGS, FLYE_INFO and FLYE_IGNORE_LIST.
# Consumed by: fix_start (contigs + ignore list), compare_hybrid_assemblies
#              (50_polish.smk, FLYE_CONTIGS) and build_replicons
#              (shared/15_replicons.smk, FLYE_INFO -> the topology column).
#
# THE IGNORE LIST: dnaapler ROTATES sequences, which only makes sense for a
# circular molecule — rotating a linear contig would move its real ends into the
# middle. IGNORE_LIST_CMD (00_common) prints every contig whose "circ." column in
# assembly_info.txt is not "Y", and dnaapler is told to leave those alone. The
# file is legitimately 0 bytes when every contig is circular; dnaapler copes.
#
# LAYOUT LANDMINE: declare directory(FLYE_DIR), never the sample directory — the
# sample directory holds other rules' output and a directory() output is deleted
# and recreated on every re-run.
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

        # Contigs Flye did NOT call circular (column "circ." != "Y").
        awk {IGNORE_LIST_CMD:q} {output.flye_info} > {output.ignore_list}
        """


# ── Rule: fix_start — rotate each circular replicon to a standard start ──────
# NOTE: identical to rules/nanopore/20_assembly.smk — keep in sync.
#
# `dnaapler all` looks for dnaA (chromosome), repA (plasmid), terL (phage
# terminase) and cog1474 in one pass and rotates each circular contig so the best
# hit starts at position 1 on the forward strand. Two things come out of that: the
# assembly becomes comparable to any other assembly of the strain, and the marker
# that was found is itself evidence of what kind of replicon each contig is —
# which is what shared/15_replicons.smk turns into Bakta's --replicons table.
#
# Takes in: FLYE_CONTIGS + FLYE_IGNORE_LIST.
# Does:     dnaapler (e-value 1e-10, fixed seed 42 for reproducibility), then
#           FASTA_HEAD_CMD (cut headers at the first whitespace token — every
#           later join in the pipeline keys on that token) and FASTA_LIN_CMD (one
#           line per sequence).
# Produces: DNAAPLER_DIR, DNAAPLER_REORIENTED, DNAAPLER_FIXED, DNAAPLER_SUMMARY.
# Consumed by: long_read_consensus (below) or, when Medaka is off, directly by
#              short_read_correction (50_polish.smk); compare_hybrid_assemblies;
#              and build_replicons (shared/15_replicons.smk).
#
# NOTE that in hybrid mode DNAAPLER_FIXED is NOT DRAFT_CONTIGS — the draft that
# gets screened is the ILLUMINA assembly. That is the one structural difference
# from nanopore mode, and it is expressed entirely through 00_common's constants.
#
# DNAAPLER_SUMMARY is NEW as a declared output (v1 wrote it but never declared or
# used it).
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


# Medaka only exists when the user did not switch it off. USE_MEDAKA is resolved
# once at parse time in 00_common section 8.
if USE_MEDAKA:

    # ── Rule: long_read_consensus — ONT consensus polishing (Medaka) ─────────
    # Same body as rules/nanopore/30_polish.smk, with TWO deliberate differences:
    #   1. input.contigs = MEDAKA_INPUT, which in HYBRID mode is DNAAPLER_FIXED —
    #      the PRE-decontamination reoriented assembly. The ONT leg has no screen
    #      of its own (see this file's banner); in nanopore mode the same constant
    #      is DECONTAM_CONTIGS instead.
    #   2. threads: capped_cpus(8) rather than capped_cpus(24). Both are the v1
    #      values for their own mode and are kept as-is for now; unify after the
    #      end-to-end gate if the timings say so.
    #
    # The model NAME comes pre-validated/resolved from check_medaka_model
    # (shared/12_medaka_check.smk), the same as nanopore. This is one v1->v2 change
    # worth noting: v1 hybrid auto-inference read the RAW ONT FASTQ, on the worry
    # that filtlong might drop a tag-carrying read. The shared check resolves from
    # FILT_LONG instead (uncompressed, headers preserved — the same reasoning the
    # Stage-4 review used to put nanopore on filtlong). In every realistic case the
    # tag is identical in both, so the resolved model is unchanged; the modes are
    # now aligned and a bad model fails before assembly rather than after.
    #
    # Takes in:
    #   reads   = FILT_LONG          the reads Flye assembled
    #   contigs = MEDAKA_INPUT       = DNAAPLER_FIXED here
    #   model   = MEDAKA_MODEL_RESOLVED  the confirmed/resolved model name
    # Produces: MEDAKA_DIR + MEDAKA_CONSENSUS.
    # Consumed by: short_read_correction and compare_hybrid_assemblies
    #              (50_polish.smk).
    #
    # (v1 message: "--- Medaka: Improve contig consensus with long reads. ---")
    rule long_read_consensus:
        input:
            reads = FILT_LONG,
            contigs = MEDAKA_INPUT,
            model = MEDAKA_MODEL_RESOLVED,
        output:
            consensus_dir = directory(MEDAKA_DIR),
            consensus_contigs = MEDAKA_CONSENSUS,
        conda:
            "../../envs/medaka.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/long_read_consensus_{sample}.log"
        priority: 9
        shell:
            """
            model=$(cat {input.model})
            echo "Polishing with Medaka model: $model" > {log}
            medaka_consensus \
              -i {input.reads} \
              -d {input.contigs} \
              -t {threads} \
              -m "$model" \
              -o {output.consensus_dir} >> {log} 2>&1
            """
