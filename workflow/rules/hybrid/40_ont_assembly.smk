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
    #   2. resources: capped_cpus(8) rather than capped_cpus(24). Both are the v1
    #      values for their own mode and are kept as-is for now; unify after the
    #      end-to-end gate if the timings say so.
    #
    # Takes in:
    #   reads    = FILT_LONG        the reads Flye assembled
    #   contigs  = MEDAKA_INPUT     = DNAAPLER_FIXED here
    #   raw_long = the RAW ONT FASTQ, used only by `medaka tools resolve_model`,
    #              because the basecaller tag lives in the original read headers
    #              and filtlong is not guaranteed to keep a read that carries it.
    # Produces: MEDAKA_DIR + MEDAKA_CONSENSUS.
    # Consumed by: short_read_correction and compare_hybrid_assemblies
    #              (50_polish.smk).
    #
    # params.model is the explicit model name, or "" when the user asked for
    # automatic inference — that empty string is what the bash if/else branches on.
    #
    # (v1 message: "--- Medaka: Improve contig consensus with long reads. ---")
    rule long_read_consensus:
        input:
            reads = FILT_LONG,
            contigs = MEDAKA_INPUT,
            raw_long = os.path.join(NANOPORE_DIR, ONT),
        output:
            consensus_dir = directory(MEDAKA_DIR),
            consensus_contigs = MEDAKA_CONSENSUS,
        params:
            model = MEDAKA_MODEL if MEDAKA_MODEL is not None else "",
        conda:
            "../../envs/medaka.yaml"
        threads: capped_cpus(8)
        log:
            LOGS + "/long_read_consensus_{sample}.log"
        priority: 9
        shell:
            # A NOTE ON BACKSLASHES BELOW: this shell block is an ordinary Python
            # string, so a backslash-n written once would become a real newline
            # before bash ever sees it. Where bash itself needs the two characters
            # \n (in printf and tr), they are written doubled: \\n.
            """
            if [ -n "{params.model}" ]; then

              # ── The user named a model explicitly ─────────────────────────
              model="{params.model}"

              # PRE-FLIGHT CHECK: a wrong model name is the most common Medaka
              # failure, and without this check it only surfaces after Medaka has
              # already loaded the reads. Skip the check when the value points at
              # a local model file on disk.
              if [ ! -e "$model" ]; then
                models=$(medaka tools list_models 2>&1) || {{
                  {{
                    echo "ERROR: Unable to query Medaka models in the current conda environment."
                    echo "  medaka_model: $model"
                    echo "  reason: 'medaka tools list_models' failed."
                    echo ""
                    echo "$models"
                  }} > {log}
                  cat {log} >&2
                  exit 1
                }}
                if ! printf '%s\\n' "$models" | sed -n 's/^Available: //p' | tr ',' '\\n' | sed 's/^ *//; s/ *$//' | grep -Fxq "$model"; then
                  {{
                    echo "ERROR: Invalid Medaka model configured for consensus polishing."
                    echo "  medaka_model: $model"
                    echo "  reason: this model is not available in the current Medaka environment."
                    echo "  next steps: choose a model from 'medaka tools list_models', set 'parameters.hybrid.medaka_model' to auto to infer it from the FASTQ headers, set it to FALSE to skip Medaka, or use a Medaka 1.x environment if this exact legacy model is required."
                  }} > {log}
                  cat {log} >&2
                  exit 1
                fi
              fi

              medaka_consensus \
                -i {input.reads} \
                -d {input.contigs} \
                -t {threads} \
                -m "$model" \
                -o {output.consensus_dir} > {log} 2>&1 || {{
                  cat {log}
                  echo "" >&2
                  echo "Medaka failed while using the explicit model setting '$model'." >&2
                  echo "If this is an older ONT model, it may be deprecated in Medaka v2 or its model weights may not be installed locally." >&2
                  echo "For Medaka v2, prefer a supported model such as 'r941_min_fast_g507' over deprecated legacy names like 'r941_min_fast_g303'." >&2
                  echo "Otherwise set 'parameters.hybrid.medaka_model' to FALSE to skip Medaka." >&2
                  exit 1
                }}

            else

              # ── Infer the model from the RAW read headers ─────────────────
              resolved_model=$(medaka tools resolve_model --auto_model consensus_bacteria {input.raw_long} 2> {log}) || {{
                cat {log}
                echo "" >&2
                echo "Medaka could not auto-infer a consensus model from {input.raw_long}." >&2
                echo "This usually means the ONT FASTQ headers do not contain exactly one basecaller model reference." >&2
                echo "Set 'parameters.hybrid.medaka_model' in the config to auto, an explicit Medaka model name, or FALSE to skip Medaka." >&2
                exit 1
              }}
              echo "Resolved Medaka model: $resolved_model" >> {log}

              medaka_consensus \
                -i {input.reads} \
                -d {input.contigs} \
                -t {threads} \
                -m "$resolved_model" \
                -o {output.consensus_dir} >> {log} 2>&1

            fi
            """
