# The ONT leg proper: assemble the filtlong-selected long reads, rotate each
# circular replicon to a conventional start, and polish the consensus. The
# delivered genome comes out of the NEXT module — Polypolish, in
# hybrid/50_polish.smk.
#
# Chain: FILT_LONG (hybrid/30_ont_reads.smk) → ont_assembly → FLYE_CONTIGS,
# FLYE_INFO and FLYE_IGNORE_LIST → fix_start → DNAAPLER_REORIENTED,
# DNAAPLER_FIXED and DNAAPLER_SUMMARY → long_read_consensus (only when Medaka is
# on) → MEDAKA_CONSENSUS → short_read_correction (hybrid/50_polish.smk) →
# FINAL_CONTIGS. Two side branches leave this module: FLYE_INFO and
# DNAAPLER_SUMMARY go to build_replicons (shared/15_replicons.smk), which turns
# them into Bakta's --replicons table; FLYE_CONTIGS, DNAAPLER_FIXED and
# MEDAKA_CONSENSUS go to compare_hybrid_assemblies (hybrid/50_polish.smk).
#
# The ONT leg is never decontaminated here. Unlike nanopore mode, Medaka polishes
# the PRE-screen reoriented assembly (DNAAPLER_FIXED), because the ONT reads were
# already filtered against the decontaminated Illumina reads back in
# hybrid/30_ont_reads.smk. There is no ONT-side screen to wait for, and adding
# one would be redundant.
#
# ont_assembly        : Flye, plus one awk pass over Flye's own summary table to
#                       build the list of contigs dnaapler must not rotate.
# fix_start           : dnaapler all — rotate each circular contig to a standard
#                       start, and record which marker gene was found where.
# long_read_consensus : Medaka. Only exists when Medaka is switched on.
#
# ont_assembly and fix_start are identical to rules/nanopore/20_assembly.smk —
# keep them in sync. long_read_consensus differs from nanopore's copy of the same
# rule (rules/nanopore/30_polish.smk) in three places, all listed in its comment.
#
# Every constant used here comes from shared/00_common.smk.


# ──────────────────────── Long-read assembly (Flye) ────────────
# Takes in: filt_long = FILT_LONG, the short-read-guided filtlong output from
#           filter_long_reads (hybrid/30_ont_reads.smk). The second input,
#           medaka_ok, carries no data — it is a gate, explained on the input
#           line itself.
# Does:     Flye in {FLYE_INPUT_MODE} — either --nano-hq or --nano-raw, resolved
#           once at parse time in shared/00_common.smk section 8 — with 5
#           internal polishing iterations, then one awk pass (IGNORE_LIST_CMD)
#           over Flye's own summary table.
# Produces:
#   flye_dir     = FLYE_DIR, Flye's own output tree
#   flye_contigs = FLYE_CONTIGS, the long-read assembly
#   flye_info    = FLYE_INFO, Flye's per-contig summary table, which carries the
#                  length, coverage and circularity of each contig
#   ignore_list  = FLYE_IGNORE_LIST, the contigs dnaapler must not rotate
# Consumed by: fix_start below (contigs plus ignore list),
#              compare_hybrid_assemblies (hybrid/50_polish.smk, FLYE_CONTIGS) and
#              build_replicons (shared/15_replicons.smk, which takes the topology
#              column out of FLYE_INFO).
#
# Identical to rules/nanopore/20_assembly.smk — keep in sync.
#
# The ignore list exists because dnaapler ROTATES sequences, and rotating only
# makes sense for a circular molecule — do it to a linear contig and its real ends
# end up in the middle. IGNORE_LIST_CMD (shared/00_common.smk) prints every contig
# whose "circ." column in assembly_info.txt is not "Y", and dnaapler is told to
# leave those alone. The file is legitimately 0 bytes when every contig came out
# circular; dnaapler copes with that.
#
# Declare directory(FLYE_DIR), never the sample directory — the sample directory
# holds other rules' output, and Snakemake deletes and recreates a directory()
# output on every re-run.
rule ont_assembly:
    input:
        filt_long = FILT_LONG,
        # Not used by the shell: this is a gate. check_medaka_model
        # (shared/12_medaka_check.smk) writes MEDAKA_MODEL_RESOLVED straight after
        # read filtering, so a bad model name fails in seconds instead of after an
        # hour of assembly. Bound to [] when Medaka is off, which removes the gate.
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


# ────────────────────── Replicon reorientation (dnaapler) ──────
# `dnaapler all` looks for dnaA (chromosome), repA (plasmid), terL (phage
# terminase) and cog1474 in one pass, and rotates each circular contig so the best
# hit starts at position 1 on the forward strand. Two useful things fall out of
# that: the assembly becomes comparable to any other assembly of the strain, and
# the marker that was found is itself evidence of what kind of replicon each
# contig is — which is what shared/15_replicons.smk turns into Bakta's
# --replicons table.
#
# Takes in: flye_contigs = FLYE_CONTIGS and ignore_list = FLYE_IGNORE_LIST, both
#           from ont_assembly above.
# Does:     dnaapler (e-value 1e-10, fixed seed 42 so a re-run gives the same
#           answer), then FASTA_HEAD_CMD cuts each header at its first whitespace
#           token — every later join in the pipeline keys on that token — and
#           FASTA_LIN_CMD puts each sequence on one line.
# Produces:
#   dnaapler_dir     = DNAAPLER_DIR, dnaapler's own output tree
#   dnaapler_contigs = DNAAPLER_REORIENTED, the rotated assembly as dnaapler
#                      wrote it
#   fixed_contigs    = DNAAPLER_FIXED, the same sequences after header trimming
#                      and linearisation
#   summary          = DNAAPLER_SUMMARY, which marker gene was found on which
#                      contig
# Consumed by: long_read_consensus below — or, when Medaka is off,
#              short_read_correction (hybrid/50_polish.smk) directly;
#              compare_hybrid_assemblies (hybrid/50_polish.smk); and
#              build_replicons (shared/15_replicons.smk).
#
# Identical to rules/nanopore/20_assembly.smk — keep in sync.
#
# In hybrid mode DNAAPLER_FIXED is NOT DRAFT_CONTIGS: the draft that gets screened
# is the ILLUMINA assembly. That is the one structural difference from nanopore
# mode, and it is expressed entirely through 00_common's constants rather than
# through this rule.
#
# DNAAPLER_SUMMARY is new as a declared output — v1 wrote the file but never
# declared it and never used it.
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


# ──────────────────────── Consensus polishing (Medaka) ─────────
# Medaka exists only when the user did not switch it off. USE_MEDAKA is resolved
# once at parse time in shared/00_common.smk section 8 from
# parameters.hybrid.medaka_model (FALSE there means skip). When it is false the
# rule below is never defined at all, and POLISH_INPUT resolves to DNAAPLER_FIXED
# instead of MEDAKA_CONSENSUS, so the stage simply drops out of the DAG rather
# than being built and ignored.
if USE_MEDAKA:

    # ── Medaka consensus on the reoriented assembly ──
    # Takes in:
    #   reads   = FILT_LONG, the same reads Flye assembled
    #             (filter_long_reads, hybrid/30_ont_reads.smk)
    #   contigs = MEDAKA_INPUT, which in hybrid mode is DNAAPLER_FIXED, from
    #             fix_start above
    #   model   = MEDAKA_MODEL_RESOLVED, a one-line file holding the confirmed
    #             model name, from check_medaka_model (shared/12_medaka_check.smk)
    # Does:     medaka_consensus, re-calling the consensus base by base from the
    #           long reads with the basecaller-matched model.
    # Produces:
    #   consensus_dir     = MEDAKA_DIR, Medaka's own output tree
    #   consensus_contigs = MEDAKA_CONSENSUS, the polished assembly
    # Consumed by: short_read_correction and compare_hybrid_assemblies
    #              (hybrid/50_polish.smk).
    #
    # Same body as rules/nanopore/30_polish.smk, with three deliberate differences
    # — two you can see here, one hidden inside a constant:
    #   1. threads: capped_cpus(8) here rather than capped_cpus(24). Both are the
    #      v1 values for their own mode and are kept as-is for now; unify after
    #      the end-to-end gate if the timings say so.
    #   2. the failure hint names parameters.hybrid.medaka_model, where nanopore
    #      names parameters.nanopore.medaka_model.
    #   3. MEDAKA_INPUT resolves to DNAAPLER_FIXED here — the PRE-decontamination
    #      reoriented assembly, since the ONT leg has no screen of its own (see
    #      this file's header) — where in nanopore mode it is DECONTAM_CONTIGS.
    #
    # The model NAME arrives pre-validated from check_medaka_model
    # (shared/12_medaka_check.smk), exactly as in nanopore mode. That is one
    # v1 → v2 change worth remembering: v1 hybrid auto-inference read the RAW ONT
    # FASTQ, on the worry that filtlong might drop a tag-carrying read. The shared
    # check resolves from FILT_LONG instead (uncompressed, headers preserved — the
    # same reasoning the Stage-4 review used to put nanopore on filtlong). In every
    # realistic case the tag is identical in both, so the resolved model is
    # unchanged; the two modes are now aligned, and a bad model fails before
    # assembly rather than after it.
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
            # check_medaka_model already validated the NAME (or resolved auto), so
            # a failure here is a RUNTIME one — model weights not installed and no
            # network to fetch them, or a legacy model broken under Medaka v2. The
            # hint keeps that case actionable instead of a bare stack trace.
            """
            model=$(cat {input.model})
            echo "Polishing with Medaka model: $model" > {log}
            medaka_consensus \
              -i {input.reads} \
              -d {input.contigs} \
              -t {threads} \
              -m "$model" \
              -o {output.consensus_dir} >> {log} 2>&1 || {{
                cat {log} >&2
                echo "" >&2
                echo "Medaka failed at polishing with the pre-validated model '$model'." >&2
                echo "The name was accepted by check_medaka_model, so this is a RUNTIME failure: the model weights may not be installed locally (no network to fetch them), or a legacy model may be broken under Medaka v2." >&2
                echo "Set 'parameters.hybrid.medaka_model' to a supported model, or to FALSE to skip Medaka." >&2
                exit 1
              }}
            """
