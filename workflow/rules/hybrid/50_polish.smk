# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — hybrid front end, short-read correction + comparison
# (rules/hybrid/50_polish.smk)
#
# The biology: an ONT assembly has the right structure but still carries residual
# base errors, most of them indels in homopolymers, which frameshift genes and
# ruin an annotation. Illumina reads have the opposite profile — they cannot span
# repeats but they get the base right. Polypolish maps the short reads to the long
# -read assembly and corrects it, and it is specifically designed to do this in
# REPEATS, which is why each mate is aligned separately with `bwa mem -a` (report
# ALL alignments, not just the best one).
#
# THIS MODULE PRODUCES THE DELIVERED GENOME. short_read_correction writes
# FINAL_CONTIGS directly; hybrid needs no `finalize_contigs` rule because, unlike
# nanopore, its last step is unconditional (Polypolish always runs, whether or not
# Medaka did).
#
# Data flow through this module:
#
#   POLISH_INPUT ──┐                  (= MEDAKA_CONSENSUS if Medaka ran,
#   (40_ont_asm)   │                     else DNAAPLER_FIXED)
#   SEL_R1/SEL_R2 ─┴──► short_read_correction (bwa mem -a + polypolish) ──┐
#   (30_ont_reads)                                                        │
#                                                                         ▼
#                                                                  FINAL_CONTIGS
#                                                          = 02.assembly/{sample}/
#                                                              contigs_final.fasta
#                                                                         │
#          ┌──────────────────────────────────────────────────────────────┤
#          ▼                                                              ▼
#   compare_hybrid_assemblies (Snippy)                    the whole shared tail
#   -> SNPS_SUMMARY                                       (QC, taxonomy, Bakta,
#                                                          AMR, plasmids, phages)
#
# THE TWO GENOMES THIS MODE DELIVERS, and how shared/20_qc.smk finds them:
#   comparator "{sample}_illumina" = DECONTAM_CONTIGS (contaminants/contigs_sel.fasta),
#                                    written by the SHARED select_contigs rule
#   primary    "{sample}_ont"      = FINAL_CONTIGS, written HERE by Polypolish
# This front end's only obligation to QC_GENOMES is that both files exist;
# stage_qc_genomes in shared/20_qc.smk does the staging, naming and labelling. No
# rule in this module may write anything under 02.assembly/{sample}/eval/.
#
# Everything referenced here comes from 00_common.smk: POLISH_INPUT, SEL_R1/R2,
# POLYPOLISH_DIR, FINAL_CONTIGS, DECONTAM_CONTIGS, FLYE_CONTIGS, DNAAPLER_FIXED,
# MEDAKA_CONSENSUS, USE_MEDAKA, SNPS_DIR, SNPS_SUMMARY, FASTA_HEAD_CMD,
# FASTA_LIN_CMD, LOGS, capped_cpus.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: short_read_correction — polish the ONT genome with Illumina reads ──
# Takes in:
#   draft_contigs = POLISH_INPUT, a PARSE-TIME constant from 00_common: the Medaka
#                   consensus when Medaka ran, otherwise the dnaapler-reoriented
#                   assembly. Being parse-time means the DAG is fixed before the
#                   run and `snakemake -n` shows the real chain (v1 used an input
#                   FUNCTION, polishing_input_contigs(wc), that re-read the config
#                   while the DAG was being built).
#   r1 / r2       = SEL_R1 / SEL_R2, the decontaminated Illumina pairs from
#                   map_sel_contigs (30_ont_reads.smk).
# Does, in order:
#   1. copy the draft into the working directory (bwa writes its index files next
#      to the reference, and we do not want those landing in another rule's dir);
#   2. `bwa mem -a` ONCE PER MATE, separately — Polypolish requires unpaired,
#      all-alignment SAMs so it can see every place a read could have come from;
#   3. `polypolish filter` removes alignments that are inconsistent with the
#      insert size, which is what makes repeat correction safe;
#   4. `polypolish polish` writes the corrected sequence;
#   5. FASTA_HEAD_CMD + FASTA_LIN_CMD normalise the headers to their first token
#      and put each sequence on one line — the same normalisation the rest of the
#      pipeline joins on (Bakta --replicons, Platon, geNomad, the BLAST screen);
#   6. delete the five bwa index files, which are not declared outputs.
# Produces: six temp() working files inside POLYPOLISH_DIR, and FINAL_CONTIGS.
# Consumed by: compare_hybrid_assemblies (below) and every shared downstream module.
#
# HOUSE RULE — the directory is NOT declared. v1 nested six temp() files inside
# directory("07.Illumina_correction/{sample}"); v2 forbids that, because a
# directory() output is deleted and recreated on a re-run, which would also remove
# files another rule may have put there. The temps are declared individually and
# the directory is created with mkdir -p in the shell.
#
# (v1 message: "--- Polypolish: Correct contigs with short reads. ---")
rule short_read_correction:
    input:
        draft_contigs = POLISH_INPUT,
        r1 = SEL_R1,
        r2 = SEL_R2,
    output:
        draft_copy       = temp(POLYPOLISH_DIR + "/draft.fasta"),
        sam_1            = temp(POLYPOLISH_DIR + "/align_1.sam"),
        sam_2            = temp(POLYPOLISH_DIR + "/align_2.sam"),
        filt_sam_1       = temp(POLYPOLISH_DIR + "/filt_align_1.sam"),
        filt_sam_2       = temp(POLYPOLISH_DIR + "/filt_align_2.sam"),
        polished_contigs = temp(POLYPOLISH_DIR + "/polished_contigs.fasta"),
        final = FINAL_CONTIGS,
    conda:
        "../../envs/polypolish.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/short_read_correction_{sample}.log"
    priority: 9
    shell:
        """
        mkdir -p $(dirname {output.draft_copy})

        cp {input.draft_contigs} {output.draft_copy}

        bwa index {output.draft_copy} > {log} 2>&1
        bwa mem -t {resources.cpus} -a {output.draft_copy} {input.r1} 2>> {log} > {output.sam_1}
        bwa mem -t {resources.cpus} -a {output.draft_copy} {input.r2} 2>> {log} > {output.sam_2}

        polypolish filter \
          --in1 {output.sam_1} \
          --in2 {output.sam_2} \
          --out1 {output.filt_sam_1} \
          --out2 {output.filt_sam_2} >> {log} 2>&1

        polypolish polish \
          {output.draft_copy} \
          {output.filt_sam_1} \
          {output.filt_sam_2} 2>> {log} > {output.polished_contigs}

        # Normalise headers (first token only) and linearise -> the delivered genome.
        awk {FASTA_HEAD_CMD:q} {output.polished_contigs} | \
        awk {FASTA_LIN_CMD:q} > {output.final}

        # bwa's index files are not declared outputs; remove them explicitly.
        rm -f \
          {output.draft_copy}.amb \
          {output.draft_copy}.ann \
          {output.draft_copy}.bwt \
          {output.draft_copy}.pac \
          {output.draft_copy}.sa
        """


# ── Rule: compare_hybrid_assemblies — how much did each step change? ─────────
# Biology: hybrid mode builds the same genome twice, by two technologies with
# opposite error profiles. Aligning each stage of the ONT genome against the
# decontaminated Illumina assembly and counting the variants shows what the ONT
# pipeline actually changed: the Flye assembly, then reorientation, then Medaka,
# then Polypolish. A well-behaved run shows the variant count falling towards zero
# as polishing proceeds.
#
# Takes in (a star comparison — one reference, four query assemblies):
#   ref            = DECONTAM_CONTIGS, the decontaminated ILLUMINA assembly
#   flye_contigs   = FLYE_CONTIGS      (01. raw long-read assembly)
#   dnaapler       = DNAAPLER_FIXED    (02. after reorientation)
#   medaka         = MEDAKA_CONSENSUS  (03. after long-read consensus) — or [] when
#                    Medaka is off, which Snakemake renders as an EMPTY STRING in
#                    the shell. The params.use_medaka guard means that empty value
#                    is never used. (v1 pointed this input at the SAME file as
#                    dnaapler_contigs purely to keep the DAG valid; needlessly
#                    confusing, so v2 binds it to [].)
#   polished       = FINAL_CONTIGS     (04. after short-read correction)
# Produces: four Snippy output directories plus SNPS_SUMMARY, a plain-text digest.
# Consumed by: nobody — it is a terminal report, which is exactly why 00_common's
#              _frontend_targets_for() has to REQUEST SNPS_SUMMARY explicitly.
#
# `|| true` ON THE FOUR grep LINES is new in v2 and is a real fix: Snakemake runs
# the shell under bash strict mode, and grep exits 1 when it matches nothing, so a
# Snippy run that produced no "Variant" line would abort the whole rule. The
# summary should simply have an empty section instead.
#
# THE INTERPRETATION CAVEAT AT THE END IS KEPT VERBATIM FROM v1 and matters: a
# variant count of zero does NOT mean two assemblies are identical. Snippy only
# reports small variants in ALIGNED regions; contig joins, structural differences
# and unaligned sequence are invisible to it.
#
# (v1 message: "--- Snippy: Differences between long-read and short-read
#  assemblies. ---")
rule compare_hybrid_assemblies:
    input:
        sel_contigs = DECONTAM_CONTIGS,
        flye_contigs = FLYE_CONTIGS,
        dnaapler_contigs = DNAAPLER_FIXED,
        medaka_contigs = MEDAKA_CONSENSUS if USE_MEDAKA else [],
        polished_contigs = FINAL_CONTIGS,
    output:
        flye_snps_dir = directory(SNPS_DIR + "/01.flye_snps_dir"),
        dnaapler_snps_dir = directory(SNPS_DIR + "/02.dnaapler_snps_dir"),
        medaka_snps_dir = directory(SNPS_DIR + "/03.medaka_snps_dir"),
        polypolish_snps_dir = directory(SNPS_DIR + "/04.polypolish_snps_dir"),
        snps_summary = SNPS_SUMMARY,
    params:
        use_medaka = "true" if USE_MEDAKA else "false",
    conda:
        "../../envs/snippy.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/compare_hybrid_assemblies_{sample}.log"
    priority: 5
    shell:
        """
        snippy \
          --prefix "01.long-read_assembly" \
          --ref {input.sel_contigs} \
          --ctgs {input.flye_contigs} \
          --cpus {resources.cpus} \
          --outdir {output.flye_snps_dir} > {log} 2>&1

        snippy \
          --prefix "02.replicon_reorientation" \
          --ref {input.sel_contigs} \
          --ctgs {input.dnaapler_contigs} \
          --cpus {resources.cpus} \
          --outdir {output.dnaapler_snps_dir} >> {log} 2>&1

        if [[ "{params.use_medaka}" == "true" ]]; then
          snippy \
            --prefix "03.long-read_correction" \
            --ref {input.sel_contigs} \
            --ctgs {input.medaka_contigs} \
            --cpus {resources.cpus} \
            --outdir {output.medaka_snps_dir} >> {log} 2>&1
        else
          mkdir -p {output.medaka_snps_dir}
          printf "%s\n" "Medaka skipped by configuration." > {output.medaka_snps_dir}/skipped.txt
        fi

        snippy \
          --prefix "04.short-read_correction" \
          --ref {input.sel_contigs} \
          --ctgs {input.polished_contigs} \
          --cpus {resources.cpus} \
          --outdir {output.polypolish_snps_dir} >> {log} 2>&1

        echo "01. Long-read assembly" > {output.snps_summary}
        grep "Variant" {output.flye_snps_dir}/*.txt >> {output.snps_summary} || true

        echo -e "\n02. Replicon reorientation" >> {output.snps_summary}
        grep "Variant" {output.dnaapler_snps_dir}/*.txt >> {output.snps_summary} || true

        echo -e "\n03. Long-read correction" >> {output.snps_summary}
        if [[ "{params.use_medaka}" == "true" ]]; then
          grep "Variant" {output.medaka_snps_dir}/*.txt >> {output.snps_summary} || true
        else
          echo "Medaka was skipped." >> {output.snps_summary}
        fi

        echo -e "\n04. Short-read correction" >> {output.snps_summary}
        grep "Variant" {output.polypolish_snps_dir}/*.txt >> {output.snps_summary} || true

        echo -e "\nShort-read assembled selected contigs were used as reference." >> {output.snps_summary}
        echo "Interpretation note:" >> {output.snps_summary}
        echo "Variant counts reported above are based on Snippy comparisons against the short-read assembled selected contigs used as reference." >> {output.snps_summary}
        echo "These values summarize only small variants detected in aligned regions and should not be interpreted as proof that two assemblies are globally identical." >> {output.snps_summary}
        echo "Differences in contig number, contig joins, genome structure, unaligned regions, or extra sequence may still be present even when VariantTotal is 0." >> {output.snps_summary}
        echo "For this reason, the SNP summary should be interpreted together with assembly size, contig count, and the broader assembly evaluation results." >> {output.snps_summary}
        """
