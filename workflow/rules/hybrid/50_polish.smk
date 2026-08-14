# Short-read correction of the ONT genome, and a stage-by-stage comparison of the
# two technologies. This module produces the delivered genome.
#
# The biology: an ONT assembly has the right structure but still carries residual
# base errors, most of them indels in homopolymers, which frameshift genes and
# ruin an annotation. Illumina reads have the opposite profile — they cannot span
# repeats but they get the base right. Polypolish maps the short reads onto the
# long-read assembly and corrects it, and it is built specifically to do this
# INSIDE repeats, which is why each mate is aligned separately with `bwa mem -a`
# (report ALL alignments, not just the best one).
#
# Chain: POLISH_INPUT (MEDAKA_CONSENSUS if Medaka ran, else DNAAPLER_FIXED, both
# from hybrid/40_ont_assembly.smk) plus SEL_R1/SEL_R2 (hybrid/30_ont_reads.smk) →
# short_read_correction → FINAL_CONTIGS = 02.assembly/{sample}/contigs_final.fasta
# → compare_hybrid_assemblies (below) and the whole shared tail (QC, taxonomy,
# Bakta, AMR, plasmids, phages).
#
# Hybrid needs no `finalize_contigs` rule, unlike nanopore, because its last step
# is unconditional: Polypolish always runs, whether or not Medaka did.
#
# The two genomes this mode delivers, and how shared/20_qc.smk finds them:
#   comparator "{sample}_illumina" = DECONTAM_CONTIGS
#                                    (contaminants/contigs_sel.fasta), written by
#                                    the shared select_contigs rule
#   primary    "{sample}_ont"      = FINAL_CONTIGS, written here by Polypolish
# This front end's only obligation to QC_GENOMES is that both files exist;
# stage_qc_genomes in shared/20_qc.smk does the staging, naming and labelling. No
# rule in this module may write anything under 02.assembly/{sample}/eval/.
#
# short_read_correction     : bwa mem -a per mate → polypolish filter → polish →
#                             header normalisation → FINAL_CONTIGS.
# compare_hybrid_assemblies : Snippy, four ONT stages against the Illumina
#                             assembly, digested into SNPS_SUMMARY.
#
# Every constant used here comes from shared/00_common.smk: POLISH_INPUT,
# SEL_R1/R2, POLYPOLISH_DIR, FINAL_CONTIGS, DECONTAM_CONTIGS, FLYE_CONTIGS,
# DNAAPLER_FIXED, MEDAKA_CONSENSUS, USE_MEDAKA, SNPS_DIR, SNPS_SUMMARY,
# FASTA_HEAD_CMD, FASTA_LIN_CMD, LOGS, capped_cpus.


# ────────────────────── Short-read correction (Polypolish) ─────
# Polish the ONT genome with the decontaminated Illumina pairs and publish the
# result as FINAL_CONTIGS — the genome every shared downstream module annotates.
#
# Takes in:
#   draft_contigs = POLISH_INPUT, a PARSE-TIME constant from
#                   shared/00_common.smk: the Medaka consensus
#                   (long_read_consensus) when Medaka ran, otherwise the
#                   dnaapler-reoriented assembly (fix_start), both in
#                   hybrid/40_ont_assembly.smk. Parse-time matters — the DAG is
#                   fixed before the run starts, so `snakemake -n` shows the real
#                   chain. v1 used an input FUNCTION, polishing_input_contigs(wc),
#                   that re-read the config while the DAG was being built.
#   r1, r2        = SEL_R1 / SEL_R2, the decontaminated Illumina pairs from
#                   map_sel_contigs (hybrid/30_ont_reads.smk).
# Does:     six things, in order:
#   1. copy the draft into the working directory, because bwa writes its index
#      files next to the reference and those must not land in another rule's
#      directory;
#   2. `bwa mem -a` once per mate, separately — Polypolish requires unpaired,
#      all-alignment SAMs so it can see every place a read could have come from;
#   3. `polypolish filter` drops alignments inconsistent with the insert size,
#      which is what makes correction inside repeats safe;
#   4. `polypolish polish` writes the corrected sequence;
#   5. FASTA_HEAD_CMD and FASTA_LIN_CMD cut each header to its first token and put
#      each sequence on one line — the same normalisation everything downstream
#      joins on (Bakta --replicons, Platon, geNomad, the BLAST screen);
#   6. delete the five bwa index files, which are not declared outputs.
# Produces:
#   draft_copy       = temp(), the working copy made in step 1
#   sam_1, sam_2     = temp(), the two per-mate `bwa mem -a` alignments
#   filt_sam_1,
#   filt_sam_2       = temp(), the same alignments after `polypolish filter`
#   polished_contigs = temp(), Polypolish's output before header normalisation
#   final            = FINAL_CONTIGS, 02.assembly/{sample}/contigs_final.fasta —
#                      the delivered genome
# Consumed by: compare_hybrid_assemblies below; stage_qc_genomes
#              (shared/20_qc.smk), which stages FINAL_CONTIGS as the "{sample}_ont"
#              primary genome and is also how taxonomic_assignment
#              (shared/30_taxonomy.smk) sees it; build_replicons
#              (shared/15_replicons.smk); blast_final_contigs
#              (shared/10_decontam.smk); annotation (shared/40_annotation.smk);
#              amr_contigs (shared/50_amr.smk); plasmid_search
#              (shared/60_plasmid.smk); the phage caller
#              (viral_identification_virsorter2, or genomad_end_to_end when
#              PHAGE_CALLER is genomad — shared/70_phage.smk); and, when the
#              mobilome module is on, amrfinderplus, isescan and six more rules
#              in shared/80_mobilome.smk.
#
# The working DIRECTORY is deliberately not declared. v1 nested six temp() files
# inside directory("07.Illumina_correction/{sample}"); v2 forbids that, because
# Snakemake deletes and recreates a directory() output on a re-run and would take
# any other rule's files in that directory with it. The temps are declared one by
# one and the directory is created with mkdir -p in the shell.
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
    threads: capped_cpus(24)
    log:
        LOGS + "/short_read_correction_{sample}.log"
    priority: 9
    shell:
        """
        mkdir -p $(dirname {output.draft_copy})

        cp {input.draft_contigs} {output.draft_copy}

        bwa index {output.draft_copy} > {log} 2>&1
        bwa mem -t {threads} -a {output.draft_copy} {input.r1} 2>> {log} > {output.sam_1}
        bwa mem -t {threads} -a {output.draft_copy} {input.r2} 2>> {log} > {output.sam_2}

        polypolish filter \
          --in1 {output.sam_1} \
          --in2 {output.sam_2} \
          --out1 {output.filt_sam_1} \
          --out2 {output.filt_sam_2} >> {log} 2>&1

        polypolish polish \
          {output.draft_copy} \
          {output.filt_sam_1} \
          {output.filt_sam_2} 2>> {log} > {output.polished_contigs}

        # Normalise headers (first token only) and linearise → the delivered genome.
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


# ────────────────────── Stage-by-stage comparison (Snippy) ─────
# Hybrid mode builds the same genome twice, with two technologies whose error
# profiles are opposites. Aligning each stage of the ONT genome against the
# decontaminated Illumina assembly and counting the variants shows what the ONT
# pipeline actually changed: the Flye assembly, then reorientation, then Medaka,
# then Polypolish. A well-behaved run shows the variant count falling towards zero
# as polishing proceeds.
#
# Takes in: a star comparison — one reference, four query assemblies:
#   sel_contigs      = DECONTAM_CONTIGS, the decontaminated ILLUMINA assembly and
#                      Snippy's --ref, from select_contigs
#                      (shared/10_decontam.smk)
#   flye_contigs     = FLYE_CONTIGS     (01. raw long-read assembly, ont_assembly)
#   dnaapler_contigs = DNAAPLER_FIXED   (02. after reorientation, fix_start)
#   medaka_contigs   = MEDAKA_CONSENSUS (03. after long-read consensus,
#                      long_read_consensus) — those three all come from
#                      hybrid/40_ont_assembly.smk
#   polished_contigs = FINAL_CONTIGS    (04. after short-read correction), from
#                      short_read_correction above
# Does:     four Snippy runs, one per ONT stage, each against the same Illumina
#           reference, then greps the "Variant" lines out of the four reports
#           into a single digest.
# Produces:
#   flye_snps_dir       = SNPS_DIR + "/01.flye_snps_dir"
#   dnaapler_snps_dir   = SNPS_DIR + "/02.dnaapler_snps_dir"
#   medaka_snps_dir     = SNPS_DIR + "/03.medaka_snps_dir"
#   polypolish_snps_dir = SNPS_DIR + "/04.polypolish_snps_dir", the four Snippy
#                         output directories, one per stage
#   snps_summary        = SNPS_SUMMARY, the plain-text digest of all four
# Consumed by: nobody — SNPS_SUMMARY is a terminal report for the user, which is
#              exactly why _frontend_targets_for() in shared/00_common.smk has to
#              request it BY NAME under `if IS_HYBRID`, or it would never be
#              built.
#
# When Medaka is off, medaka_contigs is bound to [], which Snakemake renders in
# the shell as an EMPTY STRING. The params.use_medaka guard means that empty value
# is never reached. v1 instead pointed this input at the same file as
# dnaapler_contigs purely to keep the DAG valid, which read as a real comparison
# and was not one.
#
# `|| true` on the four grep lines is new in v2 and fixes a real failure:
# Snakemake runs the shell under bash strict mode, and grep exits 1 when it
# matches nothing, so a Snippy run that produced no "Variant" line used to abort
# the whole rule. An empty section in the summary is the right outcome instead.
#
# The interpretation note printed at the end is kept verbatim from v1 and matters:
# a variant count of zero does NOT mean two assemblies are identical. Snippy
# reports only small variants in ALIGNED regions, so contig joins, structural
# differences and unaligned sequence are invisible to it.
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
    threads: capped_cpus(24)
    log:
        LOGS + "/compare_hybrid_assemblies_{sample}.log"
    priority: 5
    shell:
        """
        snippy \
          --prefix "01.long-read_assembly" \
          --ref {input.sel_contigs} \
          --ctgs {input.flye_contigs} \
          --cpus {threads} \
          --outdir {output.flye_snps_dir} > {log} 2>&1

        snippy \
          --prefix "02.replicon_reorientation" \
          --ref {input.sel_contigs} \
          --ctgs {input.dnaapler_contigs} \
          --cpus {threads} \
          --outdir {output.dnaapler_snps_dir} >> {log} 2>&1

        if [[ "{params.use_medaka}" == "true" ]]; then
          snippy \
            --prefix "03.long-read_correction" \
            --ref {input.sel_contigs} \
            --ctgs {input.medaka_contigs} \
            --cpus {threads} \
            --outdir {output.medaka_snps_dir} >> {log} 2>&1
        else
          mkdir -p {output.medaka_snps_dir}
          printf "%s\n" "Medaka skipped by configuration." > {output.medaka_snps_dir}/skipped.txt
        fi

        snippy \
          --prefix "04.short-read_correction" \
          --ref {input.sel_contigs} \
          --ctgs {input.polished_contigs} \
          --cpus {threads} \
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
