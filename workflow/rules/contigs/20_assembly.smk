# BacFlux v2.0.0 — contigs front end: make a user-supplied FASTA usable.
#
# This mode takes genomes that are ALREADY assembled — the user hands us FASTA
# files, not reads — so there is no read QC and no assembler. The whole front end
# is the one rule below, whose job is to hand the shared decontamination screen a
# FASTA it can work with.
#
# Stage chain: {sample}.fasta from config input.contigs_dir → filter_contigs →
# DRAFT_CONTIGS → map_contigs (self-mapping) and blast_contigs → BlobTools →
# select_contigs, all in shared/10_decontam.smk → DECONTAM_CONTIGS.
#
# Where this mode ends: as in illumina mode, decontamination is the last assembly
# step, so 00_common sets DECONTAM_CONTIGS = FINAL_CONTIGS (the same Python
# object) and the shared `select_contigs` rule writes the delivered genome. This
# front end must NOT declare FINAL_CONTIGS.
#
# OWNERSHIP — the self-mapping `map_contigs` (minimap2 of the contigs against
# themselves, which exists only to give BlobTools a BAM) is NOT defined here. It
# is owned by shared/10_decontam.smk, in the `else:` branch of its three-way
# parse-time split. Do not add a second copy.
#
# Everything referenced here comes from 00_common.smk: CONTIGS_DIR, CONTIGS,
# DRAFT_CONTIGS, FASTA_LIN_CMD, FASTA_SEL_CMD, FASTA_HEAD_CMD, LOGS.


# ─────────── Normalise an arbitrary input FASTA, two ways ──────
# The problem this solves: a user-supplied assembly can come from anywhere. If it
# came from SPAdes, its headers carry the length and the coverage of every contig
# and we can apply the same quality cutoffs the illumina mode uses. If it came
# from NCBI, from a different assembler, or from a colleague, those numbers do not
# exist, and inventing a filter would silently delete real sequence. So the rule
# LOOKS at the first header and branches.
#
# Takes in: contigs = the user's own genome, CONTIGS_DIR/CONTIGS — that is
#           {sample}.<ext>, for the single FASTA extension 00_common.smk accepted
#           for this run. NO rule produces it: CONTIGS_DIR is the config key
#           input.contigs_dir, and 00_common.smk discovers the sample names by
#           globbing it, so this is where a contigs run enters the workflow.
# Does:     takes one of two branches, decided by that first header:
#             Branch A, it matches 'length_<n>' or 'cov_<n>' (SPAdes style):
#               linearise (FASTA_LIN_CMD, one sequence line per record) then keep
#               contigs with coverage >= 2.0 and length >= 500 (FASTA_SEL_CMD,
#               splitting the header on "_" so field 4 is the length and field 6
#               the coverage). Headers are kept EXACTLY as they are — the full
#               SPAdes header.
#             Branch B, anything else: trim each header to its first whitespace
#               token (FASTA_HEAD_CMD) and change NOTHING else. No length filter,
#               no coverage filter, line wrapping preserved.
#           Echoes which branch ran to the log.
# Produces: contigs = DRAFT_CONTIGS, 02.assembly/{sample}/contigs_filt.fasta — the
#           draft genome the shared contamination screen reads.
# Consumed by: map_contigs, blast_contigs, blob_json and select_contigs, all in
#              shared/10_decontam.smk.
#
# BE AWARE (inherited from v1, and worth saying out loud): the rule NAME is
# misleading in branch B. Nothing is filtered there; the rule only normalises the
# headers, which every later join in the pipeline depends on (Bakta, Platon,
# geNomad and the BLAST screen all key on the first token).
#
# `|| true` after the grep is MANDATORY, not decoration: grep exits 1 when it
# matches nothing, and Snakemake runs the shell in strict mode (set -e), so the
# rule would die on every non-SPAdes input — which is precisely the case branch B
# exists to serve. v1 had it; keep it.
#
# v1 → v2 (additive): v1 echoed which branch it took to the console, where it is
# lost in a batch run. It now goes to a log file, because "SPAdes format detected"
# versus "fixing FASTA headers" decides whether any filtering happened at all, and
# this project audits that kind of decision.
#
# conda: NONE — head, grep and awk only, as in v1 FastaFlux.
#
# (v1 message: "Contig filtering.")
rule filter_contigs:
    input:
        contigs = os.path.join(CONTIGS_DIR, CONTIGS),
    output:
        contigs = DRAFT_CONTIGS,
    log:
        LOGS + "/filter_contigs_{sample}.log"
    shell:
        """
        # Look at the first header only: SPAdes writes the same style for every
        # record, so one line is enough to decide.
        header_format=$(head -n 1 {input.contigs} | grep -E 'length_[0-9]+|cov_[0-9.]+' || true)

        if [ -n "$header_format" ]; then
            echo "SPAdes-style header detected for {wildcards.sample}: removing contigs <500 bp and <2x coverage; headers kept unchanged." > {log}
            awk {FASTA_LIN_CMD:q} {input.contigs} | \
            awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}
        else
            echo "Non-SPAdes header for {wildcards.sample}: trimming headers to their first token; NO length or coverage filtering applied." > {log}
            awk {FASTA_HEAD_CMD:q} {input.contigs} > {output.contigs}
        fi
        """
