# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — contigs front end (rules/contigs/20_assembly.smk)
#
# This mode takes genomes that are ALREADY assembled — the user hands us FASTA
# files, not reads — so there is no read QC and no assembler. The whole front end
# is one rule whose job is to hand the shared decontamination screen a FASTA it
# can work with.
#
# Data flow through this module:
#
#   {sample}.fasta ──► filter_contigs ──► DRAFT_CONTIGS
#   (config input.contigs_dir)                 │
#                                              ▼
#                     shared/10_decontam.smk: map_contigs (self-mapping) +
#                     blast_contigs -> BlobTools -> select_contigs
#                                              │
#                                              ▼
#                     DECONTAM_CONTIGS, which in this mode IS FINAL_CONTIGS
#
# WHERE THIS MODE ENDS: as in illumina mode, decontamination is the last assembly
# step, so 00_common sets DECONTAM_CONTIGS = FINAL_CONTIGS (the same Python
# object) and the shared `select_contigs` rule writes the delivered genome. This
# front end must NOT declare FINAL_CONTIGS.
#
# OWNERSHIP NOTE: the self-mapping `map_contigs` (minimap2 of the contigs against
# themselves, which exists only to give BlobTools a BAM) is NOT defined here. It
# is owned by shared/10_decontam.smk, in the `else:` branch of its three-way
# parse-time split. Do not add a second copy.
#
# Everything referenced here comes from 00_common.smk: CONTIGS_DIR, CONTIGS,
# DRAFT_CONTIGS, FASTA_LIN_CMD, FASTA_SEL_CMD, FASTA_HEAD_CMD, LOGS.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: filter_contigs — make an arbitrary input FASTA usable, two ways ────
# The problem this solves: a user-supplied assembly can come from anywhere. If it
# came from SPAdes, its headers carry the length and the coverage of every contig
# and we can apply the same quality cutoffs the illumina mode uses. If it came
# from NCBI, from a different assembler, or from a colleague, those numbers do
# not exist, and inventing a filter would silently delete real sequence.
#
# So the rule LOOKS at the first header and branches:
#
#   Branch A — the header matches 'length_<n>' or 'cov_<n>' (SPAdes style):
#              linearise (FASTA_LIN_CMD) then keep contigs with coverage >= 2.0
#              and length >= 500 (FASTA_SEL_CMD, splitting the header on "_" so
#              field 4 is the length and field 6 the coverage).
#              Headers are kept EXACTLY as they are — the full SPAdes header.
#
#   Branch B — anything else: trim each header to its first whitespace token
#              (FASTA_HEAD_CMD) and change NOTHING else. No length filter, no
#              coverage filter, line wrapping preserved.
#
# BE AWARE (inherited from v1, and worth saying out loud): the rule NAME is
# misleading in branch B. Nothing is filtered there; the rule only normalises the
# headers, which every later join in the pipeline depends on (Bakta, Platon,
# geNomad and the BLAST screen all key on the first token).
#
# Takes in: the user's FASTA from config input.contigs_dir.
# Produces: DRAFT_CONTIGS.
# Consumed by: map_contigs, blast_contigs, blob_json and select_contigs, all in
#              shared/10_decontam.smk.
#
# `|| true` after the grep is MANDATORY, not decoration: grep exits 1 when it
# matches nothing, and Snakemake runs the shell in strict mode (set -e), so the
# rule would die on every non-SPAdes input — which is precisely the case branch B
# exists to serve. v1 had it; keep it.
#
# conda: NONE — head, grep and awk only, as in v1 FastaFlux.
#
# v1 -> v2 (additive): v1 echoed which branch it took to the console, where it is
# lost in a batch run. It now goes to a log file, because "SPAdes format detected"
# versus "fixing FASTA headers" decides whether any filtering happened at all, and
# this project audits that kind of decision.
#
# (v1 message: "Contig filtering.")
rule filter_contigs:
    input:
        contigs = os.path.join(CONTIGS_DIR, CONTIGS),
    output:
        contigs = DRAFT_CONTIGS,
    log:
        LOGS + "/filter_contigs_{sample}.log"
    priority: 9
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
