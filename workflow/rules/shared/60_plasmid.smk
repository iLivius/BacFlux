# BacFlux v2.0.0 — plasmid calling on the finished genome (stage 06.plasmids).
# Platon is the primary caller; geNomad is an optional second opinion. Like the
# rest of the shared tail this module reads only FINAL_CONTIGS (the D2 hand-off),
# so it behaves identically in all four modes.
#
# Which file is the terminal deliverable depends on config phage.caller, because
# geNomad calls viruses and plasmids in ONE run and only runs when it was opted in
# as the phage caller (shared/70_phage.smk):
#   virsorter2 (default) — Platon alone, and Platon's verified_plasmids.txt is the
#                          deliverable, exactly as in v1.
#   genomad    (opt-in)  — Platon and geNomad are joined into a per-contig
#                          concordance table (decision D9, see
#                          docs/unification_migration_plan.md).
#
# geNomad is opt-in and never the default because Berkeley Lab licenses it for
# academic / non-commercial use only, and BacFlux is MIT — so the concordance is
# opt-in too (00_common.smk §1b, decision D8).
# Platon runs either way, and either way v1's supplementary "does the nt BLAST hit
# call this contig a plasmid?" check rides along as one visible annotation — never
# as the decision.
#
#   contigs_final.fasta ─► plasmid_search (Platon) ─► 06.plasmids/{sample}/platon/
#                                                    verified_plasmids.txt ──┐
#   07.phages/genomad/{sample}/…_plasmid_summary.tsv ─────┐                  │
#          (only when geNomad was opted in)               ▼                  ▼
#                                   plasmid_concordance (plasmid_concordance.py)
#                         ─► 06.plasmids/{sample}/{sample}_plasmid_concordance.tsv
#
# plasmid_search      : Platon over the finished contigs, plus the kept v1
#                       BLAST-text check. Always runs.
# plasmid_concordance : joins Platon's and geNomad's calls into one
#                       confidence-tiered table. Defined only when geNomad was
#                       opted in as the phage caller.
#
# The geNomad summary the concordance reads physically lives under 07.phages/.
# Reading "backwards" from stage 06 into stage 07 is safe: the stage numbers group
# tools for the reader (D1), and Snakemake orders work by the input/output DAG,
# not by directory name.
#
# Defined once in 00_common.smk and never re-derived here: FINAL_CONTIGS,
# DIR_PLASMIDS, PLATONDB, PLASMID_BLASTOUT, PLATON_DIR, PLATON_VERIFIED,
# GENOMAD_DIR, GENOMAD_PREFIX, PLASMID_CONCORDANCE, PLASMID_CONCORDANCE_SCRIPT,
# PHAGE_CALLER, LOGS, capped_cpus.
#
# conda: env paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/x.yaml.


# ────────────────── Plasmid calling (Platon) ───────────────────
# Platon (v1.7, pinned in envs/platon.yaml) sorts every contig of the finished
# assembly into plasmid or chromosome from replicon-distribution scores — protein
# families weighted by how often they sit on a plasmid rather than a chromosome.
# Always runs: the default path and the geNomad path both need its calls.
#
# On top of that call the rule keeps v1's supplementary check: look each
# Platon-plasmid contig up in the contamination-screen BLAST hits and see whether
# the subject title contains the word "plasmid". It stays an annotation and never
# a filter, because a mobile element carried by a genuinely chromosomal contig can
# push Platon's score AND that contig's best nt hit the same wrong way — so
# agreement here is weaker evidence than it looks (D9).
#
# Platon's own output lands in a platon/ SUB-directory of the sample directory (v1
# wrote it straight into {sample}/) so the concordance TSV can sit beside it as a
# sibling. Otherwise a second rule would be writing inside this rule's directory()
# output, which Snakemake forbids.
#
# verified_plasmids.txt is the terminal plasmid product on the default path, and an
# intermediate feeding plasmid_concordance when geNomad is opted in.
# (v1 announced this stage as "--- Platon: Plasmid identification. ---".)
rule plasmid_search:
    input:
        contigs = FINAL_CONTIGS,
        # The BLAST table is searched BY CONTIG ID (`grep -m 1 -F` below), so it
        # only works if it was computed over the same contigs Platon reported on.
        # Which table that is depends on the mode, and PLASMID_BLASTOUT resolves it
        # once in 00_common.smk: in illumina and contigs it is the decontamination
        # screen's own BLASTOUT (rule blast_contigs, shared/10_decontam.smk),
        # because the screen already ran on this contig set; in nanopore and hybrid
        # it is a second blastn over FINAL_CONTIGS (rule blast_final_contigs),
        # because nothing guarantees Medaka preserves headers and the hybrid screen
        # runs on the Illumina draft, whose SPAdes NODE_… names can never match
        # Flye's contig_… names. Feed the wrong table in and every plasmid comes
        # back "not verified by BLAST search" — silently, with no error.
        # CONTRACT: BLAST outfmt 6 with the subject title (stitle) LAST, or the
        # `grep -qi "plasmid"` below has nothing to match.
        blast = PLASMID_BLASTOUT,
    output:
        platon_dir = directory(PLATON_DIR),
        plasmids = PLATON_VERIFIED,
    params:
        platon_db = PLATONDB,
        # Platon names every output file after the INPUT basename. The input is
        # FINAL_CONTIGS (contigs_final.fasta), so Platon writes contigs_final.*
        # where v1's contigs_sel.fasta gave contigs_sel.*. GENOMAD_PREFIX holds
        # that basename, and geNomad names its files the same way, so one constant
        # covers both tools and this shell cannot drift onto a different stem.
        contigs_prefix = GENOMAD_PREFIX,
    conda:
        "../../envs/platon.yaml"
    threads: capped_cpus(24)
    log:
        LOGS + "/plasmid_search_{sample}.log"
    shell:
        # The BLAST-text check is v1's, unified onto the HARDENED form the two
        # long-read v1 workflows used (BacFluxL / BacFluxL+), because the plain
        # illumina form silently mis-reports in several real situations:
        #   * awk '/^>/{{sub(/^>/,""); print $1}}' takes the contig ID as the FIRST
        #     TOKEN of the header. The illumina form (grep ">" | sed 's/^>//g')
        #     keeps the WHOLE header line; Flye/dnaapler/Medaka/Polypolish headers
        #     can carry a description, and an ID with spaces can never match a
        #     tab-separated qseqid field — so every plasmid would come back
        #     "not verified" with no error.
        #   * grep -F  — fixed-string, so a contig ID containing regex
        #     metacharacters (SPAdes cov_12.34 style) cannot false-match.
        #   * grep -qi — NCBI subject titles routinely capitalise ("… Plasmid
        #     unnamed1"), which a case-sensitive grep silently misses.
        #   * `: > {output.plasmids}` truncates first, so a re-run does not append
        #     to a stale verified_plasmids.txt.
        # Platon itself runs with errors trapped (set +e … platon_rc): a non-zero
        # exit writes one explanatory line into verified_plasmids.txt instead of
        # killing the sample (v1 BacFluxL+ behaviour). That line is also how
        # plasmid_concordance.py tells "Platon crashed" from "Platon found nothing".
        """
        set +e
        platon \
          --db {params.platon_db} \
          --output {output.platon_dir} \
          --verbose \
          --threads {threads} \
          {input.contigs} > {log} 2>&1
        platon_rc=$?
        set -e

        : > {output.plasmids}

        if [ "$platon_rc" -ne 0 ]; then
            echo "{wildcards.sample}: Platon exited with status $platon_rc; see the log." >> {output.plasmids}
        elif [[ -s {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta ]] && grep -q ">" {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta; then
            while IFS= read -r i; do
                if grep -m 1 -F "$i" {input.blast} | grep -qi "plasmid"; then
                    echo "{wildcards.sample}: $i is a plasmid." >> {output.plasmids}
                else
                    echo "{wildcards.sample}: $i was not verified by BLAST search." >> {output.plasmids}
                fi
            done < <(awk '/^>/{{sub(/^>/,""); print $1}}' {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta)
        else
            echo "Platon found no plasmid in sample {wildcards.sample}." >> {output.plasmids}
        fi
        """


# ──────────────── Platon + geNomad concordance ─────────────────
# Defined only when geNomad was opted in as the phage caller (config
# phage.caller: genomad) — the same guard the geNomad rules in shared/70_phage.smk
# sit behind, so this rule exists exactly when genomad_end_to_end has a plasmid
# summary to compare against. On the default path the rule is never defined and
# rule all asks for verified_plasmids.txt instead (see 00_common.smk).
if PHAGE_CALLER == "genomad":

    # ── Join the two callers into one per-contig table (D9) ──
    # The two callers fail differently — Platon scores protein families by how
    # plasmid-like they are, geNomad classifies from gene content with its own
    # marker set — so genuine agreement between them is far stronger evidence than
    # Platon agreeing with a screening BLAST that was run for another purpose.
    # Agreement is what drives the confidence column: both callers → high, one
    # only → medium, a real clash (Platon says chromosome, geNomad says plasmid)
    # → low. Nothing is dropped; a disagreement is flagged and kept, so the table
    # is its own audit trail, in the spirit of contig_taxonomy_decisions.tsv. The
    # kept v1 BLAST-text state rides along as one clearly supplementary column.
    #
    # Both inputs are per-sample DIRECTORIES and the shell reaches inside them for
    # the individual files, which keeps those inner files off the DAG (the same
    # arrangement shared/40_annotation.smk uses for Bakta's output directory). Four
    # files are read, and each carries one part of the picture:
    #   contigs_final.tsv               — Platon's plasmid calls and their RDS score
    #   contigs_final.chromosome.fasta  — the contig IDs Platon called chromosome
    #   verified_plasmids.txt           — the kept v1 BLAST-text state per contig
    #   …_plasmid_summary.tsv           — geNomad's calls, score and FDR, from
    #                                     genomad_end_to_end in shared/70_phage.smk
    # The join, the column contract and the tie rules live in
    # scripts/60_plasmid/plasmid_concordance.py — stdlib-only Python, so it borrows Platon's
    # pinned interpreter through the platon env instead of adding a dependency.
    #
    # Produces 06.plasmids/{sample}/{sample}_plasmid_concordance.tsv, the terminal
    # plasmid deliverable on this path — rule all asks for this file, and asking for
    # it is what pulls both callers into the run.
    rule plasmid_concordance:
        input:
            platon_dir = PLATON_DIR,
            genomad_dir = GENOMAD_DIR,
        output:
            concordance = PLASMID_CONCORDANCE,
        params:
            # Shared "contigs_final" stem (see plasmid_search) naming both Platon's
            # and geNomad's per-sample output files.
            contigs_prefix = GENOMAD_PREFIX,
            # Fixed relative path of geNomad's plasmid summary inside genomad_dir.
            # NON-templated (no {sample}), so it is safe as a plain params string.
            genomad_summary_rel = GENOMAD_PREFIX + "_summary/" + GENOMAD_PREFIX + "_plasmid_summary.tsv",
        conda:
            "../../envs/platon.yaml"
        log:
            LOGS + "/plasmid_concordance_{sample}.log"
        shell:
            # PLASMID_CONCORDANCE_SCRIPT is a plain 00_common.smk global, which
            # Snakemake substitutes into the shell string like any other name, so
            # the script path needs no params entry. The sample name comes from
            # {wildcards.sample} rather than being parsed out of a path.
            """
            python {PLASMID_CONCORDANCE_SCRIPT} \
              --sample {wildcards.sample} \
              --platon-tsv {input.platon_dir}/{params.contigs_prefix}.tsv \
              --platon-chromosome {input.platon_dir}/{params.contigs_prefix}.chromosome.fasta \
              --verified-plasmids {input.platon_dir}/verified_plasmids.txt \
              --genomad-plasmid-summary {input.genomad_dir}/{params.genomad_summary_rel} \
              --output {output.concordance} > {log} 2>&1
            """
