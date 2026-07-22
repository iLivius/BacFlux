# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Stage 06 plasmid module (rules/shared/60_plasmid.smk)  [D9]
#
# This module decides which contigs are plasmids and how confident that call is.
# Like the rest of the shared tail it consumes only FINAL_CONTIGS (D2), so it is
# identical in all four modes.
#
# The plasmid deliverable depends on whether the user opted in to geNomad
# (config.phage.caller: genomad — see PHAGE_CALLER in 00_common.smk):
#   - DEFAULT (VirSorter2 phage caller, geNomad OFF): Platon alone. The terminal
#     deliverable is Platon's verified_plasmids.txt (v1 behaviour). geNomad does
#     not run, so there is no second opinion to concord with.
#   - geNomad OPT-IN: a Platon + geNomad CONCORDANCE table (D9). geNomad — a
#     completely different, gene-content method — gives a second opinion, and the
#     two are joined into a per-contig confidence-tiered table. geNomad is run ONCE
#     in 70_phage.smk (genomad_end_to_end); this module only CONSUMES its plasmid
#     summary. (geNomad is academic/non-commercial-licensed, which is exactly why
#     the concordance is opt-in and not the default — see 00_common.smk §1b.)
#
# Platon always runs and always keeps v1's supplementary "does the nt BLAST hit say
# plasmid?" check as one visible, non-authoritative annotation (D9) — never the
# decision. When the concordance runs, that check rides along as one column.
#
# Data flow (top to bottom):
#
#   contigs_final.fasta ─► plasmid_search (Platon) ─► 06.plasmids/{sample}/platon/
#          │                                            (contigs_final.tsv,
#          │                                             contigs_final.chromosome.fasta,
#          │                                             verified_plasmids.txt ← default deliverable)
#          │                                                     │
#          │   (only if geNomad opted in)                        │
#   07.phages/genomad/{sample}/…_plasmid_summary.tsv ───┐        │
#          (from genomad_end_to_end in 70_phage.smk)     ▼        ▼
#                                       plasmid_concordance (plasmid_concordance.py)
#                                       ─► 06.plasmids/{sample}/{sample}_plasmid_concordance.tsv
#
# CROSS-STAGE NOTE (concordance path only): the geNomad plasmid summary this module
# reads physically lives under 07.phages/. That backwards-numbered 06←07 edge is
# safe — stage numbers are organisational only (D1 groups by tool family);
# Snakemake orders work by the input/output DAG, not by directory number.
#
# Everything referenced here is defined once in 00_common.smk (never re-derived):
# FINAL_CONTIGS, DIR_PLASMIDS, PLATONDB, BLASTOUT, PLATON_DIR, PLATON_VERIFIED,
# GENOMAD_DIR, GENOMAD_PREFIX, PLASMID_CONCORDANCE, PLASMID_CONCORDANCE_SCRIPT,
# PHAGE_CALLER, LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ -> workflow/envs/x.yaml.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: plasmid_search — primary plasmid calling + kept v1 check (Platon) ───
# ALWAYS runs (both the default and the geNomad-opt-in paths need Platon).
#
# Biology: Platon classifies each assembly contig as plasmid or chromosome from
# replicon-distribution scores. BacFlux then keeps v1's supplementary check: for
# each Platon-plasmid contig, grep the general contamination-screen BLAST hits for
# the word "plasmid" in the subject title. That check is a WEAK, non-authoritative
# signal (D9) — it stays as one visible annotation, never the decision.
#
# Takes in:
#   contigs = FINAL_CONTIGS — the finished, decontaminated assembly (D2). (v1 used
#             contigs_sel.fasta; the new input is contigs_final.fasta, so Platon
#             now names its outputs contigs_final.* instead of contigs_sel.*.)
#   blast   = BLASTOUT — the contamination-screen BLAST table. CROSS-STAGE edge:
#             PRODUCED by blast_contigs in the future shared/10_decontam.smk
#             (Stage 3). Both sides use the one BLASTOUT constant from 00_common so
#             they cannot drift (same pattern as COMPOSITION). CONTRACT for Stage 3:
#             the table must carry subject titles (BLAST outfmt 6 with stitle last)
#             or the `grep -q "plasmid"` check below has nothing to match.
# Does: run Platon over the contigs, then run the two-phase BLAST-text check
#       verbatim from v1 (only the input filename changes contigs_sel→contigs_final).
# Produces:
#   platon_dir = 06.plasmids/{sample}/platon/ (a DIRECTORY). v1 wrote Platon output
#                straight into {sample}/; moving it into platon/ lets the
#                concordance TSV sit as a clean SIBLING, so no second rule writes
#                inside this rule's directory() output.
#   plasmids   = verified_plasmids.txt — the kept v1 check. On the DEFAULT path this
#                IS the terminal plasmid deliverable (rule-all requests it); on the
#                geNomad path it becomes an intermediate feeding the concordance.
# Consumed by: plasmid_concordance when geNomad is opted in (reaches into platon_dir
#              for contigs_final.tsv, contigs_final.chromosome.fasta, and
#              verified_plasmids.txt); otherwise the user directly.
#
# (v1 message: "--- Platon: Plasmid identification. ---")
rule plasmid_search:
    input:
        contigs = FINAL_CONTIGS,
        blast = BLASTOUT,
    output:
        platon_dir = directory(PLATON_DIR),
        plasmids = PLATON_VERIFIED,
    params:
        platon_db = PLATONDB,
        # Platon names every output file after the INPUT basename. Our input is
        # FINAL_CONTIGS (contigs_final.fasta), so Platon writes contigs_final.*.
        # GENOMAD_PREFIX holds that shared basename ("contigs_final") — the same
        # constant single-sources the geNomad output names — so this shell and the
        # concordance rule cannot drift onto a different stem.
        contigs_prefix = GENOMAD_PREFIX,
    conda:
        "../../envs/platon.yaml"
    resources:
        cpus = capped_cpus(24)
    log:
        LOGS + "/plasmid_search_{sample}.log"
    priority: 4
    shell:
        # The if/while block is byte-for-byte v1, with only the plasmid FASTA
        # filename swapped to {params.contigs_prefix}. $i is each plasmid contig's
        # header (the contig ID, since front-end headers are trimmed to one token);
        # `grep -m 1 "$i"` finds its first BLAST line, and `grep -q "plasmid"` tests
        # whether that hit's subject title mentions a plasmid.
        """
        platon \
          --db {params.platon_db} \
          --output {output.platon_dir} \
          --verbose \
          --threads {resources.cpus} \
          {input.contigs} > {log} 2>&1

        if [[ -s {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta ]] && grep -q ">" {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta; then
            while IFS= read -r i; do
                if grep -m 1 "$i" {input.blast} | grep -q "plasmid"; then
                    echo "{wildcards.sample}: $i is a plasmid." >> {output.plasmids}
                else
                    echo "{wildcards.sample}: $i was not verified by BLAST search." >> {output.plasmids}
                fi
            done < <(grep ">" {output.platon_dir}/{params.contigs_prefix}.plasmid.fasta | sed 's/^>//g')
        else
            echo "Platon found no plasmid in sample {wildcards.sample}." > {output.plasmids}
        fi
        """


# ── geNomad concordance (only defined when PHAGE_CALLER == "genomad") ─────────
# Defined behind the same guard as the geNomad rules in 70_phage.smk, so it only
# exists when geNomad was opted in (and therefore genomad_end_to_end actually
# produced a plasmid summary to concord with). On the default path this rule does
# not exist, and rule-all's terminal plasmid target is verified_plasmids.txt above.
if PHAGE_CALLER == "genomad":

    # ── Rule: plasmid_concordance — join Platon + geNomad plasmid calls (D9) ──
    # Biology: build the per-contig concordance table. Confidence is driven by
    # whether the two INDEPENDENT callers agree (both plasmid → high; one only →
    # medium; disagreement → low, flagged not discarded). The kept v1 BLAST-text
    # state rides along as one visible, clearly-supplementary column. The file
    # itself is the audit trail, in the spirit of contig_taxonomy_decisions.tsv.
    #
    # Takes in (both {sample} DIRECTORIES; the shell reaches inside for files):
    #   platon_dir  = PLATON_DIR — for contigs_final.tsv (plasmid calls + RDS),
    #                 contigs_final.chromosome.fasta (chromosome IDs), and
    #                 verified_plasmids.txt (kept BLAST-text check).
    #   genomad_dir = GENOMAD_DIR — for the plasmid summary produced by
    #                 genomad_end_to_end in 70_phage.smk.
    # Does: run plasmid_concordance.py (stdlib-only Python — no new dependency, so
    #       it reuses Platon's pinned Python via the platon env).
    # Produces: 06.plasmids/{sample}/{sample}_plasmid_concordance.tsv.
    # Consumed by: the report / the user (the terminal plasmid product on this path;
    #              requested by rule all, transitively pulling in both callers).
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
        priority: 3
        shell:
            # {wildcards.sample} is used for --sample so the name is substituted
            # reliably. PLASMID_CONCORDANCE_SCRIPT is a 00_common global.
            """
            python {PLASMID_CONCORDANCE_SCRIPT} \
              --sample {wildcards.sample} \
              --platon-tsv {input.platon_dir}/{params.contigs_prefix}.tsv \
              --platon-chromosome {input.platon_dir}/{params.contigs_prefix}.chromosome.fasta \
              --verified-plasmids {input.platon_dir}/verified_plasmids.txt \
              --genomad-plasmid-summary {input.genomad_dir}/{params.genomad_summary_rel} \
              --output {output.concordance} > {log} 2>&1
            """
