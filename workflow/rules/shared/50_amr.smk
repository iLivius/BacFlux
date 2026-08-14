# Stage 05 — AMR and virulence screening of the finished genome.
#
# BacFlux looks for resistance genes by three methods that fail in different
# ways, which is the whole point: each one sees things the others cannot. Two of
# them live in this file.
#
#   * ABRicate on the ASSEMBLED CONTIGS, once per curated database. It consumes
#     nothing but FINAL_CONTIGS (the D2 hand-off), so it is identical in all four
#     modes.
#   * BBMap mapping the TRIMMED READS straight onto CARD's reference sequences,
#     gated on HAS_SHORT_READS. Assembly collapses repeats, so a resistance gene
#     present in several copies — or one sitting on a repeat-rich mobile element
#     that broke the assembly — can be under-represented or missing in the contigs
#     while sitting plainly in the reads. Mapping side-steps the assembly.
#
# The third leg, AMRFinderPlus, is not here. It runs inside the mobilome module
# (rule amrfinderplus in shared/80_mobilome.smk) and only when mobilome.run is
# true, because its job there is to feed structured AMR calls into the AMR ×
# mobile-element co-localisation — see docs/mobilome_module_SPEC.md §4.
#
# amr_contigs           : one ABRicate run per (sample, database), at the EFSA
#                         reporting thresholds. Fanned out over DATABASES by
#                         rule all, not by an expand() in this file.
# AMR_summary           : abricate --summary over one sample's eight tables.
# download_amr_db       : fetch and unpack CARD from links.card_link.
# download_amr_db_local : symlink an already-extracted CARD directory instead.
#                         Mutually exclusive with download_amr_db.
# map_amr_db            : BBMap the trimmed pairs onto CARD twice, at strict and
#                         relaxed read identity, plus v1's AMR_legend.
# card_mapping_report   : join those two passes with CARD's own aro_index.tsv.
#                         This is the file to open; the rest is the evidence.
#
# Data flow:
#
#   contigs_final.fasta ─► amr_contigs (once per database) ─► {db}.tsv  (×8)
#                                                                  │
#                                                    AMR_summary ◄──┘
#                                                        └─► AMR_summary.txt
#
#   [short-read modes only]
#   card_db/ (download_amr_db | download_amr_db_local) ─┐
#   {sample}_trim_R{1,2}.fastq ─────────────────────────┴─► map_amr_db
#     (rule trim_adapters, in whichever front end ran)         │
#                                             ┌────────────────┘
#                                             ├─► {sample}_AMR_legend.tsv
#                                             └─► {sample}_covstats.tsv
#                                                 {sample}_covstats_relaxed.tsv
#                                                       │
#                                   card_mapping_report ◄┘
#                                       └─► {sample}_CARD_report.tsv
#
# Inherited from 00_common.smk and never re-derived here: FINAL_CONTIGS, DIR_AMR,
# LOGS, DATABASES, DATABASE_PATTERN, and — for the CARD leg — CARDDB, CARD_LINK,
# CARD_TARBALL, CARD_DB_DIR, CARD_STRICT_ID, CARD_RELAXED_ID, CARD_MIN_COVERED,
# CARD_REPORT_SCRIPT, TRIM_R1, TRIM_R2, RAM, capped_cpus, HAS_SHORT_READS. The
# eight-database list lives ONCE in 00_common and drives (a) the per-db fan-out in
# rule all, (b) the {db} wildcard constraint, and (c) the summary input — one
# list, one source of truth (D7).
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/abricate.yaml" climbs shared/ → rules/ → workflow/ and lands on
# workflow/envs/abricate.yaml (the single shared env copy).
#
# The two ABRicate rules declare no `threads:`. They were single-threaded in v1
# and no thread request was added, so Snakemake counts them as one core each.


# ──────────────────── Contig screen (ABRicate) ─────────────────
# ABRicate BLASTs the finished contigs against one curated resistance / virulence
# database and reports every hit clearing the identity and coverage cutoffs.
# Different databases cover different gene sets (CARD, ResFinder, NCBI, VFDB, …),
# so the SAME contigs are screened against each one separately and the eight
# tables are kept side by side rather than merged.
#
# minid 80 / mincov 70 are the EFSA reporting thresholds — a hit must show ≥80%
# identity over ≥70% of the reference gene's length (see README, AMR section).
# They belong to this leg ONLY. The CARD read leg below filters on per-read
# identity instead, and AMRFinderPlus in shared/80_mobilome.smk applies its own
# curated per-gene cutoffs; a blanket 80/70 across all three would fight those.
#
# Takes in: contigs = FINAL_CONTIGS, the finished decontaminated assembly.
# Produces: 05.amr/abricate/{sample}/{db}.tsv — one hit table per (sample, db).
# Consumed by: AMR_summary below, which collates a sample's eight tables.
#
# The rule declares ONE templated output and contains NO expand(). The fan-out
# over the eight databases lives in _downstream_targets() in 00_common.smk, which
# asks rule all for one {db}.tsv per (sample × database); Snakemake then
# instantiates this rule once per pair, and params.db hands each instance its own
# database name. In v1 that expand() sat in the Snakefile's own rule all.
#
# (v1 message: "--- ABRicate: AMR detection ({wildcards.db}). ---")
rule amr_contigs:
    input:
        contigs = FINAL_CONTIGS,
    output:
        amr_tab = DIR_AMR + "/abricate/{sample}/{db}.tsv",
    wildcard_constraints:
        # Restrict {db} to the eight names in DATABASES. Snakemake's default
        # wildcard pattern is `.+`, which matches slashes too, so without this
        # both {sample} and {db} are free to split a path under abricate/ in more
        # than one way.
        db = DATABASE_PATTERN,
    params:
        db = lambda wc: wc.db,
        minid = 80,
        mincov = 70,
    conda:
        "../../envs/abricate.yaml"
    log:
        LOGS + "/amr_{db}_in_{sample}_contigs.log"
    priority: 4
    shell:
        # :q lets Snakemake shell-quote each value safely (paths / db name / log).
        """
        abricate \
          --db {params.db:q} \
          {input.contigs:q} \
          --minid {params.minid} \
          --mincov {params.mincov} \
          --nopath \
          --quiet > {output.amr_tab:q} 2> {log:q}
        """


# abricate --summary turns one sample's eight per-database hit tables into a
# single presence/absence matrix (gene × database) — the readable overview of the
# whole AMR / virulence screen for that isolate, and a terminal product.
#
# The doubled brace in expand(... "{{sample}}" ...) keeps `sample` a wildcard while
# expand fills in only the `db` list, so the input is exactly this sample's eight
# tables, in DATABASES order. `abricate --summary` takes them as a plain
# positional list and its column order follows the input-file order — reordering
# DATABASES in 00_common.smk therefore reorders this file's columns.
#
# Takes in: this sample's {db}.tsv tables from amr_contigs.
# Produces: 05.amr/abricate/{sample}/AMR_summary.txt.
# Consumed by: the user.
#
# D7 single-source: v1 hard-coded all eight database filenames TWICE (as eight
# named input: keys AND again in the shell). Here input: is driven by the shared
# DATABASES list and the shell passes {input} positionally, so adding or removing
# a database in 00_common.DATABASES now propagates by itself to rule all, the
# {db} wildcard constraint and this summary.
#
# v1's AMR_summary had NO log: and redirected only stdout, so ABRicate's
# summary-step stderr was lost. v2 adds a log and captures stderr, matching every
# other rule in this file. The ABRicate command itself is unchanged.
rule AMR_summary:
    input:
        # Unnamed on purpose: the shell passes the whole list positionally as
        # {input}, which is the form `abricate --summary` expects.
        expand(DIR_AMR + "/abricate/{{sample}}/{db}.tsv", db=DATABASES),
    output:
        amr_summary = DIR_AMR + "/abricate/{sample}/AMR_summary.txt",
    conda:
        "../../envs/abricate.yaml"
    log:
        LOGS + "/AMR_summary_{sample}.log"
    priority: 3
    shell:
        """
        abricate --summary {input} > {output.amr_summary} 2> {log}
        """


# ──────────────────── CARD read-mapping leg (BBMap) ────────────
# Gated on the CAPABILITY flag HAS_SHORT_READS (illumina + hybrid), not on the
# mode name (D7). _downstream_targets() in 00_common.smk requests these files
# behind the SAME flag, so rule all and the rule definitions can never disagree:
# either both exist or neither does. nanopore and contigs get no CARD leg at all
# — v1 had none for them, and neither mode's config even declares links.card_link.
# Because envs/bbmap.yaml is referenced only from inside this block, those two
# modes never build that env either.
if HAS_SHORT_READS:

    # ── Fetch and unpack the CARD reference sequences ──
    # wget the .tar.bz2 named by links.card_link, then extract it. Every sample
    # waits on this one download.
    #
    # Produces card_tarball and card_dir, both temp(): nothing downstream needs
    # the raw database, so Snakemake deletes them once every job that reads them
    # is done. Consumed by map_amr_db (the reference sequences) and by
    # card_mapping_report (aro_index.tsv), so the delete waits for both.
    #
    # The tarball is a SIBLING of the extraction directory, not a file inside it.
    # v1 declared it inside the same rule's directory() output; the house rule
    # (stated in shared/40_annotation.smk) is that temp scratch never nests inside
    # a directory() output. That is why the extraction directory has to be created
    # explicitly here — v1 got it for free as the tarball's parent.
    #
    # No conda: — wget and tar come from the launch environment, as in v1. Same
    # deferred decision as cazyme_db_download; see docs/README_notes.md item 3.
    #
    # Defined only when directories.card_db is empty, i.e. when BacFlux is the one
    # doing the downloading. Mutually exclusive with download_amr_db_local below,
    # same directory()-wipe safety reasoning as checkv_db / checkv_db_local in
    # shared/70_phage.smk.
    #
    # (v1 message: "--- Download AMR features from CARD repository. ---")
    if not CARDDB:

        rule download_amr_db:
            output:
                card_tarball = temp(CARD_TARBALL),
                card_dir = temp(directory(CARD_DB_DIR)),
            params:
                # Resolved and validated once in 00_common.smk §7: a short-read
                # config that omits links.card_link is rejected at parse time, by
                # key name, instead of dying inside wget here.
                link = CARD_LINK,
            log:
                LOGS + "/download_amr.log"
            priority: 9
            shell:
                """
                mkdir -p {output.card_dir}

                wget {params.link} -O {output.card_tarball} > {log} 2>&1
                tar -xjvf {output.card_tarball} -C {output.card_dir} >> {log} 2>&1
                """

    # ── Use a CARD database the user already holds ──
    # Defined only when directories.card_db is set, replacing download_amr_db.
    # Symlinks CARD's reference files (aro_index.tsv,
    # nucleotide_fasta_protein_homolog_model.fasta, …) into BacFlux's own
    # directory rather than pointing the rules at the user's directory directly:
    # card_dir is a directory() output and Snakemake WIPES a directory() output
    # before re-running its rule, so the user's own copy must never be the thing
    # declared. Same reasoning as checkv_db_local in shared/70_phage.smk.
    #
    # The simplest of the workflow's six *_local database rules, because CARD
    # needs no local rebuild at all: map_amr_db only ever READS these files, and
    # BBMap builds its own index per sample into a separate temp directory, never
    # into card_dir.
    if CARDDB:

        rule download_amr_db_local:
            input:
                src = CARDDB,
            output:
                card_dir = temp(directory(CARD_DB_DIR)),
            log:
                LOGS + "/download_amr_db_local.log"
            priority: 9
            shell:
                """
                mkdir -p {output.card_dir}
                {{
                  echo "Building a local CARD database view"
                  echo "  source (read-only): {input.src}"
                  echo "  view:               {output.card_dir}"
                }} > {log}
                for f in "{input.src}"/*; do
                    ln -sfn "$f" "{output.card_dir}/$(basename "$f")"
                done
                """

    # ── Map the trimmed reads onto CARD (BBMap) ──
    # Align this sample's quality-trimmed Illumina pairs against CARD's
    # protein-homolog-model nucleotide sequences and report how much of each
    # reference gene the reads actually covered. A gene is called present when a
    # large fraction of its LENGTH is covered — length coverage, not "some reads
    # hit it", which is what stops a short conserved domain producing a false
    # positive.
    #
    # Takes in:
    #   r1 / r2  = TRIM_R1 / TRIM_R2, the fastp-trimmed pairs from rule
    #              trim_adapters — illumina/10_reads.smk in illumina mode,
    #              hybrid/10_reads.smk in hybrid, the two kept identical on
    #              purpose. Both write to the same path, hence one shared
    #              constant. Their producer declares them temp(), so Snakemake
    #              keeps them until the assembler, map_contigs and this rule are
    #              all done.
    #   card_dir = the extracted CARD database from download_amr_db (or the
    #              symlinked view from download_amr_db_local).
    # Does: BBMap TWICE over the same reference — once at a strict read-identity
    #       filter (0.99, near-exact) and once relaxed (0.95) — then two
    #       post-processing steps:
    #       (1) re-sort each covstats by descending Covered_percent, header kept;
    #       (2) build the v1 human-readable legend from the STRICT pass: for every
    #           feature covered ≥70%, pull its row out of CARD's aro_index.tsv so
    #           the output names the drug class and mechanism, not an accession.
    #       card_mapping_report below then joins the two sorted covstats into one
    #       annotated table.
    # Produces:
    #   covstats         = 05.amr/mapping/{sample}/{sample}_covstats.tsv
    #   covstats_relaxed = 05.amr/mapping/{sample}/{sample}_covstats_relaxed.tsv
    #   amr_legend       = 05.amr/mapping/{sample}/{sample}_AMR_legend.tsv
    #   plus three temp() intermediates: BBMap's ref/ index and the two unsorted
    #   covstats. covstats and amr_legend are the paths _downstream_targets()
    #   requests by name; covstats_relaxed rides along because the same rule
    #   writes it.
    # Consumed by: the user (terminal AMR products, not fed into MultiQC) and by
    #              card_mapping_report below.
    #
    # Why two identity filters. 0.99 alone answers "is this exact reference allele
    # here?" and says nothing when the isolate carries a divergent member of the
    # same family. Mapping again at 0.95 gives that second answer, and reporting
    # the two side by side keeps the strict call's specificity while making a
    # divergent hit visible instead of simply absent.
    #
    # How much that actually buys, measured rather than assumed: on four real
    # environmental isolates (one Pseudomonad, three Bacilli) the relaxed pass
    # mapped ~2% more bases and promoted NO extra CARD sequence past the 70%
    # coverage threshold. So the honest claim is not "0.95 finds a lot more" — it
    # is that "nothing at 0.95 either" is now a recorded observation instead of an
    # untested assumption. Note also how far 0.95 reaches: identity is per READ, so
    # on 150 bp reads it tolerates ~7 mismatches and recovers genes down to roughly
    # 95% nucleotide identity to the CARD reference — within-family allelic
    # variation, not a distant homolog. Anything more divergent needs a
    # protein-level search, which is what the ABRicate and AMRFinderPlus legs on
    # the contigs are for.
    #
    # The second pass reuses the index the first one built (no ref= on the second
    # invocation, so BBMap loads it from path= rather than rebuilding it). Measured
    # cost of the whole rule: 9.7 s for the first pass including indexing, 8.7 s
    # for the second, on a 5 Mb genome at 16 threads.
    #
    # PRESERVED FROM v1, ALL KNOWN WARTS, DELIBERATELY NOT FIXED:
    #   * path={output.bbmap_temp} where bbmap_temp already ends in /ref, so BBMap
    #     creates ref/ref/. Odd, harmless, kept.
    #   * the legend's `grep $i {card_dir}/aro_index.tsv` is unquoted, unanchored
    #     and not -F, so an ARO accession that is a substring of another one can
    #     pull in extra rows. Kept as-is.
    #   * NO `set -euo pipefail` in this rule. The legend's `for ... grep ... done`
    #     loop legitimately exits non-zero when nothing clears 70%, and `set -e`
    #     would turn "this isolate has no strongly covered AMR gene" — a perfectly
    #     normal result — into a hard pipeline failure.
    #
    # One deliberate v1→v2 change: -Xmx32g was hard-coded, which fails outright on
    # a machine with less RAM than that. It now asks for min(RAM, 32) GB, the same
    # class of fix as Qualimap's --java-mem-size in shared/20_qc.smk (both taken
    # together, for consistency). Like Qualimap's, this is a GB figure passed
    # through to the JVM, not something the scheduler reserves.
    #
    # (v1 message: "--- Map trimmed reads against CARD db. ---")
    rule map_amr_db:
        input:
            card_dir = CARD_DB_DIR,
            r1 = TRIM_R1,
            r2 = TRIM_R2,
        output:
            bbmap_temp = temp(directory(DIR_AMR + "/mapping/{sample}/ref")),
            covstats_temp = temp(DIR_AMR + "/mapping/{sample}/{sample}_covstats_temp.tsv"),
            covstats_relaxed_temp = temp(DIR_AMR + "/mapping/{sample}/{sample}_covstats_relaxed_temp.tsv"),
            covstats = DIR_AMR + "/mapping/{sample}/{sample}_covstats.tsv",
            covstats_relaxed = DIR_AMR + "/mapping/{sample}/{sample}_covstats_relaxed.tsv",
            amr_legend = DIR_AMR + "/mapping/{sample}/{sample}_AMR_legend.tsv",
        params:
            # CARD ships several models; the protein homolog model is the one that
            # holds acquired resistance genes (not the mutation-based models).
            card_target = "nucleotide_fasta_protein_homolog_model.fasta",
            # Read-identity filters for the two passes. Kept as plain numbers here,
            # not config keys, for the same reason the 0.99 always was: they are a
            # methods decision, not a per-run knob. CARD_STRICT_ID / CARD_RELAXED_ID
            # are 00_common globals, shared with card_mapping_report so the two
            # rules cannot drift apart on which filter produced which file.
            min_id = CARD_STRICT_ID,
            min_id_relaxed = CARD_RELAXED_ID,
            max_ram = min(RAM, 32),
        conda:
            "../../envs/bbmap.yaml"
        threads: capped_cpus(24)
        log:
            LOGS + "/map_amr_{sample}.log"
        priority: 5
        shell:
            # Column 5 of BBMap's covstats is Covered_percent; the >=70 threshold
            # and the legend header text below are v1's, kept verbatim. The awk
            # braces are doubled because this awk is written directly in the shell.
            """
            bbmap.sh \
              -in={input.r1} \
              -in2={input.r2} \
              ref={input.card_dir}/{params.card_target} \
              path={output.bbmap_temp} \
              idfilter={params.min_id} \
              idtag \
              -Xmx{params.max_ram}g \
              threads={threads} \
              ambiguous=best \
              secondary=f \
              covstats={output.covstats_temp} > {log} 2>&1

            (head -n 1 {output.covstats_temp} > {output.covstats}) && \
            tail -n +2 {output.covstats_temp} | awk -F'\t' '{{print $5 "\t" $0}}' | sort -t$'\t' -k1,1nr | cut -f2- >> {output.covstats}

            # Second, relaxed pass. The missing ref= is deliberate: with only path=
            # given, BBMap loads the index the strict pass just wrote instead of
            # rebuilding it, so this costs one extra alignment pass and no extra
            # indexing. Everything else is identical to the strict pass, because the
            # two coverage figures are only comparable if nothing but the identity
            # filter changed between them.
            bbmap.sh \
              -in={input.r1} \
              -in2={input.r2} \
              path={output.bbmap_temp} \
              idfilter={params.min_id_relaxed} \
              idtag \
              -Xmx{params.max_ram}g \
              threads={threads} \
              ambiguous=best \
              secondary=f \
              covstats={output.covstats_relaxed_temp} >> {log} 2>&1

            (head -n 1 {output.covstats_relaxed_temp} > {output.covstats_relaxed}) && \
            tail -n +2 {output.covstats_relaxed_temp} | awk -F'\t' '{{print $5 "\t" $0}}' | sort -t$'\t' -k1,1nr | cut -f2- >> {output.covstats_relaxed}

            echo "#AMR features with a covered length of at least 70%" > {output.amr_legend}
            (head -n 1 {input.card_dir}/aro_index.tsv >> {output.amr_legend}) && \
            for i in $(tail -n +2 {output.covstats} | awk -F'\t' '$5 >=70' | cut -f1 | awk -F'|' '{{print $5}}'); do \
                grep $i {input.card_dir}/aro_index.tsv; \
            done >> {output.amr_legend}
            """

    # ── Join the two CARD passes into one readable table ──
    # The AMR_legend written above answers "which CARD sequences were covered?".
    # It does not say what KIND of thing each one is, and CARD's protein homolog
    # model is not a list of acquired resistance genes — it also holds the parts of
    # multi-subunit efflux pumps, the transcriptional regulators of those pumps, and
    # a handful of entries where resistance comes from the gene being ABSENT. On an
    # environmental Gram-negative, whose chromosome encodes whole RND efflux
    # repertoires as normal core machinery, a flat count of "AMR genes found" is
    # therefore badly inflated. CARD already classifies every entry in aro_index.tsv;
    # this rule joins that classification onto the hits so the inflation is visible
    # instead of silent, and pairs each hit's strict and relaxed coverage so a
    # divergent allele shows up as `divergent` rather than as nothing at all.
    #
    # Takes in:
    #   covstats         = the strict (0.99) sorted coverage table from map_amr_db.
    #   covstats_relaxed = the relaxed (0.95) one, same rule, same reference.
    #   card_dir         = the CARD database, for aro_index.tsv.
    # Does: runs card_mapping_report.py, which joins the two tables on the ARO
    #       accession carried in the BBMap defline, attaches CARD's AMR Gene Family
    #       / Drug Class / Resistance Mechanism, and labels every row with one of
    #       five categories (resistance_determinant, efflux_other, efflux_component,
    #       regulator, presence_indicates_susceptibility). The script's docstring
    #       carries the full reasoning, including what it deliberately does NOT do:
    #       it never calls a gene intrinsic or acquired, because that is a
    #       population-level judgement no single genome can support — see
    #       docs/methods_amr_intrinsic_acquired.md.
    # Produces:
    #   card_report = 05.amr/mapping/{sample}/{sample}_CARD_report.tsv
    # Consumed by: the user. This is the file to read; covstats and the legend are
    #              kept as the raw evidence behind it.
    #
    # A separate rule rather than more shell inside map_amr_db, for two reasons: the
    # bbmap conda env is Java-only and ships no Python interpreter, and keeping them
    # apart means editing the report script re-runs the report, not the mapping.
    # The platon env is reused purely because it already ships a Python and is
    # already built in every run (rule plasmid_search) — the script is stdlib-only
    # and uses nothing from Platon. Same trick, same reason, as plasmid_concordance
    # in shared/60_plasmid.smk.
    rule card_mapping_report:
        input:
            card_dir = CARD_DB_DIR,
            covstats = DIR_AMR + "/mapping/{sample}/{sample}_covstats.tsv",
            covstats_relaxed = DIR_AMR + "/mapping/{sample}/{sample}_covstats_relaxed.tsv",
        output:
            card_report = DIR_AMR + "/mapping/{sample}/{sample}_CARD_report.tsv",
        params:
            report_script = CARD_REPORT_SCRIPT,
            # Bare integers, used only to name the report's two coverage columns
            # (covered_percent_id99 / covered_percent_id95) so each figure says
            # which filter produced it. Derived from the same globals map_amr_db
            # filters on, so renaming a filter renames its column with it.
            id_label = round(CARD_STRICT_ID * 100),
            id_label_relaxed = round(CARD_RELAXED_ID * 100),
            # Minimum covered length for a CARD sequence to reach the report. The
            # same 70% the v1 legend uses, so the two files agree on what counts
            # as present.
            min_covered = CARD_MIN_COVERED,
        conda:
            "../../envs/platon.yaml"
        log:
            LOGS + "/card_mapping_report_{sample}.log"
        priority: 5
        shell:
            """
            python {params.report_script} \
              --sample {wildcards.sample} \
              --covstats-strict {input.covstats} \
              --covstats-relaxed {input.covstats_relaxed} \
              --aro-index {input.card_dir}/aro_index.tsv \
              --strict-id {params.id_label} \
              --relaxed-id {params.id_label_relaxed} \
              --min-covered {params.min_covered} \
              --out {output.card_report} > {log} 2>&1
            """
