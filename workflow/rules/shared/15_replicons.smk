# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — Bakta replicon table (rules/shared/15_replicons.smk)
#
# WHAT THIS IS FOR. Bakta assumes every sequence it is given is a LINEAR CONTIG
# unless told otherwise, and that assumption costs genes. Pyrodigal (Bakta's gene
# caller) branches on topology: on a sequence declared circular it may call a gene
# that runs ACROSS THE ORIGIN and it merges the first and last partial CDS into
# one feature. On a closed bacterial chromosome that is typically a handful of
# genes at position 1 — very often including dnaA itself, which is precisely the
# gene dnaapler rotated to the start.
#
# Long-read assemblies are the ONLY ones where we actually know the topology:
# Flye reports, per contig, whether it closed the sequence into a circle. So this
# module exists in the long-read modes only, and it is gated accordingly.
#
# Data flow:
#
#   FINAL_CONTIGS ────┐   (the join anchor: one row per record, in file order)
#   FLYE_INFO ────────┼──► build_replicons ──┬──► BAKTA_REPLICONS ─► rule annotation
#   (topology)        │   (helper script)    │                       (--replicons)
#   DNAAPLER_SUMMARY ─┘                      └──► BAKTA_REPLICONS_AUDIT
#   (replicon type)                               (the per-contig reasons)
#
# WHY IT IS SHARED-AND-GATED rather than one copy in rules/nanopore/ and one in
# rules/hybrid/: 00_common pins FLYE_INFO, DNAAPLER_SUMMARY, FINAL_CONTIGS and the
# two output paths to the SAME strings in both long-read modes, so the rule would
# be byte-identical in both — exactly the duplication this migration exists to
# remove. `if HAS_LONG_READS:` is already the house pattern (`if HAS_READS:` in
# 20_qc, `if NEEDS_FINAL_BLAST:` in 10_decontam, `if HAS_SHORT_READS:` in 50_amr),
# and Snakefile globs rules/shared/*.smk, so nothing else needs changing.
#
# THE NUMBER 15 reads as "after decontamination, before QC", which is where this
# sits in the DAG: in nanopore mode FINAL_CONTIGS only exists after stage 10.
#
# WHAT GOES IN THE TABLE (decided design; the reasoning lives in the helper
# script's docstring, and every number is in the audit file):
#   topology — ALWAYS from Flye's "circ." column. A measurement, not a guess.
#   type     — from dnaapler ALONE, and only on a strong marker hit
#              (Coverage >= 80% AND Identity >= 40%): dnaA/cog1474 -> chromosome,
#              repA -> plasmid. Everything else, including terL, stays the
#              neutral "contig". Platon is deliberately NOT consulted: it runs at
#              stage 06 and annotation at stage 04, so asking for its opinion
#              would invert the DAG and serialise annotation behind plasmid
#              calling.
#
# ⚠ PIN-SENSITIVE — RE-CHECK ON EVERY BAKTA VERSION BUMP. Writing the conservative
# type="contig" while keeping topology="circular" only works because Bakta
# 1.12.0's parser (bakta/utils.py) has the line that would force a "contig" row
# back to linear written as a COMPARISON instead of an assignment
# (`topology == TOPOLOGY_LINEAR`), so it does nothing. If a future Bakta fixes
# that typo, every "contig" row becomes linear and the entire benefit of this
# table disappears silently, with no error. This is on the end-to-end gate
# checklist.
#
# Everything referenced here comes from 00_common.smk: HAS_LONG_READS,
# FINAL_CONTIGS, FLYE_INFO, DNAAPLER_SUMMARY, BAKTA_REPLICONS,
# BAKTA_REPLICONS_AUDIT, REPLICONS_SCRIPT, LOGS.
# ─────────────────────────────────────────────────────────────────────────────

if HAS_LONG_READS:

    # ── Rule: build_replicons — write the Bakta --replicons table + its audit ─
    # Takes in:
    #   contigs   = FINAL_CONTIGS    the assembly Bakta will annotate. It is the
    #                                JOIN ANCHOR, so the table always describes
    #                                exactly the sequences Bakta sees.
    #   flye_info = FLYE_INFO        Flye's assembly_info.txt -> the topology.
    #   dnaapler  = DNAAPLER_SUMMARY dnaapler's reorientation summary -> the type.
    # Produces:
    #   BAKTA_REPLICONS       5 tab-separated fields per contig, NO header. Bakta
    #                         unpacks each row into exactly five variables inside
    #                         a bare try/except, so a wrong field count exits with
    #                         one unhelpful line and nothing else.
    #   BAKTA_REPLICONS_AUDIT one row per contig with the raw dnaapler numbers and
    #                         a reason column, per the project convention that
    #                         every filtering decision is auditable.
    # Consumed by: rule annotation in shared/40_annotation.smk (the replicons
    #              file); the audit file is terminal, which is why
    #              _frontend_targets_for() in 00_common has to request it
    #              explicitly or it would never be built.
    #
    # conda: NONE — the helper is stdlib-only Python and runs in the environment
    # Snakemake was launched from, exactly like select_contigs_by_taxonomy.py and
    # plasmid_concordance.py. No new dependency.
    #
    # The script never crashes on a bad join: if the contig IDs do not line up it
    # warns loudly (into {log}) and still writes an all-neutral table. That
    # warning is the detector for the one silent failure mode here — Medaka
    # renaming contigs in nanopore mode would otherwise produce a
    # perfectly-formatted table that says nothing.
    rule build_replicons:
        input:
            contigs = FINAL_CONTIGS,
            flye_info = FLYE_INFO,
            dnaapler = DNAAPLER_SUMMARY,
        output:
            replicons = BAKTA_REPLICONS,
            audit = BAKTA_REPLICONS_AUDIT,
        params:
            script = REPLICONS_SCRIPT,
            # The "strong hit" thresholds. Coverage is the alignment length as a
            # percentage of the reference protein and is what actually separates a
            # genuine full-length replication initiator (observed ~100%) from a
            # marginal fragment (observed 14-61%). The identity floor only excludes
            # the protein-alignment twilight zone (~20-35%); it is deliberately not
            # set near the observed 68-79%, because dnaapler's reference proteins
            # routinely come from another genus.
            min_cov = 80.0,
            min_id = 40.0,
        log:
            LOGS + "/build_replicons_{sample}.log"
        priority: 6
        shell:
            """
            python {params.script:q} \
              --sample {wildcards.sample:q} \
              --contigs {input.contigs:q} \
              --flye-info {input.flye_info:q} \
              --dnaapler-summary {input.dnaapler:q} \
              --out-replicons {output.replicons:q} \
              --out-audit {output.audit:q} \
              --min-coverage {params.min_cov} \
              --min-identity {params.min_id} > {log} 2>&1
            """
