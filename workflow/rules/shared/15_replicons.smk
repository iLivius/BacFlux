# BacFlux v2.0.0 — the Bakta replicon table (long-read modes only).
#
# build_replicons — write the small table that tells Bakta which contigs are
# circular, plus the audit trail behind every call.
#
# Bakta treats every sequence it is given as a LINEAR CONTIG unless told
# otherwise, and that assumption costs genes. Pyrodigal, Bakta's gene caller,
# branches on topology: on a sequence declared circular it may call a gene that
# runs ACROSS THE ORIGIN, merging the first and last partial CDS into one feature.
# On a closed bacterial chromosome that is typically a handful of genes at
# position 1 — very often including dnaA itself, which is precisely the gene
# dnaapler rotated to the start.
#
# Long-read assemblies are the only ones where the topology is actually known:
# Flye reports, per contig, whether it closed the sequence into a circle. Hence
# the `if HAS_LONG_READS:` gate — in illumina and contigs mode no rule here is
# defined at all, and BAKTA_REPLICON_INPUT (00_common.smk) is an empty list, so
# rule annotation simply adds no --replicons flag.
#
# Takes in:
#   contigs   = FINAL_CONTIGS    the assembly Bakta will annotate. It is the JOIN
#                                ANCHOR — one output row per record, in file
#                                order — so the table always describes exactly
#                                the sequences Bakta sees.
#   flye_info = FLYE_INFO        Flye's assembly_info.txt → the topology.
#   dnaapler  = DNAAPLER_SUMMARY dnaapler's reorientation summary → the type.
# Produces:
#   BAKTA_REPLICONS       5 tab-separated fields per contig, NO header. Bakta
#                         unpacks each row into exactly five variables inside a
#                         bare try/except, so a wrong field count exits with one
#                         unhelpful line and nothing else.
#   BAKTA_REPLICONS_AUDIT one row per contig with the raw dnaapler numbers and a
#                         reason column, per the project convention that every
#                         filtering decision is auditable.
# Consumed by: rule annotation in shared/40_annotation.smk (the replicons file).
#              The audit file is terminal — nothing reads it — which is why
#              _frontend_targets_for() in 00_common.smk has to request it by name
#              or it would never be built.
#
# What goes in the table (the reasoning lives in the helper script's docstring,
# and every number behind a call is in the audit file):
#   topology — ALWAYS from Flye's "circ." column. A measurement, not a guess.
#   type     — from dnaapler ALONE, and only on a strong marker hit (Coverage
#              ≥ 80% AND Identity ≥ 40%): dnaA/cog1474 → chromosome, repA →
#              plasmid. Everything else, including terL, stays the neutral
#              "contig". Platon is deliberately NOT consulted: it runs at stage
#              06 and annotation at stage 04, so asking for its opinion would
#              invert the DAG and serialise annotation behind plasmid calling.
#
# PIN-SENSITIVE — RE-CHECK ON EVERY BAKTA VERSION BUMP. Writing the conservative
# type="contig" while keeping topology="circular" only works because Bakta
# 1.12.0's parser (bakta/utils.py) has the line that would force a "contig" row
# back to linear written as a COMPARISON instead of an assignment
# (`topology == TOPOLOGY_LINEAR`), so it does nothing. If a future Bakta fixes
# that typo, every "contig" row becomes linear and the entire benefit of this
# table disappears silently, with no error. This is on the end-to-end gate
# checklist.
#
# The script never crashes on a bad join: if the contig IDs do not line up it
# warns loudly into {log} and still writes an all-neutral table. That warning is
# the detector for the one silent failure mode here — Medaka renaming contigs in
# nanopore mode would otherwise produce a perfectly-formatted table that says
# nothing.
#
# conda: NONE — the helper is stdlib-only Python and runs in the environment
# Snakemake was launched from, exactly like select_contigs_by_taxonomy.py and
# plasmid_concordance.py. No new dependency.
#
# Why this file is shared-and-gated rather than one copy in rules/nanopore/ and
# one in rules/hybrid/: 00_common.smk pins FLYE_INFO, DNAAPLER_SUMMARY,
# FINAL_CONTIGS and the two output paths to the SAME strings in both long-read
# modes, so the rule would be byte-identical in both — exactly the duplication the
# v2 unification exists to remove. `if HAS_LONG_READS:` is already the house
# pattern (`if HAS_READS:` in shared/20_qc.smk, `if NEEDS_FINAL_BLAST:` in
# shared/10_decontam.smk, `if HAS_SHORT_READS:` in shared/50_amr.smk), and the
# Snakefile globs rules/shared/*.smk, so nothing else needs changing.
#
# THE NUMBER 15 reads as "after decontamination, before QC", which is where this
# sits in the DAG: in nanopore mode FINAL_CONTIGS only exists after stage 10.
#
# Everything referenced here comes from 00_common.smk: HAS_LONG_READS,
# FINAL_CONTIGS, FLYE_INFO, DNAAPLER_SUMMARY, BAKTA_REPLICONS,
# BAKTA_REPLICONS_AUDIT, REPLICONS_SCRIPT, LOGS.

if HAS_LONG_READS:

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
            # The "strong hit" thresholds — this project's convention, not a
            # dnaapler or Bakta default. Coverage is the alignment length as a
            # percentage of the reference protein, and it is the number that
            # actually separates a genuine full-length replication initiator
            # (observed ~100%) from a marginal fragment (observed 14–61%).
            # Identity does not separate them on its own: an 18%-coverage
            # terminase fragment scored 57.5% while a real DnaA scored 68.3%. So
            # the identity floor is there only to exclude the protein-alignment
            # twilight zone (~20–35%), and is deliberately not set near the
            # observed 68–79% — dnaapler's reference proteins routinely come from
            # another genus, and a high cutoff would demote real chromosomes.
            # The same two numbers are the defaults of classify_type() in
            # scripts/15_replicons/build_bakta_replicons.py.
            min_cov = 80.0,
            min_id = 40.0,
        log:
            LOGS + "/build_replicons_{sample}.log"
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
