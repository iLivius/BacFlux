# ─────────────────────────────────────────────────────────────────────────────
# BacFlux v2.0.0 — nanopore front end, assembly (rules/nanopore/20_assembly.smk)
#
# The biology: ONT reads are long enough to span the repeats that break a
# short-read assembly, so Flye usually closes a bacterial chromosome into a
# single circular contig. A circular sequence has no natural start, though, and
# Flye's arbitrary start point makes two assemblies of the same isolate look
# different. dnaapler rotates each replicon so it begins at a conventional
# landmark — dnaA for a chromosome, repA for a plasmid, terL for a phage — which
# makes assemblies comparable and puts the origin of replication where a reader
# expects it.
#
# Data flow through this module:
#
#   FILT_LONG ──► ont_assembly (Flye) ──► FLYE_CONTIGS
#   (10_reads)          │                      │
#                       ├──► FLYE_INFO ────────┼──► (awk) FLYE_IGNORE_LIST
#                       │    (per-contig            (contigs Flye did NOT call
#                       │     circularity)           circular; dnaapler skips them)
#                       │                      │
#                       │                      ▼
#                       │                  fix_start (dnaapler)
#                       │                      │
#                       │        ┌─────────────┼──────────────┐
#                       │        ▼             ▼              ▼
#                       │  DNAAPLER_      DNAAPLER_FIXED   DNAAPLER_SUMMARY
#                       │  REORIENTED     = DRAFT_CONTIGS  (which marker was
#                       │  (raw dnaapler)       │           found per contig)
#                       │                       ▼                │
#                       │        shared/10_decontam.smk screen   │
#                       │        -> DECONTAM_CONTIGS -> Medaka   │
#                       │                                        │
#                       └────────────────────────────────────────┴──►
#                                        shared/15_replicons.smk (build_replicons)
#                                        -> the Bakta --replicons table
#
# D3 ORDERING, preserved from v1: dnaapler runs BEFORE the contamination screen
# and Medaka runs AFTER it. That is why DRAFT_CONTIGS (what the screen reads) is
# the dnaapler output, and why 30_polish.smk polishes DECONTAM_CONTIGS.
#
# Everything referenced here comes from 00_common.smk: FILT_LONG, FLYE_DIR,
# FLYE_CONTIGS, FLYE_INFO, FLYE_IGNORE_LIST, FLYE_INPUT_MODE, DNAAPLER_DIR,
# DNAAPLER_REORIENTED, DNAAPLER_FIXED (== DRAFT_CONTIGS), DNAAPLER_SUMMARY,
# IGNORE_LIST_CMD, FASTA_HEAD_CMD, FASTA_LIN_CMD, LOGS, CPUS, capped_cpus.
# ─────────────────────────────────────────────────────────────────────────────


# ── Rule: ont_assembly — long-read de novo assembly (Flye) ───────────────────
# Takes in: FILT_LONG, the filtlong-selected reads.
# Does:     one Flye run, then one awk over Flye's own summary table.
#           --iterations 5 is Flye's internal long-read polishing, run five times
#           (v1 value). {FLYE_INPUT_MODE} is --nano-hq or --nano-raw, resolved
#           once in 00_common section 8 from the basecalling quality.
# Produces: FLYE_DIR (the whole Flye working directory) plus three named files:
#             FLYE_CONTIGS      the assembly itself
#             FLYE_INFO         assembly_info.txt: per contig length, coverage,
#                               circularity, repeat status
#             FLYE_IGNORE_LIST  the contig names Flye did NOT flag circular
# Consumed by: fix_start (contigs + ignore list), shared/15_replicons.smk
#              (FLYE_INFO -> the topology column of the Bakta replicon table),
#              and compare_hybrid_assemblies in hybrid mode.
#
# WHY THE IGNORE LIST: dnaapler rotates a sequence, which only makes sense for a
# circular molecule. Rotating a linear contig would move its true ends into the
# middle. IGNORE_LIST_CMD (00_common) reads assembly_info.txt and prints the name
# of every contig whose "circ." column is not "Y"; dnaapler is told to leave those
# alone. The file is legitimately EMPTY (0 bytes) when every contig is circular —
# dnaapler handles that, verified in v1 runs.
#
# LAYOUT LANDMINE (same as SPAdes in illumina mode): declare directory(FLYE_DIR),
# never directory of the sample folder. The sample folder also holds
# contaminants/, eval/, fix_start/ and medaka/ written by other rules, and a
# directory() output is wiped and recreated on a re-run.
#
# FLYE_INFO is a DECLARED output (v1 declared it too) and is now also a REQUESTED
# one, because shared/15_replicons.smk consumes it.
#
# (v1 rule name: `assembly` in BacFluxL, `ONT_assembly` in BacFluxL+; D5 unifies
#  on `ont_assembly`. v1 message: "--- Flye: Genome assembly with long reads. ---")
rule ont_assembly:
    input:
        filt_long = FILT_LONG,
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

        # Contigs Flye did NOT call circular (column "circ." != "Y") -> dnaapler
        # must not rotate them. Positional awk, unchanged from v1.
        awk {IGNORE_LIST_CMD:q} {output.flye_info} > {output.ignore_list}
        """


# ── Rule: fix_start — rotate each circular replicon to a standard start ──────
# Biology: `dnaapler all` searches every contig for the three canonical start
# genes at once — dnaA (chromosome), repA (plasmid), terL (phage terminase) —
# plus cog1474 for archaeal-type origins, and rotates the sequence so the best
# hit begins at position 1 on the forward strand. Two consequences we use later:
# assemblies of the same strain become directly comparable, and the marker that
# was found is itself evidence of what kind of replicon the contig is.
#
# Takes in: FLYE_CONTIGS and FLYE_IGNORE_LIST from ont_assembly.
# Does:     dnaapler (e-value 1e-10, fixed seed 42 so a re-run is reproducible),
#           then two awk passes over dnaapler's FASTA:
#             FASTA_HEAD_CMD — cut each header at the first whitespace token.
#                              This matters far beyond tidiness: Bakta's
#                              --replicons table, Platon, geNomad and the
#                              BlobTools/BLAST screen all join on that first
#                              token, so trimming here is what makes every later
#                              join work.
#             FASTA_LIN_CMD  — one line per sequence.
# Produces: DNAAPLER_DIR, DNAAPLER_REORIENTED (dnaapler's own output),
#           DNAAPLER_FIXED (header-trimmed + linearised; this IS DRAFT_CONTIGS)
#           and DNAAPLER_SUMMARY.
# Consumed by: the contamination screen in shared/10_decontam.smk (via
#              DRAFT_CONTIGS) and build_replicons in shared/15_replicons.smk
#              (via DNAAPLER_SUMMARY).
#
# DNAAPLER_SUMMARY is NEW as a declared output. v1 produced
# {sample}_all_reorientation_summary.tsv and never declared or surfaced it, even
# though its Gene_Reoriented / Coverage / Identity_Percentage columns already
# say, per contig, which start gene was found and how convincingly. That is
# exactly the input shared/15_replicons.smk needs to tell Bakta which contigs are
# chromosomes and which are plasmids, so we now declare it.
#
# (v1 message: "--- dnaapler: Re-orient replicons. ---")
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
