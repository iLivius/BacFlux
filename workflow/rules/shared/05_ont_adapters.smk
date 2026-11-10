# BacFlux — ONT adapter and barcode trimming (rules/shared/05_ont_adapters.smk).
#
# Removes the Oxford Nanopore sequencing adapter, and the barcode block that follows it
# in a barcoded library, from the ends of the ONT reads before Filtlong and Flye see
# them. Long-read modes only (nanopore and hybrid).
#
# Why. An ONT read starts with the adapter that brought the DNA molecule to the pore;
# in a barcoded library the barcode and its flanks follow it. Dorado removes both when
# it is run by hand with the kit named, but reads basecalled live or by older
# basecallers often keep them, and the FASTQ file does not say which. Left in, they end
# up in the assembly wherever many reads end at the same place: at the ends of linear
# contigs. In the hybrid assembly of K. pneumoniae TUM24772 (ONT reads SRR30870266,
# Rapid Barcoding Kit SQK-RBK004), the linear contig_2 ended in about 110 bp of adapter,
# rapid barcode RB11 and its flanks, carried unchanged through Flye, dnaapler, Medaka
# and Polypolish. A circular contig has no ends, so a complete genome is rarely
# affected; an incomplete one, or a linear replicon, is. Trimming also lets fastplong
# split chimeric reads (two molecules joined into one read) at the adapter between them.
#
# Three rules, once per sample:
#
#   raw ONT FASTQ (os.path.join(NANOPORE_DIR, ONT), the user's file)
#     --find_ont_adapter-->  01.reads/{sample}/ont/adapters/
#                              {sample}_adapter_scan.tsv       (ONT_ADAPTER_SCAN)
#                              {sample}_adapter_choice.tsv     (ONT_ADAPTER_CHOICE)
#                              {sample}_fastplong_options.txt  (ONT_FASTPLONG_OPTIONS)
#                              {sample}_adapter_sequence.fasta (ONT_ADAPTER_SEQUENCE)
#     --trim_ont_adapters--> 01.reads/{sample}/ont/{sample}_trimmed.fastq.gz (TRIMMED_LONG,
#                            temporary) + fastplong's own JSON and HTML reports (kept)
#       --> rule filter_long_reads (nanopore/10_reads.smk, hybrid/30_ont_reads.smk) reads
#           TRIMMED_LONG instead of the raw file (FILTLONG_INPUT in 00_common.smk) and
#           writes FILT_LONG, the reads Flye assembles (kept)
#         --check_ont_trimming--> 01.reads/{sample}/ont/adapters/{sample}_adapter_check.tsv
#   The sequence FASTA also goes to rule adapter_search_input (shared/20_qc.smk), which
#   adds it to the adapters BLAST searches the genome for, so the screen looks for exactly
#   what was trimmed; the choice and check files go to rule screen_adapters, which puts
#   their numbers into the MultiQC report.
#
# What decides the sequence (parameters.<mode>.ont_adapter, read in 00_common.smk
# section 8 as ONT_ADAPTER_SETTING):
#   detect (the default) : find_ont_adapter.py works the sequence out from all the
#                          reads: the adapter and the barcode block after it, or a
#                          barcode block left without its adapter, or nothing when the
#                          reads are already trimmed. It names the barcode from ONT's
#                          published list (ONT_BARCODE_LIST) and, in a file that holds
#                          several barcodes, trims each of them. The script's opening
#                          text explains the method and every threshold.
#   a DNA sequence       : that sequence is trimmed as written (the scan still runs and
#                          reports what the reads carry).
#   off                  : none of these rules exists, and Filtlong reads the raw reads.
# fastplong's own adapter guess is never used: it looks only at the first reads of the
# file, which can be a batch of already-trimmed reads, and takes the most common
# sequence at the read ends, which can be genome.
#
# The Snakefile includes this file with the other shared files. Its rules exist only
# when ONT_ADAPTER_TRIM is true (00_common.smk section 8: a long-read mode, and
# ont_adapter not off), so in illumina and contigs mode, and with trimming off, this
# file defines nothing.
#
# Everything referenced here comes from 00_common.smk: ONT_ADAPTER_TRIM, NANOPORE_DIR,
# ONT, FILT_LONG, ONT_ADAPTER_SETTING, ONT_BARCODE_LIST, ONT_ADAPTER_SCAN,
# ONT_ADAPTER_CHOICE, ONT_FASTPLONG_OPTIONS, ONT_ADAPTER_SEQUENCE, TRIMMED_LONG,
# FASTPLONG_JSON, FASTPLONG_HTML, ONT_ADAPTER_CHECK, ONT_ADAPTER_SCRIPT,
# LONG_RUNTIME_MINUTES, LOGS, capped_cpus.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/x.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/x.yaml. The
# two script rules run in envs/python.yaml (a pinned Python): the script uses only the
# standard library.

if ONT_ADAPTER_TRIM:

    # ──────────────── What the reads carry (find_ont_adapter) ────────────────
    # Takes in: reads = the sample's RAW ONT FASTQ (plain or .gz), the same file
    #           raw_long_read_qc profiles; no rule makes it, NANOPORE_DIR and ONT come
    #           from the config. ONT_BARCODE_LIST, ONT's published native and rapid
    #           barcodes (workflow/data/adapters/ont_chtd.fasta), goes in as a param.
    # Does:     python find_ont_adapter.py scan with the setting (detect, or the
    #           sequence given). The script reads every read once; for a typical
    #           isolate (0.2-1 Gb of reads) that takes seconds to a minute.
    #           The setting is a param, so a changed ont_adapter re-runs the scan and,
    #           through its files, the trimming and everything after it.
    #           The barcode list is a param and not an input on purpose: a file in the
    #           repository gets a new modification time at every checkout or update,
    #           and as an input that alone would re-run the scan, the trimming and the
    #           whole assembly of every long-read sample in an existing output folder.
    #           The cost: a change to the list does not re-run anything by itself.
    #           The last line of the command copies the script's [BacFlux] summary
    #           line, and any NOTE, from the log to the console.
    # Produces: ONT_ADAPTER_SCAN       one row per 10,000 reads: how many start with
    #                                  the adapter, with a barcode block and no adapter
    #                                  core, and with the Mu end and no adapter core;
    #                                  a batch of trimmed reads in an untrimmed file
    #                                  shows up here.
    #           ONT_ADAPTER_CHOICE     the decision and the numbers behind it, among
    #                                  them the barcodes found and the second most
    #                                  common barcode (bleed-through from another
    #                                  sample of the same flow cell).
    #           ONT_FASTPLONG_OPTIONS  one line: "-s <the most common sequence>", with
    #                                  --disable_chimera_splitting when a sequence to
    #                                  trim is shorter than 50 bp, or "-A" when there
    #                                  is nothing to trim.
    #           ONT_ADAPTER_SEQUENCE   every sequence to trim, each in both
    #                                  orientations (FASTA), or an empty file.
    #           All kept: they are small, and they are the audit of what was done to
    #           the reads.
    # Consumed by: trim_ont_adapters (options and sequences), check_ont_trimming (the
    #              choice), and, in shared/20_qc.smk, adapter_search_input (the
    #              sequences) and screen_adapters (the choice).
    rule find_ont_adapter:
        input:
            reads = os.path.join(NANOPORE_DIR, ONT),
        output:
            scan = ONT_ADAPTER_SCAN,
            choice = ONT_ADAPTER_CHOICE,
            options = ONT_FASTPLONG_OPTIONS,
            sequence = ONT_ADAPTER_SEQUENCE,
        params:
            script = ONT_ADAPTER_SCRIPT,
            setting = ONT_ADAPTER_SETTING,
            barcode_list = ONT_BARCODE_LIST,
        resources:
            # Cluster requests (00_common.smk §3a). The script is single-threaded.
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/find_ont_adapter_{sample}.log"
        shell:
            """
            python {params.script:q} scan \
              --reads {input.reads:q} \
              --setting {params.setting:q} \
              --sample {wildcards.sample:q} \
              --barcode-list {params.barcode_list:q} \
              --out-scan {output.scan:q} \
              --out-choice {output.choice:q} \
              --out-options {output.options:q} \
              --out-sequence {output.sequence:q} > {log} 2>&1

            grep -F '[BacFlux]' {log} >&2 || true
            """


    # ──────────────── Trimming the reads (fastplong) ────────────────
    # Takes in: reads = the sample's RAW ONT FASTQ, as above; options and sequence =
    #           ONT_FASTPLONG_OPTIONS and ONT_ADAPTER_SEQUENCE from find_ont_adapter.
    # Does:     fastplong with the options the scan wrote (read from the file with
    #           $(cat ...)), the sequence file with -a, and these fixed options:
    #             -Q -L                  fastplong's quality and length filters off:
    #                                    Filtlong, the next step, is the one read filter
    #                                    of the workflow.
    #             --trimming_extension 0 cut exactly the sequence found. fastplong's
    #                                    default (10) also cuts the next 10 bases, which
    #                                    after a whole adapter and barcode are genome.
    #             -z 1                   fast compression: the file is temporary.
    #           With -s, fastplong searches the start of every read for the sequence and
    #           the end for its reverse complement, allows a quarter of its bases to
    #           differ, and removes it together with everything on the outer side of it
    #           (so an adapter the scan did not recognise, in front of a barcode block,
    #           goes too). The -a file adds every sequence to trim, searched at both
    #           ends as written; it holds each one in both orientations, because -a adds
    #           no reverse complement itself. This is how a file holding several
    #           barcodes gets all of them trimmed, and how a rapid block copied to the
    #           read end in its forward orientation is found. A sequence found in the
    #           middle of a read splits it in two, named chimera-segment-1-<read> and
    #           chimera-segment-2-<read>; the scan switches splitting off for sequences
    #           shorter than 50 bp, which would match genome by chance.
    #           With "-A", fastplong copies the reads unchanged (checked byte for byte)
    #           and the empty sequence file is not read. Read headers are kept as they
    #           are, including the basecaller tags check_medaka_model reads.
    #           The command first makes sure that the options file is not empty:
    #           without -s or -A, fastplong would guess the adapter itself, which must
    #           never happen.
    # Produces: reads = TRIMMED_LONG, {sample}_trimmed.fastq.gz, TEMPORARY: Snakemake
    #           deletes it once Filtlong has read it. Filtlong is its only reader, so
    #           only a change to the trimming itself can make Snakemake write it, and
    #           re-run Filtlong and the assembly, again.
    #           json, html = FASTPLONG_JSON and FASTPLONG_HTML, fastplong's own reports
    #           (reads and bases in and out, reads trimmed, reads split), kept for
    #           people. They are NOT given to MultiQC: MultiQC 1.33 has no fastplong
    #           module, and its fastp module would read the JSON as a fastp run. They
    #           sit in adapters/, a folder MultiQC is never pointed at.
    # Consumed by: filter_long_reads (through FILTLONG_INPUT in 00_common.smk).
    #
    # Threads: fastplong is asked for at most 16 (capped_cpus); it is quick at that.
    rule trim_ont_adapters:
        input:
            reads = os.path.join(NANOPORE_DIR, ONT),
            options = ONT_FASTPLONG_OPTIONS,
            sequence = ONT_ADAPTER_SEQUENCE,
        output:
            reads = temp(TRIMMED_LONG),
            json = FASTPLONG_JSON,
            html = FASTPLONG_HTML,
        conda:
            "../../envs/fastplong.yaml"
        threads: capped_cpus(16)
        resources:
            # Cluster requests (00_common.smk §3a).
            runtime = LONG_RUNTIME_MINUTES,
        log:
            LOGS + "/trim_ont_adapters_{sample}.log"
        shell:
            """
            if [ ! -s {input.options:q} ]; then
                echo "[BacFlux] trim_ont_adapters {wildcards.sample}: the options file {input.options} is empty; fastplong is not started, because without them it would guess the adapter itself." | tee {log} >&2
                exit 1
            fi

            fastplong \
              -i {input.reads:q} \
              -o {output.reads:q} \
              -Q -L \
              --trimming_extension 0 \
              -z 1 \
              -w {threads} \
              -j {output.json:q} \
              -h {output.html:q} \
              -a {input.sequence:q} \
              $(cat {input.options:q}) > {log} 2>&1
            """


    # ──────────────── What is left after trimming (check_ont_trimming) ────────────────
    # Takes in: reads = FILT_LONG, the reads Filtlong kept from the trimmed reads (rule
    #           filter_long_reads), and choice = ONT_ADAPTER_CHOICE from
    #           find_ont_adapter (to know which sequence was trimmed). FILT_LONG is
    #           read, not the temporary trimmed reads: a re-run of this small check
    #           (after a BacFlux update, or when its file was deleted) then never makes
    #           Snakemake trim, filter and assemble again. These are also exactly the
    #           reads Flye assembles.
    # Does:     python find_ont_adapter.py check: counts the reads that still start
    #           with the adapter core and, when a barcode block was trimmed, with the
    #           last 20 bp of the trimmed sequence. More than 1 % of reads still
    #           starting with either prints a NOTE on the console; it does not stop the
    #           run. The second count matters for older kits (such as SQK-RBK004),
    #           whose adapter lacks the core the first count looks for.
    # Produces: check = ONT_ADAPTER_CHECK, "name<TAB>value" lines, kept.
    # Consumed by: screen_adapters (shared/20_qc.smk), which puts the numbers into the
    #              Adapters table of the MultiQC report; _frontend_targets_for() in
    #              00_common.smk asks for it too.
    rule check_ont_trimming:
        input:
            reads = FILT_LONG,
            choice = ONT_ADAPTER_CHOICE,
        output:
            check = ONT_ADAPTER_CHECK,
        params:
            script = ONT_ADAPTER_SCRIPT,
        resources:
            # Cluster requests (00_common.smk §3a). The script is single-threaded.
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/check_ont_trimming_{sample}.log"
        shell:
            """
            python {params.script:q} check \
              --reads {input.reads:q} \
              --choice {input.choice:q} \
              --sample {wildcards.sample:q} \
              --out {output.check:q} > {log} 2>&1

            grep -F '[BacFlux]' {log} >&2 || true
            """
