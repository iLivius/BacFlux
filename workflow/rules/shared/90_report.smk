# BacFlux v2.0.0 — Stage 90 report module (rules/shared/90_report.smk)
#
# One rule: gather every QC artefact the run produced into a single MultiQC HTML
# page, with sample names a human can read. This module is included LAST by
# Snakefile (step 4 of the include order) so it sees every upstream rule.
#
# What gets aggregated depends on what the mode actually produced, so the input
# set is assembled from the capability flags (HAS_SHORT_READS / HAS_LONG_READS /
# HAS_READS, defined in 00_common.smk) rather than from the mode name. That is
# decision D7 of the unification — gate on what a mode produces, not on what it is
# called — see docs/unification_migration_plan.md:
#
#   input                        illumina  nanopore  hybrid  contigs
#   ───────────────────────────  ────────  ────────  ──────  ───────
#   fastp JSON (read QC)             y         -        y        -
#   NanoPlot raw + filtered          -         y        y        -
#   Qualimap (mapping QC)            y         y        y        -
#   QUAST (assembly QC)              y         y        y (2)    y
#   CheckM (staged, relabelled)      y         y        y (2)    y
#   GTDB-Tk (staged, relabelled)     y         y        y (2)    y
#   Bakta (annotation)               y         y        y        y
#                                                    (2) = two genomes per sample
#
# Two things in here are easy to break, and both fail SILENTLY rather than loudly,
# so each is explained in full at its point of use:
#   1. the `cd $OUT` + relative-path trick that makes MultiQC's -d prefixes short
#      and machine-independent (see the params block);
#   2. the single-vs-double brace rule for awk (see the staging block).
#
# Everything referenced here is defined once in 00_common.smk and never
# re-derived: OUT, DIR_REPORT, DIR_ANNOTATION, LOGS, SAMPLES, CHECKM_STATS,
# GTDBTK_DIR, QUAST_DIR, QUALIMAP_DIR, FASTP_JSON, NANOPLOT_RAW_DIR,
# NANOPLOT_FILT_DIR, CHECKM_RELABEL_AWK, GTDBTK_RELABEL_AWK, QC_GENOMES,
# IS_HYBRID and the capability flags.
#
# "the D1 layout", which the rename map refers to, is v2's unified stage
# numbering: all technology-specific work is grouped under 01.reads and
# 02.assembly, so every shared stage lands on the same number in every mode
# (03.taxonomy, 04.annotation, …). That is why v1's rename regexes needed only
# their path components changed, not their form — see
# docs/unification_migration_plan.md.
#
# conda: paths resolve relative to THIS file (workflow/rules/shared/), so
# "../../envs/multiqc.yaml" climbs shared/ → rules/ → workflow/ → workflow/envs/multiqc.yaml.


# ── Which QC artefacts exist in THIS mode ──
def multiqc_qc_inputs():
    """Every QC artefact MultiQC aggregates, for the ACTIVE mode only.

    Built from the capability flags rather than the mode name (D7), so each line
    reads as a plain statement about what the mode produces. CheckM and GTDB-Tk
    are NOT in here — they are staged and relabelled first, so the rule names them
    separately (see the rule's input: block for why that split matters).
    """
    paths = []
    if HAS_SHORT_READS:
        paths += expand(FASTP_JSON, sample=SAMPLES)
    if HAS_LONG_READS:
        paths += expand(NANOPLOT_RAW_DIR, sample=SAMPLES)
        paths += expand(NANOPLOT_FILT_DIR, sample=SAMPLES)
    if HAS_READS:
        # No Qualimap in contigs mode: its BAM is a self-alignment made only to
        # feed BlobTools, so a mapping-quality report on it says nothing.
        paths += expand(QUALIMAP_DIR, sample=SAMPLES)
    paths += expand(QUAST_DIR, sample=SAMPLES)
    paths += expand(DIR_ANNOTATION + "/bakta/{sample}", sample=SAMPLES)
    return paths


# Built once at parse time (while Snakemake reads the workflow, before any job
# runs): the rule needs this SAME list in two shapes — absolute paths in input:,
# where Snakemake needs real ones to build the DAG, and their output-root-relative
# twins in params.qc_rel, which is what MultiQC is handed after the cd. Both are
# derived from this one list, so they cannot disagree.
MULTIQC_QC_INPUTS = multiqc_qc_inputs()


# ── The sample-name rewrite map ──
def multiqc_replace_block():
    """The body of MultiQC's `sample_names_replace:` mapping, as YAML text.

    MultiQC's -d flag builds each sample's display name out of the directory
    components of the file it came from, joined with " | ". Left alone that gives
    unreadable names like "02.assembly | S1 | eval | quast | S1". These regexes
    turn each one into "<what it is> | <sample>", e.g. "assembly QC | S1".

    Patterns are raw strings so the \\1 backreferences reach MultiQC intact. Both
    the regex form and the backreferences are inherited from v1, so they are known
    to work with this MultiQC version — only the path components changed (the D1
    layout). Hybrid needs two QUAST lines because a hybrid sample is QC'd as two
    genomes, {sample}_illumina and {sample}_ont.
    """
    lines = []
    if HAS_SHORT_READS:
        lines.append("  # read QC (fastp)")
        lines.append(r"  '^01\.reads \| ([^|]+) \| illumina \| .+$': 'read QC | \1'")
    if HAS_LONG_READS:
        lines.append("  # read QC before and after filtering (NanoPlot)")
        lines.append(r"  '^01\.reads \| ([^|]+) \| ont \| raw_qc \| \1$': 'raw read QC | \1'")
        lines.append(r"  '^01\.reads \| ([^|]+) \| ont \| filt_qc \| \1$': 'filtered read QC | \1'")
    if HAS_READS:
        # One catch-all tail (.+) covers both files Qualimap names after the BAM
        # and the ones under raw_data_qualimapReport/; v1 needed two lines.
        #
        # In hybrid the label is tagged "Illumina" because that is literally what
        # the panel shows: DECONTAM_BAM is Illumina reads mapped to the
        # PRE-decontamination SPAdes draft (10_decontam takes the short-read branch
        # in hybrid, so the ONT reads are never mapped for QC). Left untagged, a
        # reader would naturally take this depth/evenness panel as coverage of the
        # delivered ONT genome — which it is not. The other hybrid panels (QUAST,
        # CheckM, GTDB-Tk) are already technology-tagged from QC_GENOMES.
        mapping_label = "mapping Illumina QC" if IS_HYBRID else "mapping QC"
        lines.append("  # mapping QC (Qualimap)")
        lines.append(
            r"  '^02\.assembly \| ([^|]+) \| eval \| qualimap \| .+$': '"
            + mapping_label + r" | \1'"
        )
    # QUAST — generated by looping over QC_GENOMES, the SAME list that drives the
    # staged FASTA names, the CheckM/GTDB-Tk bin ids and the relabel awk bodies.
    # (It used to branch on IS_HYBRID with the suffixes typed in by hand, so adding
    # or renaming a QC genome updated CheckM/GTDB-Tk but silently NOT QUAST. The
    # real condition was never "is this hybrid" but "is there more than one genome".)
    lines.append("  # assembly QC (QUAST)")
    for genome in QC_GENOMES:
        # "assembly QC" when a mode has one genome; "assembly Illumina QC" /
        # "assembly ONT QC" when it has several (hybrid).
        label = " ".join(part for part in ("assembly", genome.label, "QC") if part)
        lines.append(
            r"  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1" + genome.suffix
            + r"$': '" + label + r" | \1'"
        )
    lines.append("  # annotation (Bakta)")
    lines.append(r"  '^04\.annotation \| bakta \| ([^|]+) \| \1$': 'annotation | \1'")

    # CheckM / GTDB-Tk are relabelled by awk into a staging directory BEFORE
    # MultiQC sees them, but -d then prepends that staging path to the name the awk
    # produced ("09.report | multiqc_inputs | checkm | S1 | completeness Illumina |
    # S1"). Without these two rules the carefully relabelled rows — the only place
    # the hybrid Illumina-vs-ONT distinction reaches the report — display as long
    # staging paths. Strip the prefix, keeping the awk-produced tail. The literal is
    # derived from the real staging location so it cannot drift from the rule.
    staging_pattern = os.path.relpath(DIR_REPORT + "/multiqc_inputs", OUT).replace(".", r"\.").replace("/", r" \| ")
    lines.append("  # staged CheckM / GTDB-Tk tables: drop the staging prefix -d adds")
    lines.append(r"  '^" + staging_pattern + r" \| checkm \| [^|]+ \| (.+)$': '\1'")
    lines.append(r"  '^" + staging_pattern + r" \| gtdbtk \| [^|]+ \| (.+)$': '\1'")
    return "\n".join(lines)


MULTIQC_REPLACE_BLOCK = multiqc_replace_block()


# ── NanoPlot's duplicate-stats guard (long-read modes only) ──
# Given a filtering threshold, NanoPlot writes NanoStats.txt AND
# NanoStats_post_filtering.txt into the same directory. MultiQC reads both and
# reports the sample twice, so v1 told it to ignore the second, and v2 keeps that.
# The find is only a readability guard: --ignore on a pattern that matches nothing
# is harmless, so this exists to make the intent visible, not to prevent an error.
#
# Preserved verbatim from v1, INCLUDING the unquoted $IGNORE_ARG at the use site.
# That is load-bearing, not an oversight: unquoted, the value word-splits into the
# two arguments MultiQC expects (--ignore, then the pattern), and an empty
# IGNORE_ARG vanishes from the command line altogether. Quoted, it would arrive as
# a single argument — and as a single EMPTY one in every mode where the guard
# never fires.
#
# The whole guard is baked into a parse-time STRING because a Snakemake shell:
# block is static text. In illumina or contigs mode there are no NanoPlot
# directories to name, so the guard must not be in that text at all — hence the
# else branch, which substitutes a shell comment and nothing else.
if HAS_LONG_READS:
    _nanoplot_dirs = " ".join(
        expand(NANOPLOT_RAW_DIR, sample=SAMPLES) + expand(NANOPLOT_FILT_DIR, sample=SAMPLES)
    )
    MULTIQC_IGNORE_GUARD = (
        "if find " + _nanoplot_dirs + " \\\n"
        "  -type f -name '*NanoStats_post_filtering.txt' -print -quit | grep -q .; then\n"
        "  IGNORE_ARG=\"--ignore '*NanoStats_post_filtering.txt'\"\n"
        "fi"
    )
else:
    MULTIQC_IGNORE_GUARD = "# (this mode has no long reads: nothing for MultiQC to ignore)"


# ── multiqc — one HTML report for the whole run ──
# Takes in:
#   checkm_stats / gtdbtk_dir — NAMED, because the shell loops over them to
#       rewrite bin ids before MultiQC sees them. Present in ALL four modes, so
#       naming them is safe.
#   qc_inputs — an UNNAMED flat list, because its contents are mode-dependent. A
#       Snakemake shell: string is static, so writing {input.fastp_json} in it
#       would be a parse error in nanopore mode where that key does not exist.
#       Keeping the mode-varying items anonymous and passing them positionally
#       (via a params twin, see below) avoids that whole class of problem.
# Does: writes a MultiQC config, stages + relabels the CheckM and GTDB-Tk tables,
#       then runs MultiQC from the output root.
# Produces:
#   report_html  = 09.report/multiqc_report.html — the terminal deliverable that
#                  _downstream_targets() requests
#   multiqc_yaml = the generated config, kept so the rename rules are inspectable
#   staging      = a temp() directory holding the relabelled CheckM/GTDB-Tk copies
# Consumed by: the user. Runs in every mode — nothing switches this rule off.
#
# The report directory itself is deliberately NOT a declared output. v1 declared
# `multiqc_dir = directory("09.report")`, which cannot work here: a directory()
# output may not also contain the staging directory as a separate declared output,
# and Snakemake wipes a directory() output before re-running its rule, taking
# anything else the user kept in 09.report with it. See the --force paragraph in
# the shell body — it is the other half of the same decision.
rule multiqc:
    input:
        checkm_stats = expand(CHECKM_STATS, sample=SAMPLES),
        gtdbtk_dir = expand(GTDBTK_DIR, sample=SAMPLES),
        qc_inputs = MULTIQC_QC_INPUTS,
    output:
        report_html = DIR_REPORT + "/multiqc_report.html",
        multiqc_yaml = DIR_REPORT + "/multiqc_config.yaml",
        staging = temp(directory(DIR_REPORT + "/multiqc_inputs")),
    params:
        # ── Why the relative twins below exist ──
        # v1 got short, stable MultiQC sample names for free: `workdir:` made the
        # process CWD equal to output_dir, so every path handed to MultiQC was
        # relative and -d produced prefixes like "02.assembly | S1 | eval | quast".
        # v2 deliberately removed workdir: (see 00_common §2) and every constant is
        # ABSOLUTE — which would make -d prepend the user's entire installation
        # path, and NO sample_names_replace regex can be written against that,
        # because it differs on every machine.
        #
        # Fix: keep the rule's input: paths absolute (Snakemake needs real paths to
        # build the DAG) but cd into the output root and hand MultiQC the
        # output-root-relative twins computed here. The prefixes are then exactly
        # the ones the regexes above are written for, on any machine.
        #
        # Rejected alternatives: --dirs-depth N (one N cannot normalise tools
        # sitting at different depths — fastp is 3 levels down, QUAST 4); staging
        # every input into one flat directory (would copy whole Qualimap, QUAST
        # and Bakta trees).
        out_root = OUT,
        report_dir = DIR_REPORT,
        # qc_rel is a single space-joined STRING, not a list, so it lands in the
        # shell command as separate words. That is safe here because none of
        # these paths can contain a space: the stage and sub-directory names are
        # fixed literals, and 00_common.smk's BAD_CHARS check rejects any sample
        # name with a space in it before the run starts.
        qc_rel = " ".join(os.path.relpath(path, OUT) for path in MULTIQC_QC_INPUTS),
        staging_rel = os.path.relpath(DIR_REPORT + "/multiqc_inputs", OUT),
        # Parse-time text blocks (see the helpers at the top of this file).
        replace_block = MULTIQC_REPLACE_BLOCK,
        ignore_guard = MULTIQC_IGNORE_GUARD,
        checkm_relabel = CHECKM_RELABEL_AWK,
        gtdbtk_relabel = GTDBTK_RELABEL_AWK,
    conda:
        "../../envs/multiqc.yaml"
    log:
        LOGS + "/multiqc.log"
    priority: 2
    shell:
        # r""" so backslashes in the regex block survive to the YAML file.
        #
        # BRACE RULE — the single most confusable thing in this file:
        #   * awk written DIRECTLY in this shell string needs DOUBLED braces
        #     ({{ print }}), because Snakemake formats this template and would
        #     otherwise read { print } as a placeholder;
        #   * awk arriving through a params VALUE ({params.checkm_relabel}) uses
        #     SINGLE braces, because Snakemake does NOT re-scan substituted values.
        # Both appear a few lines apart below. Getting it backwards corrupts the
        # relabelling silently instead of raising an error.
        r"""
        set -euo pipefail

        # 1. Write the MultiQC config. printf lays down the two display options,
        #    then a QUOTED heredoc ('EOF' — no shell expansion) appends the rename
        #    rules exactly as written.
        #    The heredoc body and its EOF terminator sit at column 0 on purpose and
        #    must stay there: those lines ARE the file content, and MultiQC needs
        #    sample_names_replace and friends as top-level YAML keys. There is no
        #    <<- to strip leading whitespace, so indenting them to match the block
        #    around them writes an invalid config.
        mkdir -p "{params.report_dir}"

        printf "%s\n" "show_analysis_paths: False" "show_analysis_time: False" > "{output.multiqc_yaml}"

        cat >> "{output.multiqc_yaml}" << 'EOF'
sample_names_replace_regex: true
sample_names_replace_exact: true

sample_names_replace:
{params.replace_block}
EOF

        # 2. Stage CheckM and GTDB-Tk with rewritten bin ids.
        #    Neither tool's table is keyed by a filename — the sample name is the
        #    first COLUMN — so a copy is made with that column rewritten into the
        #    report label ("completeness | S1", or "completeness Illumina | S1" and
        #    "completeness ONT | S1" for a hybrid sample's two genomes). The awk
        #    body comes from QC_GENOMES, the same list that named the staged FASTAs
        #    in the first place.
        #    The staging tree is rebuilt from scratch each time: MultiQC is pointed
        #    at the whole directory further down, so a per-sample table left behind
        #    by an interrupted run would otherwise be aggregated into this report.
        REPORT_INPUT_DIR="{output.staging}"
        rm -rf "$REPORT_INPUT_DIR"
        mkdir -p "$REPORT_INPUT_DIR/checkm" "$REPORT_INPUT_DIR/gtdbtk"

        # The sample name is recovered from the FILE name here, not the directory:
        # under the D1 layout every sample's CheckM directory is called "checkm".
        for f in {input.checkm_stats}; do
          sample=$(basename "$f" _checkm_stats.tsv)
          mkdir -p "$REPORT_INPUT_DIR/checkm/$sample"
          awk -F '\t' -v OFS='\t' -v sample="$sample" '
            NR == 1 {{ print; next }}
            {params.checkm_relabel}
            {{ print }}
          ' "$f" > "$REPORT_INPUT_DIR/checkm/$sample/checkm_stats.tsv"
        done

        # GTDB-Tk keeps one directory per sample, so basename still gives the
        # sample. The summary lookup is defensive (kept verbatim from v1): try
        # classify/ first, then the directory itself, and skip the sample silently
        # if neither holds a summary — a genome GTDB-Tk could not place should not
        # take the whole report down.
        for gtdbtk_dir in {input.gtdbtk_dir}; do
          sample=$(basename "$gtdbtk_dir")
          summary=""
          if [[ -d "$gtdbtk_dir/classify" ]]; then
            summary=$(find "$gtdbtk_dir/classify" -maxdepth 1 -name 'gtdbtk.*.summary.tsv' -print -quit)
          fi
          if [[ -z "$summary" ]]; then
            summary=$(find "$gtdbtk_dir" -maxdepth 1 -name 'gtdbtk.*.summary.tsv' -print -quit)
          fi
          if [[ -n "$summary" ]]; then
            mkdir -p "$REPORT_INPUT_DIR/gtdbtk/$sample"
            awk -F '\t' -v OFS='\t' -v sample="$sample" '
              NR == 1 {{ print; next }}
              {params.gtdbtk_relabel}
              {{ print }}
            ' "$summary" > "$REPORT_INPUT_DIR/gtdbtk/$sample/$(basename "$summary")"
          fi
        done

        # 3. Run MultiQC from the output root, on relative paths (see params).
        #    Only the INPUT paths below are relative — that is the entire point of
        #    the cd. The config, the outdir and the log redirect stay absolute
        #    (00_common.smk builds all three off OUT), so they still resolve after
        #    the working directory moves.
        cd "{params.out_root}"

        IGNORE_ARG=""
        {params.ignore_guard}

        # --force is REQUIRED for re-runs, not a convenience. Without it MultiQC
        # refuses to overwrite an existing multiqc_report.html and silently writes
        # multiqc_report_1.html instead ("Existing reports found, adding suffix to
        # filenames"). The rule then fails on a missing declared output even though
        # MultiQC exited 0 — so the workflow works on a clean directory and breaks
        # on the first re-run, which is the normal case when a sample is added.
        # It stayed invisible through the whole v2 validation for exactly that
        # reason: every earlier run went into an empty output directory.
        #
        # v1 never hit this only because it declared the whole 09.report/ as a
        # directory() output, which Snakemake wipes before each re-run — the same
        # data-loss hazard fixed elsewhere in v2 (it also deleted anything else the
        # user kept in 09.report). v2 declares the real files instead and lets
        # --force do the overwriting.
        multiqc \
          $IGNORE_ARG \
          --config "{output.multiqc_yaml}" \
          --force \
          -d \
          {params.qc_rel} \
          "{params.staging_rel}/checkm" \
          "{params.staging_rel}/gtdbtk" \
          --outdir "{params.report_dir}" \
          --filename "multiqc_report.html" > "{log}" 2>&1
        """
