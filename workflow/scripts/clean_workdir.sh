#!/usr/bin/env bash
# clean_workdir.sh
#
# Strip a finished bacterial BioFlux output directory down to results, reports
# and audit trails, by deleting the bulky intermediates and every database the
# workflow downloaded into the output tree.
#
# Dry-run by default: with no flags it prints what it would delete and touches
# nothing. Only --run deletes.
#
# When to run this: after you have checked the results and are ready to archive
# or share the output directory. It is deliberately aggressive, so rerunning the
# workflow on a cleaned directory recomputes a lot (assembly and polishing
# included). That is the intended trade-off, not an oversight: this is an
# archive/hand-over step, not a between-runs tidy.
#
# Raw input reads are never touched: they live in input.illumina_dir /
# input.nanopore_dir and are read in place, never copied into the output tree.
#
# Supported layouts:
#   BacFlux v2 (all four modes: illumina, nanopore, hybrid, contigs — one
#   directory numbering scheme, so no per-mode branching is needed)
#   plus the retired v1 family kept for any output dirs still around from
#   before the v2 cutover: BacFlux, FastaFlux, BacFluxL, BacFluxL+
#
# Usage:
#   clean_workdir.sh [OUTPUT_DIR]
#   clean_workdir.sh --target OUTPUT_DIR
#   clean_workdir.sh --run [OUTPUT_DIR]
#   clean_workdir.sh --run --target OUTPUT_DIR
#   clean_workdir.sh --run --include-snakemake OUTPUT_DIR
#
# If OUTPUT_DIR is omitted, the current directory is used.
#
# --include-snakemake (opt-in, OFF by default)
#   Also delete a .snakemake/ working directory if one is found INSIDE the
#   output directory. Whether there is one depends entirely on where the
#   workflow was launched from: Snakemake creates .snakemake/ in the CURRENT
#   WORKING DIRECTORY, not in output_dir. Launch from the repo (the documented
#   way: `snakemake --sdm conda --configfile config/config.yaml`) and it lands
#   next to the Snakefile, out of this script's reach. Launch from inside the
#   output folder and it lands there, where this flag can reach it.
#
#   It is opt-in because it is categorically heavier than everything else here.
#   The rest of this script deletes regenerable intermediates; .snakemake/ holds
#   the per-rule CONDA ENVIRONMENTS (tens of GB — routinely the largest single
#   item in a finished project), the provenance metadata that drives Snakemake's
#   rerun decisions, and the run logs. Deleting it means the next run on that
#   directory rebuilds every environment from scratch, which needs network
#   access and time.
#
#   Use it when archiving or handing over a finished project. Do not use it on
#   a directory you still intend to rerun soon.

set -euo pipefail

DO_RUN=0
TARGET_DIR=""
INCLUDE_SNAKEMAKE=0

usage() {
  # --help IS the header comment block above: everything from the title line
  # down to the last line before `set -euo pipefail`. The end of the range comes
  # from that marker, not from a line number — v1.3.0 hard-coded `sed -n '1,25p'`
  # and leaked the first lines of shell code into the help text. Whatever goes in
  # the header is what users see with --help.
  sed -n '2,/^set -euo pipefail/p' "$0" | sed '$d'
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Options may come in any order, and the output directory may be given either
# positionally or after -t/--target/--output-dir. Last flag wins, which is why
# --dry-run exists at all despite being the default: it cancels an earlier --run
# on the same command line.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DO_RUN=1
      shift
      ;;
    --dry-run)
      DO_RUN=0
      shift
      ;;
    --include-snakemake)
      INCLUDE_SNAKEMAKE=1
      shift
      ;;
    -t|--target|--output-dir)
      [[ $# -ge 2 ]] || die "$1 requires a directory argument."
      TARGET_DIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      [[ -z "$TARGET_DIR" ]] || die "Only one output directory can be specified."
      TARGET_DIR="$1"
      shift
      ;;
  esac
done

[[ $# -eq 0 ]] || die "Unexpected extra argument(s): $*"

# The guard rails, in the order they have to happen. `pwd -P` resolves every
# symlink, so the two refusals below cannot be walked around by pointing at a
# link that happens to lead to / or to the home directory — checking the string
# the user typed would let exactly that through. Everything after this point is
# relative globbing from inside the target, which is why the resolve-and-cd has
# to come before any pattern is expanded.
TARGET_DIR="${TARGET_DIR:-.}"
[[ -d "$TARGET_DIR" ]] || die "Target is not a directory: $TARGET_DIR"

TARGET_DIR="$(cd "$TARGET_DIR" && pwd -P)"
[[ "$TARGET_DIR" != "/" ]] || die "Refusing to operate on /"
[[ "$TARGET_DIR" != "$HOME" ]] || die "Refusing to operate on HOME: $TARGET_DIR"

cd "$TARGET_DIR"

# True only when EVERY listed directory exists — the layout fingerprints below
# are all-or-nothing, so one missing marker rules a layout out.
has_dirs() {
  local d
  for d in "$@"; do
    [[ -d "$d" ]] || return 1
  done
  return 0
}

# Work out which workflow produced this directory by fingerprinting the stage
# folders it contains, and echo the name — that name selects the target list.
#
# v2 unified all four modes (illumina/nanopore/hybrid/contigs) onto one
# directory numbering scheme, so a single marker set covers every mode — unlike
# v1, where each repo in the family renumbered its own stages, which is why the
# chain below needs one test per repo. 08.mobilome is deliberately NOT part of
# the v2 check: the module defaults to off, so its absence does not mean "not v2
# output".
#
# The order of the chain is load-bearing. v1 BacFlux's marker set is a strict
# SUBSET of v1 BacFluxL's — the only difference is 09.report — so a BacFluxL
# directory satisfies the BacFlux test too. Testing BacFluxL first is what stops
# a BacFluxL output being cleaned with the shorter BacFlux list, which would
# silently leave the Flye working directories behind. Do not reorder these
# branches, and put any new layout above the one whose markers it contains.
detect_workflow() {
  if has_dirs 01.reads 02.assembly 04.annotation 09.report; then
    echo "BacFlux-v2"
  elif has_dirs 02.Illumina_assembly 03.post-processing 04.ONT_assembly 10.annotation 13.phages 14.report; then
    echo "BacFluxL+"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.taxonomy 05.annotation 08.phages 09.report; then
    echo "BacFluxL"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.taxonomy 05.annotation 08.phages; then
    echo "BacFlux"
  elif has_dirs 01.pre-processing 02.post-processing 03.taxonomy 04.annotation 07.phages 08.report; then
    echo "FastaFlux"
  else
    die "Target does not look like a supported BacFlux/FastaFlux/BacFluxL/BacFluxL+ output directory."
  fi
}

WORKFLOW="$(detect_workflow)"
TARGETS=()

# One glob list per layout. Patterns are relative to TARGET_DIR (we cd'd there
# above) and are expanded, not passed to rm as text — see the nullglob block
# further down.
case "$WORKFLOW" in
  BacFlux-v2)
    TARGETS+=(
      # ── Intermediate read files (01.reads) ─────────────────────────────────
      # EVERY sequence file under 01.reads is an intermediate: trimmed
      # (fastp), decontaminated (bowtie2-selected pairs) or length/quality
      # filtered (filtlong) reads. The RAW reads are never copied into the
      # output tree — they are read in place from input.illumina_dir /
      # input.nanopore_dir — so nothing here is anyone's only copy.
      #
      # Deliberately a blanket sequence-file glob rather than a list of the
      # specific names (*_sel_R1, *_filt, ...): it needs no maintenance when a
      # mode gains a new read-processing step, and it cannot silently miss one.
      # It is also, in practice, the single biggest win in a hybrid run —
      # decontaminated Illumina pairs alone routinely outweigh every other
      # remaining file combined.
      #
      # This DOES retrigger assembly and polishing if the workflow is rerun
      # later, which is accepted: this script is an archive/share step run once
      # the results are known good, not a between-runs tidy. Note what that costs
      # in hybrid mode: the decontaminated pairs (*_sel_R1/R2) are deliberately
      # NOT temp() there, because re-running Polypolish alone would otherwise
      # drag SPAdes and the whole contamination screen back through the DAG to
      # recreate them (hybrid/30_ont_reads.smk). Deleting them here accepts that.
      #
      # NOT touched here: fastp JSON/HTML and the NanoPlot raw_qc/filt_qc
      # directories in the same folders. They are QC REPORTS (and MultiQC
      # inputs), not sequence data.
      "01.reads/*/*/*.fastq"
      "01.reads/*/*/*.fq"
      "01.reads/*/*/*.fastq.gz"
      "01.reads/*/*/*.fq.gz"

      # PhiX spike-in reference plus its bowtie2 index, downloaded per project
      # (links.phix_link) for the short-read decontamination step. A database,
      # not a result — re-fetched automatically if the workflow is rerun.
      "01.reads/phix"

      # ── Per-sample assembly scaffolding (02.assembly) ──────────────────────
      # Every mode writes its raw assembler's working directory alongside the
      # final, decontaminated, polished genome (contigs_final.fasta). Once that
      # file exists the working directories below are redundant — they were
      # only ever inputs to later steps, never something a colleague reads
      # directly.
      "02.assembly/*/spades"           # illumina/hybrid — SPAdes K-mer dirs, misc, tmp, assembly graphs
      "02.assembly/*/flye"             # nanopore/hybrid — Flye 00-assembly .. 40-polishing
      "02.assembly/*/medaka"           # nanopore/hybrid — polishing BAM + consensus, already in contigs_final.fasta
      "02.assembly/*/snps/*_snps_dir"  # hybrid — per-stage Snippy BAM/VCF/HTML dirs; keeps snps/SNPs_summary.txt

      # dnaapler's reorientation working files. NANOPORE AND HYBRID, not hybrid
      # only: rule fix_start exists in both (nanopore/20_assembly.smk and
      # hybrid/40_ont_assembly.smk). What survives is
      # *_all_reorientation_summary.tsv, the audit record of which replicon was
      # rotated, to which start gene and how well it matched. That file is a
      # result, not scratch: shared/15_replicons.smk turns it into the replicon
      # table Bakta annotates against.
      #
      # In NANOPORE mode *_fixed.fasta is the assembly the contamination screen
      # ran on (it is DRAFT_CONTIGS there), so deleting it makes a later rerun
      # redo dnaapler, the screen, and the Medaka polishing that follows both.
      # Same accepted trade-off as the reads above, just a longer one.
      #
      # The .fai and .map-ont.mmi are not dnaapler's. Medaka's mini_align writes
      # a samtools index and a minimap2 index next to whatever assembly it
      # polishes, and in hybrid mode that is *_fixed.fasta. (In nanopore mode
      # Medaka polishes contaminants/assembly_decontam.fasta instead, and its
      # two index files are caught by the contaminants patterns below.)
      "02.assembly/*/fix_start/*_fixed.fasta"
      "02.assembly/*/fix_start/*_fixed.fasta.fai"
      "02.assembly/*/fix_start/*_fixed.fasta.map-ont.mmi"
      "02.assembly/*/fix_start/*_reoriented.fasta"
      "02.assembly/*/fix_start/*_MMseqs2_output.txt"
      "02.assembly/*/fix_start/dnaapler_*.log"
      "02.assembly/*/fix_start/logs"

      # The same pair of Medaka indexes on the nanopore-mode reference,
      # contaminants/assembly_decontam.fasta. These are NOT declared Snakemake
      # outputs, so nothing invalidates them between reruns automatically — and
      # they are, concretely, the exact files that let Medaka silently reuse a
      # stale alignment index after the decontamination step is rerun, polishing
      # the OLD contig set while reporting success (observed in practice, not
      # hypothetical). Routine cleanup removes the trap, not just the disk usage.
      "02.assembly/*/contaminants/*.fai"
      "02.assembly/*/contaminants/*.mmi"

      # CheckM's own scratch (bins/storage/lineage.ms) — the kept result is
      # *_checkm_stats.tsv.
      "02.assembly/*/eval/checkm/bins"
      "02.assembly/*/eval/checkm/storage"
      "02.assembly/*/eval/checkm/*.ms"

      # GTDB-Tk (03.taxonomy/{sample}) has NO cleanup target here. Unlike v1,
      # where align/identify were separate scratch dirs distinct from wherever
      # the final summary landed, v2's classify_wf (GTDB-Tk >=2.7, per
      # workflow/rules/shared/30_taxonomy.smk's own rule comment) nests the
      # real result INSIDE classify/: classify/gtdbtk.bac120.summary.tsv is
      # the file MultiQC reads, reached from the top level only via a
      # convenience symlink. A "03.taxonomy/*/classify" target here would
      # delete that real result and leave the symlink dangling — confirmed
      # the hard way, and recoverable only by rerunning taxonomic_assignment
      # for every affected sample. Do not re-add a classify/ target without
      # checking where classify_wf's own output actually lives for the GTDB-Tk
      # version pinned in workflow/envs/gtdbtk.yaml.

      # 05.amr also has no targets, for a different reason: in v2 the CARD
      # database and BBMap's per-sample ref/ index are declared temp() outputs
      # (shared/50_amr.smk), so Snakemake removes them itself once the mapping
      # rule finishes and there is nothing left here to clean. v1 kept them, which
      # is why the retired BacFluxL+ list further down still carries 11.AMR/AMR_db
      # and 11.AMR/AMR_mapping/*/ref.

      # antiSMASH's bundled reference database, fetched into the output tree
      # per project rather than pointed at a shared path. Reliably the single
      # largest item here (multi-GB).
      #
      # Both this and the dbCAN directory below have two possible shapes. With
      # directories.antismash_db / directories.dbcan_db left empty (the default)
      # BacFlux downloads its own copy and this deletes the real thing. With
      # either key set, the rule builds a directory of SYMLINKS pointing at the
      # copy you already hold, and rm -rf then removes only the links — the
      # shared database itself is never inside the output tree and is never at
      # risk. Both cases are safe to delete; only the first costs a re-download.
      "04.annotation/antismash/databases"

      # dbCAN, same two shapes as antiSMASH above. The version in this path is
      # NOT a constant: 00_common.smk derives the folder name from the archive
      # named in links.dbcan_link, so this literal matches the shipped default
      # link and would silently miss the database of a run configured against a
      # different dbCAN release.
      "04.annotation/dbcan/dbcan_db_v5.1.2"

      # eggNOG-mapper's temp working directory (emptied by the tool itself,
      # left behind as an empty shell).
      "04.annotation/eggnog/*_tmp"

      # Phage-calling stage (07.phages): the shared databases, whichever caller
      # is in use — genomad_db, vs2_db and checkv_db all end in _db, so one glob
      # covers the caller's database and CheckV's. Each is either downloaded
      # here or built as a view of directories.genomad_db / vs2_db / checkv_db,
      # as described for antiSMASH above; CheckV's view also holds a DIAMOND
      # index rebuilt locally, so deleting it costs a re-index as well.
      #
      # The two virsorter patterns below are carried over from the v1 target
      # list and have not been re-verified against a v2 VirSorter2 run.
      "07.phages/*_db"
      "07.phages/virsorter/*/iter-*"                        # unverified against a v2 VirSorter2 run; harmless no-op via nullglob if absent
      "07.phages/virsorter/*/log"

      # Mobilome module (08.mobilome): shared model and database directories,
      # fetched once per project rather than per sample. All five are
      # third-party and none is redistributed by BacFlux; the two model
      # packages are CC BY-NC-SA (academic / non-commercial only), which is why
      # the workflow fetches them at the user's request rather than shipping
      # them.
      #
      # conjscan_models is the one worth pausing over, because it is the only
      # entry here with no local-copy escape hatch in the config. Rule
      # conjscan_models (shared/80_mobilome.smk) always installs the CONJScan
      # package with `macsydata`/`msf_data install`, which pulls it from the
      # macsy-models organisation through the GitHub API — 60 requests per hour
      # for an unauthenticated client, counted PER HOST and shared with
      # everything else running on that machine. Delete this and the next
      # mobilome run has to fetch it again; a host that has already run several
      # mobilome jobs can exhaust the hourly budget and fail at the fetch.
      # Consider keeping this directory if another run is due soon.
      #
      # The other four can all be re-supplied from a copy you already hold,
      # without a download: mobilome.icescan.dir, mobilome.tncentral.dir,
      # mobilome.iceberg.dir and mobilome.isosdb.dir. icescan_models comes from
      # EBI's https mirror of the ICEfinder2 bundle (~61 MB), not from GitHub.
      "08.mobilome/conjscan_models"
      "08.mobilome/icescan_models"
      "08.mobilome/iceberg_db"
      "08.mobilome/tncentral_db"
      # ISOSDB: the IS nucleotide database behind the read-mapping copy-number
      # leg, so it only appears in illumina/hybrid runs (nanopore and contigs
      # have no reads to map).
      "08.mobilome/isosdb_db"
      # Per-sample mobilome tool directories (conjscan/, icescan/, isescan/
      # under 08.mobilome/{sample}/) are deliberately NOT listed: like
      # 04.annotation/bakta/{sample} or 04.annotation/dbcan/{sample}, they are
      # per-tool RESULTS, not scratch, and the same restraint v1 already
      # applied to Bakta/dbCAN output applies here too.
    )
    ;;

  # ── The four retired v1 layouts ────────────────────────────────────────────
  # Frozen, and kept only so output directories produced before the v2 cutover
  # can still be cleaned. Each repo in the v1 family numbered its own stages, so
  # the same content sits at a different number in each list — that repetition
  # is the reason v2 unified the numbering, not an oversight here. Read these
  # top to bottom against the v2 list above if you need to know where a stage
  # moved to. Nothing new should be added to them.
  BacFlux)
    TARGETS+=(
      "02.assembly/*/K*"
      "02.assembly/*/misc"
      "02.assembly/*/pipeline_state"
      "02.assembly/*/tmp"
      "03.post-processing/completeness_evaluation/*/bins"
      "03.post-processing/completeness_evaluation/*/storage"
      "03.post-processing/completeness_evaluation/*/*.ms"
      "04.taxonomy/*/align"
      "04.taxonomy/*/identify"
      "05.annotation/antismash/databases"
      "05.annotation/dbcan/dbcan_db_v5.1.2"
      "08.phages/*_db"
      "08.phages/checkv/*/tmp"
      "08.phages/virsorter/*/iter-*"
      "08.phages/virsorter/*/log"
    )
    ;;
  FastaFlux)
    TARGETS+=(
      "02.post-processing/completeness_evaluation/*/bins"
      "02.post-processing/completeness_evaluation/*/storage"
      "02.post-processing/completeness_evaluation/*/*.ms"
      "03.taxonomy/*/align"
      "03.taxonomy/*/identify"
      "04.annotation/antismash/databases"
      "04.annotation/dbcan/dbcan_db_v5.1.2"
      "07.phages/*_db"
      "07.phages/checkv/*/tmp"
      "07.phages/virsorter/*/iter-*"
      "07.phages/virsorter/*/log"
    )
    ;;
  BacFluxL)
    TARGETS+=(
      "01.pre-processing/*.fastq"
      "01.pre-processing/*.fq"
      "01.pre-processing/*.fastq.gz"
      "01.pre-processing/*.fq.gz"
      "02.assembly/*/00-assembly"
      "02.assembly/*/10-consensus"
      "02.assembly/*/20-repeat"
      "02.assembly/*/30-contigger"
      "02.assembly/*/40-polishing"
      "03.post-processing/completeness_evaluation/*/bins"
      "03.post-processing/completeness_evaluation/*/storage"
      "03.post-processing/completeness_evaluation/*/*.ms"
      "03.post-processing/consensus/*/*.bam*"
      "03.post-processing/consensus/*/*.bed"
      "03.post-processing/consensus/*/*.hdf"
      "03.post-processing/*/*.fai"
      "03.post-processing/*/*.mmi"
      "03.post-processing/contaminants/*/*.fai"
      "03.post-processing/contaminants/*/*.mmi"
      "04.taxonomy/*/align"
      "04.taxonomy/*/identify"
      "05.annotation/antismash/databases"
      "05.annotation/dbcan/dbcan_db_v5.1.2"
      "08.phages/*_db"
      "08.phages/checkv/*/tmp"
      "08.phages/virsorter/*/iter-*"
      "08.phages/virsorter/*/log"
    )
    ;;
  BacFluxL+)
    TARGETS+=(
      "02.Illumina_assembly/*/K*"
      "02.Illumina_assembly/*/misc"
      "02.Illumina_assembly/*/pipeline_state"
      "02.Illumina_assembly/*/tmp"
      "03.post-processing/*_sel_R1.fastq"
      "03.post-processing/*_sel_R2.fastq"
      "03.post-processing/*_ont_filt.fastq"
      "04.ONT_assembly/*/00-assembly"
      "04.ONT_assembly/*/10-consensus"
      "04.ONT_assembly/*/20-repeat"
      "04.ONT_assembly/*/30-contigger"
      "04.ONT_assembly/*/40-polishing"
      "05.ONT_consensus/*/*.bam*"
      "05.ONT_consensus/*/*.bed"
      "05.ONT_consensus/*/*.hdf"
      "08.SNPs/*/*_snps_dir"
      "09.taxonomy/*/align"
      "09.taxonomy/*/identify"
      "10.annotation/antismash/databases"
      "10.annotation/dbcan/dbcan_db_v5.1.2"
      "11.AMR/AMR_db"
      "11.AMR/AMR_mapping/*/ref"
      "13.phages/*_db"
      "13.phages/checkv/*/tmp"
      "13.phages/virsorter/*/iter-*"
      "13.phages/virsorter/*/log"
    )
    ;;
esac

# nullglob makes a wildcard that matches nothing expand to nothing instead of to
# itself (the -e test below covers the case it does NOT reach). dotglob makes "*"
# match names beginning with a dot as well, so a hidden file cannot slip through
# a pattern meant to catch everything in a directory.
shopt -s nullglob dotglob

# Expand every pattern into the concrete paths that actually exist, so the list
# printed below is exactly the list rm is handed.
to_delete=()
for pattern in "${TARGETS[@]}"; do
  matches=( $pattern )
  for match in "${matches[@]}"; do
    # The -e test is required, not belt-and-braces: `nullglob` only drops
    # patterns that CONTAIN a wildcard. A literal path with no glob character
    # (e.g. "01.reads/phix", or "08.mobilome/iceberg_db" for an optional database
    # that was never enabled) is not subject to pathname expansion at all, so it
    # passes through as itself even when nothing is there. Without this check the
    # dry run lists paths that do not exist and will not be deleted — which is
    # the wrong way round for the preview people rely on before running --run.
    # -L as well as -e so a broken symlink still counts as something to remove.
    [[ -e "$match" || -L "$match" ]] || continue
    to_delete+=( "$match" )
  done
done

# ── The zero-byte .snakemake_timestamp markers ───────────────────────────────
# Snakemake drops one inside every directory() output to track whether that
# directory is up to date. They are bookkeeping, not results: meaningless to
# anyone opening the archive, and they survive the patterns above whenever they
# sit in a directory this script KEEPS (Bakta, isescan, platon, ...). Collected
# with `find` rather than a glob because they occur at several different depths.
#
# `dotglob` above already makes the "*" patterns match hidden entries, and
# anything inside a deleted directory goes with it — this step exists only for
# the markers left behind in directories that are kept.
#
# Removing them makes Snakemake treat those directory outputs as out of date on
# a later rerun. That is the same accepted trade-off as the rest of this script.
#
# Do not confuse these with the .snakemake/ WORKING DIRECTORY handled separately
# below — same prefix, entirely different thing.
#
# Kept in their own array so the listing below can report them as a single
# count. Printing all of them (86 in a routine two-sample run) would bury the
# entries that actually matter in the preview people read before --run.
timestamp_markers=()
while IFS= read -r marker; do
  timestamp_markers+=( "$marker" )
done < <(find . -name ".snakemake_timestamp" -type f 2>/dev/null | sed 's|^\./||')

# ── The .snakemake/ working directory (--include-snakemake only) ─────────────
# Whether one is here at all depends on where the workflow was launched from,
# NOT on any config setting: Snakemake creates .snakemake/ in the current
# working directory. Launching from the repo (the documented way) puts it next
# to the Snakefile, where this script never looks. Launching from inside the
# output folder puts it here.
#
# Only the top level is checked, deliberately. A `find` for .snakemake at any
# depth could match one belonging to a DIFFERENT project that happens to live
# under this tree, and deleting tens of GB of someone else's conda environments
# is not something this script should ever do by inference. If a nested one
# really does need to go, delete it by hand.
snakemake_workdir=()
if (( INCLUDE_SNAKEMAKE )) && [[ -d ".snakemake" ]]; then
  snakemake_workdir=( ".snakemake" )
fi

# ── Report what was found, and delete it only if --run was given ─────────────
# The listing below is the whole safety mechanism: it is printed identically in
# both modes, so what a dry run shows is exactly what a --run would remove.
echo "== Bacterial BioFlux output cleanup =="
echo "Workflow: $WORKFLOW"
echo "Target: $TARGET_DIR"
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Mode: RUN (will delete)"
else
  echo "Mode: DRY-RUN (no deletion)"
fi
echo

if (( ${#to_delete[@]} == 0 && ${#timestamp_markers[@]} == 0 && ${#snakemake_workdir[@]} == 0 )); then
  echo "Nothing to clean; no matching files or directories were found."
  # If the only thing here was a .snakemake/ the user did not opt into, say so
  # rather than reporting a clean directory and leaving tens of GB in place.
  if [[ -d ".snakemake" ]]; then
    echo
    echo "NOTE: a .snakemake/ working directory ($(du -sh .snakemake 2>/dev/null | cut -f1)) is present"
    echo "      but was NOT touched. Re-run with --include-snakemake to remove it."
  fi
  exit 0
fi

echo "Targets:"
for path in "${to_delete[@]}"; do
  echo "  - $path"
done
if (( ${#timestamp_markers[@]} )); then
  echo "  - ${#timestamp_markers[@]} x .snakemake_timestamp (zero-byte Snakemake directory markers)"
fi
if (( ${#snakemake_workdir[@]} )); then
  echo "  - .snakemake/ working directory ($(du -sh .snakemake 2>/dev/null | cut -f1)) -- conda envs, metadata and logs"
fi
echo

# Flag the untouched .snakemake/ whenever there IS one and it was not opted
# into, so its size is visible at exactly the moment someone is deciding what
# to archive — otherwise it is the largest thing left and the easiest to miss.
if (( ! INCLUDE_SNAKEMAKE )) && [[ -d ".snakemake" ]]; then
  echo "NOTE: .snakemake/ ($(du -sh .snakemake 2>/dev/null | cut -f1)) is present and will NOT be removed."
  echo "      It holds this project's conda environments, provenance metadata and logs."
  echo "      Add --include-snakemake to delete it too (rebuilt from scratch on the next run)."
  echo
fi

if [[ "$DO_RUN" -eq 1 ]]; then
  # All three arrays in one call, in the order they were listed. Overlap is
  # harmless: a timestamp marker inside a directory that is itself being deleted
  # no longer exists by the time rm reaches it, and -f makes that a no-op rather
  # than an error. The `--` stops any path that begins with a hyphen from being
  # read as an rm option.
  rm -rf -- "${to_delete[@]}" "${timestamp_markers[@]}" "${snakemake_workdir[@]}"
  echo "Cleanup finished."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
