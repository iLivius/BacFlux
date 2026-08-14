#!/usr/bin/env python3
"""Keep the assembly contigs that belong to the isolate, drop the rest, and record
a reason for every contig either way.

Runs as rule select_contigs (shared/10_decontam.smk), the last step of the
contamination screen. Two legs feed that screen — this sample's own reads mapped
back to its own draft (per-contig depth) and megablast against NCBI nt (per-contig
taxonomy) — and BlobTools joins them into one table. Only the taxonomy half
matters here: the genus is blobtools' bestsum call, i.e. the taxon with the
highest summed BLAST bitscore. Depth sits in the same table and this script never
looks at it. The step runs BEFORE annotation because a contaminant contig left in
inflates CheckM's contamination estimate, pulls GTDB-Tk off the right lineage, and
pollutes every gene, AMR and plasmid call after it.

What it reads
-------------
  --bestscore  BLOB_TABLE, written by rule blob_table — one row per contig with
               its best-scoring taxon at every rank. The genus is column 22.
  --contigs    DRAFT_CONTIGS, the same contigs BlobTools judged.
  the policy   seven flags carrying parameters.decontamination from
               config/config.yaml, resolved once by _decontam_settings() in
               shared/00_common.smk so a whole batch follows one policy: mode
               (auto | include | exclude | off), discard_no_hit, and the five
               genus-list and override paths.

What it writes
--------------
  --output-fasta  DECONTAM_CONTIGS — the kept sequences, wrapped at 80 columns.
                  Which rule consumes them depends on the mode; the per-mode
                  table in shared/10_decontam.smk's header says which.
  --output-list   the kept contig IDs, one per line. No other rule reads it.
  --composition   COMPOSITION — each genus's share of the assembly. Read back by
                  rule annotation in shared/40_annotation.smk to pick Bakta's
                  --genus hint, so its shape is a contract; see write_composition.
  --decisions     CONTIG_DECISIONS — the audit trail the project convention asks
                  for, four columns:
                    contig          first whitespace token of the contig name
                                    source: column 1 of the BlobTools table
                    assigned_genus  the genus BlobTools settled on, or the
                                    literal "no-hit"
                                    source: column 22 of the BlobTools table
                    action          keep | remove
                                    source: decide()
                    reason          one token saying why — decide() carries the
                                    full vocabulary
                                    source: decide()

Before trusting auto mode
-------------------------
Auto mode votes on the NUMBER OF CONTIGS per genus, not on how much DNA each
genus holds, so a genus can win on many short contigs while another genus holds
most of the genome. That matters because BLAST regularly spreads ONE organism
across several related genus names whenever the isolate is thinly represented in
nt: the vote then lands on the wrong half and the other half leaves as
contamination. Three things push back, none of them a cure — the curated aliases
in GENUS_EQUIVALENCE_ALIASES fold the commonest splits back into one target,
while the two advisory warnings in warn_if_selection_looks_wrong() and the DNA
column of the composition report only make the damage visible.

In HYBRID mode a dropped contig can take its ONT reads with it and disappear from
the assembly, not merely from the taxonomy table. The worked case and the fixes
are in the decontamination block of config/config.yaml.

Stdlib only and no conda environment: it runs in the environment Snakemake was
launched from, like build_bakta_replicons.py and plasmid_concordance.py.
"""

import argparse
import csv
import sys
from collections import Counter, OrderedDict


# The four values parameters.decontamination.mode accepts:
#   auto     keep the most abundant genus BlobTools reports for this sample
#   include  keep only the genera the user listed
#   exclude  drop only the genera the user listed
#   off      keep everything (the audit files are still written)
VALID_MODES = {"auto", "include", "exclude", "off"}

# Booleans reach this script as text — a YAML value Snakemake stringified, a TSV
# cell, or a word on the command line. The accepted spellings are listed here
# rather than guessed at; parse_bool() says what happens to anything else.
FALSE_VALUES = {"false", "0", "no", "off"}
TRUE_VALUES = {"true", "1", "yes", "on"}

# Genus names that BLAST and BlobTools routinely split ONE organism across,
# treated here as a single target so a genuine isolate contig is not thrown out
# as a contaminant.
#
# Peribacillus and Priestia are 2020 splits of Bacillus sensu lato, so hits for a
# single genome come back under either name. One isolate lost its entire 4.5 Mb
# chromosome to exactly this: the contig count tied 2-2 between Bacillus and
# Peribacillus, the alphabetical tiebreak in choose_auto_genus() handed the vote
# to Bacillus, and the chromosome left as contamination. A later reassembly of the
# same reads called that small contig Priestia instead — which is why both names
# are here; patching only one would have hit the same bug from the other side.
#
# Pseudoarthrobacter (with the 'o') is deliberately absent: it is not a validly
# published name. LPSN lists only Pseudarthrobacter (Busse 2016) — the extra 'o'
# is what you get from concatenating pseudo- and arthrobacter without the Latin
# elision. It was removed on purpose and should not be added back.
#
# This is a safeguard against false removal, not taxonomic reconciliation.
GENUS_EQUIVALENCE_ALIASES = {
    "paenibacillus": ("bacillus",),
    "peribacillus": ("bacillus",),
    "priestia": ("bacillus",),
    "pseudarthrobacter": ("arthrobacter",),
    "paenarthrobacter": ("arthrobacter",),
    "paraburkholderia": ("burkholderia",),
}

# The same idea, blunter: strip a leading prefix to reach the parent genus.
# Bradyrhizobium/Mesorhizobium/Neorhizobium/Sinorhizobium collapse to rhizobium,
# Aeribacillus/Caldibacillus/Geobacillus to bacillus. Carried over from the
# empirical regex the original BacFlux selector used, and kept because it covers
# splits the curated table above has no entry for.
GENUS_EQUIVALENCE_PREFIXES = (
    "brady",
    "meso",
    "neo",
    "sino",
    "aeri",
    "caldi",
    "geo",
)


# ── Normalise genus names and fold known genus splits into one ───────────────

def normalize(value):
    """Return value as a stripped string; None and empty TSV cells become ""."""
    return str(value or "").strip()


def normalize_genus(value):
    """Lower-case a genus name so "Pseudomonas" as BlobTools writes it and
    "pseudomonas" as someone typed it into the config compare equal."""
    return normalize(value).lower()


def genus_aliases(value):
    """Return every name a genus is allowed to match.

    Always the genus itself, lower-cased, plus whatever the curated table maps it
    to, plus the parent left after stripping a known prefix:

    - Peribacillus -> peribacillus, bacillus
    - Paraburkholderia -> paraburkholderia, burkholderia
    - Bradyrhizobium -> bradyrhizobium, rhizobium

    A prefix is stripped only when four or more characters survive it, so a short
    name is never cut down to a stub that would match half the database.

    Used by auto and include matching only. Exclude stays exact, because a
    too-broad alias there removes contigs rather than rescuing them.
    """
    genus = normalize_genus(value)
    aliases = {genus} if genus else set()
    aliases.update(GENUS_EQUIVALENCE_ALIASES.get(genus, ()))
    for prefix in GENUS_EQUIVALENCE_PREFIXES:
        if genus.startswith(prefix) and len(genus) > len(prefix) + 3:
            aliases.add(genus[len(prefix):])
    return aliases


def genus_matches(genus, targets):
    """Return True when a genus matches any target directly or by alias."""
    genus_keys = genus_aliases(genus)
    return any(genus_keys & genus_aliases(target) for target in targets)


# ── Read the keep/drop policy from config, genus files and overrides ─────────

def parse_bool(value, default=False):
    """Parse a config-style boolean: true/false, yes/no, on/off, 1/0.

    An empty value falls back to default. Anything unrecognised raises, because
    a misspelling would otherwise become False and silently flip whether
    unclassified contigs are discarded.
    """
    text = normalize(value).lower()
    if not text:
        return default
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"Cannot parse boolean value '{value}'.")


def split_genera(value):
    """Split one genus list into individual genus names.

    The same string may arrive in any shape the config and the override tables
    allow — a single genus, "GenusA;GenusB", "GenusA,GenusB", or a tab-separated
    TSV cell — so every separator is turned into a tab first, then empty items
    are dropped.
    """
    text = normalize(value)
    if not text:
        return []
    for separator in [",", ";"]:
        text = text.replace(separator, "\t")
    return [item.strip() for item in text.split("\t") if item.strip()]


def read_lines_file(path):
    """Read the one-genus-per-line file behind --exclude-genera-file.

    Blank lines and lines starting with # are skipped so the file can be
    commented; a line may still hold several genera in any of split_genera's
    shapes. An empty path means the option was not configured, and returns [].
    """
    if not path:
        return []
    genera = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genera.extend(split_genera(line))
    return genera


def read_overrides(path):
    """Read per-sample settings from the --sample-overrides TSV, keyed by sample.

    Optional. Columns sample and mode are required; include_genera,
    exclude_genera and discard_no_hit may follow. A row replaces the global
    setting for that one sample and leaves every other sample alone — the usual
    case is one run-wide policy plus one awkward isolate.

    The trap: an EMPTY include/exclude cell means "keep the global list", NOT
    "match nothing". A row must therefore still carry the full set of tabs even
    where fields are blank, or the columns shift and the wrong value is read.
    """
    if not path:
        return {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "mode"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Override file '{path}' must contain at least columns: sample, mode."
            )
        overrides = {}
        for row in reader:
            sample = normalize(row.get("sample"))
            if not sample or sample.startswith("#"):
                continue
            mode = normalize(row.get("mode")).lower()
            if mode and mode not in VALID_MODES:
                raise ValueError(
                    f"Invalid mode '{mode}' for sample '{sample}' in '{path}'."
                )
            # The two genus columns are parsed into lists here, but discard_no_hit
            # stays raw TEXT. Parsing it now would need a default, and the natural
            # False would make a blank cell look like an explicit "false" and
            # quietly override the run-wide setting. Keeping the raw string lets
            # main() tell "not specified" from "specified as false".
            overrides[sample] = {
                "mode": mode,
                "include": split_genera(row.get("include_genera")),
                "exclude": split_genera(row.get("exclude_genera")),
                "discard_no_hit": row.get("discard_no_hit"),
            }
    return overrides


def read_include_by_sample(path, sample):
    """Read this sample's include genera from the --include-genera-by-sample TSV.

    Two columns, sample and genus, and a sample may appear on several rows. The
    compact alternative to a full override table when many samples each need
    their own include list but nothing else differs.
    """
    if not path:
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "genus"}.issubset(reader.fieldnames):
            raise ValueError(
                f"Include file '{path}' must contain columns: sample, genus."
            )
        genera = []
        for row in reader:
            if normalize(row.get("sample")) == sample:
                genera.extend(split_genera(row.get("genus")))
        return genera


# ── Read the BlobTools genus table and the assembly FASTA ────────────────────

def read_bestscore(path):
    """Read the BlobTools table into [(contig, genus), ...] plus a per-genus count.

    Row order is preserved so the decisions file lists contigs in the order
    BlobTools reported them. The counts are what auto mode votes on and what the
    composition report's contig column is built from.

    Raises rather than returning an empty list, and raises HERE rather than
    letting main() fail later on "no contigs were kept" — that message sends the
    reader to the decontamination policy, when the real fault is an unreadable or
    truncated BlobTools table.
    """
    records = []
    counts = Counter()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 22:
                continue
            contig = fields[0].strip()

            # The genus sits in column 22 of `blobtools view --rank all`, i.e.
            # fields[21] counting from zero — which is also why a row with fewer
            # than 22 fields is skipped above rather than half-read. An empty
            # cell means BLAST placed nothing, and is folded into "no-hit" so the
            # rest of the script has one spelling to test for.
            genus = fields[21].strip() or "no-hit"
            if not contig:
                continue

            records.append((contig, genus))
            counts[genus] += 1
    if not records:
        raise ValueError(f"No contig taxonomy records could be parsed from '{path}'.")
    return records, counts


def read_fasta(path):
    """Read the assembly into {full header: sequence}, in the input record order.

    Order is kept so the filtered FASTA comes out in the same layout as the
    assembly it came from. Sequence lines are accumulated and joined at the end,
    which is what lets this read WRAPPED FASTA directly — v1 FastaFlux
    re-linearised the file into contigs_filt_lin.fasta before calling the
    selector, and that step was dropped in v2 as redundant.
    """
    records = OrderedDict()
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # A new header closes the previous record.
                if header is not None:
                    records[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        # The last record has no following header to close it, so store it here.
        records[header] = "".join(seq_chunks)
    if not records:
        raise ValueError(f"No FASTA records could be parsed from '{path}'.")
    return records


def fasta_key(header):
    """Return the first whitespace token of a FASTA header, which is the contig ID.

    A header can carry description text after the ID — ">contig_1 length=1234" —
    while the BlobTools table names only "contig_1". The first token is what the
    two are joined on, here and everywhere else in BacFlux.
    """
    return header.split()[0]


# ── Decide keep or remove, one contig at a time ──────────────────────────────

def is_no_hit(genus):
    """Return True for a contig BLAST could not place.

    A substring test, not an equality test: read_bestscore() folds an empty genus
    cell into the bare "no-hit", but BlobTools' own label is not guaranteed to be
    exactly that string, and every variant carrying it means the same thing here.
    """
    return "no-hit" in normalize_genus(genus)


def choose_auto_genus(counts, discard_no_hit):
    """Pick the genus auto mode keeps: the one carried by the most contigs.

    Counts CONTIGS, not bases. A genus can win here on many short contigs while
    another genus holds most of the assembly's DNA — write_composition() prints
    both figures side by side precisely so that disagreement is visible.

    With discard_no_hit true, "no-hit" is dropped from the vote so unclassified
    contigs can never become the target.

    Ties are broken alphabetically so the answer is the same on every run, and
    that tiebreak is not harmless: a 2-2 tie between Bacillus and Peribacillus
    once handed the vote to Bacillus and sent a 4.5 Mb chromosome out as
    contamination. The aliases in GENUS_EQUIVALENCE_ALIASES exist so that tie
    never has to be broken in the first place.
    """
    candidates = [
        (genus, count)
        for genus, count in counts.items()
        if not discard_no_hit or not is_no_hit(genus)
    ]
    if not candidates:
        return None

    candidates.sort(key=lambda item: (-item[1], normalize_genus(item[0])))
    return candidates[0][0]


def decide(contig, genus, mode, include, exclude, discard_no_hit, auto_genus):
    """Return ("keep"|"remove", reason) for one contig.

    reason becomes the fourth column of contig_taxonomy_decisions.tsv and is the
    only record of why a contig went. The full vocabulary, so an audit file can
    be read without opening this function:

      mode_off               mode is off; nothing is filtered
      discarded_no_hit       unplaced by BLAST, and discard_no_hit is true
      auto_genus:<G>         auto mode, this is the winning genus
      auto_genus_alias:<G>   auto mode, an alias of the winning genus
      not_auto_genus:<G>     auto mode, some other genus
      auto_no_genus          auto mode found no genus at all to keep
      included_genus         named in include_genera
      included_genus_alias   an alias of something in include_genera
      not_included_genus     absent from include_genera
      excluded_genus         named in exclude_genera (exact match only)
      not_excluded_genus     survived exclude mode
    """
    genus_key = normalize_genus(genus)
    include_keys = {normalize_genus(item) for item in include}
    exclude_keys = {normalize_genus(item) for item in exclude}

    # Mode off is an explicit bypass: keep everything, including no-hit contigs.
    if mode == "off":
        return "keep", "mode_off"

    # discard_no_hit is applied BEFORE the mode logic, so it removes unclassified
    # contigs in auto, include and exclude alike — the only escape is mode off.
    # An unplaced contig is often short and low-coverage, but it can equally be a
    # small plasmid nt has no near neighbour for, which is why the option exists
    # rather than being hard-coded.
    if is_no_hit(genus) and discard_no_hit:
        return "remove", "discarded_no_hit"

    # Auto mode keeps only the winning genus chosen by choose_auto_genus(), and
    # is the convenience setting for a clean single-organism culture. Aliases are
    # allowed here: a genome split across Bacillus and Peribacillus should come
    # through whole, not half.
    if mode == "auto":
        if auto_genus is None:
            return "remove", "auto_no_genus"
        if genus_key == normalize_genus(auto_genus):
            return "keep", f"auto_genus:{auto_genus}"
        if genus_matches(genus, [auto_genus]):
            return "keep", f"auto_genus_alias:{auto_genus}"
        return "remove", f"not_auto_genus:{auto_genus}"

    # Include mode is strict: the genus list is mandatory and everything outside
    # it goes. Empty here means the user asked for include and gave no genus, so
    # raise rather than delete the whole assembly.
    if mode == "include":
        if not include_keys:
            raise ValueError("Mode 'include' requires at least one genus.")
        if genus_key in include_keys:
            return "keep", "included_genus"
        if genus_matches(genus, include):
            return "keep", "included_genus_alias"
        return "remove", "not_included_genus"

    # Exclude mode is the inverse: only the named contaminant genera go, and
    # anything else stays. Matching here is EXACT, with no alias expansion — a
    # broad alias in include mode rescues a contig, but the same alias here would
    # delete contigs the user never named.
    if mode == "exclude":
        if genus_key in exclude_keys:
            return "remove", "excluded_genus"
        return "keep", "not_excluded_genus"
    raise ValueError(f"Unsupported decontamination mode '{mode}'.")


# ── Write the composition report and warn when the call looks wrong ──────────

def write_composition(path, counts, total, bases=None, total_bases=0):
    """Write genus composition as relative frequencies for quick inspection.

    The output is not used for filtering. It is a human-readable summary that
    helps decide whether auto/include/exclude settings make biological sense.

    Two figures are given per genus, and the difference between them matters:

      bases   - the share of the assembly's DNA carried by that genus
      contigs - the share of the contig COUNT assigned to that genus

    Auto mode picks its target genus on the contig count, so a genus can win the
    vote on many short contigs while another genus holds far more of the actual
    genome. When that happens the two columns disagree and the file says so at a
    glance. A real case: one isolate listed Bacillus first on contig share (0.30
    vs 0.26) while Aneurinibacillus held more DNA (0.40 vs 0.29) — the wrong half
    of a single genome was kept, and the count-only report gave no hint of it.

    Sorted by DNA, because that is the more honest ranking of "what is this
    sample mostly made of".

    The line SHAPE is a contract, not just a report: rule annotation in
    shared/40_annotation.smk parses this file to pick Bakta's --genus hint, and
    adding the second figure changed what that parse returns. Read the comment on
    its genus= line before changing the format again.
    """
    with open(path, "w", encoding="utf-8") as handle:
        if bases and total_bases:
            for genus, bp in sorted(bases.items(), key=lambda item: (-item[1], item[0])):
                handle.write(
                    f"{genus}: bases {bp / total_bases:.2f}  "
                    f"contigs {counts.get(genus, 0) / total:.2f}\n"
                )
        else:
            # The v1 one-figure form, kept as the fallback for when no contig in
            # the BlobTools table was found in the FASTA and there are no lengths
            # to divide by. Same "Genus: 0.87" shape, ranked by contig count.
            for genus, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
                handle.write(f"{genus}: {count / total:.2f}\n")


# A genus must hold at least this share of the assembly's DNA before it counts
# as a "major" assignment, and at least this many major genera make the sample
# worth a second look. Both are REPORTING thresholds — they change no keep/remove
# decision, they only decide whether a warning is printed.
#
# Two is deliberate, and measured. A clean isolate has ONE genus holding
# essentially all of its DNA: 54 of 56 isolates in the batch these numbers came
# from looked exactly like that, so the second major genus is already abnormal.
# The two exceptions were the only two problem samples in the batch, and they
# failed in opposite directions — which is why the warning names both causes:
#   - one genome split across related genera (Aneurinibacillus 40% / Bacillus 29%
#     / Paenibacillus 16% / Brevibacillus 11%; ended up 28% complete, 0% contaminated)
#   - a genuine two-organism culture (Priestia 61% / Bacillus 39%; 100% complete
#     and 104% CONTAMINATED, i.e. two genomes in one assembly)
# Zero false positives against the other 54. Still only one batch, so treat it as
# a well-supported starting point rather than a universal constant.
MAJOR_GENUS_BASE_FRACTION = 0.05
CONFUSED_GENUS_COUNT = 2

# Removing more than this share of the assembly is worth announcing. A genuinely
# contaminated culture can legitimately exceed it; so can a filter that has just
# thrown away the target genome. The warning does not distinguish them — it asks
# a human to look.
LARGE_REMOVAL_BASE_FRACTION = 0.20


def warn_if_selection_looks_wrong(bases, total_bases, removed_bases):
    """Print a warning when the genus assignment, or the amount being discarded,
    suggests the keep/remove call should be checked by a human.

    Why this exists. Auto mode assumes one genus dominates and everything else is
    contamination. That assumption breaks when BLAST spreads ONE genome across
    several related genera — which happens when the organism is thinly
    represented in the nucleotide database, so different contigs match different
    relatives. The selector cannot tell that apart from real contamination, and
    without this warning it proceeds silently either way.

    Measured on a 56-isolate batch: 55 samples had a single genus holding
    essentially all the DNA and discarded almost nothing (<3%). The one failure
    had FOUR genera above 5% and discarded 52% of the assembly — half of a single
    Aneurinibacillus genome that BLAST had scattered across 14 genus labels.
    The separation was total, but it is one bad sample against 55 good ones, so
    treat these numbers as a first cut rather than a calibrated cutoff.

    Input:  per-genus base counts, the assembly total, and how much is being
            removed (all in bp).
    Output: nothing; prints to stdout, which the Snakemake rule captures into
            logs/select_contigs_{sample}.log.
    """
    if not total_bases:
        return

    major = [g for g, bp in bases.items() if bp / total_bases >= MAJOR_GENUS_BASE_FRACTION]
    removed_fraction = removed_bases / total_bases

    if len(major) >= CONFUSED_GENUS_COUNT:
        listed = ", ".join(
            f"{g} {bases[g] / total_bases:.0%}"
            for g in sorted(major, key=lambda g: -bases[g])
        )
        print(
            f"WARNING: {len(major)} genera each hold at least "
            f"{MAJOR_GENUS_BASE_FRACTION:.0%} of this assembly ({listed}). "
            "Either the sample is a mixed culture, or one genome is being split "
            "across related genera by the BLAST assignment. Check the kept/removed "
            "split in contig_taxonomy_decisions.tsv before trusting this assembly."
        )

    if removed_fraction > LARGE_REMOVAL_BASE_FRACTION:
        print(
            f"WARNING: decontamination is discarding {removed_fraction:.0%} of the "
            f"assembly ({removed_bases:,} of {total_bases:,} bp). If the kept "
            "assembly then looks incomplete but NOT contaminated, the filter has "
            "most likely removed genome rather than contamination."
        )


# ── Apply the policy and write the four output files ─────────────────────────

def main():
    # Every flag is filled in by rule select_contigs in shared/10_decontam.smk
    # from parameters.decontamination.
    #
    # The six optional flags use nargs="?" with a const value (empty string, or
    # "false" for --discard-no-hit): an unset YAML key reaches the shell as an
    # empty word, and without const argparse would stop the run with "expected one
    # argument" instead of reading it as "not configured". This is why leaving
    # include_genera blank in the config is safe.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bestscore", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=sorted(VALID_MODES), default="auto")
    parser.add_argument("--include-genera", nargs="?", const="", default="")
    parser.add_argument("--include-genera-by-sample", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera-file", nargs="?", const="", default="")
    parser.add_argument("--sample-overrides", nargs="?", const="", default="")
    parser.add_argument("--discard-no-hit", nargs="?", const="false", default="false")
    parser.add_argument("--output-list", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--composition", required=True)
    parser.add_argument("--decisions", required=True)
    args = parser.parse_args()

    # Start from the run-wide policy. Each genus list is the inline config value
    # PLUS whatever the optional file adds, so a long shared exclude list can live
    # in a file while one or two genera stay visible in config.yaml.
    mode = args.mode
    include = split_genera(args.include_genera)
    include.extend(read_include_by_sample(args.include_genera_by_sample, args.sample))
    exclude = split_genera(args.exclude_genera)
    exclude.extend(read_lines_file(args.exclude_genera_file))
    discard_no_hit = parse_bool(args.discard_no_hit, default=False)

    # Then let this sample's override row REPLACE what it names, and only that.
    # A blank cell falls through and the run-wide value stands — which is how one
    # awkward isolate gets its own genus list without disturbing the batch.
    overrides = read_overrides(args.sample_overrides)
    if args.sample in overrides:
        override = overrides[args.sample]
        mode = override["mode"] or mode
        if override["include"]:
            include = override["include"]
        if override["exclude"]:
            exclude = override["exclude"]
        if override["discard_no_hit"] is not None and normalize(override["discard_no_hit"]):
            discard_no_hit = parse_bool(override["discard_no_hit"], default=discard_no_hit)

    # The two real inputs. They are joined on the contig ID, so a rename anywhere
    # upstream shows up as the missing-contig warning at the end of this function.
    records, counts = read_bestscore(args.bestscore)
    fasta_records = read_fasta(args.contigs)

    # Index the sequences by the short contig ID BlobTools uses, but carry the
    # FULL original header alongside so the filtered assembly keeps the
    # length/coverage metadata the assembler wrote into it.
    fasta_by_id = {fasta_key(header): (header, sequence) for header, sequence in fasta_records.items()}
    total = sum(counts.values())
    auto_genus = choose_auto_genus(counts, discard_no_hit) if mode == "auto" else None

    # One pass over the BlobTools rows: decide, write the audit line, tally the
    # sequence, and collect the survivors.
    kept = []
    seen = set()
    # How much SEQUENCE each genus holds, and how much of it is leaving. These
    # totals feed the composition report and the two warnings only — the
    # keep/remove decision above is still made on contig counts. The lengths come
    # from the FASTA already in memory, so nothing extra is read from disk.
    bases_by_genus = Counter()
    total_bases = 0
    removed_bases = 0
    with open(args.decisions, "w", encoding="utf-8") as decisions:
        decisions.write("contig\tassigned_genus\taction\treason\n")
        for contig, genus in records:
            action, reason = decide(contig, genus, mode, include, exclude, discard_no_hit, auto_genus)
            decisions.write(f"{contig}\t{genus}\t{action}\t{reason}\n")

            # A contig with no FASTA record contributes no sequence, so it is
            # skipped here; the missing-contig warning further down reports those.
            if contig in fasta_by_id:
                length = len(fasta_by_id[contig][1])
                bases_by_genus[genus] += length
                total_bases += length
                if action != "keep":
                    removed_bases += length

            # Kept, present in the FASTA, and not already collected. The `seen`
            # test guards against a BlobTools table that lists the same contig on
            # more than one row: writing it twice would hand Bakta and CheckM a
            # duplicated sequence.
            if action == "keep" and contig in fasta_by_id and contig not in seen:
                kept.append(contig)
                seen.add(contig)

    # The composition report describes the assembly BlobTools saw, kept and
    # dropped contigs alike — that is what makes it useful for judging whether
    # the policy was right. Written here, after the loop, so it exists even when
    # almost nothing survived.
    write_composition(args.composition, counts, total, bases_by_genus, total_bases)

    # Announce a taxonomically confused sample, or an unusually large removal,
    # before the run moves on. Purely advisory: no decision above depends on it.
    warn_if_selection_looks_wrong(bases_by_genus, total_bases, removed_bases)

    # The kept IDs as plain text. It is a declared output of rule select_contigs,
    # but no other rule reads it — the filtered FASTA written next is what
    # everything downstream consumes.
    with open(args.output_list, "w", encoding="utf-8") as handle:
        for contig in kept:
            handle.write(f"{contig}\n")

    # The filtered assembly itself, wrapped at 80 columns. Where it goes next
    # depends on the mode — the per-mode table in shared/10_decontam.smk's header
    # says whether this file is already the delivered genome or an intermediate
    # that Medaka or Polypolish still has to polish.
    with open(args.output_fasta, "w", encoding="utf-8") as handle:
        for contig in kept:
            header, sequence = fasta_by_id[contig]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    # Contigs BlobTools judged that are not in the FASTA at all. The join is on
    # the contig ID, so this is the detector for a mismatched pair of inputs or
    # for an upstream step having rewritten the headers. It goes to stderr, which
    # rule select_contigs folds into logs/select_contigs_{sample}.log.
    missing = sorted({contig for contig, _ in records if contig not in fasta_by_id})
    if missing:
        print(
            f"WARNING: {len(missing)} contigs from BlobTools table were not found in FASTA.",
            file=sys.stderr,
        )
    if not kept:
        # Nothing survived. Raise rather than write an empty FASTA: Snakemake
        # would otherwise mark the rule complete and every stage from annotation
        # onwards would run on zero sequence. The decisions file written above is
        # the first place to look — every contig there carries the reason it went,
        # which separates a policy that names the wrong genus from a sample the
        # taxonomy screen could not place at all.
        raise ValueError(
            f"No contigs were kept for sample '{args.sample}' with mode '{mode}'. "
            "Check taxonomy assignments and decontamination settings."
        )
    print(f"Sample {args.sample} decontamination mode: {mode}")
    print(f"Kept {len(kept)} of {len(records)} contig assignments.")


if __name__ == "__main__":
    main()
