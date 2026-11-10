#!/usr/bin/env python3
r"""Search the finished genome for sequencing adapters and barcodes, and report what is found.

Why this is needed
------------------
Sequencing adapters and barcodes that stay in the reads can end up in the assembly, most often
at the ends of linear contigs, where many reads end at the same place. NCBI screens genome
submissions for contamination, and provides FCS-adaptor to find adaptor before a submission to
GenBank; it trims adaptor found within 100 bp of a sequence end. BacFlux
trims ONT reads before assembly (rules/shared/05_ont_adapters.smk); this screen checks the
delivered genome, contigs_final.fasta, in every mode, including assemblies made elsewhere
(contigs mode) and Illumina assemblies. It reports and removes nothing.

In the BacFlux hybrid assembly of Klebsiella pneumoniae TUM24772 made without trimming, the
linear contig_2 ended in 112 bp of rapid-kit construct (ONT adapter, rapid barcode RB11 and its
flanks); this screen reports it as one end locus.

What goes in (subcommand prepare, rule adapter_search_input)
-----------------------------------------------------------
  * the adapters searched for: the ONT adapters and barcodes and the PacBio and MGI adapters
    shipped in workflow/data/adapters/, the Illumina records of NCBI UniVec (downloaded by rule
    univec_illumina) and, for an ONT sample whose reads were trimmed, the very sequences that
    were trimmed (rule find_ont_adapter), which catches a barcode no published list holds.
    They are joined into one query FASTA, kept, so the report knows what each query is.
  * the genome, with every CIRCULAR contig extended by its own first 200 bp. A circular
    contig's two ends are an arbitrary cut (dnaapler moves it to the start gene), so an adapter
    lying across the cut would otherwise be seen as two short pieces at the two ends, or not at
    all. Topology comes from the Bakta replicon table of the long-read modes (Flye's circular
    flag, column 4); in illumina and contigs mode it is not known and every contig is treated
    as linear.
Rule blast_adapters then runs blastn -task blastn-short (BLAST's setting for short queries)
with -evalue 1, adapters as queries, the extended genome as subject.

The rules (subcommand report, rule screen_adapters; fixed, no setting changes them)
-----------------------------------------------------------------------------------
A hit counts when at least 90 % of its aligned bases are identical and:
  * at the end of a LINEAR contig: it lies within 100 bp of the end (the distance NCBI's
    FCS-adaptor uses to call a hit terminal) and is at least 20 bp long. Nothing is asked of how
    much of the adapter record it covers: an adapter is often cut short at the very end of a
    contig, and a barcode left behind after its adapter was trimmed must still count;
  * inside a contig: at least 28 bp long, the length of a whole ONT adapter. Shorter internal
    matches turn up by chance. A circular contig has no ends, so every hit on it is judged by
    this rule, the junction included.
Hits on one contig that overlap, or lie within 10 bp of each other, form one locus, as several
adapter records often match the same stretch.

One sequence needs care in bacteria. The 50 bp that follow every rapid barcode are the right
end of phage Mu (identical to the last 50 bp of the Mu genome, NC_000929.1), and a genome that
carries a Mu prophage holds them naturally. A locus whose every hit lies inside those 50 bp of
its query is therefore marked as a Mu end:
  * inside a contig it is listed but not counted (reason mu_end_inside_contig): a Mu prophage
    end, since a rapid-kit remnant inside a contig would bring its barcode along;
  * at the end of a linear contig it counts (reason mu_end_at_contig_end): a rapid-kit remnant
    whose barcode part was cut off with the contig end, or a prophage end at a contig break;
    the reason names both so a reader can check.

What comes out
--------------
  --out      one row per locus, with a counted (yes/no) and a reason column, as every filtering
             decision in BacFlux: contig, contig length, topology, start and end (1-based,
             inclusive; on a circular contig a locus across the cut has an end smaller than its
             start), distance to the nearer contig end (NA on a circular contig), position (end
             or internal), the longest hit and its identity, and the adapters that matched,
             best hit first. A header line alone when nothing was found. The BLAST table of rule
             blast_adapters keeps every hit, those that pass these rules and those that do not.
  --out-mqc  one row per sample for the Adapters table of the MultiQC report: for an ONT sample
             whose reads BacFlux trimmed, the share of reads that started with the adapter or a
             barcode before trimming, the barcodes found, the second barcode (bleed-through from
             another sample of the same flow cell), what was trimmed and what was left after
             trimming (from --read-choice and --read-check); for every sample, the contig ends
             with adapter, the counted loci inside contigs, and a verdict.
The verdict never says that a genome is ready for submission: the lists cover the common
adapters only. NCBI's own screen, FCS-adaptor, is the check before a submission.

Usage (rules adapter_search_input and screen_adapters, workflow/rules/shared/20_qc.smk):
    python adapter_screen.py prepare --contigs FASTA [--replicons TSV] --lists FASTA... \
        [--trimmed FASTA] --out-query FASTA --out-subject FASTA
    python adapter_screen.py report --blast TSV --contigs FASTA [--replicons TSV] \
        --query FASTA --sample NAME [--read-choice TSV] [--read-check TSV] --out TSV --out-mqc TSV

  --contigs      the finished genome, contigs_final.fasta
  --replicons    optional: the Bakta replicon table of the long-read modes (contig, new name,
                 type, topology, name; no header); without it every contig counts as linear
  --lists        the adapter FASTA files (the shipped lists and the UniVec Illumina records)
  --trimmed      optional: the sequences fastplong trimmed from this sample's reads
  --blast        the BLAST table of rule blast_adapters, columns: qseqid sseqid pident length
                 qstart qend sstart send evalue bitscore qlen slen
  --query        the query FASTA written by prepare
  --read-choice  optional: the choice file of rule find_ont_adapter (ONT samples)
  --read-check   optional: the check file of rule check_ont_trimming (ONT samples)
An optional file given as an empty value counts as not given (the rules pass an empty list in
the modes and runs that have no such file).

Lines starting with [BacFlux] are printed for the console (the rule copies them from the log).
The report part was first written for FunFlux, BacFlux's sister workflow for fungal genomes;
the circular contigs, the Mu end, the reason column and the read-side columns were added here.
Uses only the Python standard library. Tests: test_adapter_screen.py (same folder).
"""

import argparse
import sys

END_DISTANCE = 100        # bp: a hit this close to a linear contig end is an end hit
MIN_END_LENGTH = 20       # bp aligned, for an end hit
MIN_INTERNAL_LENGTH = 28  # bp aligned, for a hit inside a contig
MIN_IDENTITY = 90.0       # percent identical bases
MERGE_GAP = 10            # bp: hits closer than this on one contig form one locus
MAX_NAMES = 5             # adapter names listed per locus; more are counted
CIRCULAR_EXTENSION = 200  # bp of a circular contig's start appended to its end for the search
MU_SLACK = 2              # bp a Mu-end hit may reach past the 50 bp Mu end on its query

# The right end of phage Mu: the 50 bp that follow every ONT rapid barcode, identical to the
# last 50 bp of the Mu genome (NC_000929.1). Also RAPID_TAIL in
# workflow/scripts/05_ont_adapters/find_ont_adapter.py.
MU_RIGHT_END = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

LOCUS_COLUMNS = [
    "contig", "contig_length", "topology", "start", "end", "distance_to_end", "position",
    "longest_hit", "identity", "adapters", "counted", "reason",
]


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Read a FASTA file as an ordered dict: first word of the header -> (full header, sequence).

    Takes in:  path, a FASTA file (one or many lines per sequence).
    Used by:   prepare() (the genome and the adapter lists) and report() (the genome and the
               query file).
    """
    records = {}
    name = None
    header = ""
    parts = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = (header, "".join(parts).upper())
                header = line[1:]
                name = header.split()[0] if header.split() else ""
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        records[name] = (header, "".join(parts).upper())
    return records


def read_topology(path):
    """Contig name -> "circular" or "linear", from the Bakta replicon table; {} without one.

    Takes in:  path, BAKTA_REPLICONS (5 tab-separated fields, no header: original locus id, new
               locus id, type, topology, name), or "" in the modes that have none.
    Produces:  the topology of every contig the table lists. The long-read modes write one row
               per contig of contigs_final.fasta, so every contig is covered there.
    Used by:   prepare() and report().
    """
    topology = {}
    if not path:
        return topology
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[0]:
                topology[fields[0]] = "circular" if fields[3].strip().lower() == "circular" else "linear"
    return topology


# ======================================================================
# SUBCOMMAND prepare: the files blastn reads
# ======================================================================

def prepare(args):
    """Write the query FASTA (all adapters) and the subject FASTA (the genome, circular contigs extended).

    Takes in:  args.contigs, args.replicons, args.lists, args.trimmed, args.out_query,
               args.out_subject.
    Does:      copies every adapter record into one query file, under the first word of its
               header (BLAST names queries by that word), and every contig into the subject
               file; a circular contig gets its own first CIRCULAR_EXTENSION bp (the whole
               contig when it is shorter) appended to its end.
    Produces:  the two FASTA files, and one [BacFlux] line with the counts.
    """
    queries = 0
    with open(args.out_query, "w") as out:
        for path in list(args.lists) + ([args.trimmed] if args.trimmed else []):
            for name, (header, sequence) in read_fasta(path).items():
                if sequence:
                    out.write(f">{header}\n{sequence}\n")
                    queries += 1
    topology = read_topology(args.replicons)
    circular = 0
    with open(args.out_subject, "w") as out:
        for name, (header, sequence) in read_fasta(args.contigs).items():
            if topology.get(name) == "circular":
                sequence = sequence + sequence[:CIRCULAR_EXTENSION]
                circular += 1
            out.write(f">{name}\n{sequence}\n")
    print(f"[BacFlux] adapter_search_input: {queries} adapter sequences; {circular} circular "
          "contig(s) extended across their start for the search.")


# ======================================================================
# SUBCOMMAND report: hits -> loci -> tables
# ======================================================================

def mu_end_intervals(query_records):
    """Where the Mu right end sits on each query: name -> list of (first, last) positions, 1-based.

    Takes in:  the query FASTA as read_fasta() returns it.
    Does:      looks for the 50 bp Mu right end, and for its reverse complement, in every query
               sequence (the rapid barcode records hold it at positions 41-90, the record of the
               rapid 3' flank is nothing else, a trimmed rapid block holds it too).
    Used by:   report(), to tell a Mu-end hit from an adapter hit.
    """
    intervals = {}
    for name, (header, sequence) in query_records.items():
        for target in (MU_RIGHT_END, reverse_complement(MU_RIGHT_END)):
            position = sequence.find(target)
            while position != -1:
                intervals.setdefault(name, []).append((position + 1, position + len(target)))
                position = sequence.find(target, position + 1)
    return intervals


def is_mu_end_hit(query, qstart, qend, intervals):
    """True if the aligned part of the query lies inside one of its Mu-end intervals."""
    first, last = min(qstart, qend), max(qstart, qend)
    for low, high in intervals.get(query, []):
        if first >= low - MU_SLACK and last <= high + MU_SLACK:
            return True
    return False


def read_hits(path, contig_lengths, topology, intervals):
    """Read the BLAST table and keep the hits that pass the rules.

    Takes in:  the BLAST table, the contig lengths (of the genome as delivered, not extended),
               the topology and the Mu-end intervals of the queries.
    Does:      for every hit, works out its start and end on the extended contig (BLAST writes
               them in reverse order for a hit on the minus strand). On a circular contig, a hit
               that lies wholly in the appended copy of the start is dropped (it repeats a hit
               at the start); a hit across the cut is kept. On a linear contig the distance to
               the nearer end decides between an end hit and an internal one; on a circular
               contig every hit is internal. A hit too short or too different for its position
               is dropped.
    Produces:  a list of dicts: contig, start, end (on the extended contig), distance (or None),
               position, length, identity, bitscore, name, mu (True for a Mu-end hit).
    Used by:   report().
    """
    hits = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            name, contig, identity, length = fields[0], fields[1], float(fields[2]), int(fields[3])
            qstart, qend = int(fields[4]), int(fields[5])
            start = min(int(fields[6]), int(fields[7]))
            end = max(int(fields[6]), int(fields[7]))
            bitscore = float(fields[9])
            contig_length = contig_lengths.get(contig, int(fields[11]))
            if topology.get(contig) == "circular":
                if start > contig_length:
                    continue
                distance = None
                position = "internal"
                long_enough = length >= MIN_INTERNAL_LENGTH
            else:
                distance = min(start - 1, contig_length - end)
                if distance <= END_DISTANCE:
                    position = "end"
                    long_enough = length >= MIN_END_LENGTH
                else:
                    position = "internal"
                    long_enough = length >= MIN_INTERNAL_LENGTH
            if identity < MIN_IDENTITY or not long_enough:
                continue
            hits.append({
                "contig": contig, "start": start, "end": end, "distance": distance,
                "position": position, "length": length, "identity": identity,
                "bitscore": bitscore, "name": name,
                "mu": is_mu_end_hit(name, qstart, qend, intervals),
            })
    return drop_start_copies_of_cross_cut_hits(hits, contig_lengths, topology)


def drop_start_copies_of_cross_cut_hits(hits, contig_lengths, topology):
    """On a circular contig, drop the part of an adapter across the cut that BLAST reports twice.

    Takes in:  the kept hits, the contig lengths and the topology.
    Does:      an adapter that lies across the cut of a circular contig is found whole on the
               extended copy (from before the cut to past it), and BLAST also reports its second
               part where it really sits, at the very start of the contig. That partial hit is
               dropped when, moved by one contig length, it overlaps (or lies within MERGE_GAP bp
               of) a hit across the cut on the same contig; a hit at the start that no cross-cut
               hit covers is kept.
    Produces:  the hits without those copies.
    Used by:   read_hits().
    """
    cross_cut = {}
    for hit in hits:
        contig_length = contig_lengths.get(hit["contig"])
        if topology.get(hit["contig"]) == "circular" and contig_length and hit["start"] <= contig_length < hit["end"]:
            cross_cut.setdefault(hit["contig"], []).append(hit)
    kept = []
    for hit in hits:
        covered = False
        if hit["contig"] in cross_cut and hit["end"] <= CIRCULAR_EXTENSION:
            contig_length = contig_lengths[hit["contig"]]
            moved_start, moved_end = hit["start"] + contig_length, hit["end"] + contig_length
            covered = any(other["start"] - MERGE_GAP <= moved_end and moved_start <= other["end"] + MERGE_GAP
                          for other in cross_cut[hit["contig"]])
        if not covered:
            kept.append(hit)
    return kept


def merge_loci(hits, contig_lengths, topology):
    """Merge overlapping hits on one contig into loci, and decide what each locus counts as.

    Takes in:  the kept hits, the contig lengths and the topology.
    Does:      sorts the hits by contig (in genome order) and start, and joins a hit to the
               locus before it when it starts at most MERGE_GAP bp after that locus ends. A
               locus is an end locus when any of its hits is an end hit, and a Mu-end locus when
               all of its hits are Mu-end hits. Then it sets counted and reason (see the opening
               text) and maps a circular locus's end back onto the contig.
    Produces:  a list of loci (dicts with the LOCUS_COLUMNS, plus "side" for an end locus).
    Used by:   report().
    """
    order = {name: index for index, name in enumerate(contig_lengths)}
    hits = sorted(hits, key=lambda hit: (order.get(hit["contig"], len(order)), hit["start"]))
    loci = []
    for hit in hits:
        if loci and loci[-1]["contig"] == hit["contig"] and hit["start"] <= loci[-1]["end"] + MERGE_GAP:
            locus = loci[-1]
            locus["end"] = max(locus["end"], hit["end"])
            if hit["distance"] is not None:
                locus["distance"] = min(locus["distance"], hit["distance"])
            if hit["position"] == "end":
                locus["position"] = "end"
            if hit["length"] > locus["longest_hit"]:
                locus["longest_hit"] = hit["length"]
                locus["identity"] = hit["identity"]
            locus["mu"] = locus["mu"] and hit["mu"]
            locus["scores"][hit["name"]] = max(locus["scores"].get(hit["name"], 0.0), hit["bitscore"])
        else:
            loci.append({
                "contig": hit["contig"], "start": hit["start"], "end": hit["end"],
                "distance": hit["distance"], "position": hit["position"],
                "longest_hit": hit["length"], "identity": hit["identity"], "mu": hit["mu"],
                "scores": {hit["name"]: hit["bitscore"]},
            })

    for locus in loci:
        contig_length = contig_lengths.get(locus["contig"], "NA")
        locus["contig_length"] = contig_length
        locus["topology"] = topology.get(locus["contig"], "not known")
        # A circular locus across the cut ends in the appended copy: map the end back.
        if locus["topology"] == "circular" and contig_length != "NA" and locus["end"] > contig_length:
            locus["end"] = locus["end"] - contig_length
        locus["distance_to_end"] = "NA" if locus["distance"] is None else locus["distance"]
        # The adapters that matched, best hit first (highest bit score), at most MAX_NAMES.
        names = sorted(locus.pop("scores").items(), key=lambda item: (-item[1], item[0]))
        listed = ",".join(name for name, score in names[:MAX_NAMES])
        if len(names) > MAX_NAMES:
            listed += f" and {len(names) - MAX_NAMES} more"
        locus["adapters"] = listed
        locus["identity"] = f"{locus['identity']:.1f}"
        if locus["position"] == "end":
            from_start = locus["start"] - 1
            from_end = contig_length - locus["end"] if contig_length != "NA" else from_start + 1
            locus["side"] = "start" if from_start <= from_end else "end"
            locus["counted"] = "yes"
            locus["reason"] = "mu_end_at_contig_end" if locus["mu"] else "adapter_at_contig_end"
        else:
            locus["side"] = ""
            if locus["mu"]:
                locus["counted"] = "no"
                locus["reason"] = "mu_end_inside_contig"
            else:
                locus["counted"] = "yes"
                locus["reason"] = "adapter_inside_contig"
    return loci


def contig_ends_with_adapter(loci):
    """Count the contig ends that carry adapter: each (contig, start or end) pair once.

    Two separate loci at one end, or a short contig with a locus near each end, are then counted
    correctly. Used by write_mqc() and report().
    """
    return len({(locus["contig"], locus["side"]) for locus in loci
                if locus["position"] == "end" and locus["counted"] == "yes"})


def counts(loci):
    """(contig ends with adapter, counted loci inside contigs, Mu ends inside contigs)."""
    internal = sum(1 for locus in loci if locus["position"] == "internal" and locus["counted"] == "yes")
    mu_inside = sum(1 for locus in loci if locus["reason"] == "mu_end_inside_contig")
    return contig_ends_with_adapter(loci), internal, mu_inside


def verdict_text(end_loci, internal_loci, mu_inside):
    """The verdict column of the MultiQC table and the console line."""
    parts = []
    if end_loci:
        parts.append(f"adapter at {end_loci} contig end{'s' if end_loci != 1 else ''}")
    if internal_loci:
        parts.append(f"{internal_loci} adapter match{'es' if internal_loci != 1 else ''} inside contigs to check")
    if not parts:
        parts.append("none of the listed adapters found")
    if mu_inside:
        parts.append(f"{mu_inside} Mu end{'s' if mu_inside != 1 else ''} inside contigs, as a Mu prophage holds (not counted)")
    return "; ".join(parts)


def read_key_values(path):
    """Read "name<TAB>value" lines (the read-level files of 05_ont_adapters.smk) as a dict."""
    values = {}
    if not path:
        return values
    with open(path) as handle:
        for line in handle:
            name, _, value = line.rstrip("\n").partition("\t")
            values[name] = value
    return values


def write_loci(path, loci):
    with open(path, "w") as out:
        out.write("\t".join(LOCUS_COLUMNS) + "\n")
        for locus in loci:
            out.write("\t".join(str(locus[column]) for column in LOCUS_COLUMNS) + "\n")


def write_mqc(path, sample, loci, choice, check):
    """Write the sample's row of the Adapters table, behind MultiQC's custom-content header lines.

    Rule multiqc lists the file among its inputs; the rows of all samples merge into one table
    because they share the id below, and MultiQC shows each row under the sample name written in
    its first column.
    """
    end_loci, internal_loci, mu_inside = counts(loci)

    def from_file(values, key):
        # A read-level number, or NA for a sample whose reads BacFlux did not trim.
        return values.get(key, "NA") if values else "NA"

    trimmed = choice.get("summary", "see the adapter choice file") if choice else "reads not trimmed by BacFlux"
    with open(path, "w") as out:
        out.write("# id: 'bacflux_adapters'\n")
        out.write("# section_name: 'Adapters'\n")
        out.write(
            "# description: 'ONT adapter and barcode trimming of the reads (ONT samples, when "
            "BacFlux trimmed them) and the search of the finished genome for ONT, PacBio, MGI and "
            "Illumina adapters and barcodes. An end hit lies within 100 bp of the end of a linear "
            "contig; a circular contig has no ends. The lists cover the common adapters only: before "
            "a GenBank submission, also run NCBI FCS-adaptor.'\n"
        )
        out.write("# plot_type: 'table'\n")
        out.write("# pconfig: {namespace: 'Adapters'}\n")
        out.write("\t".join([
            "Sample", "Reads with adapter (%)", "Reads with barcode, no adapter (%)",
            "ONT barcodes in the reads", "Second barcode", "Trimmed from the reads",
            "Adapter left after trimming (%)", "Barcode end left after trimming (%)",
            "Contig ends with adapter", "Adapter matches inside contigs", "Verdict",
        ]) + "\n")
        out.write("\t".join([
            sample,
            from_file(choice, "percent_with_adapter"),
            from_file(choice, "percent_with_barcode_without_adapter"),
            from_file(choice, "published_barcodes"),
            from_file(choice, "second_barcode"),
            trimmed,
            from_file(check, "percent_with_adapter"),
            from_file(check, "percent_with_end_of_trimmed_sequence"),
            str(end_loci), str(internal_loci),
            verdict_text(end_loci, internal_loci, mu_inside),
        ]) + "\n")


def report(args):
    """Turn the BLAST hits into the locus table and the MultiQC row.

    Takes in:  args.blast, args.contigs, args.replicons, args.query, args.sample,
               args.read_choice, args.read_check, args.out, args.out_mqc.
    Produces:  the two tables and one [BacFlux] line with the verdict.
    """
    contig_lengths = {name: len(sequence) for name, (header, sequence) in read_fasta(args.contigs).items()}
    topology = read_topology(args.replicons)
    intervals = mu_end_intervals(read_fasta(args.query))
    hits = read_hits(args.blast, contig_lengths, topology, intervals)
    choice = read_key_values(args.read_choice)
    check = read_key_values(args.read_check)
    loci = merge_loci(hits, contig_lengths, topology)
    write_loci(args.out, loci)
    write_mqc(args.out_mqc, args.sample, loci, choice, check)
    end_loci, internal_loci, mu_inside = counts(loci)
    print(
        f"[BacFlux] screen_adapters {args.sample}: {verdict_text(end_loci, internal_loci, mu_inside)} "
        f"({len(contig_lengths)} contigs searched)."
    )


def build_parser():
    parser = argparse.ArgumentParser(description="Search the finished genome for adapters, and report them.")
    subcommands = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subcommands.add_parser("prepare", help="write the files blastn reads")
    prepare_parser.add_argument("--contigs", required=True)
    prepare_parser.add_argument("--replicons", nargs="?", const="", default="")
    prepare_parser.add_argument("--lists", nargs="+", required=True)
    prepare_parser.add_argument("--trimmed", nargs="?", const="", default="")
    prepare_parser.add_argument("--out-query", required=True)
    prepare_parser.add_argument("--out-subject", required=True)

    report_parser = subcommands.add_parser("report", help="turn the BLAST hits into the tables")
    report_parser.add_argument("--blast", required=True)
    report_parser.add_argument("--contigs", required=True)
    report_parser.add_argument("--replicons", nargs="?", const="", default="")
    report_parser.add_argument("--query", required=True)
    report_parser.add_argument("--sample", required=True)
    report_parser.add_argument("--read-choice", nargs="?", const="", default="")
    report_parser.add_argument("--read-check", nargs="?", const="", default="")
    report_parser.add_argument("--out", required=True)
    report_parser.add_argument("--out-mqc", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "prepare":
            prepare(args)
        else:
            report(args)
    except (OSError, ValueError, IndexError) as error:
        sys.exit(f"[BacFlux] adapter_screen {args.command}: cannot read an input file: {error}")


if __name__ == "__main__":
    main()
